"""
Mapping-driven transition-state builder.

Everything here is derived from the reaction SMILES + Step-2 atom mapping,
so the user does not have to place molecules by hand:

1. ReactionModel   - one shared atom index space (reactant complex order),
                     completed R->P bijection (incl. transferred H),
                     formed / broken bonds, total charge, multiplicity.
2. attack_mode()   - chemically motivated approach direction per reactive
                     atom (backside of a leaving group, lone-pair side,
                     normal of an sp2 plane, D-H axis for proton transfer).
3. dock()          - rigid-body placement of reactant fragments that
                     satisfies all inter-fragment forming bonds + sterics.
4. bond_drive()    - synchronous restrained xTB scan R -> P in the SAME
                     frame (the product is reached continuously, so there
                     is no R/P orientation/permutation mismatch).
5. refine_ts()     - Sella saddle search from the scan maximum and a
                     frequency check that the imaginary mode is the
                     expected bond reorganisation.

Note on product 3D geometries: independently embedded products are NOT used
as a template - e.g. for SN2 the CH3 umbrella is inverted in the product, so
superimposing on it would put the nucleophile on the leaving-group side.
The product *graph* is what carries the information.
"""
import os
import tempfile

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Geometry import Point3D

_PT = Chem.GetPeriodicTable()
EV_TO_KCAL = 23.060548

DOCK_EXTRA = 1.0      # pre-reaction complex: forming bond = d0 + 1.0 A
BREAK_EXTRA = 1.4     # end of drive: breaking bond = d0 + 1.4 A


def rcov(sym):
    return _PT.GetRcovalent(_PT.GetAtomicNumber(sym))


def rvdw(sym):
    return _PT.GetRvdw(_PT.GetAtomicNumber(sym))


def _unit(v):
    n = np.linalg.norm(v)
    return v / n if n > 1e-9 else v * 0.0


def dir_pen(mode, R, u):
    """Penalty (0 = ideal) for approaching along unit vector u."""
    if mode is None:
        return 0.0
    kind, v = mode
    c = float(np.dot(R @ v, u))
    if kind == "vec":
        return 1.0 - c
    if kind == "axis":
        return 1.0 - c * c
    return c * c  # perp


# ---------------------------------------------------------------------------
# 1. Reaction model
# ---------------------------------------------------------------------------
class ReactionModel:
    def __init__(self, raw_smiles, mapping_json):
        r_part, p_part = raw_smiles.split(">>")
        self.r_smiles = [s.strip() for s in r_part.split(".") if s.strip()]
        self.p_smiles = [s.strip() for s in p_part.split(".") if s.strip()]

        # Same construction as Step 2 / Step 3 so atom indices agree
        self.r_mols = [Chem.AddHs(Chem.MolFromSmiles(s)) for s in self.r_smiles]
        self.p_mols = [Chem.AddHs(Chem.MolFromSmiles(s)) for s in self.p_smiles]
        if any(m is None for m in self.r_mols + self.p_mols):
            raise ValueError("Could not parse reaction SMILES")

        self.offsets, off = [], 0
        for m in self.r_mols:
            self.offsets.append(off)
            off += m.GetNumAtoms()
        self.n_atoms = off
        self.frag_of, self.symbols = [], []
        for fi, m in enumerate(self.r_mols):
            for a in m.GetAtoms():
                self.frag_of.append(fi)
                self.symbols.append(a.GetSymbol())

        # reactant graph in complex indices
        self.r_bonds = set()
        self.neighbors = [[] for _ in range(self.n_atoms)]
        for fi, m in enumerate(self.r_mols):
            o = self.offsets[fi]
            for b in m.GetBonds():
                i, j = o + b.GetBeginAtomIdx(), o + b.GetEndAtomIdx()
                self.r_bonds.add((min(i, j), max(i, j)))
                self.neighbors[i].append(j)
                self.neighbors[j].append(i)

        # mapping
        r2p = {}
        for k, v in (mapping_json or {}).items():
            try:
                _, rfi, rai = k.split("_")
                _, pfi, pai = v.split("_")
                ci = self.offsets[int(rfi)] + int(rai)
                r2p[ci] = (int(pfi), int(pai))
            except (ValueError, IndexError):
                continue
        # drop duplicate targets (keep first)
        seen, clean = set(), {}
        for ci in sorted(r2p):
            if r2p[ci] not in seen:
                seen.add(r2p[ci])
                clean[ci] = r2p[ci]
        self.n_explicit = len(clean)
        self.r2p = self._complete_mapping(clean)
        self.p2r = {v: k for k, v in self.r2p.items()}

        # product graph pulled back to complex indices
        self.p_bonds = set()
        for pfi, m in enumerate(self.p_mols):
            for b in m.GetBonds():
                i = self.p2r.get((pfi, b.GetBeginAtomIdx()))
                j = self.p2r.get((pfi, b.GetEndAtomIdx()))
                if i is not None and j is not None:
                    self.p_bonds.add((min(i, j), max(i, j)))

        mapped = set(self.r2p)
        self.formed = sorted(b for b in self.p_bonds - self.r_bonds)
        self.broken = sorted(b for b in self.r_bonds - self.p_bonds
                             if b[0] in mapped and b[1] in mapped)

        self.charge = int(sum(a.GetFormalCharge() for m in self.r_mols for a in m.GetAtoms()))
        n_rad = int(sum(a.GetNumRadicalElectrons() for m in self.r_mols for a in m.GetAtoms()))
        self.multiplicity = n_rad + 1

    def _complete_mapping(self, r2p):
        """Extend a partial mapping to a bijection (transferred H etc.)."""
        r2p = dict(r2p)
        used = set(r2p.values())
        # pass 1: an unmapped atom stays on the image of its own neighbour
        for ci in range(self.n_atoms):
            if ci in r2p:
                continue
            for nb in self.neighbors[ci]:
                if nb not in r2p:
                    continue
                pfi, pai = r2p[nb]
                patom = self.p_mols[pfi].GetAtomWithIdx(pai)
                for pn in patom.GetNeighbors():
                    key = (pfi, pn.GetIdx())
                    if pn.GetSymbol() == self.symbols[ci] and key not in used:
                        r2p[ci] = key
                        used.add(key)
                        break
                if ci in r2p:
                    break
        # pass 2: any free product atom of the same element
        free = [(pfi, a.GetIdx(), a.GetSymbol())
                for pfi, m in enumerate(self.p_mols) for a in m.GetAtoms()
                if (pfi, a.GetIdx()) not in used]
        for ci in range(self.n_atoms):
            if ci in r2p:
                continue
            for k, (pfi, pai, sym) in enumerate(free):
                if sym == self.symbols[ci]:
                    r2p[ci] = (pfi, pai)
                    free.pop(k)
                    break
        return r2p

    # helpers -------------------------------------------------------------
    def local(self, ci):
        fi = self.frag_of[ci]
        return fi, ci - self.offsets[fi]

    def label(self, ci):
        fi, li = self.local(ci)
        return f"{self.symbols[ci]}(R{fi + 1}:{li})"

    def d0(self, i, j):
        return rcov(self.symbols[i]) + rcov(self.symbols[j])

    def summary_lines(self):
        lines = [f"Charge {self.charge:+d}, multiplicity {self.multiplicity}",
                 f"Mapped atoms: {len(self.r2p)}/{self.n_atoms} "
                 f"({self.n_explicit} from Step 2, rest completed automatically)"]
        for i, j in self.formed:
            lines.append(f"  FORMING  {self.label(i)} ... {self.label(j)}")
        for i, j in self.broken:
            lines.append(f"  BREAKING {self.label(i)} -- {self.label(j)}")
        if not self.formed and not self.broken:
            lines.append("  (no bond changes detected - check the mapping)")
        return lines

    def reaction_info(self):
        return {
            "forming": [[list(self.local(i)), list(self.local(j))] for i, j in self.formed],
            "breaking": [[list(self.local(i)), list(self.local(j))] for i, j in self.broken],
            "summary": self.summary_lines(),
            "charge": self.charge,
            "multiplicity": self.multiplicity,
        }

    def apply_override(self, override):
        """Replace formed/broken bonds by user-picked pairs.

        override = {"forming": [[[fi, ai], [fj, aj]], ...], "breaking": [...]}
        """
        if not override:
            return False

        def conv(pairs):
            out = []
            for (fi, ai), (fj, aj) in pairs:
                i = self.offsets[int(fi)] + int(ai)
                j = self.offsets[int(fj)] + int(aj)
                if i != j and 0 <= i < self.n_atoms and 0 <= j < self.n_atoms:
                    out.append((min(i, j), max(i, j)))
            return sorted(set(out))

        forming = conv(override.get("forming", []))
        breaking = conv(override.get("breaking", []))
        if not forming and not breaking:
            return False
        self.formed, self.broken = forming, breaking
        self.overridden = True
        return True

    def pose_problems(self, positions):
        """Is this a usable pre-reaction arrangement? Returns list of reasons if not."""
        reasons = []
        coords = lambda c: positions[c]
        for i, j in self.formed:
            if self.frag_of[i] == self.frag_of[j]:
                continue
            d = float(np.linalg.norm(positions[j] - positions[i]))
            if d > self.d0(i, j) + DOCK_EXTRA + 1.0:
                reasons.append(f"{self.label(i)}...{self.label(j)} are {d:.1f} A apart")
                continue
            u = _unit(positions[j] - positions[i])
            for a, b, uu in ((i, j, u), (j, i, -u)):
                pen = dir_pen(self.attack_mode(a, coords), np.eye(3), uu)
                if pen > 0.5:  # > 60 deg off the ideal approach
                    reasons.append(f"{self.label(b)} approaches {self.label(a)} from a bad direction")
        return reasons

    def split(self, positions):
        return [positions[o:o + m.GetNumAtoms()] for o, m in zip(self.offsets, self.r_mols)]

    # -------------------------------------------------------------------
    # 2. Attack directions
    # -------------------------------------------------------------------
    def attack_mode(self, ci, coords_of):
        """Return ('vec'|'axis'|'perp', unit vector) or None.

        coords_of(cj) -> xyz of complex atom cj in the frame of ci's fragment.
        """
        x = coords_of(ci)
        nbrs = self.neighbors[ci]
        leaving = [j for j in nbrs if (min(ci, j), max(ci, j)) in self.broken]
        if leaving:  # substitution / transfer: come in anti to what leaves
            v = -sum(_unit(coords_of(j) - x) for j in leaving)
            if np.linalg.norm(v) > 1e-3:
                return ("vec", _unit(v))
        if not nbrs:
            return None
        units = [_unit(coords_of(j) - x) for j in nbrs]
        s = np.sum(units, axis=0)
        if len(nbrs) == 1 or np.linalg.norm(s) > 0.35:  # lone pair side
            return ("vec", _unit(-s))
        if len(nbrs) == 2:  # linear sp centre: attack perpendicular
            return ("perp", units[0])
        if len(nbrs) == 3:  # planar sp2 centre: attack along the normal
            n = np.cross(units[1] - units[0], units[2] - units[0])
            if np.linalg.norm(n) > 1e-6:
                return ("axis", _unit(n))
        return None

    # -------------------------------------------------------------------
    # Fragment embedding (conformer chosen for intramolecular forming bonds)
    # -------------------------------------------------------------------
    def embed_fragments(self, seed=42):
        out = []
        for fi, m0 in enumerate(self.r_mols):
            mol = Chem.Mol(m0)
            o = self.offsets[fi]
            intra = [(i - o, j - o) for i, j in self.formed
                     if self.frag_of[i] == fi and self.frag_of[j] == fi]
            params = AllChem.ETKDGv3()
            params.randomSeed = seed
            nconf = 30 if intra else 1
            cids = list(AllChem.EmbedMultipleConfs(mol, nconf, params))
            if not cids:
                AllChem.Compute2DCoords(mol)
                out.append(mol)
                continue
            try:
                AllChem.MMFFOptimizeMoleculeConfs(mol)
            except Exception:
                pass
            best = cids[0]
            if intra:
                def score(cid):
                    p = mol.GetConformer(cid).GetPositions()
                    return sum(np.linalg.norm(p[i] - p[j]) for i, j in intra)
                best = min(cids, key=score)
            out.append(Chem.Mol(mol, False, best))
        return out

    # -------------------------------------------------------------------
    # 3. Rigid-body docking
    # -------------------------------------------------------------------
    def dock(self, frag_coords, variant=0, n_starts=32, seed=7):
        """frag_coords: list of (n_i,3). Returns (list of (n_i,3), n_variants, score)."""
        from scipy.optimize import minimize
        from scipy.spatial.transform import Rotation

        nf = len(frag_coords)
        if nf < 2:
            return [np.array(c, float) for c in frag_coords], 1, 0.0

        local = [np.array(c, float) - np.mean(c, axis=0) for c in frag_coords]
        anchor = int(np.argmax([len(c) for c in local]))

        def coords_in_frag(fi):
            o = self.offsets[fi]
            return lambda cj: local[fi][cj - o]

        pairs = []
        for i, j in self.formed:
            fi, fj = self.frag_of[i], self.frag_of[j]
            if fi == fj:
                continue
            d = self.d0(i, j) + DOCK_EXTRA
            mi = self.attack_mode(i, coords_in_frag(fi))
            mj = self.attack_mode(j, coords_in_frag(fj))
            pairs.append((fi, i - self.offsets[fi], fj, j - self.offsets[fj], d, mi, mj))

        reactive = {(self.frag_of[a], a - self.offsets[self.frag_of[a]])
                    for b in self.formed for a in b}
        connected = {p[0] for p in pairs} | {p[2] for p in pairs}

        # steric thresholds between fragment pairs
        thr = {}
        for f in range(nf):
            for g in range(f + 1, nf):
                sf = [self.symbols[self.offsets[f] + k] for k in range(len(local[f]))]
                sg = [self.symbols[self.offsets[g] + k] for k in range(len(local[g]))]
                T = 0.75 * (np.array([rvdw(s) for s in sf])[:, None] +
                            np.array([rvdw(s) for s in sg])[None, :])
                for k in range(len(sf)):
                    for l in range(len(sg)):
                        if (f, k) in reactive and (g, l) in reactive:
                            T[k, l] = 0.0
                thr[(f, g)] = T

        movers = [f for f in range(nf) if f != anchor]

        def transforms(p):
            T = {anchor: (np.eye(3), np.zeros(3))}
            for n, f in enumerate(movers):
                R = Rotation.from_rotvec(p[6 * n:6 * n + 3]).as_matrix()
                T[f] = (R, p[6 * n + 3:6 * n + 6])
            return T

        def objective(p):
            T = transforms(p)
            X = [local[f] @ T[f][0].T + T[f][1] for f in range(nf)]
            E = 0.0
            for fi, li, fj, lj, d, mi, mj in pairs:
                r = X[fj][lj] - X[fi][li]
                dist = np.linalg.norm(r) + 1e-9
                u = r / dist
                E += (dist - d) ** 2
                E += d * d * (dir_pen(mi, T[fi][0], u) + dir_pen(mj, T[fj][0], -u))
            for (f, g), Tm in thr.items():
                D = np.linalg.norm(X[f][:, None, :] - X[g][None, :, :], axis=2)
                E += 3.0 * np.sum(np.clip(Tm - D, 0.0, None) ** 2)
            for f in movers:  # keep spectators close but not bonded
                if f not in connected:
                    E += 0.1 * max(0.0, np.linalg.norm(T[f][1]) - 6.0) ** 2
            return E

        rng = np.random.default_rng(seed)
        size = [np.max(np.linalg.norm(c, axis=1)) if len(c) > 1 else 0.0 for c in local]
        results = []
        for _ in range(n_starts):
            p0 = []
            for f in movers:
                rv = Rotation.random(random_state=rng).as_rotvec()
                direction = _unit(rng.normal(size=3))
                t = direction * (size[anchor] + size[f] + 2.5)
                p0.extend(list(rv) + list(t))
            res = minimize(objective, np.array(p0), method="L-BFGS-B",
                           options={"maxiter": 400})
            results.append((res.fun, res.x))
        results.sort(key=lambda r: r[0])

        # distinct poses (for the "next pose" button)
        unique = []
        for fun, x in results:
            T = transforms(x)
            X = np.vstack([local[f] @ T[f][0].T + T[f][1] for f in range(nf)])
            if all(np.sqrt(np.mean(np.sum((X - U) ** 2, axis=1))) > 0.5 for _, _, U in unique):
                unique.append((fun, x, X))
            if len(unique) >= 6:
                break
        fun, x, X = unique[variant % len(unique)]
        X = X - X.mean(axis=0)
        out, k = [], 0
        for f in range(nf):
            n = len(local[f])
            out.append(X[k:k + n])
            k += n
        return out, len(unique), float(fun)


# ---------------------------------------------------------------------------
# Mol / coordinate helpers
# ---------------------------------------------------------------------------
def set_positions(mol, xyz):
    mol = Chem.Mol(mol)
    if mol.GetNumConformers() == 0:
        mol.AddConformer(Chem.Conformer(mol.GetNumAtoms()), assignId=True)
    conf = mol.GetConformer()
    for i, p in enumerate(xyz):
        conf.SetAtomPosition(i, Point3D(float(p[0]), float(p[1]), float(p[2])))
    return mol


def molblocks_to_positions(molblocks):
    mols = [Chem.MolFromMolBlock(mb, removeHs=False) for mb in molblocks]
    if any(m is None for m in mols):
        raise ValueError("Failed to parse molblocks")
    return mols, [m.GetConformer().GetPositions() for m in mols]


def positions_to_molblocks(mols, positions):
    out, k = [], 0
    for m in mols:
        n = m.GetNumAtoms()
        out.append(Chem.MolToMolBlock(set_positions(m, positions[k:k + n])))
        k += n
    return out


def auto_dock_molblocks(model, variant=0):
    mols = model.embed_fragments()
    coords = [m.GetConformer().GetPositions() for m in mols]
    docked, n_var, score = model.dock(coords, variant=variant)
    blocks = [Chem.MolToMolBlock(set_positions(m, x)) for m, x in zip(mols, docked)]
    return blocks, n_var, score


# ---------------------------------------------------------------------------
# 4. Restrained xTB bond drive
# ---------------------------------------------------------------------------
def _make_restrained_calc(base_atoms, restraints, k):
    from ase.calculators.calculator import Calculator, all_changes

    class RestrainedCalculator(Calculator):
        implemented_properties = ["energy", "forces"]

        def __init__(self):
            super().__init__()
            self.e_base = None

        def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
            super().calculate(atoms, properties, system_changes)
            pos = self.atoms.get_positions()
            base_atoms.set_positions(pos)
            e = base_atoms.get_potential_energy()
            f = base_atoms.get_forces().copy()
            er = 0.0
            for i, j, t in restraints:
                r = pos[j] - pos[i]
                d = np.linalg.norm(r)
                u = r / d
                dd = d - t
                er += 0.5 * k * dd * dd
                f[j] -= k * dd * u
                f[i] += k * dd * u
            self.e_base = e
            self.results = {"energy": e + er, "forces": f}

    return RestrainedCalculator()


def make_xtb(charge, multiplicity):
    from tblite.ase import TBLite
    return TBLite(method="GFN2-xTB", charge=charge, multiplicity=multiplicity,
                  verbosity=0)


def bond_drive(atoms, model, log, n_points=13, k=20.0):
    """Synchronous restrained scan of all forming/breaking bonds.

    Returns (frames, energies_eV, lambdas).
    """
    from ase.optimize import LBFGS

    pairs = model.formed + model.broken
    pos = atoms.get_positions()
    d_cur = np.array([np.linalg.norm(pos[i] - pos[j]) for i, j in pairs])
    # Start from reactant-like distances whatever the input geometry is
    # (auto-docked complex, a hand-made TS-like guess, or a previous TS):
    # forming bonds open, breaking bonds intact.  The restraint at lambda=0
    # moves the geometry there first.
    d_start, d_end = [], []
    for (i, j), dc in zip(pairs, d_cur):
        d0 = model.d0(i, j)
        if (i, j) in model.formed:
            d_start.append(min(max(dc, d0 + DOCK_EXTRA), d0 + 2.5))
            d_end.append(d0)
        else:
            d_start.append(min(dc, d0 + 0.1))
            d_end.append(d0 + BREAK_EXTRA)
    d_start, d_end = np.array(d_start), np.array(d_end)

    base = atoms.copy()
    base.calc = make_xtb(model.charge, model.multiplicity)

    frames, energies, lams = [], [], []
    i_max = 0
    for n, lam in enumerate(np.linspace(0.0, 1.0, n_points)):
        targets = d_start + lam * (d_end - d_start)
        calc = _make_restrained_calc(base, [(i, j, t) for (i, j), t in zip(pairs, targets)], k)
        atoms.calc = calc
        opt = LBFGS(atoms, logfile=None, maxstep=0.2)
        opt.run(fmax=0.1, steps=120)
        atoms.get_potential_energy()
        e = calc.e_base
        frames.append(atoms.get_positions().copy())
        energies.append(e)
        lams.append(float(lam))
        if e > energies[i_max]:
            i_max = n
        rel = (e - energies[0]) * EV_TO_KCAL
        bar = "#" * max(0, min(40, int(round(rel / 2.0))))
        dist_txt = ", ".join(f"{np.linalg.norm(atoms.positions[i] - atoms.positions[j]):.2f}"
                             for i, j in pairs)
        log(f"  lambda={lam:4.2f}  dE={rel:+7.2f} kcal/mol  d=[{dist_txt}]  {bar}")
        # stop once we are clearly past the barrier
        if n >= i_max + 2 and (energies[i_max] - e) * EV_TO_KCAL > 4.0:
            log("  (past the barrier - stopping scan early)")
            break
    atoms.calc = None
    return frames, energies, lams


def refine_guess(atoms, model, log, k=20.0):
    """Manual mode: the given geometry IS the TS guess.

    Freeze the reacting bonds at their current lengths (harmonic restraints)
    and relax everything else, so rigid-body docking artefacts (e.g. a CH3
    umbrella that has not inverted yet) are fixed before the saddle search.
    """
    from ase.optimize import LBFGS

    pairs = model.formed + model.broken
    if not pairs:
        log("  no reacting bonds defined - skipping constrained pre-relaxation")
        return
    pos = atoms.get_positions()
    restraints = [(i, j, float(np.linalg.norm(pos[i] - pos[j]))) for i, j in pairs]
    for i, j, d in restraints:
        log(f"  freezing {model.label(i)}-{model.label(j)} at {d:.2f} A")
    base = atoms.copy()
    base.calc = make_xtb(model.charge, model.multiplicity)
    atoms.calc = _make_restrained_calc(base, restraints, k)
    opt = LBFGS(atoms, logfile=None, maxstep=0.2)
    opt.run(fmax=0.08, steps=200)
    log(f"  pre-relaxation done in {opt.nsteps} steps, "
        f"E(xTB) = {atoms.calc.e_base:.4f} eV")
    atoms.calc = None


# ---------------------------------------------------------------------------
# 5. TS refinement + validation
# ---------------------------------------------------------------------------
def analyse_mode(atoms, model, log, max_atoms=60):
    """Finite-difference Hessian; report imaginary modes and their character."""
    if len(atoms) > max_atoms:
        log(f"  (skipping frequency check: {len(atoms)} atoms)")
        return None
    from ase.vibrations import Vibrations

    with tempfile.TemporaryDirectory() as tmp:
        vib = Vibrations(atoms, name=os.path.join(tmp, "vib"), delta=0.005)
        vib.run()
        data = vib.get_vibrations()
        freqs = data.get_frequencies()
        modes = data.get_modes()
    imag = [(k, f.imag) for k, f in enumerate(freqs) if f.imag > 50.0]
    report = {"n_imag": len(imag), "imag_freqs": [round(f, 1) for _, f in imag]}
    if not imag:
        log("  No imaginary frequency > 50i cm-1: this is a MINIMUM, not a TS.")
        return report
    log("  Imaginary frequencies: " + ", ".join(f"{f:.0f}i" for _, f in imag) + " cm-1")
    k0 = imag[int(np.argmax([f for _, f in imag]))][0]
    m = modes[k0]
    pos = atoms.get_positions()

    def stretch(i, j):
        return float(np.dot(_unit(pos[j] - pos[i]), m[j] - m[i]))

    all_bonds = model.r_bonds | model.p_bonds
    total = sum(abs(stretch(i, j)) for i, j in all_bonds) + 1e-12
    ch_f = [stretch(i, j) for i, j in model.formed]
    ch_b = [stretch(i, j) for i, j in model.broken]
    share = (sum(map(abs, ch_f)) + sum(map(abs, ch_b))) / total
    signs_ok = True
    if ch_f and ch_b:
        sf, sb = np.sign(ch_f), np.sign(ch_b)
        signs_ok = bool(np.all(sf == sf[0]) and np.all(sb == -sf[0]))
    report.update({"mode_share": round(share, 2), "mode_consistent": signs_ok})
    log(f"  Main imaginary mode: {share * 100:.0f}% of bond-stretch amplitude on the "
        f"reacting bonds; forming/breaking in antiphase: {'yes' if signs_ok else 'NO'}")
    return report


def refine_ts(atoms, model, log_fh, log, fmax=0.01, steps=250):
    from sella import Sella

    atoms.calc = make_xtb(model.charge, model.multiplicity)
    opt = Sella(atoms, internal=True, order=1, logfile=log_fh)
    converged = opt.run(fmax=fmax, steps=steps)
    return bool(converged)
