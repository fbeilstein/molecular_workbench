"""
Mapping-driven transition-state builder & Multi-strategy TS optimizers.

1. ReactionModel:
   - Tracks shared atom index space, R->P mapping.
   - Categorizes bonds into:
       * Formed bonds (new covalent connections)
       * Broken bonds (severed covalent connections)
       * Order-changing bonds (e.g. C=C -> C-C, pi-bond breakage in electrophilic additions)
       * Spectator bonds (inviolable bonds whose connectivity and order do NOT change)
   - Computes chemical attack vectors, charge, and spin multiplicity.

2. Fragment Shape Refinement (Molecular Explorer methods):
   - MMFF94 force-field clean-up
   - GFN2-xTB quantum geometry optimization
   - In-place Kabsch alignment so docked orientation and center-of-mass are preserved.

3. TS Optimization Strategies:
   - Strategy 'dimer': Dimer method (ase.mep.dimer) for crude manual starts (force-based mode search).
   - Strategy 'drive': Constrained coordinate drive across reaction coordinate -> Sella polish.
   - Strategy 'sella': Constrained pre-relaxation of spectators -> Seeded Sella (v0 along reaction vector).
   - Strategy 'neb': Climbing Image Nudged Elastic Band (CI-NEB) with IDPP.

4. Safety & Geometry Sanity:
   - Spectator bond guards prevent unphysical C-H or C-C tearing.
   - Pre- and post-optimization sanity checks detect divergence without corrupting the user scene.
"""
import os
import tempfile
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Geometry import Point3D
from scipy.spatial.transform import Rotation

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


def align_to_reference(p_orig, p_new):
    """Rigid Kabsch alignment: rotate and translate p_new to best fit p_orig."""
    c_orig = np.mean(p_orig, axis=0)
    c_new = np.mean(p_new, axis=0)
    p0 = p_orig - c_orig
    p1 = p_new - c_new
    if len(p0) < 2:
        return p_new - c_new + c_orig
    try:
        rot, _ = Rotation.align_vectors(p0, p1)
        return rot.apply(p1) + c_orig
    except Exception:
        return p_new - c_new + c_orig


# ---------------------------------------------------------------------------
# 1. Reaction model with Bond Order Tracking
# ---------------------------------------------------------------------------
class ReactionModel:
    def __init__(self, raw_smiles, mapping_json):
        r_part, p_part = raw_smiles.split(">>")
        self.r_smiles = [s.strip() for s in r_part.split(".") if s.strip()]
        self.p_smiles = [s.strip() for s in p_part.split(".") if s.strip()]

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

        # Reactant graph with bond orders
        self.r_bonds = set()
        self.r_bond_orders = {}
        self.neighbors = [[] for _ in range(self.n_atoms)]
        for fi, m in enumerate(self.r_mols):
            o = self.offsets[fi]
            for b in m.GetBonds():
                i, j = o + b.GetBeginAtomIdx(), o + b.GetEndAtomIdx()
                pair = (min(i, j), max(i, j))
                self.r_bonds.add(pair)
                self.r_bond_orders[pair] = float(b.GetBondTypeAsDouble())
                self.neighbors[i].append(j)
                self.neighbors[j].append(i)

        # Mapping dictionary
        r2p = {}
        for k, v in (mapping_json or {}).items():
            try:
                _, rfi, rai = k.split("_")
                _, pfi, pai = v.split("_")
                ci = self.offsets[int(rfi)] + int(rai)
                r2p[ci] = (int(pfi), int(pai))
            except (ValueError, IndexError):
                continue

        seen, clean = set(), {}
        for ci in sorted(r2p):
            if r2p[ci] not in seen:
                seen.add(r2p[ci])
                clean[ci] = r2p[ci]
        self.n_explicit = len(clean)
        self.r2p = self._complete_mapping(clean)
        self.p2r = {v: k for k, v in self.r2p.items()}

        # Product graph with bond orders pulled back to reactant atom indices
        self.p_bonds = set()
        self.p_bond_orders = {}
        for pfi, m in enumerate(self.p_mols):
            for b in m.GetBonds():
                i = self.p2r.get((pfi, b.GetBeginAtomIdx()))
                j = self.p2r.get((pfi, b.GetEndAtomIdx()))
                if i is not None and j is not None:
                    pair = (min(i, j), max(i, j))
                    self.p_bonds.add(pair)
                    self.p_bond_orders[pair] = float(b.GetBondTypeAsDouble())

        mapped = set(self.r2p)
        self.formed = sorted(b for b in self.p_bonds - self.r_bonds)
        self.broken = sorted(b for b in self.r_bonds - self.p_bonds
                             if b[0] in mapped and b[1] in mapped)

        # Detect order-changing bonds (e.g. C=C -> C-C, pi-bond breaking)
        self.order_changed = []
        for pair in self.r_bonds.intersection(self.p_bonds):
            r_ord = self.r_bond_orders.get(pair, 1.0)
            p_ord = self.p_bond_orders.get(pair, 1.0)
            if abs(r_ord - p_ord) > 0.4:
                self.order_changed.append((pair[0], pair[1], r_ord, p_ord))

        # Spectator bonds: bonds that exist in both and DO NOT change bond order
        changing_pairs = set(self.formed).union(set(self.broken)).union({(i, j) for i, j, _, _ in self.order_changed})
        self.spectator_bonds = sorted([b for b in self.r_bonds if b not in changing_pairs])

        self.charge = int(sum(a.GetFormalCharge() for m in self.r_mols for a in m.GetAtoms()))
        n_rad = int(sum(a.GetNumRadicalElectrons() for m in self.r_mols for a in m.GetAtoms()))
        self.multiplicity = n_rad + 1
        self.overridden = False

    def _complete_mapping(self, r2p):
        r2p = dict(r2p)
        used = set(r2p.values())
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
        for i, j, ro, po in self.order_changed:
            lines.append(f"  PI/ORDER {self.label(i)}--{self.label(j)}: order {ro:.1f} -> {po:.1f}")
        if not self.formed and not self.broken and not self.order_changed:
            lines.append("  (no bond changes detected - check the mapping)")
        return lines

    def reaction_info(self):
        return {
            "forming": [[list(self.local(i)), list(self.local(j))] for i, j in self.formed],
            "breaking": [[list(self.local(i)), list(self.local(j))] for i, j in self.broken],
            "order_changed": [[[list(self.local(i)), list(self.local(j))], ro, po]
                              for i, j, ro, po in self.order_changed],
            "summary": self.summary_lines(),
            "charge": self.charge,
            "multiplicity": self.multiplicity,
        }

    def apply_override(self, override):
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
        changing = set(self.formed).union(set(self.broken)).union({(i, j) for i, j, _, _ in self.order_changed})
        self.spectator_bonds = sorted([b for b in self.r_bonds if b not in changing])
        self.overridden = True
        return True

    def attack_mode(self, ci, coords_of):
        x = coords_of(ci)
        nbrs = self.neighbors[ci]
        leaving = [j for j in nbrs if (min(ci, j), max(ci, j)) in self.broken]
        if leaving:
            v = -sum(_unit(coords_of(j) - x) for j in leaving)
            if np.linalg.norm(v) > 1e-3:
                return ("vec", _unit(v))
        if not nbrs:
            return None
        units = [_unit(coords_of(j) - x) for j in nbrs]
        s = np.sum(units, axis=0)
        if len(nbrs) == 1 or np.linalg.norm(s) > 0.35:
            return ("vec", _unit(-s))
        if len(nbrs) == 2:
            return ("perp", units[0])
        if len(nbrs) == 3:
            n = np.cross(units[1] - units[0], units[2] - units[0])
            if np.linalg.norm(n) > 1e-6:
                return ("axis", _unit(n))
        return None

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

    def dock(self, frag_coords, variant=0, n_starts=32, seed=7):
        from scipy.optimize import minimize
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
            for f in movers:
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
# 2. Fragment In-Place Refinement (Molecular Explorer methods)
# ---------------------------------------------------------------------------
def set_positions(mol, xyz):
    mol = Chem.Mol(mol)
    if mol.GetNumConformers() == 0:
        mol.AddConformer(Chem.Conformer(mol.GetNumAtoms()), assignId=True)
    conf = mol.GetConformer()
    for i, p in enumerate(xyz):
        conf.SetAtomPosition(i, Point3D(float(p[0]), float(p[1]), float(p[2])))
    return mol


def refine_fragment_shape(mol_block, method="mmff", charge=0, mult=1):
    """Refine isolated molecule geometry in-place using MMFF94 or GFN2-xTB."""
    mol = Chem.MolFromMolBlock(mol_block, removeHs=False)
    if mol is None:
        raise ValueError("Invalid molblock")
    p_orig = mol.GetConformer().GetPositions().copy()

    if method.lower() == "mmff":
        AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
        p_new = mol.GetConformer().GetPositions()
    elif method.lower() == "xtb":
        from ase import Atoms
        from ase.optimize import LBFGS
        from tblite.ase import TBLite
        symbols = [a.GetSymbol() for a in mol.GetAtoms()]
        frag_atoms = Atoms(symbols=symbols, positions=p_orig.copy())
        frag_atoms.calc = TBLite(method="GFN2-xTB", charge=charge, multiplicity=mult, verbosity=0)
        opt = LBFGS(frag_atoms, logfile=None, maxstep=0.2)
        opt.run(fmax=0.04, steps=100)
        p_new = frag_atoms.get_positions()
    else:
        raise ValueError(f"Unknown refinement method: {method}")

    # Align back to current pose so docking orientation is preserved
    p_aligned = align_to_reference(p_orig, p_new)
    return Chem.MolToMolBlock(set_positions(mol, p_aligned))


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
# 3. Calculators, Spectator Restraints & Sanity Checks
# ---------------------------------------------------------------------------
def make_xtb(charge, multiplicity):
    from tblite.ase import TBLite
    return TBLite(method="GFN2-xTB", charge=charge, multiplicity=multiplicity, verbosity=0,
                  electronic_temperature=1000.0, max_iterations=600, mixer_damping=0.2)


def _make_restrained_calc(base_atoms, restraints, k, spectator_restraints=None, k_spec=50.0):
    """Calculator adding harmonic restraints for reaction coordinates and spectator protection."""
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

            # Active reaction coordinate restraints
            for i, j, t in restraints:
                r = pos[j] - pos[i]
                d = np.linalg.norm(r)
                u = r / d
                dd = d - t
                er += 0.5 * k * dd * dd
                f[j] -= k * dd * u
                f[i] += k * dd * u

            # Inviolable spectator bond protection
            if spectator_restraints:
                for i, j, t0 in spectator_restraints:
                    r = pos[j] - pos[i]
                    d = np.linalg.norm(r)
                    u = r / d
                    # Soft half-harmonic: only penalize elongation > +0.2 A
                    if d > t0 + 0.15:
                        dd = d - (t0 + 0.15)
                        er += 0.5 * k_spec * dd * dd
                        f[j] -= k_spec * dd * u
                        f[i] += k_spec * dd * u

            self.e_base = e
            self.results = {"energy": e + er, "forces": f}

    return RestrainedCalculator()


def check_geometry_sanity(atoms, model, log, tolerance=0.35):
    """Verify that non-reacting spectator bonds were not stretched into unphysical states."""
    pos = atoms.get_positions()
    bad_bonds = []
    for i, j in model.spectator_bonds:
        d = float(np.linalg.norm(pos[i] - pos[j]))
        d_ideal = model.d0(i, j)
        if d > d_ideal + tolerance:
            bad_bonds.append((i, j, d, d_ideal))
    if bad_bonds:
        log("\n⚠️ SANITY CHECK FAILED: The following spectator bond(s) were severely distorted:")
        for i, j, d, di in bad_bonds:
            log(f"   {model.label(i)}-{model.label(j)}: length {d:.2f} A (expected ~{di:.2f} A)")
        return False
    return True


def build_reaction_eigenvector(atoms, model):
    """Construct an explicit reaction vector (compressing formed, stretching broken/order-changed)."""
    pos = atoms.get_positions()
    v = np.zeros_like(pos)
    for i, j in model.formed:
        u = _unit(pos[j] - pos[i])
        v[i] += u
        v[j] -= u
    for i, j in model.broken:
        u = _unit(pos[j] - pos[i])
        v[i] -= u
        v[j] += u
    for i, j, _, _ in model.order_changed:
        u = _unit(pos[j] - pos[i])
        v[i] -= u * 0.5
        v[j] += u * 0.5
    norm = np.linalg.norm(v)
    return v / norm if norm > 1e-6 else None


# ---------------------------------------------------------------------------
# 4. Multi-Strategy Optimizers for Manual Guesses
# ---------------------------------------------------------------------------

def run_dimer_ts(atoms, model, log, steps=100, fmax=0.04):
    """Strategy 1: Dimer Method (Robust saddle search from crude manual guess)."""
    from ase.mep.dimer import DimerControl, MinModeAtoms, MinModeTranslate
    from ase.optimize import LBFGS
    from sella import Sella
    log("\n=== Strategy: Dimer Method (ase.mep.dimer) ===")
    log("Dynamically searching and following lowest-curvature reaction mode using force replicas...")

    # If the user placed fragments far apart (> d0 + 1.0 A), pre-contract forming bonds
    # into the interaction zone (~d0 + 0.65 A) so Dimer does not wander into the dissociation limit
    pos = atoms.get_positions()
    far_forming = [pair for pair in model.formed if np.linalg.norm(pos[pair[0]] - pos[pair[1]]) > model.d0(pair[0], pair[1]) + 1.0]
    if far_forming:
        log("  Initial manual placement has large forming distance; pre-contracting into the TS zone (~2.2 A)...")
        base = atoms.copy()
        base.calc = make_xtb(model.charge, model.multiplicity)
        targets = [(i, j, model.d0(i, j) + 0.65) for i, j in model.formed]
        spec = [(i, j, float(np.linalg.norm(pos[i] - pos[j]))) for i, j in model.spectator_bonds]
        calc = _make_restrained_calc(base, targets, 35.0, spectator_restraints=spec)
        atoms.calc = calc
        pre_opt = LBFGS(atoms, logfile=None, maxstep=0.2)
        pre_opt.run(fmax=0.1, steps=25)
        atoms.set_positions(atoms.get_positions())
        log(f"  Pre-contraction complete. Initializing Dimer saddle search...")

    atoms.calc = make_xtb(model.charge, model.multiplicity)
    d_control = DimerControl(
        initial_eigenmode_method='displacement',
        displacement_method='vector',
        dimer_separation=0.005,
        maximum_translation=0.08,
        trial_trans_step=0.005,
        f_rot_min=0.05,
        f_rot_max=0.5,
        max_num_rot=3
    )
    d_atoms = MinModeAtoms(atoms, d_control)

    v0 = build_reaction_eigenvector(atoms, model)
    if v0 is not None:
        try:
            d_atoms.initialize_eigenmodes(eigenmodes=[v0])
            log("  Seeded dimer eigenmode along forming/breaking bonds.")
        except Exception:
            pass

    opt = MinModeTranslate(d_atoms, logfile=None)

    def _log_dimer():
        try:
            f = float(np.max(np.linalg.norm(d_atoms.get_forces(), axis=1)))
            e = atoms.get_potential_energy()
            if opt.nsteps % 5 == 0 or f < fmax:
                log(f"  Dimer step {opt.nsteps:3d}: E = {e:.4f} eV, max|F| = {f:.4f} eV/A")
        except Exception:
            pass

    opt.attach(_log_dimer)
    try:
        opt.run(fmax=fmax, steps=steps)
        log("  Dimer stage finished.")
    except Exception as e:
        log(f"  Dimer search note: {e}")

    # Quick polish with Sella
    try:
        log("  Polishing Dimer saddle point with Sella...")
        s_opt = Sella(atoms, internal=True, order=1, logfile=None)
        s_opt.run(fmax=0.03, steps=30)
    except Exception as e:
        log(f"  Sella polish note: {e} (using Dimer saddle geometry)")
    return True


def run_drive_sella_ts(atoms, model, log, n_points=11, k=25.0):
    """Strategy 2: Constrained Coordinate Drive across barrier peak -> Sella."""
    from ase.optimize import LBFGS
    from sella import Sella
    log("\n=== Strategy: Constrained Coordinate Drive -> Sella ===")
    pairs = model.formed + model.broken
    pos = atoms.get_positions()
    d_cur = np.array([np.linalg.norm(pos[i] - pos[j]) for i, j in pairs])

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

    spec_restraints = [(i, j, float(np.linalg.norm(pos[i] - pos[j]))) for i, j in model.spectator_bonds]

    base = atoms.copy()
    base.calc = make_xtb(model.charge, model.multiplicity)

    frames, energies, lams = [], [], []
    i_max = 0
    for n, lam in enumerate(np.linspace(0.0, 1.0, n_points)):
        targets = d_start + lam * (d_end - d_start)
        calc = _make_restrained_calc(base, [(i, j, t) for (i, j), t in zip(pairs, targets)], k,
                                     spectator_restraints=spec_restraints)
        atoms.calc = calc
        opt = LBFGS(atoms, logfile=None, maxstep=0.15)
        opt.run(fmax=0.08, steps=100)
        e = calc.e_base
        frames.append(atoms.get_positions().copy())
        energies.append(e)
        lams.append(float(lam))
        if e > energies[i_max]:
            i_max = n
        rel = (e - energies[0]) * EV_TO_KCAL
        log(f"  lambda={lam:4.2f}  dE={rel:+7.2f} kcal/mol")
        # Only stop early if we have truly climbed a barrier and descended past it
        if i_max >= 2 and (energies[i_max] - energies[0]) * EV_TO_KCAL > 3.0:
            if (energies[i_max] - e) * EV_TO_KCAL > 6.0:
                log("  (passed barrier maximum - stopping scan early)")
                break

    atoms.set_positions(frames[i_max])
    log(f"  Barrier peak selected at lambda={lams[i_max]:.2f}. Refining with Sella...")
    atoms.calc = make_xtb(model.charge, model.multiplicity)
    try:
        opt = Sella(atoms, internal=True, order=1, logfile=None)
        opt.run(fmax=0.03, steps=120)
    except Exception as e:
        log(f"  Sella polish note: {e}")
    profile = [[l, (e - energies[0]) * EV_TO_KCAL] for l, e in zip(lams, energies)]
    return profile


def run_seeded_sella_ts(atoms, model, log):
    """Strategy 3: Constrained Pre-Relaxation -> Sella."""
    from ase.optimize import LBFGS
    from sella import Sella
    log("\n=== Strategy: Pre-Relaxation -> Sella ===")
    pairs = model.formed + model.broken
    pos = atoms.get_positions()
    restraints = [(i, j, float(np.linalg.norm(pos[i] - pos[j]))) for i, j in pairs]
    spec_restraints = [(i, j, float(np.linalg.norm(pos[i] - pos[j]))) for i, j in model.spectator_bonds]

    base = atoms.copy()
    base.calc = make_xtb(model.charge, model.multiplicity)
    atoms.calc = _make_restrained_calc(base, restraints, 25.0, spectator_restraints=spec_restraints)
    opt = LBFGS(atoms, logfile=None, maxstep=0.15)
    opt.run(fmax=0.06, steps=120)
    log(f"  Pre-relaxation complete. Relaxed energy = {atoms.calc.e_base:.4f} eV")

    atoms.calc = make_xtb(model.charge, model.multiplicity)
    try:
        opt = Sella(atoms, internal=True, order=1, logfile=None)
        opt.run(fmax=0.03, steps=150)
    except Exception as e:
        log(f"  Sella optimization note: {e}")
    return []


def run_neb_ts(atoms, model, log, n_images=7):
    """Strategy 4: Climbing Image Nudged Elastic Band (CI-NEB with IDPP)."""
    from ase.optimize import LBFGS, MDMin
    from ase.mep import NEB
    from ase.mep.neb import idpp_interpolate
    from sella import Sella
    log("\n=== Strategy: Climbing Image NEB (CI-NEB) ===")

    # 1. Product endpoint P generated cleanly via restrained relaxation
    p_atoms = atoms.copy()
    p_targets = [(i, j, model.d0(i, j)) for i, j in model.formed]
    p_targets += [(i, j, model.d0(i, j) + 1.25) for i, j in model.broken]
    spec_restraints = [(i, j, float(np.linalg.norm(atoms.positions[i] - atoms.positions[j]))) for i, j in model.spectator_bonds]

    p_base = p_atoms.copy()
    p_base.calc = make_xtb(model.charge, model.multiplicity)
    p_atoms.calc = _make_restrained_calc(p_base, p_targets, k=30.0, spectator_restraints=spec_restraints)
    p_opt = LBFGS(p_atoms, logfile=None, maxstep=0.2)
    p_opt.run(fmax=0.08, steps=50)

    # 2. Build NEB band between current manual pose (R) and relaxed product (P)
    images = [atoms.copy() for _ in range(n_images)]
    images[0].set_positions(atoms.get_positions())
    images[-1].set_positions(p_atoms.get_positions())
    idpp_interpolate(images, mic=False)

    for img in images:
        img.calc = make_xtb(model.charge, model.multiplicity)

    neb = NEB(images, climb=True, method='improvedtangent')
    opt = MDMin(neb, logfile=None)
    opt.run(fmax=0.10, steps=80)

    # 3. Intermediate image with highest energy is the TS guess
    energies = [img.get_potential_energy() for img in images]
    # Restrict saddle search to intermediate images (1 to n_images-2)
    mid_idx = int(np.argmax(energies[1:-1])) + 1
    atoms.set_positions(images[mid_idx].get_positions())
    log(f"  CI-NEB finished. Image {mid_idx} selected as barrier top (E = {energies[mid_idx]:.4f} eV). Refining with Sella...")

    atoms.calc = make_xtb(model.charge, model.multiplicity)
    try:
        s_opt = Sella(atoms, internal=True, order=1, logfile=None)
        s_opt.run(fmax=0.03, steps=100)
    except Exception as e:
        log(f"  Sella polish note: {e}")
    profile = [[float(idx), float((e - energies[0]) * EV_TO_KCAL)] for idx, e in enumerate(energies)]
    return profile


# ---------------------------------------------------------------------------
# 5. Vibrational Frequency Verification
# ---------------------------------------------------------------------------
def analyse_mode(atoms, model, log, max_atoms=60):
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
