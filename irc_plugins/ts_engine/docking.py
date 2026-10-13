import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from scipy.spatial.transform import Rotation
from scipy.optimize import minimize
from .common import rvdw, _unit, dir_pen, set_positions, DOCK_EXTRA


def _select_lowest_energy_conf(mol, cids, solvent="water"):
    if not cids or len(cids) <= 1:
        return cids[0] if cids else 0
    chg = int(sum(a.GetFormalCharge() for a in mol.GetAtoms()))
    if chg == 0:
        try:
            res = AllChem.MMFFOptimizeMoleculeConfs(mol)
            energies = [e for _, e in res]
            if energies and max(energies) - min(energies) > 1e-4:
                return int(np.argmin(energies))
        except Exception:
            pass
    try:
        from ase import Atoms
        from .calculators import make_xtb
        mult = int(sum(a.GetNumRadicalElectrons() for a in mol.GetAtoms())) + 1
        calc = make_xtb(chg, mult, solvent=solvent or "water")
        syms = [a.GetSymbol() for a in mol.GetAtoms()]
        xtb_energies = []
        for cid in cids:
            p = mol.GetConformer(cid).GetPositions()
            at = Atoms(symbols=syms, positions=p)
            at.calc = calc
            xtb_energies.append(float(at.get_potential_energy()))
        return int(np.argmin(xtb_energies))
    except Exception:
        return cids[0]


def _select_product_aligned_conf(mol, fi, model, cids, seed=42):
    if not cids:
        return 0
    if len(cids) == 1:
        return cids[0]
    o = model.offsets[fi]
    n_atoms = len(mol.GetAtoms())
    p_counts = {}
    for loc_i in range(n_atoms):
        pfi, _ = model.r2p.get(o + loc_i, (-1, -1))
        if pfi >= 0:
            p_counts[pfi] = p_counts.get(pfi, 0) + 1
    if p_counts and getattr(model, "p_mols", None):
        best_pfi = max(p_counts, key=p_counts.get)
        if best_pfi < len(model.p_mols):
            pmol = Chem.Mol(model.p_mols[best_pfi])
            AllChem.EmbedMolecule(pmol, randomSeed=seed)
            try:
                AllChem.MMFFOptimizeMolecule(pmol)
            except Exception:
                pass
            p_pos = pmol.GetConformer().GetPositions()
            conserved_r, conserved_p = [], []
            for loc_i in range(n_atoms):
                pfi, p_idx = model.r2p.get(o + loc_i, (-1, -1))
                if pfi == best_pfi:
                    is_breaking = any(o + loc_i in (b1, b2) for b1, b2 in model.broken)
                    if not is_breaking or p_counts[best_pfi] <= 3:
                        conserved_r.append(loc_i)
                        conserved_p.append(p_idx)
            if len(conserved_r) >= 3:
                P_sub = p_pos[conserved_p]
                cP = np.mean(P_sub, axis=0)
                def rmsd_score(cid):
                    r_pos = mol.GetConformer(cid).GetPositions()
                    R_sub = r_pos[conserved_r]
                    cR = np.mean(R_sub, axis=0)
                    H = (P_sub - cP).T @ (R_sub - cR)
                    U, S, Vt = np.linalg.svd(H)
                    d = np.linalg.det(Vt.T @ U.T)
                    V = Vt.T
                    if d < 0:
                        V[:, -1] *= -1
                    alg_P = (P_sub - cP) @ (V @ U.T).T + cR
                    return float(np.sqrt(np.mean(np.sum((alg_P - R_sub)**2, axis=1))))
                return min(cids, key=rmsd_score)
    return _select_lowest_energy_conf(mol, cids, solvent=getattr(model, "solvent", "water"))


def generate_lowest_energy_molblock(smiles: str, seed=42) -> str:
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return ""
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    cids = list(AllChem.EmbedMultipleConfs(mol, 16, params))
    if not cids:
        AllChem.Compute2DCoords(mol)
        return Chem.MolToMolBlock(mol)
    best = _select_lowest_energy_conf(mol, cids)
    return Chem.MolToMolBlock(mol, confId=best)


def embed_fragments(model, seed=42):
    """Generate 3D conformers for each reactant fragment, selecting product-aligned rotamers."""
    out = []
    for fi, m0 in enumerate(model.r_mols):
        mol = Chem.Mol(m0)
        o = model.offsets[fi]
        intra = [(i - o, j - o) for i, j in model.formed
                 if model.frag_of[i] == fi and model.frag_of[j] == fi]

        params = AllChem.ETKDGv3()
        params.randomSeed = seed
        nconf = 30 if intra else 24
        cids = list(AllChem.EmbedMultipleConfs(mol, nconf, params))

        if not cids:
            AllChem.Compute2DCoords(mol)
            out.append(mol)
            continue

        if intra:
            def score(cid):
                p = mol.GetConformer(cid).GetPositions()
                return sum(np.linalg.norm(p[i] - p[j]) for i, j in intra)
            best = min(cids, key=score)
        else:
            best = _select_product_aligned_conf(mol, fi, model, cids, seed=seed)

        out.append(Chem.Mol(mol, False, best))
    return out


def _build_clash_thresholds(model, nf, local, reactive):
    thr = {}
    for f in range(nf):
        for g in range(f + 1, nf):
            sf = [model.symbols[model.offsets[f] + k] for k in range(len(local[f]))]
            sg = [model.symbols[model.offsets[g] + k] for k in range(len(local[g]))]
            T = 0.75 * (np.array([rvdw(s) for s in sf])[:, None] +
                        np.array([rvdw(s) for s in sg])[None, :])
            f_reac = np.array([(f, k) in reactive for k in range(len(sf))])
            g_reac = np.array([(g, l) in reactive for l in range(len(sg))])
            T[f_reac[:, None] & g_reac[None, :]] = 0.0
            thr[(f, g)] = T
    return thr


def dock_fragments(model, frag_coords, variant=0, n_starts=32, seed=7):
    """Grid/stochastic docking of reactant fragments into a reactive encounter complex."""
    nf = len(frag_coords)
    if nf < 2:
        return [np.array(c, float) for c in frag_coords], 1, 0.0

    local = [np.array(c, float) - np.mean(c, axis=0) for c in frag_coords]
    anchor = int(np.argmax([len(c) for c in local]))

    def coords_in_frag(fi):
        o = model.offsets[fi]
        return lambda cj: local[fi][cj - o]

    pairs = []
    for i, j in model.formed:
        fi, fj = model.frag_of[i], model.frag_of[j]
        if fi == fj:
            continue
        d = model.d0(i, j) + DOCK_EXTRA
        mi = model.attack_mode(i, coords_in_frag(fi))
        mj = model.attack_mode(j, coords_in_frag(fj))
        pairs.append((fi, i - model.offsets[fi], fj, j - model.offsets[fj], d, mi, mj))

    reactive = {(model.frag_of[a], a - model.offsets[model.frag_of[a]])
                for b in model.formed for a in b}
    connected = {p[0] for p in pairs} | {p[2] for p in pairs}
    thr = _build_clash_thresholds(model, nf, local, reactive)
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
        res = minimize(objective, np.array(p0), method="L-BFGS-B", options={"maxiter": 400})
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


def auto_dock_molblocks(model, variant=0):
    mols = embed_fragments(model)
    coords = [m.GetConformer().GetPositions() for m in mols]
    docked, n_var, score = dock_fragments(model, coords, variant=variant)
    blocks = [Chem.MolToMolBlock(set_positions(m, x)) for m, x in zip(mols, docked)]
    return blocks, n_var, score
