import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from scipy.spatial.transform import Rotation
from scipy.optimize import minimize
from .common import rvdw, _unit, dir_pen, set_positions, DOCK_EXTRA


def embed_fragments(model, seed=42):
    """Generate 3D conformers for each reactant fragment."""
    out = []
    for fi, m0 in enumerate(model.r_mols):
        mol = Chem.Mol(m0)
        o = model.offsets[fi]
        intra = [(i - o, j - o) for i, j in model.formed
                 if model.frag_of[i] == fi and model.frag_of[j] == fi]

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
