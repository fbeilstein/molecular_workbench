import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from scipy.optimize import minimize
from ase.optimize import LBFGS
from sella import Sella
from .calculators import make_xtb, _make_restrained_calc


def _kabsch_align(P_pts, R_pts, indices):
    """Rigidly align points P onto R using the specified indices via Kabsch SVD."""
    if len(indices) < 3:
        return P_pts - np.mean(P_pts, axis=0) + np.mean(R_pts, axis=0)
    P_sub, R_sub = P_pts[indices], R_pts[indices]
    c_P, c_R = np.mean(P_sub, axis=0), np.mean(R_sub, axis=0)
    H = (P_sub - c_P).T @ (R_sub - c_R)
    U, S, Vt = np.linalg.svd(H)
    d = np.linalg.det(Vt.T @ U.T)
    V = Vt.T
    if d < 0:
        V[:, -1] *= -1
    R_mat = V @ U.T
    return (P_pts - c_P) @ R_mat.T + c_R


def _embed_aligned_products(model, pos_R):
    """Generate aligned, fully relaxed product 3D geometry preserving reactant scaffold & dihedrals."""
    from ase import Atoms
    from ase.optimize import BFGS
    from .fragment_relax import _rigid_align

    pos_P = pos_R.copy()
    if getattr(model, "broken", None):
        for bi, bj in model.broken:
            pfi_i = model.r2p.get(bi, (-1, -1))[0]
            pfi_j = model.r2p.get(bj, (-1, -1))[0]
            if pfi_i >= 0 and pfi_j >= 0 and pfi_i != pfi_j:
                vec = pos_R[bj] - pos_R[bi]
                u = vec / (np.linalg.norm(vec) + 1e-6)
                j_atoms = [k for k in range(model.n_atoms) if model.r2p.get(k, (-1, -1))[0] == pfi_j]
                pos_P[j_atoms] += u * 4.0

    if getattr(model, "formed", None):
        for fi, fj in model.formed:
            rfi_i = model.frag_of[fi]
            rfi_j = model.frag_of[fj]
            if rfi_i != rfi_j:
                vec = pos_P[fj] - pos_P[fi]
                dist = float(np.linalg.norm(vec))
                target_d = model.d0(fi, fj)
                if dist > target_d:
                    shift = (dist - target_d) * (vec / (dist + 1e-6))
                    j_atoms = [k for k in range(model.n_atoms) if model.frag_of[k] == rfi_j]
                    pos_P[j_atoms] -= shift

    solvent = getattr(model, "solvent", "water")
    for pfi, pmol in enumerate(getattr(model, "p_mols", [])):
        r_atoms = [ci for ci in range(model.n_atoms) if model.r2p.get(ci, (-1, -1))[0] == pfi]
        if len(r_atoms) > 1:
            try:
                chg = int(sum(a.GetFormalCharge() for a in pmol.GetAtoms()))
                mult = int(sum(a.GetNumRadicalElectrons() for a in pmol.GetAtoms())) + 1
                sub_syms = [model.symbols[ci] for ci in r_atoms]
                at = Atoms(symbols=sub_syms, positions=pos_P[r_atoms].copy())
                at.calc = make_xtb(chg, mult, solvent=solvent)
                opt = BFGS(at, logfile=None)
                opt.run(fmax=0.05, steps=30)
                pos_P[r_atoms] = _rigid_align(at.get_positions(), pos_P[r_atoms])
            except Exception:
                pass

    return pos_P


def build_lego_ts_guess(pos_R, model):
    """Universal Dual-Endpoint reaction midpoint interpolation via IDPP pair potential."""
    pos_P = _embed_aligned_products(model, pos_R)
    D_R = np.linalg.norm(pos_R[:, None, :] - pos_R[None, :, :], axis=2)
    D_P = np.linalg.norm(pos_P[:, None, :] - pos_P[None, :, :], axis=2)
    D_target = 0.5 * (D_R + D_P)

    # Universal TS bond target metrics across all reaction classes
    for i, j in model.formed:
        D_target[i, j] = D_target[j, i] = model.d0(i, j) + 0.65
    for i, j in model.broken:
        D_target[i, j] = D_target[j, i] = model.d0(i, j) + 0.35

    def idpp_obj(x):
        X = x.reshape(-1, 3)
        D = np.linalg.norm(X[:, None, :] - X[None, :, :], axis=2)
        diff = D - D_target
        np.fill_diagonal(diff, 0.0)
        w = 1.0 / (D_target**4 + 1e-4)
        np.fill_diagonal(w, 0.0)
        val = np.sum(w * (diff**2))
        grad = np.zeros_like(X)
        for a in range(len(X)):
            for b in range(len(X)):
                if a == b:
                    continue
                dist = D[a, b] + 1e-7
                grad[a] += 4.0 * w[a, b] * diff[a, b] * (X[a] - X[b]) / dist
        return val, grad.ravel()

    x0 = 0.5 * (pos_R + pos_P).ravel()
    res = minimize(idpp_obj, x0, jac=True, method="L-BFGS-B", options={"maxiter": 250})
    return res.x.reshape(-1, 3)


def run_lego_ts(atoms, model, log):
    """Universal Dual-Endpoint IDPP interpolation -> Valley Relax -> Sella saddle polish."""
    log("=== Strategy: Dual-Endpoint IDPP Reaction Interpolation ===")
    log("  Step 1: Generating aligned Product conformers and reaction distance metric...")
    pos_guess = build_lego_ts_guess(atoms.positions.copy(), model)
    atoms.set_positions(pos_guess)
    log(f"  Step 2: IDPP geometric midpoint synthesized for {len(model.formed)} forming & {len(model.broken)} breaking bonds.")

    # Stage 2: Quantum Valley Internal Relaxation
    solvent = getattr(model, "solvent", "water")
    log(f"  Step 3: Quantum valley relaxation of internal coordinates in {solvent.upper()}...")
    base = atoms.copy()
    base.calc = make_xtb(model.charge, model.multiplicity, solvent=solvent)
    targets = [(i, j, model.d0(i, j) + 0.65) for i, j in model.formed]
    targets += [(i, j, model.d0(i, j) + 0.35) for i, j in model.broken]
    pos_cur = atoms.positions.copy()
    spec = [(i, j, float(np.linalg.norm(pos_cur[i] - pos_cur[j]))) for i, j in model.spectator_bonds]
    calc = _make_restrained_calc(base, targets, k=30.0, spectator_restraints=spec, model=model)
    atoms.calc = calc
    opt_valley = LBFGS(atoms, logfile=None, maxstep=0.12)
    opt_valley.run(fmax=0.08, steps=35)
    log(f"  Internal valley relaxation finished in {opt_valley.nsteps} steps.")

    # Stage 3: Sella DLC Saddle Point Polish (guarded)
    log("  Step 4: Polishing saddle point with Sella in internal coordinates...")
    from .calculators import guard_saddle_polish
    from .vibrations import analyse_mode
    atoms.calc = make_xtb(model.charge, model.multiplicity, solvent=solvent)
    ts_rep = None
    try:
        ts_rep = analyse_mode(atoms, model, log)
    except Exception:
        pass

    if ts_rep and ts_rep.get("n_imag") == 1 and ts_rep.get("mode_consistent", True):
        log(f"  ✓ Saddle point confirmed by vibrational analysis ({ts_rep['imag_freqs'][0]:.0f}i cm-1).")
        return [{"step": 1, "energy_ev": float(atoms.get_potential_energy()), "fmax": float(np.max(np.linalg.norm(atoms.get_forces(), axis=1)))}]

    if ts_rep and ts_rep.get("n_imag") == 0:
        log("  ℹ️ Local minimum detected in scan: steering along reaction coordinate into saddle region...")

    return guard_saddle_polish(atoms, model, log, max_steps=50, max_drift=0.35)

