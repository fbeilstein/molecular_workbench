import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from ase import Atoms
from sella import Sella
from .calculators import make_xtb
from .common import rcov

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

def _embed_product_positions(model):
    """Generate aligned 3D conformers for product molecules."""
    p_coords = []
    for pmol in model.p_mols:
        m = Chem.Mol(pmol)
        AllChem.EmbedMolecule(m, randomSeed=42)
        try:
            AllChem.MMFFOptimizeMolecule(m)
        except Exception:
            pass
        p_coords.append(m.GetConformer().GetPositions())
    pos_P = np.zeros((model.n_atoms, 3))
    for ci in range(model.n_atoms):
        pfi, pai = model.r2p[ci]
        pos_P[ci] = p_coords[pfi][pai]
    return pos_P

def build_lego_ts_guess(pos_R, model):
    """Synthesize intermediate TS guess from Reactant (R) and Product (P) geometries."""
    pos_P = _embed_product_positions(model)
    pos_ts = pos_R.copy()

    if model.formed and len(model.r_mols) > 1:
        i_center, j_reag = model.formed[0]
        fi_center = model.frag_of[i_center]
        fj_reag = model.frag_of[j_reag]

        # Primary substrate fragment containing inverting/reaction center
        frag_A_atoms = [ci for ci in range(model.n_atoms) if model.frag_of[ci] == fi_center]
        leaving_atoms = {k for i, k in model.broken if model.frag_of[k] == fi_center}
        scaffold_A = [ci for ci in frag_A_atoms if ci not in leaving_atoms]
        if len(scaffold_A) < 3:
            scaffold_A = frag_A_atoms

        # Align Product strictly onto Substrate Fragment A scaffold
        pos_P_aligned = _kabsch_align(pos_P, pos_R, scaffold_A)

        # 1. Substrate inverting center and spectator ligands:
        # Interpolate between R and P (flattens umbrella from sp3 to planar sp2)
        for ci in scaffold_A:
            pos_ts[ci] = 0.5 * (pos_R[ci] + pos_P_aligned[ci])

        # 2. Leaving group(s) on substrate: elongate along breaking vector
        for i, k in model.broken:
            u_bk = pos_R[k] - pos_R[i]
            norm_bk = np.linalg.norm(u_bk)
            dir_bk = u_bk / norm_bk if norm_bk > 1e-4 else np.array([1.0, 0.0, 0.0])
            d_target = model.d0(i, k) + 0.35
            pos_ts[k] = pos_ts[i] + dir_bk * d_target

        # 3. Attacking reagent fragment:
        # Take orientation and internal coordinates from aligned Product!
        # Shift the entire reagent fragment as a rigid body along attack vector
        frag_B_atoms = [ci for ci in range(model.n_atoms) if model.frag_of[ci] == fj_reag]
        u_att = pos_P_aligned[j_reag] - pos_P_aligned[i_center]
        norm_att = np.linalg.norm(u_att)
        dir_att = u_att / norm_att if norm_att > 1e-4 else np.array([1.0, 0.0, 0.0])
        d_form_target = model.d0(i_center, j_reag) + 0.65

        current_j_pos = pos_P_aligned[j_reag]
        target_j_pos = pos_ts[i_center] + dir_att * d_form_target
        shift_B = target_j_pos - current_j_pos

        for cj in frag_B_atoms:
            pos_ts[cj] = pos_P_aligned[cj] + shift_B
    else:
        # Unimolecular reaction or elimination: elongate leaving group
        for i, k in model.broken:
            u_bk = pos_R[k] - pos_R[i]
            norm_bk = np.linalg.norm(u_bk)
            dir_bk = u_bk / norm_bk if norm_bk > 1e-4 else np.array([1.0, 0.0, 0.0])
            extra = 0.55 if not model.formed else 0.35
            pos_ts[k] = pos_ts[i] + dir_bk * (model.d0(i, k) + extra)

    return pos_ts


def run_lego_ts(atoms, model, log):
    """Run Lego Dual-Endpoint TS guess synthesis followed by Sella DLC saddle polish."""
    log("=== Strategy: Lego Dual-Endpoint Superposition ===")
    log("  Step 1: Superimposing Product onto Reactant substrate scaffold...")
    pos_guess = build_lego_ts_guess(atoms.positions.copy(), model)
    atoms.set_positions(pos_guess)
    log("  Step 2: Synthesized intermediate TS geometry (sp2 flattened, rigid reagent attack).")

    log("  Step 3: Polishing saddle point with Sella in internal coordinates...")
    atoms.calc = make_xtb(model.charge, model.multiplicity)
    opt_sella = Sella(atoms, internal=True, order=1, logfile=None)
    profile = []

    def on_step():
        try:
            profile.append({
                "step": opt_sella.nsteps,
                "energy_ev": float(atoms.get_potential_energy()),
                "fmax": float(np.max(np.linalg.norm(atoms.get_forces(), axis=1)))
            })
        except Exception:
            pass

    opt_sella.attach(on_step)
    converged = opt_sella.run(fmax=0.03, steps=60)
    log(f"  Sella Polish {'Converged' if converged else 'Finished'} in {opt_sella.nsteps} steps.")
    return profile
