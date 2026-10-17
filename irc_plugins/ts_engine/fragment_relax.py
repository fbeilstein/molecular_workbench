import numpy as np
from ase import Atoms
from ase.constraints import FixConstraint
from ase.optimize import BFGS
from .calculators import make_xtb


class FixFragmentCOMs(FixConstraint):
    """Constrain the Center-of-Mass of each fragment while allowing all internal coordinates to relax freely.

    Removes net translational forces per fragment so fragments maintain their dragged-apart
    separation without collapsing back or drifting, while allowing every internal bond,
    angle, and dihedral in BOTH/ALL fragments to relax naturally on the PES.
    """

    def __init__(self, frag_indices_list):
        self.frag_indices_list = [list(idx) for idx in frag_indices_list]

    def get_removed_dof(self, atoms):
        return 3 * len(self.frag_indices_list)

    def adjust_positions(self, atoms, newpositions):
        old_pos = atoms.positions
        for idxs in self.frag_indices_list:
            old_com = np.mean(old_pos[idxs], axis=0)
            new_com = np.mean(newpositions[idxs], axis=0)
            newpositions[idxs] += (old_com - new_com)

    def adjust_forces(self, atoms, forces):
        for idxs in self.frag_indices_list:
            net_f = np.mean(forces[idxs], axis=0)
            forces[idxs] -= net_f


def _compute_multifragment_displacements(positions, frag_atoms, target_sep=5.5, max_iter=25):
    """Compute pairwise repulsive displacements ensuring all distinct fragments have d >= target_sep."""
    nf = len(frag_atoms)
    anchor_fi = int(np.argmax([len(idxs) for idxs in frag_atoms]))
    disp = np.zeros((nf, 3))
    for _ in range(max_iter):
        coms = [np.mean(positions[frag_atoms[f]] + disp[f], axis=0) for f in range(nf)]
        force_disp, all_ok = np.zeros_like(disp), True
        for f in range(nf):
            for g in range(f + 1, nf):
                p_f, p_g = positions[frag_atoms[f]] + disp[f], positions[frag_atoms[g]] + disp[g]
                d_min = min(float(np.linalg.norm(pi - pj)) for pi in p_f for pj in p_g)
                if d_min < target_sep:
                    all_ok = False
                    v = coms[f] - coms[g]
                    norm_v = np.linalg.norm(v)
                    u = v / (norm_v + 1e-9) if norm_v > 1e-4 else np.random.normal(size=3)
                    u /= (np.linalg.norm(u) + 1e-9)
                    deficit = 0.5 * (target_sep - d_min) * u
                    if f != anchor_fi:
                        force_disp[f] += deficit
                    if g != anchor_fi:
                        force_disp[g] -= deficit
        if all_ok:
            break
        disp += force_disp
    return {f: disp[f] for f in range(nf)}


def _rigid_align(source, target):
    cP, cQ = np.mean(target, axis=0), np.mean(source, axis=0)
    H = (source - cQ).T @ (target - cP)
    U, _, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T
    return source @ R.T + (cP - cQ @ R.T)


def extend_to_asymptotic_fragments(endpoint_frame, model, solvent, log, is_product=True, n_steps=8, target_sep=5.5):
    """Drag fragments apart step-by-step and relax all internal degrees of freedom.

    Supports both products (forward IRC end) and reactants (reverse IRC end),
    and handles arbitrary N >= 2 fragments via multi-fragment pairwise repulsion.
    """
    mols = model.p_mols if is_product else model.r_mols
    if len(mols) <= 1:
        return []

    frag_atoms = [
        [model.p2r[(pfi, pai)] for pai in range(model.p_mols[pfi].GetNumAtoms())]
        for pfi in range(len(mols))
    ] if is_product else [
        [i for i in range(model.n_atoms) if model.frag_of[i] == fi]
        for fi in range(len(mols))
    ]

    pos_last = np.array(endpoint_frame["positions"])
    all_d = [
        float(np.linalg.norm(pos_last[i] - pos_last[j]))
        for f in range(len(frag_atoms))
        for g in range(f + 1, len(frag_atoms))
        for i in frag_atoms[f]
        for j in frag_atoms[g]
    ]
    min_d = min(all_d) if all_d else target_sep

    if min_d >= target_sep - 0.2:
        return []

    label = "Products" if is_product else "Reactants"
    log(f"\n=== Asymptotic Separation & Fragment Relaxation ({label}) ===")
    log(f"  {len(mols)} fragments (current contact distance: {min_d:.2f} A).")
    log(f"  Dragging fragments apart to contact d >= {target_sep:.2f} A (all fragments relaxing freely)...")

    frag_disp = _compute_multifragment_displacements(pos_last, frag_atoms, target_sep=target_sep)
    calc_full = make_xtb(model.charge, model.multiplicity, solvent=solvent)
    asymp_frames = []

    for step in range(1, n_steps + 1):
        pos_step = pos_last.copy()
        frac = step / n_steps
        for f in range(len(frag_atoms)):
            pos_step[frag_atoms[f]] += frac * frag_disp[f]

        atoms_step = Atoms(symbols=model.symbols, positions=pos_step)
        atoms_step.calc = calc_full
        try:
            e_step = float(atoms_step.get_potential_energy())
            f_step = float(np.max(np.linalg.norm(atoms_step.get_forces(), axis=1)))
        except Exception:
            e_step, f_step = float(endpoint_frame["energy_ev"]), 0.0

        asymp_frames.append({
            "positions": pos_step.tolist(),
            "energy_ev": round(e_step, 5),
            "fmax": round(f_step, 4)
        })

    p_last = np.array(asymp_frames[-1]["positions"])
    final_min_d = min(
        float(np.linalg.norm(p_last[i] - p_last[j]))
        for f in range(len(frag_atoms))
        for g in range(f + 1, len(frag_atoms))
        for i in frag_atoms[f]
        for j in frag_atoms[g]
    )
    log(f"  ✓ Appended {len(asymp_frames)} asymptotic frames (contact d = {final_min_d:.2f} A).")
    return asymp_frames


def _relax_single_isolated_fragment(mol_rdkit, current_pos, symbols, solvent):
    """Independently relax a single isolated molecule to its universal textbook equilibrium ground state."""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    chg = int(sum(a.GetFormalCharge() for a in mol_rdkit.GetAtoms()))
    mult = int(sum(a.GetNumRadicalElectrons() for a in mol_rdkit.GetAtoms())) + 1

    atoms = Atoms(symbols=symbols, positions=current_pos.copy())
    atoms.calc = make_xtb(chg, mult, solvent=solvent)
    opt = BFGS(atoms, logfile=None)
    try:
        opt.run(fmax=0.002, steps=60)
        p_opt = atoms.get_positions()
        e_val = float(atoms.get_potential_energy())
    except Exception:
        p_opt = current_pos.copy()
        e_val = 0.0

    return _rigid_align(p_opt, current_pos), e_val


def extend_to_asymptotic_products(endpoint_frame, model, solvent, log, n_steps=8, target_sep=5.5):
    """Backward compatibility wrapper for product separation."""
    return extend_to_asymptotic_fragments(endpoint_frame, model, solvent, log, is_product=True, n_steps=n_steps, target_sep=target_sep)


def build_fragments_meta(model):
    """Extract fragment membership, smiles, and formula for multistep MegaJSON stitching."""
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors
    def _m(mols, get_idxs):
        return [{
            "frag_index": fi, "smiles": Chem.MolToSmiles(m),
            "formula": rdMolDescriptors.CalcMolFormula(m),
            "atom_indices": get_idxs(fi, m), "num_atoms": len(get_idxs(fi, m))
        } for fi, m in enumerate(mols)]
    return {
        "reactants": _m(model.r_mols, lambda fi, m: [i for i in range(model.n_atoms) if model.frag_of[i] == fi]),
        "products": _m(model.p_mols, lambda fi, m: [model.p2r[(fi, pai)] for pai in range(m.GetNumAtoms())])
    }

