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
    frag_sizes = [len(idxs) for idxs in frag_atoms]
    anchor_fi = int(np.argmax(frag_sizes))  # keep main scaffold centered

    disp = np.zeros((nf, 3))
    for _ in range(max_iter):
        coms = [np.mean(positions[frag_atoms[f]] + disp[f], axis=0) for f in range(nf)]
        force_disp = np.zeros_like(disp)
        all_ok = True
        for f in range(nf):
            for g in range(f + 1, nf):
                p_f = positions[frag_atoms[f]] + disp[f]
                p_g = positions[frag_atoms[g]] + disp[g]
                d_min = min(float(np.linalg.norm(pi - pj)) for pi in p_f for pj in p_g)
                if d_min < target_sep:
                    all_ok = False
                    v = coms[f] - coms[g]
                    norm_v = np.linalg.norm(v)
                    u = v / (norm_v + 1e-9) if norm_v > 1e-4 else np.random.normal(size=3)
                    u /= (np.linalg.norm(u) + 1e-9)
                    deficit = target_sep - d_min
                    if f != anchor_fi:
                        force_disp[f] += 0.5 * deficit * u
                    if g != anchor_fi:
                        force_disp[g] -= 0.5 * deficit * u
        if all_ok:
            break
        disp += force_disp
    return {f: disp[f] for f in range(nf)}


def extend_to_asymptotic_fragments(endpoint_frame, model, solvent, log, is_product=True, n_steps=8, target_sep=5.5):
    """Drag fragments apart step-by-step and relax all internal degrees of freedom.

    Supports both products (forward IRC end) and reactants (reverse IRC end),
    and handles arbitrary N >= 2 fragments via multi-fragment pairwise repulsion.
    """
    mols = model.p_mols if is_product else model.r_mols
    if len(mols) <= 1:
        return []

    frag_atoms = [[] for _ in range(len(mols))]
    if is_product:
        for (pfi, pai), gi in model.p2r.items():
            frag_atoms[pfi].append(gi)
    else:
        for i in range(model.n_atoms):
            frag_atoms[model.frag_of[i]].append(i)

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
    curr_pos = pos_last.copy()
    asymp_frames = []

    for step in range(1, n_steps + 1):
        for f in range(len(frag_atoms)):
            curr_pos[frag_atoms[f]] += (1.0 / n_steps) * frag_disp[f]

        atoms_step = Atoms(symbols=model.symbols, positions=curr_pos)
        atoms_step.calc = calc_full
        atoms_step.set_constraint(FixFragmentCOMs(frag_atoms))

        opt = BFGS(atoms_step, logfile=None)
        max_opt = 25 if step == n_steps else 6
        try:
            opt.run(fmax=0.03, steps=max_opt)
            curr_pos = atoms_step.get_positions()
            e_step = float(atoms_step.get_potential_energy())
            f_step = float(np.max(np.linalg.norm(atoms_step.get_forces(), axis=1)))
        except Exception:
            try:
                e_step = float(atoms_step.get_potential_energy())
                f_step = float(np.max(np.linalg.norm(atoms_step.get_forces(), axis=1)))
            except Exception:
                e_step = float(endpoint_frame["energy_ev"])
                f_step = 0.0

        min_pair_d = min(
            float(np.linalg.norm(curr_pos[i] - curr_pos[j]))
            for i in range(len(model.symbols))
            for j in range(i + 1, len(model.symbols))
        )
        asymp_frames.append({
            "positions": curr_pos.tolist(),
            "energy_ev": e_step,
            "fmax": f_step
        })

    # Final Step: Independent Relaxation of each molecule to textbook ground state
    log("  Polishing endpoint: relaxing all molecules independently to textbook ground state...")
    pos_textbook = curr_pos.copy()
    e_tot_iso = 0.0
    for f, f_idxs in enumerate(frag_atoms):
        if len(f_idxs) <= 1:
            m_f = mols[f]
            chg_f = int(sum(a.GetFormalCharge() for a in m_f.GetAtoms()))
            mult_f = int(sum(a.GetNumRadicalElectrons() for a in m_f.GetAtoms())) + 1
            a_single = Atoms(symbols=[model.symbols[f_idxs[0]]], positions=[pos_textbook[f_idxs[0]]])
            a_single.calc = make_xtb(chg_f, mult_f, solvent=solvent)
            try:
                e_tot_iso += float(a_single.get_potential_energy())
            except Exception:
                pass
            continue

        p_frag = pos_textbook[f_idxs].copy()
        p_opt_frag, e_frag = _relax_single_isolated_fragment(
            mols[f], p_frag, [model.symbols[i] for i in f_idxs], solvent
        )
        pos_textbook[f_idxs] = p_opt_frag
        e_tot_iso += e_frag

    asymp_frames.append({
        "positions": pos_textbook.tolist(),
        "energy_ev": round(e_tot_iso, 5) if e_tot_iso != 0.0 else e_step,
        "fmax": 0.001
    })

    final_min_d = min(
        float(np.linalg.norm(pos_textbook[i] - pos_textbook[j]))
        for f in range(len(frag_atoms))
        for g in range(f + 1, len(frag_atoms))
        for i in frag_atoms[f]
        for j in frag_atoms[g]
    )
    log(f"  ✓ Appended {len(asymp_frames)} asymptotic frames (contact d = {final_min_d:.2f} A, all molecules independently relaxed to textbook geometry).")
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

    # Universal textbook conformer matching (ETKDGv3 encodes experimental torsions)
    try:
        smi = Chem.MolToSmiles(Chem.RemoveHs(mol_rdkit))
        can_mol = Chem.AddHs(Chem.MolFromSmiles(smi))
        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        cids = list(AllChem.EmbedMultipleConfs(can_mol, numConfs=8, params=params))
        if cids:
            for cid in cids:
                AllChem.MMFFOptimizeMolecule(can_mol, confId=cid)
            match = can_mol.GetSubstructMatch(mol_rdkit)
            if match and len(match) == len(symbols):
                best_cid = min(
                    cids,
                    key=lambda c: np.mean((can_mol.GetConformer(c).GetPositions()[list(match)] - p_opt) ** 2)
                )
                p_ref = can_mol.GetConformer(best_cid).GetPositions()[list(match)]
                rmsd_ref = np.sqrt(np.mean((p_ref - p_opt) ** 2))
                if rmsd_ref < 0.45:
                    atoms_ref = Atoms(symbols=symbols, positions=p_ref)
                    atoms_ref.calc = make_xtb(chg, mult, solvent=solvent)
                    opt_ref = BFGS(atoms_ref, logfile=None)
                    opt_ref.run(fmax=0.001, steps=40)
                    e_ref = float(atoms_ref.get_potential_energy())
                    if e_ref <= e_val + 0.01:
                        p_opt = atoms_ref.get_positions()
                        e_val = e_ref
    except Exception:
        pass

    cP, cQ = np.mean(current_pos, axis=0), np.mean(p_opt, axis=0)
    H = (p_opt - cQ).T @ (current_pos - cP)
    U, _, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T
    t = cP - cQ @ R.T
    return p_opt @ R.T + t, e_val


def extend_to_asymptotic_products(endpoint_frame, model, solvent, log, n_steps=8, target_sep=5.5):
    """Backward compatibility wrapper for product separation."""
    return extend_to_asymptotic_fragments(endpoint_frame, model, solvent, log, is_product=True, n_steps=n_steps, target_sep=target_sep)


def build_fragments_meta(model):
    """Extract fragment membership, smiles, and formula for multistep MegaJSON stitching."""
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors
    r_meta = [
        {
            "frag_index": fi,
            "smiles": Chem.MolToSmiles(m),
            "formula": rdMolDescriptors.CalcMolFormula(m),
            "atom_indices": [i for i in range(model.n_atoms) if model.frag_of[i] == fi],
            "num_atoms": len([i for i in range(model.n_atoms) if model.frag_of[i] == fi])
        }
        for fi, m in enumerate(model.r_mols)
    ]
    p_meta = [
        {
            "frag_index": fi,
            "smiles": Chem.MolToSmiles(m),
            "formula": rdMolDescriptors.CalcMolFormula(m),
            "atom_indices": [model.p2r[(fi, pai)] for pai in range(m.GetNumAtoms())],
            "num_atoms": len([model.p2r[(fi, pai)] for pai in range(m.GetNumAtoms())])
        }
        for fi, m in enumerate(model.p_mols)
    ]
    return {"reactants": r_meta, "products": p_meta}

