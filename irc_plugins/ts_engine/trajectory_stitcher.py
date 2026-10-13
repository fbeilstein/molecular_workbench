import json
from typing import Any, Dict, List, Optional, Tuple
import numpy as np
from rdkit import Chem


def kabsch_fit(P_target: np.ndarray, Q_source: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Compute optimal rotation matrix R and translation vector t minimizing RMSD(P, Q @ R.T + t)."""
    centroid_P = np.mean(P_target, axis=0)
    centroid_Q = np.mean(Q_source, axis=0)
    P_c = P_target - centroid_P
    Q_c = Q_source - centroid_Q

    H = Q_c.T @ P_c
    U, _, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T

    t = centroid_P - centroid_Q @ R.T
    return R, t


def find_hero_overlap(
    step_prev: Dict[str, Any], step_next: Dict[str, Any]
) -> Optional[Tuple[List[int], List[int], str]]:
    """Identify the common "Hero" fragment passing from Step K (products) to Step K+1 (reactants).

    Returns (prev_atom_indices, next_atom_indices, hero_smiles) or None.
    """
    prods_prev = step_prev.get("fragments", {}).get("products", [])
    reacs_next = step_next.get("fragments", {}).get("reactants", [])

    if not prods_prev or not reacs_next:
        return None

    # Sort candidates by size descending to prioritize the larger main scaffold
    prods_sorted = sorted(prods_prev, key=lambda f: f.get("num_atoms", 0), reverse=True)
    reacs_sorted = sorted(reacs_next, key=lambda f: f.get("num_atoms", 0), reverse=True)

    p_opts = Chem.SmilesParserParams()
    p_opts.removeHs = False

    for fp in prods_sorted:
        s_p = fp.get("smiles", "")
        mol_p = Chem.MolFromSmiles(s_p, p_opts)
        if not mol_p:
            mol_p = Chem.MolFromSmiles(s_p)
        if not mol_p:
            continue
        can_p = Chem.CanonSmiles(s_p)

        for fr in reacs_sorted:
            s_r = fr.get("smiles", "")
            mol_r = Chem.MolFromSmiles(s_r, p_opts)
            if not mol_r:
                mol_r = Chem.MolFromSmiles(s_r)
            if not mol_r:
                continue
            can_r = Chem.CanonSmiles(s_r)

            if can_p == can_r:
                match = mol_p.GetSubstructMatch(mol_r)
                if match and len(match) == len(fr["atom_indices"]):
                    idxs_prev = [fp["atom_indices"][m_i] for m_i in match]
                    idxs_next = list(fr["atom_indices"])
                    return idxs_prev, idxs_next, can_p

    return None


def align_trajectory_to_previous(
    step_next: Dict[str, Any],
    step_prev: Dict[str, Any],
    pair_indices: Optional[Tuple[List[int], List[int]]] = None
) -> Dict[str, Any]:
    """Rigidly rotate and translate all frames of step_next to perfectly overlay the Hero onto step_prev.

    Ensures zero camera/molecule jump (RMSD ~ 0 A) at step boundaries.
    """
    out_step = json.loads(json.dumps(step_next))
    if not step_prev.get("frames") or not out_step.get("frames"):
        return out_step

    if pair_indices is None:
        hero = find_hero_overlap(step_prev, out_step)
        if hero is None:
            return out_step
        idxs_prev, idxs_next, _ = hero
    else:
        idxs_prev, idxs_next = pair_indices

    # Target: final product frame of previous step
    p_last_prev = np.array(step_prev["frames"][-1]["positions"])
    # Source: first reactant frame of next step
    p_first_next = np.array(out_step["frames"][0]["positions"])

    P_target = p_last_prev[idxs_prev]
    Q_source = p_first_next[idxs_next]

    R, t = kabsch_fit(P_target, Q_source)

    # Apply transformation to ALL atoms in ALL frames of step_next
    for fr in out_step["frames"]:
        pos = np.array(fr["positions"])
        pos_aligned = pos @ R.T + t
        fr["positions"] = pos_aligned.round(5).tolist()

    return out_step


def stitch_multistep_trajectories(trajectories: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Stitch a sequence of elementary IRC reaction trajectories into 1 continuous Master MegaJSON.

    Sequentially aligns each step onto the preceding step's Hero coordinates with 0 RMSD.
    """
    if not trajectories:
        return {}
    if len(trajectories) == 1:
        return trajectories[0]

    aligned_steps = [json.loads(json.dumps(trajectories[0]))]
    step_reports = []

    for k in range(1, len(trajectories)):
        prev_step = aligned_steps[k - 1]
        next_step = trajectories[k]

        hero = find_hero_overlap(prev_step, next_step)
        if hero:
            idxs_p, idxs_n, hero_smiles = hero
            aligned_next = align_trajectory_to_previous(next_step, prev_step, (idxs_p, idxs_n))
            # Verify alignment RMSD on boundary
            p_end = np.array(prev_step["frames"][-1]["positions"])[idxs_p]
            p_start = np.array(aligned_next["frames"][0]["positions"])[idxs_n]
            b_rmsd = float(np.sqrt(np.mean((p_end - p_start) ** 2)))
            step_reports.append({
                "step_index": k + 1,
                "hero_smiles": hero_smiles,
                "n_matched_atoms": len(idxs_p),
                "joint_rmsd_angstrom": round(b_rmsd, 6)
            })
        else:
            aligned_next = json.loads(json.dumps(next_step))
            step_reports.append({
                "step_index": k + 1,
                "hero_smiles": None,
                "n_matched_atoms": 0,
                "joint_rmsd_angstrom": None
            })

        aligned_steps.append(aligned_next)

    # Compile into Master MegaJSON
    master_frames = []
    step_boundaries = []
    current_frame_idx = 0

    for k, step in enumerate(aligned_steps):
        n_fr = len(step.get("frames", []))
        step_boundaries.append({
            "step_index": k + 1,
            "reaction_smiles": step.get("reaction_smiles", ""),
            "start_frame": current_frame_idx,
            "end_frame": current_frame_idx + n_fr - 1,
            "ts_frame": current_frame_idx + step.get("ts_frame_index", 0)
        })
        master_frames.extend(step.get("frames", []))
        current_frame_idx += n_fr

    megajson = {
        "title": "Master Multistep IRC Trajectory",
        "n_steps": len(aligned_steps),
        "total_frames": len(master_frames),
        "step_boundaries": step_boundaries,
        "stitching_alignment_report": step_reports,
        "steps": aligned_steps,
        "combined_frames": master_frames
    }

    return megajson
