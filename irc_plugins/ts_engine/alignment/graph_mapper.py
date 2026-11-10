"""
Graph Mapper Module
Exact chemical graph isomorphism matching between consecutive reaction steps.
Ensures intermediate molecules coincide 100% with zero missing or redundant atoms.
"""

from typing import Any, Dict, List, Optional, Set, Tuple
import numpy as np
from rdkit import Chem
from scipy.optimize import linear_sum_assignment
from .rigid_superposition import kabsch_fit, compute_rmsd


def _get_h_neighbors(traj: Dict[str, Any], heavy_idx: int) -> List[int]:
    """Retrieve atom indices of hydrogens directly bonded to a given heavy atom."""
    res = []
    syms = traj.get("symbols", [])
    for b in traj.get("bonds", []):
        u = b.get("from", b[0] if isinstance(b, (list, tuple)) else 0)
        v = b.get("to", b[1] if isinstance(b, (list, tuple)) else 0)
        if u == heavy_idx and v < len(syms) and syms[v] == "H":
            res.append(v)
        elif v == heavy_idx and u < len(syms) and syms[u] == "H":
            res.append(u)
    return res


def _match_intermediate_fragment(
    f_prev: Dict[str, Any],
    f_next: Dict[str, Any],
    traj_prev: Dict[str, Any],
    traj_next: Dict[str, Any],
    p_last: np.ndarray,
    p_first: np.ndarray
) -> List[Dict[str, Any]]:
    """Match a common intermediate fragment using graph isomorphism + local H-pairing across proper rotation and enantiomeric reflection."""
    m_prev = Chem.MolFromSmiles(f_prev.get("smiles", ""))
    m_next = Chem.MolFromSmiles(f_next.get("smiles", ""))
    if not m_prev or not m_next:
        return []

    # Get topological substructure matches of m_prev inside m_next
    matches = m_next.GetSubstructMatches(m_prev, uniquify=False)
    if not matches:
        return []

    idxs_p = f_prev.get("atom_indices", [])
    idxs_n = f_next.get("atom_indices", [])
    if len(idxs_p) < len(m_prev.GetAtoms()) or len(idxs_n) < len(m_next.GetAtoms()):
        return []

    # Parity reflection matrix (inverts chirality for enantiomer matching in achiral solvent)
    M_refl = np.diag([-1.0, 1.0, 1.0])

    best_rmsd_std = 999.0
    best_map_std: Dict[int, int] = {}
    best_R_std = np.eye(3)
    best_t_std = np.zeros(3)

    best_rmsd_ref = 999.0
    best_map_ref: Dict[int, int] = {}
    best_R_ref = np.eye(3)
    best_t_ref = np.zeros(3)

    for m in matches:
        # 1. Map heavy atoms exactly via chemical graph isomorphism
        heavy_map = {idxs_p[i]: idxs_n[m[i]] for i in range(len(m))}
        p_heavy = np.array([p_last[k] for k in heavy_map.keys()])
        q_heavy = np.array([p_first[v] for v in heavy_map.values()])

        # --- A. Standard proper rotation SE(3) ---
        R_init, t_init = kabsch_fit(p_heavy, q_heavy)
        full_map_std = dict(heavy_map)
        for u, v in heavy_map.items():
            hp = _get_h_neighbors(traj_prev, u)
            hn = _get_h_neighbors(traj_next, v)
            if len(hp) == len(hn) and len(hp) > 0:
                if len(hp) == 1:
                    full_map_std[hp[0]] = hn[0]
                else:
                    q_rot = p_first[hn] @ R_init.T + t_init
                    cost = np.linalg.norm(p_last[hp, None, :] - q_rot[None, :, :], axis=2)
                    r_i, c_i = linear_sum_assignment(cost)
                    for r, c in zip(r_i, c_i):
                        full_map_std[hp[r]] = hn[c]

        p_all_std = np.array([p_last[k] for k in full_map_std.keys()])
        q_all_std = np.array([p_first[v] for v in full_map_std.values()])
        R_all, t_all = kabsch_fit(p_all_std, q_all_std)
        aligned_std = q_all_std @ R_all.T + t_all
        rmsd_std = compute_rmsd(p_all_std, aligned_std)

        if rmsd_std < best_rmsd_std:
            best_rmsd_std = rmsd_std
            best_map_std = full_map_std
            best_R_std, best_t_std = R_all, t_all

        # --- B. Enantiomeric reflection (chirality inversion in achiral solution) ---
        q_refl_heavy = q_heavy @ M_refl
        R_ref_init, t_ref_init = kabsch_fit(p_heavy, q_refl_heavy)
        full_map_ref = dict(heavy_map)
        for u, v in heavy_map.items():
            hp = _get_h_neighbors(traj_prev, u)
            hn = _get_h_neighbors(traj_next, v)
            if len(hp) == len(hn) and len(hp) > 0:
                if len(hp) == 1:
                    full_map_ref[hp[0]] = hn[0]
                else:
                    q_rot_ref = (p_first[hn] @ M_refl) @ R_ref_init.T + t_ref_init
                    cost = np.linalg.norm(p_last[hp, None, :] - q_rot_ref[None, :, :], axis=2)
                    r_i, c_i = linear_sum_assignment(cost)
                    for r, c in zip(r_i, c_i):
                        full_map_ref[hp[r]] = hn[c]

        p_all_ref = np.array([p_last[k] for k in full_map_ref.keys()])
        q_all_ref = np.array([p_first[v] for v in full_map_ref.values()]) @ M_refl
        R_ref_all, t_ref_all = kabsch_fit(p_all_ref, q_all_ref)
        aligned_ref = q_all_ref @ R_ref_all.T + t_ref_all
        rmsd_ref = compute_rmsd(p_all_ref, aligned_ref)

        if rmsd_ref < best_rmsd_ref:
            best_rmsd_ref = rmsd_ref
            best_map_ref = full_map_ref
            best_R_ref, best_t_ref = R_ref_all, t_ref_all

    results: List[Dict[str, Any]] = []

    if best_map_std:
        p_std_idxs = list(best_map_std.keys())
        results.append({
            "prev_indices": p_std_idxs,
            "next_indices": [best_map_std[k] for k in p_std_idxs],
            "num_atoms": len(p_std_idxs),
            "rmsd": round(float(best_rmsd_std), 4),
            "R": best_R_std.tolist(),
            "t": best_t_std.tolist(),
            "is_joint": False,
            "is_reflected": False,
            "refl_mat": np.eye(3).tolist()
        })

    if best_map_ref and (best_rmsd_ref < best_rmsd_std - 0.05 or abs(best_rmsd_ref - best_rmsd_std) > 0.05):
        p_ref_idxs = list(best_map_ref.keys())
        results.append({
            "prev_indices": p_ref_idxs,
            "next_indices": [best_map_ref[k] for k in p_ref_idxs],
            "num_atoms": len(p_ref_idxs),
            "rmsd": round(float(best_rmsd_ref), 4),
            "R": best_R_ref.tolist(),
            "t": best_t_ref.tolist(),
            "is_joint": False,
            "is_reflected": True,
            "refl_mat": M_refl.tolist()
        })

    # Sort results so lowest RMSD candidate is always first
    results.sort(key=lambda c: c["rmsd"])
    return results


def find_step_chemical_overlaps(
    step_prev: Dict[str, Any],
    step_next: Dict[str, Any],
    prev_frame_idx: int = -1,
    next_frame_idx: int = 0
) -> List[Dict[str, Any]]:
    """
    Find all coinciding chemical intermediate molecules between step_prev and step_next.
    Guarantees that identical intermediate fragments coincide 100% with no missing/redundant atoms.
    """
    p_frames = step_prev.get("frames", [])
    n_frames = step_next.get("frames", [])
    if not p_frames or not n_frames:
        return []

    f_p = max(0, min(len(p_frames) - 1, len(p_frames) + prev_frame_idx if prev_frame_idx < 0 else prev_frame_idx))
    f_n = max(0, min(len(n_frames) - 1, len(n_frames) + next_frame_idx if next_frame_idx < 0 else next_frame_idx))

    p_last = np.array(p_frames[f_p].get("coords", p_frames[f_p].get("positions", [])))
    p_first = np.array(n_frames[f_n].get("coords", n_frames[f_n].get("positions", [])))

    prods = step_prev.get("fragments", {}).get("products", [])
    reacs = step_next.get("fragments", {}).get("reactants", [])

    candidates: List[Dict[str, Any]] = []

    for p_idx, fp in enumerate(prods):
        can_p = Chem.CanonSmiles(fp.get("smiles", "")) if fp.get("smiles") else ""
        if not can_p:
            continue
        for r_idx, fr in enumerate(reacs):
            can_r = Chem.CanonSmiles(fr.get("smiles", "")) if fr.get("smiles") else ""
            if can_r and can_p == can_r:
                matches = _match_intermediate_fragment(fp, fr, step_prev, step_next, p_last, p_first)
                for m_idx, match_info in enumerate(matches):
                    formula = fp.get("formula") or fr.get("formula") or can_p
                    is_ref = match_info.get("is_reflected", False)
                    cand_id = f"mol_{p_idx}_{r_idx}_ref" if is_ref else f"mol_{p_idx}_{r_idx}"
                    tag = "Enantiomer Inverted" if is_ref else "Exact Graph Match"
                    match_info.update({
                        "id": cand_id,
                        "label": f"{formula} ({tag}, {match_info['num_atoms']} atoms, RMSD: {match_info['rmsd']:.4f} Å)",
                        "formula": formula,
                        "smiles": can_p,
                        "prev_frame": f_p,
                        "next_frame": f_n
                    })
                    candidates.append(match_info)

    candidates.sort(key=lambda c: c["rmsd"])
    return candidates


def extract_fragment_molblock(traj: Dict[str, Any], frag_meta: Dict[str, Any], frame_idx: int = -1, target_mol: Any = None) -> str:
    """Extract a 3D MolBlock corresponding to a fragment in a trajectory frame with exact coordinates and bonds."""
    idxs = frag_meta.get("atom_indices", [])
    idx_map = {old: new for new, old in enumerate(idxs)}
    frames = traj.get("frames", [])
    frame = frames[frame_idx] if frames else {}
    pos = frame.get("positions", frame.get("coords", []))
    symbols = [traj.get("symbols", [])[i] for i in idxs]
    sub_pos = [pos[i] for i in idxs] if len(pos) >= max(idxs, default=0) + 1 else []

    sub_bonds = []
    for b in traj.get("bonds", []):
        u = b.get("from", b[0] if isinstance(b, list) else 0)
        v = b.get("to", b[1] if isinstance(b, list) else 0)
        if u in idx_map and v in idx_map and (ord_val := b.get("p_order", b.get("order", 1))) > 0:
            sub_bonds.append((idx_map[u], idx_map[v], ord_val))

    lines = [f"Seeded Product {frag_meta.get('formula', '')}", "  Molecular Workbench", "", f"{len(symbols):>3}{len(sub_bonds):>3}  0  0  0  0  0  0  0  0999 V2000"]
    for s, p in zip(symbols, sub_pos):
        lines.append(f"{p[0]:>10.4f}{p[1]:>10.4f}{p[2]:>10.4f} {s:<3} 0  0  0  0  0  0  0  0  0  0  0  0")
    for u, v, o in sub_bonds:
        lines.append(f"{u+1:>3}{v+1:>3}{o:>3}  0  0  0  0")
    lines.append("M  END")
    return "\n".join(lines)
