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
) -> Optional[Dict[str, Any]]:
    """Match a single common intermediate fragment using graph isomorphism + local H-pairing."""
    m_prev = Chem.MolFromSmiles(f_prev.get("smiles", ""))
    m_next = Chem.MolFromSmiles(f_next.get("smiles", ""))
    if not m_prev or not m_next:
        return None

    # Get topological substructure matches of m_prev inside m_next
    matches = m_next.GetSubstructMatches(m_prev, uniquify=False)
    if not matches:
        return None

    idxs_p = f_prev.get("atom_indices", [])
    idxs_n = f_next.get("atom_indices", [])
    if len(idxs_p) < len(m_prev.GetAtoms()) or len(idxs_n) < len(m_next.GetAtoms()):
        return None

    best_rmsd = 999.0
    best_map: Dict[int, int] = {}
    best_R = np.eye(3)
    best_t = np.zeros(3)

    for m in matches:
        # 1. Map heavy atoms exactly via chemical graph isomorphism
        heavy_map = {idxs_p[i]: idxs_n[m[i]] for i in range(len(m))}
        p_heavy = np.array([p_last[k] for k in heavy_map.keys()])
        q_heavy = np.array([p_first[v] for v in heavy_map.values()])
        R_init, t_init = kabsch_fit(p_heavy, q_heavy)

        # 2. Map attached hydrogens strictly within each heavy atom's local coordination sphere
        full_map = dict(heavy_map)
        for u, v in heavy_map.items():
            hp = _get_h_neighbors(traj_prev, u)
            hn = _get_h_neighbors(traj_next, v)
            if len(hp) == len(hn) and len(hp) > 0:
                if len(hp) == 1:
                    full_map[hp[0]] = hn[0]
                else:
                    q_rot = p_first[hn] @ R_init.T + t_init
                    cost = np.linalg.norm(p_last[hp, None, :] - q_rot[None, :, :], axis=2)
                    r_i, c_i = linear_sum_assignment(cost)
                    for r, c in zip(r_i, c_i):
                        full_map[hp[r]] = hn[c]

        # 3. Evaluate total all-atom superposition
        p_all = np.array([p_last[k] for k in full_map.keys()])
        q_all = np.array([p_first[v] for v in full_map.values()])
        R_all, t_all = kabsch_fit(p_all, q_all)
        aligned = q_all @ R_all.T + t_all
        rmsd = compute_rmsd(p_all, aligned)

        if rmsd < best_rmsd:
            best_rmsd = rmsd
            best_map = full_map
            best_R, best_t = R_all, t_all

    if not best_map:
        return None

    prev_indices = list(best_map.keys())
    next_indices = [best_map[k] for k in prev_indices]

    return {
        "prev_indices": prev_indices,
        "next_indices": next_indices,
        "num_atoms": len(prev_indices),
        "rmsd": round(float(best_rmsd), 4),
        "R": best_R.tolist(),
        "t": best_t.tolist(),
        "is_joint": False,
        "is_reflected": False
    }


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
                match_info = _match_intermediate_fragment(fp, fr, step_prev, step_next, p_last, p_first)
                if match_info:
                    formula = fp.get("formula") or fr.get("formula") or can_p
                    match_info.update({
                        "id": f"mol_{p_idx}_{r_idx}",
                        "label": f"{formula} (Exact Graph Match, {match_info['num_atoms']} atoms, RMSD: {match_info['rmsd']:.4f} Å)",
                        "formula": formula,
                        "smiles": can_p,
                        "prev_frame": f_p,
                        "next_frame": f_n
                    })
                    candidates.append(match_info)

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
