"""
Orbital Matching and Smoothing Module — Hungarian Assignment & 2D Degenerate Alignment.
"""

from collections import Counter
from typing import Any, Callable, Dict, List, Optional, Set, Tuple
import numpy as np
from scipy.optimize import linear_sum_assignment

from .common import SUBSCRIPTS


def _add_subscript(label: str, sub: str) -> str:
    for prefix in ("π*_deloc(", "π_deloc(", "π*(", "π(", "σ*(", "σ(", "LP("):
        if prefix in label:
            pos = label.index(prefix) + len(prefix) - 1
            return label[:pos] + sub + label[pos:]
    return f"{label} {sub}"


def _align_degenerate_2d_pair(
    prev_pair: np.ndarray, curr_pair: np.ndarray, s_cross: np.ndarray
) -> np.ndarray:
    """Compute optimal 2x2 Procrustes rotation aligning curr_pair to prev_pair."""
    W = prev_pair.T @ s_cross @ curr_pair
    U_w, s_w, Vt_w = np.linalg.svd(W)
    if s_w[0] > 0.50 and s_w[1] > 0.50:
        d = np.linalg.det(Vt_w.T @ U_w.T)
        if d < 0:
            Vt_w_mod = Vt_w.copy()
            Vt_w_mod[-1, :] *= -1.0
            R_w = Vt_w_mod.T @ U_w.T
        else:
            R_w = Vt_w.T @ U_w.T
        return curr_pair @ R_w
    return curr_pair


def _align_degenerate_subspace_pairs(
    active_items: List[Dict[str, Any]],
    curr_items: List[Dict[str, Any]],
    prev_mat: np.ndarray,
    curr_mat: np.ndarray,
    s_cross: np.ndarray,
    pi_type_keyword: str,
) -> None:
    """Detect and align pairs of degenerate pi/pi* orbitals across frames."""
    Q = len(curr_items)
    for i1 in range(Q):
        for i2 in range(i1 + 1, Q):
            m1, m2 = curr_items[i1]["meta"], curr_items[i2]["meta"]
            if pi_type_keyword in m1.get("type", "") and pi_type_keyword in m2.get("type", ""):
                ats1 = tuple(sorted(m1.get("atoms", [])))
                ats2 = tuple(sorted(m2.get("atoms", [])))
                if len(ats1) == 2 and ats1 == ats2:
                    prev_matches = [
                        pi_idx for pi_idx, av in enumerate(active_items)
                        if pi_type_keyword in av.get("type", "") and tuple(sorted(av.get("atoms", []))) == ats1
                    ]
                    if len(prev_matches) == 2:
                        p_pair = prev_mat[:, prev_matches]
                        c_pair = curr_mat[:, [i1, i2]]
                        c_pair_rot = _align_degenerate_2d_pair(p_pair, c_pair, s_cross)
                        curr_mat[:, [i1, i2]] = c_pair_rot
                        curr_items[i1]["vec"] = c_pair_rot[:, 0]
                        curr_items[i2]["vec"] = c_pair_rot[:, 1]


def match_subspace_orbitals(
    active_items: List[Dict[str, Any]],
    curr_items: List[Dict[str, Any]],
    s_cross: Optional[np.ndarray],
    t: int,
    subspace: str,
    threshold: float,
    next_seq: int,
    record_fn: Callable[..., None],
    prefix: str = "trk_",
    pad_zeros: int = 4,
    pi_type_keyword: Optional[str] = None,
) -> Tuple[List[Dict[str, Any]], int]:
    """
    Subspace-specific Hungarian cross-overlap assignment and track propagation.
    Returns (next_active_items, updated_next_seq).
    """
    if not curr_items:
        return [], next_seq

    def _fmt_tid(num: int) -> str:
        return f"{prefix}{num:0{pad_zeros}d}"

    next_active: List[Dict[str, Any]] = []

    if t == 0 or s_cross is None or len(active_items) == 0:
        for item in curr_items:
            tid = _fmt_tid(next_seq)
            next_seq += 1
            meta = item["meta"]
            record_fn(
                tid, item["vec"], item["pop"], meta.get("energy_ev", 0.0), meta.get("label", ""),
                meta.get("type", subspace), meta.get("category", subspace),
                meta.get("orbital_class", "bonding"), meta.get("atoms", []), subspace
            )
            next_active.append({
                "tid": tid, "vec": item["vec"],
                "atoms": meta.get("atoms", []), "type": meta.get("type", subspace)
            })
        return next_active, next_seq

    prev_mat = np.column_stack([item["vec"] for item in active_items])
    curr_mat = np.column_stack([item["vec"] for item in curr_items])

    if pi_type_keyword is not None:
        _align_degenerate_subspace_pairs(active_items, curr_items, prev_mat, curr_mat, s_cross, pi_type_keyword)

    Omega = prev_mat.T @ s_cross @ curr_mat
    r_ind, c_ind = linear_sum_assignment(-np.abs(Omega))
    matched_curr: Set[int] = set()

    for r, c in zip(r_ind, c_ind):
        if abs(Omega[r, c]) >= threshold:
            matched_curr.add(c)
            tid = active_items[r]["tid"]
            v = curr_mat[:, c].copy()
            if Omega[r, c] < 0.0:
                v = -v
            meta = curr_items[c]["meta"]
            record_fn(
                tid, v, curr_items[c]["pop"], meta.get("energy_ev", 0.0), meta.get("label", ""),
                meta.get("type", subspace), meta.get("category", subspace),
                meta.get("orbital_class", "bonding"), meta.get("atoms", []), subspace
            )
            next_active.append({
                "tid": tid, "vec": v,
                "atoms": meta.get("atoms", []), "type": meta.get("type", subspace)
            })

    for c in range(len(curr_items)):
        if c not in matched_curr:
            tid = _fmt_tid(next_seq)
            next_seq += 1
            meta = curr_items[c]["meta"]
            record_fn(
                tid, curr_items[c]["vec"], curr_items[c]["pop"], meta.get("energy_ev", 0.0),
                meta.get("label", ""), meta.get("type", subspace), meta.get("category", subspace),
                meta.get("orbital_class", "bonding"), meta.get("atoms", []), subspace
            )
            next_active.append({
                "tid": tid, "vec": curr_items[c]["vec"],
                "atoms": meta.get("atoms", []), "type": meta.get("type", subspace)
            })

    return next_active, next_seq


def _role_prefix(lbl: Optional[str]) -> str:
    if not lbl:
        return ""
    for pfx in ("π*_deloc", "π_deloc", "π*", "σ*", "π", "σ", "LP", "Core"):
        if lbl.startswith(pfx):
            return pfx
    return lbl.split("(")[0]


def smooth_and_label_tracks(track_records: Dict[str, Dict[str, Any]], N: int) -> Dict[str, Dict[str, Any]]:
    """Construct smoothed role transitions and disambiguated labels for all tracks."""
    registry_tracks: Dict[str, Dict[str, Any]] = {}

    for tid, rec in track_records.items():
        present = rec["frame_present"]
        if not any(present):
            continue
        s_f = min(i for i, p in enumerate(present) if p)
        e_f = max(i for i, p in enumerate(present) if p)

        smoothed = list(rec["frame_roles"])
        smoothed_types = list(rec["frame_types"])
        smoothed_cats = list(rec["frame_categories"])
        smoothed_atoms = list(rec["frame_atoms"])
        for i in range(s_f + 1, e_f):
            if smoothed[i - 1] is not None and smoothed[i + 1] is not None and smoothed[i] is not None:
                if smoothed[i - 1] == smoothed[i + 1] and smoothed[i] != smoothed[i - 1]:
                    smoothed[i] = smoothed[i - 1]
                    smoothed_types[i] = smoothed_types[i - 1]
                    smoothed_cats[i] = smoothed_cats[i - 1]
                    smoothed_atoms[i] = smoothed_atoms[i - 1]

        # Collect distinct role runs persisting for >= 3 frames
        filtered_runs: List[Tuple[str, int]] = []
        cur_role = None
        cur_len = 0
        for i in range(s_f, e_f + 1):
            r = smoothed[i]
            if r == cur_role:
                cur_len += 1
            else:
                if cur_role is not None and cur_len >= 3:
                    if filtered_runs and filtered_runs[-1][0] == cur_role:
                        filtered_runs[-1] = (cur_role, filtered_runs[-1][1] + cur_len)
                    else:
                        filtered_runs.append((cur_role, cur_len))
                cur_role = r
                cur_len = 1
        if cur_role is not None and cur_len >= 3:
            if filtered_runs and filtered_runs[-1][0] == cur_role:
                filtered_runs[-1] = (cur_role, filtered_runs[-1][1] + cur_len)
            else:
                filtered_runs.append((cur_role, cur_len))

        if not filtered_runs:
            valid = [r for r in smoothed[s_f:e_f + 1] if r is not None]
            if valid:
                top_r, top_c = Counter(valid).most_common(1)[0]
                filtered_runs = [(top_r, top_c)]
            else:
                filtered_runs = [(rec["type"], 1)]

        r_first, r_last = filtered_runs[0][0], filtered_runs[-1][0]
        if rec["subspace"] == "deloc_pi" or len(filtered_runs) == 1:
            summary_label = r_first if (r_first == r_last or len(filtered_runs) == 1) else f"{r_first} → {r_last}"
        elif len(filtered_runs) == 2:
            summary_label = f"{r_first} → {r_last}"
        else:
            mid_candidates = [
                (r_m, l_m) for (r_m, l_m) in filtered_runs[1:-1]
                if r_m != r_first and r_m != r_last and (
                    _role_prefix(r_m) != _role_prefix(r_first) or l_m >= 5
                )
            ]
            if mid_candidates:
                best_mid = max(mid_candidates, key=lambda x: x[1])[0]
                summary_label = f"{r_first} → {best_mid} → {r_last}"
            elif r_first == r_last:
                summary_label = r_first
            else:
                summary_label = f"{r_first} → {r_last}"

        role_trans = summary_label

        trk_type = rec["type"]
        trk_cat = rec["category"]
        if rec["subspace"] == "valence_ibo":
            t_counts = Counter(ft for ft in smoothed_types[s_f:e_f + 1] if ft is not None)
            total_t = sum(t_counts.values())
            if t_counts.get("delocalized_pi", 0) == total_t and total_t > 0:
                trk_type, trk_cat = "delocalized_pi", "conjugated_pi"
            elif (t_counts.get("pi", 0) + t_counts.get("delocalized_pi", 0)) >= 3:
                trk_type, trk_cat = "pi", "bonding"
            elif t_counts.get("sigma", 0) >= 3:
                trk_type, trk_cat = "sigma", "bonding"
            elif t_counts:
                trk_type = t_counts.most_common(1)[0][0]
                trk_cat = "non_bonding" if trk_type == "lone_pair" else "bonding"
        elif rec["subspace"] == "antibond":
            t_counts = Counter(ft for ft in smoothed_types[s_f:e_f + 1] if ft is not None)
            total_t = sum(t_counts.values())
            if t_counts.get("delocalized_pistar", 0) == total_t and total_t > 0:
                trk_type, trk_cat = "delocalized_pistar", "conjugated_pi_star"
            elif (t_counts.get("pi_star", 0) + t_counts.get("delocalized_pistar", 0)) >= 3:
                trk_type, trk_cat = "pi_star", "antibonding"
            elif t_counts.get("sigma_star", 0) >= 3:
                trk_type, trk_cat = "sigma_star", "antibonding"
            elif t_counts:
                trk_type = t_counts.most_common(1)[0][0]
                trk_cat = "antibonding"

        atoms_list = sorted(list(rec["atoms_set"]))
        registry_tracks[tid] = {
            "track_id": tid,
            "id": tid,
            "orbital_key": tid,
            "label": summary_label,
            "role_transition": role_trans,
            "short_label": role_trans,
            "category": trk_cat,
            "orbital_class": rec["orbital_class"],
            "type": trk_type,
            "atoms": atoms_list,
            "start_frame": s_f,
            "end_frame": e_f,
            "frame_map": [tid if present[f] else None for f in range(N)],
            "frame_labels": smoothed,
            "frame_types": smoothed_types,
            "frame_categories": smoothed_cats,
            "frame_atoms": smoothed_atoms,
        }

    # Disambiguate identical summary labels within same subspace
    label_groups: Dict[str, List[str]] = {}
    for tid, trk in registry_tracks.items():
        subsp = track_records[tid]["subspace"]
        grp_k = f"{subsp}::{trk['label']}"
        label_groups.setdefault(grp_k, []).append(tid)

    for grp_k, tids in label_groups.items():
        if len(tids) >= 2:
            tids.sort(key=lambda x: (registry_tracks[x]["start_frame"], registry_tracks[x]["track_id"]))
            for idx, tid in enumerate(tids):
                sub = SUBSCRIPTS[idx % len(SUBSCRIPTS)]
                trk = registry_tracks[tid]
                orig_lbl = trk["label"]
                new_lbl = _add_subscript(orig_lbl, sub)
                trk["label"] = new_lbl
                if "→" not in trk["role_transition"]:
                    trk["role_transition"] = new_lbl
                    trk["short_label"] = new_lbl

    return registry_tracks

