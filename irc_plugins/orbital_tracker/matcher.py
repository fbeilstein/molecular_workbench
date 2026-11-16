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
    curr_atom_labels: Optional[List[str]] = None,
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

    if s_cross is not None and len(active_items) > 0:
        expected_len = s_cross.shape[0]
        active_items = [item for item in active_items if len(item["vec"]) == expected_len]

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

    Omega = prev_mat.T @ s_cross @ curr_mat
    cost_matrix = -np.abs(Omega).copy()

    # Apply symmetry bonus to Hungarian cost matrix:
    # Favor matching pi with pi and sigma/lp with sigma/lp to prevent cross-symmetry channel swaps
    def _is_pi_type(typ: Any) -> bool:
        return "pi" in str(typ).lower()

    def _is_sigma_or_lp(typ: Any) -> bool:
        t_str = str(typ).lower()
        return "sigma" in t_str or "lp" in t_str or "lone_pair" in t_str

    for r in range(len(active_items)):
        t_typ = active_items[r].get("type", "")
        for c in range(len(curr_items)):
            c_typ = curr_items[c]["meta"].get("type", "")
            if (_is_pi_type(t_typ) and _is_pi_type(c_typ)) or (_is_sigma_or_lp(t_typ) and _is_sigma_or_lp(c_typ)):
                cost_matrix[r, c] *= 1.25

    r_ind, c_ind = linear_sum_assignment(cost_matrix)
    matched_curr: Set[int] = set()
    matched_active: Set[int] = set()

    for r, c in zip(r_ind, c_ind):
        # Slightly relaxed threshold for coasting tracks that missed a frame
        eff_thresh = threshold * 0.85 if active_items[r].get("missed", 0) > 0 else threshold
        if abs(Omega[r, c]) >= eff_thresh:
            matched_curr.add(c)
            matched_active.add(r)
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
                "atoms": meta.get("atoms", []), "type": meta.get("type", subspace),
                "missed": 0,
            })

    # Coasting buffer: retain unassigned active tracks for up to 2 frames
    MAX_MISSED = 2
    for r in range(len(active_items)):
        if r not in matched_active:
            cur_missed = active_items[r].get("missed", 0)
            if cur_missed < MAX_MISSED:
                if len(active_items[r]["vec"]) != curr_mat.shape[0]:
                    continue
                if curr_atom_labels is not None and active_items[r].get("atoms"):
                    if not set(active_items[r]["atoms"]).issubset(set(curr_atom_labels)):
                        continue
                next_active.append({
                    "tid": active_items[r]["tid"],
                    "vec": active_items[r]["vec"],
                    "atoms": active_items[r].get("atoms", []),
                    "type": active_items[r].get("type", subspace),
                    "missed": cur_missed + 1,
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
                "atoms": meta.get("atoms", []), "type": meta.get("type", subspace),
                "missed": 0,
            })

    return next_active, next_seq


def _role_prefix(lbl: Optional[str]) -> str:
    if not lbl:
        return ""
    for pfx in ("π*_deloc", "π_deloc", "π*", "σ*", "π", "σ", "LP", "Core"):
        if lbl.startswith(pfx):
            return pfx
    return lbl.split("(")[0]


def smooth_and_label_tracks(
    track_records: Dict[str, Dict[str, Any]],
    N: int,
    frame_outputs: Optional[List[Dict[str, Any]]] = None,
) -> Dict[str, Dict[str, Any]]:
    """Construct smoothed role transitions, heal fragmented tracks, and generate disambiguated labels."""
    # Pass 1: Reconnect fragmented and reacting tracks within the same subspace sharing reacting centers
    merged_any = True
    while merged_any:
        merged_any = False
        tids = sorted(list(track_records.keys()))
        for i in range(len(tids)):
            tid_a = tids[i]
            rec_a = track_records.get(tid_a)
            if not rec_a or not any(rec_a["frame_present"]):
                continue
            e_a = max(idx for idx, p in enumerate(rec_a["frame_present"]) if p)
            s_a = min(idx for idx, p in enumerate(rec_a["frame_present"]) if p)
            subsp_a = rec_a["subspace"]
            atoms_a = rec_a["atoms_set"]

            for j in range(len(tids)):
                if i == j:
                    continue
                tid_b = tids[j]
                rec_b = track_records.get(tid_b)
                if not rec_b or not any(rec_b["frame_present"]):
                    continue
                s_b = min(idx for idx, p in enumerate(rec_b["frame_present"]) if p)
                e_b = max(idx for idx, p in enumerate(rec_b["frame_present"]) if p)
                subsp_b = rec_b["subspace"]
                atoms_b = rec_b["atoms_set"]

                if subsp_a != subsp_b:
                    continue
                if e_b <= e_a or s_b <= s_a:
                    continue

                heavy_a = {a for a in atoms_a if not a.startswith("H")}
                heavy_b = {a for a in atoms_b if not a.startswith("H")}
                if not (heavy_a & heavy_b or (len(heavy_a) > 2 and len(heavy_b) > 2)):
                    continue

                gap = s_b - e_a
                if not (-10 <= gap <= 6):
                    continue

                can_merge = False
                flip_b = False
                if frame_outputs is not None:
                    trans_f = max(s_a, min(e_a, s_b))
                    fo = frame_outputs[trans_f]
                    tids_arr = list(fo["tids"])
                    if tid_a in tids_arr and tid_b in tids_arr:
                        va = fo["vecs"][:, tids_arr.index(tid_a)]
                        vb = fo["vecs"][:, tids_arr.index(tid_b)]
                        ov_raw = float(va.T @ vb)
                        if abs(ov_raw) > 0.45:
                            can_merge = True
                            if ov_raw < 0.0:
                                flip_b = True
                    elif gap > 0:
                        fo_a = frame_outputs[e_a]
                        fo_b = frame_outputs[s_b]
                        tids_a = list(fo_a["tids"])
                        tids_b = list(fo_b["tids"])
                        if tid_a in tids_a and tid_b in tids_b:
                            va = fo_a["vecs"][:, tids_a.index(tid_a)]
                            vb = fo_b["vecs"][:, tids_b.index(tid_b)]
                            if len(va) == len(vb):
                                ov_raw = float(va.T @ vb)
                                if ov_raw < 0.0:
                                    flip_b = True
                        can_merge = True
                else:
                    can_merge = True

                if not can_merge:
                    continue

                # Merge rec_b into rec_a starting at switch frame
                switch_f = max(s_b, e_a + 1) if gap > 0 else s_b
                for f in range(switch_f, N):
                    if rec_b["frame_present"][f]:
                        rec_a["frame_present"][f] = True
                        rec_a["frame_roles"][f] = rec_b["frame_roles"][f]
                        rec_a["frame_types"][f] = rec_b["frame_types"][f]
                        rec_a["frame_categories"][f] = rec_b["frame_categories"][f]
                        rec_a["frame_atoms"][f] = rec_b["frame_atoms"][f]
                        rec_a["frame_energies"][f] = rec_b["frame_energies"][f]
                        rec_a["frame_pops"][f] = rec_b["frame_pops"][f]

                rec_a["atoms_set"].update(rec_b["atoms_set"])

                # Remap in frame_outputs if provided (and invert sign if needed for continuity)
                if frame_outputs is not None:
                    for f in range(switch_f, N):
                        fo = frame_outputs[f]
                        tids_arr = fo.get("tids")
                        if tids_arr is not None:
                            mask = (tids_arr == tid_b)
                            if np.any(mask):
                                if flip_b:
                                    fo["vecs"][:, mask] = -fo["vecs"][:, mask]
                                tids_arr[mask] = tid_a

                del track_records[tid_b]
                merged_any = True
                break
            if merged_any:
                break

    # Pass 2: Fill internal coasting gaps (e.g. 1-2 frames) within continuous tracks
    for tid, rec in track_records.items():
        present = rec["frame_present"]
        if not any(present):
            continue
        s_f = min(idx for idx, p in enumerate(present) if p)
        e_f = max(idx for idx, p in enumerate(present) if p)
        for f in range(s_f + 1, e_f):
            if not present[f]:
                p_prev = max(k for k in range(s_f, f) if present[k])
                p_next = min(k for k in range(f + 1, e_f + 1) if present[k])
                if p_next - p_prev <= 4:
                    present[f] = True
                    rec["frame_roles"][f] = rec["frame_roles"][p_prev]
                    rec["frame_types"][f] = rec["frame_types"][p_prev]
                    rec["frame_categories"][f] = rec["frame_categories"][p_prev]
                    rec["frame_atoms"][f] = list(rec["frame_atoms"][p_prev])
                    alpha = (f - p_prev) / (p_next - p_prev)
                    e0 = rec["frame_energies"][p_prev] or 0.0
                    e1 = rec["frame_energies"][p_next] or 0.0
                    rec["frame_energies"][f] = round((1 - alpha) * e0 + alpha * e1, 2)
                    pop0 = np.array(rec["frame_pops"][p_prev] or [])
                    pop1 = np.array(rec["frame_pops"][p_next] or [])
                    if len(pop0) > 0 and len(pop0) == len(pop1):
                        rec["frame_pops"][f] = ((1 - alpha) * pop0 + alpha * pop1).round(2).tolist()

                    # In frame_outputs, interpolate vector
                    if frame_outputs is not None:
                        fo_f = frame_outputs[f]
                        fo_prev = frame_outputs[p_prev]
                        fo_next = frame_outputs[p_next]
                        idx_prev = list(fo_prev["tids"]).index(tid) if tid in fo_prev["tids"] else -1
                        idx_next = list(fo_next["tids"]).index(tid) if tid in fo_next["tids"] else -1
                        if idx_prev >= 0 and idx_next >= 0:
                            v0 = fo_prev["vecs"][:, idx_prev]
                            v1 = fo_next["vecs"][:, idx_next]
                            v1_aligned = v1.copy()
                            if len(v0) == len(v1_aligned):
                                if float(v0.T @ v1_aligned) < 0.0:
                                    v1_aligned = -v1_aligned
                            v_interp = (1 - alpha) * v0 + alpha * v1_aligned
                            norm_v = np.linalg.norm(v_interp)
                            if norm_v > 1e-6:
                                v_interp /= norm_v
                            fo_f["tids"] = np.append(fo_f["tids"], tid)
                            fo_f["vecs"] = np.column_stack([fo_f["vecs"], v_interp.astype(np.float32)])
                            p_interp = np.array(rec["frame_pops"][f], dtype=np.float32)
                            fo_f["pops"] = np.vstack([fo_f["pops"], p_interp]) if len(fo_f["pops"]) > 0 else p_interp.reshape(1, -1)
                            fo_f["energies"] = np.append(fo_f["energies"], rec["frame_energies"][f])

    # Pass 3: Strict forward phase continuity across all adjacent frames for each track
    if frame_outputs is not None:
        for tid, rec in track_records.items():
            pres = [f for f in range(N) if rec["frame_present"][f]]
            for k in range(len(pres) - 1):
                f_prev, f_curr = pres[k], pres[k + 1]
                fo_p = frame_outputs[f_prev]
                fo_c = frame_outputs[f_curr]
                tids_p = list(fo_p["tids"])
                tids_c = list(fo_c["tids"])
                if tid in tids_p and tid in tids_c:
                    idx_p = tids_p.index(tid)
                    idx_c = tids_c.index(tid)
                    vp = fo_p["vecs"][:, idx_p]
                    vc = fo_c["vecs"][:, idx_c]
                    if len(vp) == len(vc) and float(vp.T @ vc) < 0.0:
                        # Invert current frame and all subsequent frames of this track
                        for f_sub in pres[k + 1:]:
                            fo_s = frame_outputs[f_sub]
                            tids_s = list(fo_s["tids"])
                            if tid in tids_s:
                                s_idx = tids_s.index(tid)
                                fo_s["vecs"][:, s_idx] = -fo_s["vecs"][:, s_idx]

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

