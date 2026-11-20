"""
Orbital Matching and Smoothing Module — Conserved Bijective Assignment,
Atom-Indexed Löwdin Cross-Overlap (with Kabsch Rotation Guard), and O(k) Degenerate Alignment.
"""

from collections import Counter, deque
from typing import Any, Callable, Dict, List, Optional, Set, Tuple
import numpy as np
from scipy.optimize import linear_sum_assignment

from .common import SUBSCRIPTS

BOHR_TO_ANG = 0.529177210903


def _add_subscript(label: str, sub: str) -> str:
    for prefix in ("π*_deloc(", "π_deloc(", "π*(", "π(", "σ*(", "σ(", "p*(", "LP("):
        if prefix in label:
            pos = label.index(prefix) + len(prefix) - 1
            return label[:pos] + sub + label[pos:]
    return f"{label} {sub}"


def _kabsch_rotation(P: np.ndarray, Q: np.ndarray) -> np.ndarray:
    """Compute proper 3x3 rotation matrix R (det=+1) minimizing ||P @ R.T - Q||_F."""
    P_c = P - P.mean(axis=0, keepdims=True)
    Q_c = Q - Q.mean(axis=0, keepdims=True)
    H = P_c.T @ Q_c
    U, _, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T
    if np.linalg.det(R) < 0.0:
        Vt[-1, :] *= -1.0
        R = Vt.T @ U.T
    return R


def _build_fragment_p_rotations(
    prev_coords_bohr: Optional[np.ndarray],
    curr_coords_bohr: Optional[np.ndarray],
    prev_atom_labels: Optional[List[str]],
    curr_atom_labels: Optional[List[str]],
    min_rot_deg: float = 25.0,
) -> Dict[str, np.ndarray]:
    """
    Detect disconnected common fragments that rigidly rotated > min_rot_deg between
    frames t-1 and t, returning atom_label -> 3x3 rotation matrix R for p-shell alignment.
    Stationary / slowly-moving fragments (<25°) receive no rotation so genuine 90° orbital
    flips retain exact 0.0 overlap.
    """
    if (
        prev_coords_bohr is None
        or curr_coords_bohr is None
        or not prev_atom_labels
        or not curr_atom_labels
    ):
        return {}

    prev_idx = {lbl: i for i, lbl in enumerate(prev_atom_labels)}
    common = [(lbl, prev_idx[lbl], j) for j, lbl in enumerate(curr_atom_labels) if lbl in prev_idx]
    if len(common) < 3:
        return {}

    curr_ang = np.asarray(curr_coords_bohr, dtype=np.float64) * BOHR_TO_ANG
    prev_ang = np.asarray(prev_coords_bohr, dtype=np.float64) * BOHR_TO_ANG

    n_com = len(common)
    adj: List[List[int]] = [[] for _ in range(n_com)]
    for a in range(n_com):
        ja = common[a][2]
        for b in range(a + 1, n_com):
            jb = common[b][2]
            if float(np.linalg.norm(curr_ang[ja] - curr_ang[jb])) < 2.15:
                adj[a].append(b)
                adj[b].append(a)

    visited = [False] * n_com
    atom_rot: Dict[str, np.ndarray] = {}
    min_rot_rad = np.deg2rad(min_rot_deg)

    for start in range(n_com):
        if visited[start]:
            continue
        comp = []
        q = deque([start])
        visited[start] = True
        while q:
            u = q.popleft()
            comp.append(u)
            for v in adj[u]:
                if not visited[v]:
                    visited[v] = True
                    q.append(v)

        if len(comp) < 3:
            continue
        P = np.array([prev_ang[common[u][1]] for u in comp])
        Q = np.array([curr_ang[common[u][2]] for u in comp])
        # Skip nearly collinear fragments where rotation around axis is ill-posed
        P_c = P - P.mean(axis=0, keepdims=True)
        sv = np.linalg.svd(P_c, compute_uv=False)
        if len(sv) < 2 or sv[1] < 0.15:
            continue
        R = _kabsch_rotation(P, Q)
        cos_theta = np.clip((float(np.trace(R)) - 1.0) * 0.5, -1.0, 1.0)
        theta = float(np.arccos(cos_theta))
        if theta > min_rot_rad:
            for u in comp:
                atom_rot[common[u][0]] = R

    return atom_rot


def prepare_lowdin_common_rows(
    prev_ao_keys: Optional[List[Tuple[str, str]]],
    curr_ao_keys: Optional[List[Tuple[str, str]]],
    prev_coords_bohr: Optional[np.ndarray] = None,
    curr_coords_bohr: Optional[np.ndarray] = None,
    prev_atom_labels: Optional[List[str]] = None,
    curr_atom_labels: Optional[List[str]] = None,
) -> Optional[Tuple[List[int], List[int], Dict[str, np.ndarray]]]:
    """Map shared (atom_label, ao_type) rows between frames t-1 and t for Löwdin overlap."""
    if not prev_ao_keys or not curr_ao_keys:
        return None
    curr_key_to_mu = {k: mu for mu, k in enumerate(curr_ao_keys)}
    common_prev: List[int] = []
    common_curr: List[int] = []
    for mu_p, k in enumerate(prev_ao_keys):
        if k in curr_key_to_mu:
            common_prev.append(mu_p)
            common_curr.append(curr_key_to_mu[k])
    atom_rot = _build_fragment_p_rotations(
        prev_coords_bohr, curr_coords_bohr, prev_atom_labels, curr_atom_labels
    )
    return common_prev, common_curr, atom_rot


def _rotate_prev_lowdin_p_shells(
    P_low_full: np.ndarray,
    prev_ao_keys: List[Tuple[str, str]],
    atom_rot: Dict[str, np.ndarray],
) -> np.ndarray:
    """Apply per-fragment 3x3 Kabsch rotations to p-shell triplets on fast-tumbling fragments."""
    if not atom_rot:
        return P_low_full
    P_rot = P_low_full.copy()
    p_triplets: Dict[Tuple[str, str], Dict[str, int]] = {}
    for mu, (at_lbl, ao_str) in enumerate(prev_ao_keys):
        if at_lbl in atom_rot and len(ao_str) >= 3 and ao_str[-2] == "p" and ao_str[-1] in ("x", "y", "z"):
            shell = ao_str[:-1]
            p_triplets.setdefault((at_lbl, shell), {})[ao_str[-1]] = mu
    for (at_lbl, _shell), xyz_map in p_triplets.items():
        if len(xyz_map) == 3:
            idxs = [xyz_map["x"], xyz_map["y"], xyz_map["z"]]
            P_rot[idxs, :] = atom_rot[at_lbl] @ P_low_full[idxs, :]
    return P_rot


def _align_degenerate_ok(W: np.ndarray) -> Optional[np.ndarray]:
    """Compute optimal O(k) Procrustes rotation R_w = Vt.T @ U.T if overlap block is well-conditioned."""
    U_w, s_w, Vt_w = np.linalg.svd(W)
    if np.all(s_w > 0.50):
        return Vt_w.T @ U_w.T
    return None


def _align_1to2_degenerate_pair(w: np.ndarray) -> Optional[np.ndarray]:
    """Align a splitting 1->2 cylindrical degenerate pair so col 0 follows prev_vec."""
    w = w.ravel()
    norm_w = float(np.linalg.norm(w))
    if norm_w >= 0.75:
        u0 = w / norm_w
        u1 = np.array([-u0[1], u0[0]])
        return np.column_stack([u0, u1])
    return None


def _align_degenerate_subspace_clusters(
    prev_items: List[Dict[str, Any]],
    curr_items: List[Dict[str, Any]],
    prev_mat: np.ndarray,
    curr_mat: np.ndarray,
    prev_ov_mat: np.ndarray,
    curr_ov_mat: np.ndarray,
    overlap_fn: Callable[[np.ndarray, np.ndarray], np.ndarray],
    pi_type_keyword: Optional[str],
    subspace: str,
) -> None:
    """
    Detect and O(k)-Procrustes-align degenerate orbital clusters across frames:
      - 6-atom aromatic delocalized π / π* pairs (k=2, |ΔE| < 0.25 eV)
      - 2-atom linear/cylindrical π / π* pairs (2->2 if |ΔE| < 0.08 eV, 1->2 if |ΔE| < 0.18 eV)
      - 1-atom degenerate lone-pair (e.g. Cl⁻ 3p) or core (e.g. Cl 2p) clusters (k=2,3, |ΔE| < 0.08 / 0.15 eV)
    Rotates BOTH curr_mat (AO basis) and curr_ov_mat (Löwdin basis) identically.
    """
    Q = len(curr_items)
    if Q < 2 or not prev_items:
        return

    # 1. Group current items by (family, sorted_atoms_tuple)
    groups: Dict[Tuple[str, Tuple[str, ...]], List[int]] = {}
    for idx, it in enumerate(curr_items):
        meta = it["meta"]
        typ = meta.get("type", "")
        ats = tuple(sorted(meta.get("atoms", [])))
        if len(ats) == 1 and (typ == "lone_pair" or typ == "core"):
            groups.setdefault((typ, ats), []).append(idx)
        elif pi_type_keyword is not None and len(ats) in (2, 6):
            # Support unified virtual pool containing both localized pi_star and delocalized_pistar
            is_pi_match = (
                pi_type_keyword in typ
                or (subspace == "antibond" and typ in ("pi_star", "delocalized_pistar"))
                or (subspace == "valence_ibo" and typ in ("pi", "delocalized_pi"))
            )
            if is_pi_match:
                fam = "deloc_pi" if len(ats) == 6 else "local_pi"
                groups.setdefault((fam, ats), []).append(idx)

    for (fam, ats), idxs in groups.items():
        if len(idxs) < 2:
            continue
        # Cluster by energy degeneracy
        idxs_sorted = sorted(idxs, key=lambda c: float(curr_items[c]["meta"].get("energy_ev", 0.0)))
        if fam == "deloc_pi":
            max_dE = 0.25
        elif fam == "local_pi":
            max_dE = 0.18
        elif fam == "core":
            max_dE = 0.15
        else:
            max_dE = 0.08

        clusters: List[List[int]] = [[idxs_sorted[0]]]
        for c in idxs_sorted[1:]:
            e_cur = float(curr_items[c]["meta"].get("energy_ev", 0.0))
            e_prev = float(curr_items[clusters[-1][-1]]["meta"].get("energy_ev", 0.0))
            if abs(e_cur - e_prev) < max_dE:
                clusters[-1].append(c)
            else:
                clusters.append([c])

        for cl in clusters:
            k = len(cl)
            if k < 2:
                continue
            energies = [float(curr_items[c]["meta"].get("energy_ev", 0.0)) for c in cl]
            e_span = max(energies) - min(energies)
            e_mid = float(np.mean(energies))

            prev_matches = []
            for r, av in enumerate(prev_items):
                av_typ = av.get("type", "")
                av_ats = tuple(sorted(av.get("atoms", [])))
                if av_ats != ats:
                    continue
                if fam in ("lone_pair", "core") and av_typ != fam:
                    continue
                if fam == "deloc_pi" and av_typ not in ("delocalized_pi", "delocalized_pistar"):
                    continue
                if fam == "local_pi" and av_typ not in ("pi", "pi_star"):
                    continue
                if fam in ("deloc_pi", "lone_pair", "core"):
                    if abs(float(av.get("energy_ev", 0.0)) - e_mid) >= 1.0:
                        continue
                prev_matches.append(r)

            R_w = None
            if len(prev_matches) == k and (fam != "local_pi" or e_span < 0.08):
                W = overlap_fn(prev_ov_mat[:, prev_matches], curr_ov_mat[:, cl])
                R_w = _align_degenerate_ok(W)
            elif k == 2 and len(prev_matches) == 1 and fam == "local_pi":
                w = overlap_fn(prev_ov_mat[:, [prev_matches[0]]], curr_ov_mat[:, cl])
                R_w = _align_1to2_degenerate_pair(w)

            if R_w is not None:
                curr_mat[:, cl] = curr_mat[:, cl] @ R_w
                curr_ov_mat[:, cl] = curr_ov_mat[:, cl] @ R_w
                for j, c_idx in enumerate(cl):
                    curr_items[c_idx]["vec"] = curr_mat[:, c_idx].copy()
                    if "vec_low" in curr_items[c_idx]:
                        curr_items[c_idx]["vec_low"] = curr_ov_mat[:, c_idx].copy()


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
    conserved: bool = True,
    prev_atom_labels: Optional[List[str]] = None,
    curr_atom_labels: Optional[List[str]] = None,
    lowdin_map: Optional[Tuple[List[int], List[int], Dict[str, np.ndarray]]] = None,
    prev_ao_keys: Optional[List[Tuple[str, str]]] = None,
) -> Tuple[List[Dict[str, Any]], int]:
    """
    Subspace-specific orbital assignment and track propagation.

    When `conserved=True` and Löwdin vectors (`vec_low`) are provided:
      1. Partitions departed-molecule tracks and arrived-molecule orbitals cleanly by
         atom ownership (>70% Löwdin population on departed/arrived atoms).
      2. Matches all remaining common-atom orbitals via atom-indexed Löwdin cross-overlap
         (translation-invariant, rotation-sensitive, with >25° Kabsch guard).
      3. Enforces a strict 1-to-1 bijection whenever persisting count == continuing count
         (exact orbital conservation with zero spurious births or deaths).
    """
    if not curr_items:
        return [], next_seq

    next_active: List[Dict[str, Any]] = []

    def _emit(tid: str, vec: np.ndarray, item: Dict[str, Any], vec_low: Optional[np.ndarray] = None) -> None:
        meta = item["meta"]
        record_fn(
            tid, vec, item["pop"], meta.get("energy_ev", 0.0), meta.get("label", ""),
            meta.get("type", subspace), meta.get("category", subspace),
            meta.get("orbital_class", "bonding"), meta.get("atoms", []), subspace
        )
        entry: Dict[str, Any] = {
            "tid": tid,
            "vec": vec,
            "atoms": meta.get("atoms", []),
            "type": meta.get("type", subspace),
            "label": meta.get("label", ""),
            "energy_ev": float(meta.get("energy_ev", 0.0)),
        }
        if vec_low is not None:
            entry["vec_low"] = vec_low
        elif "vec_low" in item:
            entry["vec_low"] = item["vec_low"].copy()
        if "pop_dict" in item:
            entry["pop_dict"] = item["pop_dict"]
        next_active.append(entry)

    if t == 0 or len(active_items) == 0 or (s_cross is None and lowdin_map is None):
        for item in curr_items:
            tid = f"{prefix}{next_seq:0{pad_zeros}d}"
            next_seq += 1
            _emit(tid, item["vec"], item)
        return next_active, next_seq

    # Step 1: Partition departing / arriving molecule orbitals at gluing boundaries
    departed: Set[str] = set()
    arrived: Set[str] = set()
    if prev_atom_labels is not None and curr_atom_labels is not None:
        prev_set, curr_set = set(prev_atom_labels), set(curr_atom_labels)
        departed = prev_set - curr_set
        arrived = curr_set - prev_set

    persisting_idx: List[int] = []
    for r, av in enumerate(active_items):
        pdict = av.get("pop_dict")
        if departed and pdict is not None:
            dep_pop = sum(pdict.get(at, 0.0) for at in departed)
            if dep_pop > 0.70:
                continue  # Track cleanly dies with departed molecule
        persisting_idx.append(r)

    continuing_idx: List[int] = []
    newborn_idx: List[int] = []
    for c, it in enumerate(curr_items):
        pdict = it.get("pop_dict")
        if arrived and pdict is not None:
            arr_pop = sum(pdict.get(at, 0.0) for at in arrived)
            if arr_pop > 0.70:
                newborn_idx.append(c)
                continue  # Track cleanly born with arrived molecule
        continuing_idx.append(c)

    # Step 2: Build overlap operator (Löwdin AO-index preferred; fallback to s_cross)
    use_lowdin = (
        lowdin_map is not None
        and prev_ao_keys is not None
        and all("vec_low" in active_items[r] for r in persisting_idx)
        and all("vec_low" in curr_items[c] for c in continuing_idx)
    )

    if persisting_idx and continuing_idx:
        persisting_items = [active_items[r] for r in persisting_idx]
        continuing_items = [curr_items[c] for c in continuing_idx]

        prev_mat = np.column_stack([it["vec"] for it in persisting_items])
        curr_mat = np.column_stack([it["vec"] for it in continuing_items])

        if use_lowdin and lowdin_map is not None and prev_ao_keys is not None:
            common_prev, common_curr, atom_rot = lowdin_map
            prev_low_full = np.column_stack([it["vec_low"] for it in persisting_items])
            prev_low_rot = _rotate_prev_lowdin_p_shells(prev_low_full, prev_ao_keys, atom_rot)
            curr_low_full = np.column_stack([it["vec_low"] for it in continuing_items])

            def _ov_fn(P_full: np.ndarray, C_full: np.ndarray) -> np.ndarray:
                return P_full[common_prev, :].T @ C_full[common_curr, :]

            _align_degenerate_subspace_clusters(
                persisting_items, continuing_items, prev_mat, curr_mat,
                prev_low_rot, curr_low_full, _ov_fn, pi_type_keyword, subspace
            )
            Omega = _ov_fn(prev_low_rot, curr_low_full)
        else:
            assert s_cross is not None
            def _ov_ao(P_ao: np.ndarray, C_ao: np.ndarray) -> np.ndarray:
                return P_ao.T @ s_cross @ C_ao

            _align_degenerate_subspace_clusters(
                persisting_items, continuing_items, prev_mat, curr_mat,
                prev_mat, curr_mat, _ov_ao, pi_type_keyword, subspace
            )
            Omega = _ov_ao(prev_mat, curr_mat)
            curr_low_full = None

        r_ind, c_ind = linear_sum_assignment(-np.abs(Omega))
        bijective = conserved and (len(persisting_idx) == len(continuing_idx))
        matched_local: Set[int] = set()

        for ri, ci in zip(r_ind, c_ind):
            ov_signed = float(Omega[ri, ci])
            if bijective or abs(ov_signed) >= threshold:
                matched_local.add(ci)
                tid = persisting_items[ri]["tid"]
                v_ao = curr_mat[:, ci].copy()
                v_low = curr_low_full[:, ci].copy() if curr_low_full is not None else None
                if ov_signed < 0.0:
                    v_ao = -v_ao
                    if v_low is not None:
                        v_low = -v_low
                _emit(tid, v_ao, continuing_items[ci], vec_low=v_low)

        for ci in range(len(continuing_items)):
            if ci not in matched_local:
                tid = f"{prefix}{next_seq:0{pad_zeros}d}"
                next_seq += 1
                v_low = curr_low_full[:, ci].copy() if curr_low_full is not None else None
                _emit(tid, curr_mat[:, ci].copy(), continuing_items[ci], vec_low=v_low)
    else:
        for c in continuing_idx:
            tid = f"{prefix}{next_seq:0{pad_zeros}d}"
            next_seq += 1
            _emit(tid, curr_items[c]["vec"], curr_items[c])

    # Step 3: Spawn tracks cleanly for arrived molecules
    for c in newborn_idx:
        tid = f"{prefix}{next_seq:0{pad_zeros}d}"
        next_seq += 1
        _emit(tid, curr_items[c]["vec"], curr_items[c])

    return next_active, next_seq


def _role_prefix(lbl: Optional[str]) -> str:
    if not lbl:
        return ""
    for pfx in ("π*_deloc", "π_deloc", "π*", "σ*", "p*", "π", "σ", "LP", "Core"):
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
        if len(filtered_runs) == 1:
            summary_label = r_first
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
