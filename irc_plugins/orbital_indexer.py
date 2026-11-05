"""
Orbital Indexer and Multi-Step Pathway Tracker
Maps localized orbitals continuously across multi-step IRC reaction coordinates
and stitched intermediate junctions using atom correspondence maps and global track IDs.
"""

import re
from typing import Any, Callable, Dict, List, Optional, Tuple
import numpy as np
from .orbital_selection import clean_role_label


def get_atom_p_vectors(fd: Optional[Dict[str, Any]], k: str) -> Dict[int, np.ndarray]:
    """Extract atomic 3D p-vector coefficients from localized molecular orbital."""
    if not fd or "mol" not in fd or "vecs" not in fd or k not in fd["vecs"]:
        return {}
    mol = fd["mol"]
    v = fd["vecs"][k]
    atom_vecs: Dict[int, np.ndarray] = {}
    for i, aol in enumerate(mol.ao_labels()):
        parts = aol.split()
        at_idx = int(parts[0])
        orb = parts[2]
        if "2p" in orb or "3p" in orb or "4p" in orb or "p" in orb:
            if at_idx not in atom_vecs:
                atom_vecs[at_idx] = np.zeros(3)
            if "px" in orb: atom_vecs[at_idx][0] += v[i]
            elif "py" in orb: atom_vecs[at_idx][1] += v[i]
            elif "pz" in orb: atom_vecs[at_idx][2] += v[i]
    return atom_vecs


def compute_orbital_spatial_overlap(
    fd_prev: Optional[Dict[str, Any]],
    k_prev: str,
    fd_next: Optional[Dict[str, Any]],
    k_next: str,
    p2n: Dict[int, int],
    R: Optional[np.ndarray] = None
) -> float:
    """Calculate 3D spatial directional overlap between localized orbitals across aligned junction frames."""
    v_prev = get_atom_p_vectors(fd_prev, k_prev)
    v_next = get_atom_p_vectors(fd_next, k_next)
    if not v_prev or not v_next:
        return 0.0
    dots = []
    for at_prev, vec_prev in v_prev.items():
        at_next = p2n.get(at_prev)
        if at_next in v_next:
            vec_n = v_next[at_next]
            norm_p = np.linalg.norm(vec_prev)
            norm_n = np.linalg.norm(vec_n)
            if norm_p > 0.05 and norm_n > 0.05:
                vec_n_rot = (R @ vec_n) if R is not None else vec_n
                dot = abs(float(np.dot(vec_prev, vec_n_rot))) / (norm_p * norm_n)
                dots.append(dot)
    return float(np.mean(dots)) if dots else 0.0


def get_frame_groups(megajson: Dict[str, Any]) -> List[Tuple[int, int]]:
    """Determine contiguous frame groups where molecular atom counts are uniform."""
    frames = megajson.get("frames", [])
    if not frames:
        return []
    groups, curr_start = [], 0
    for i in range(1, len(frames)):
        if len(frames[i].get("symbols", [])) != len(frames[curr_start].get("symbols", [])):
            groups.append((curr_start, i - 1))
            curr_start = i
    groups.append((curr_start, len(frames) - 1))
    return groups


def _get_frame_dict(fi: int, megajson: Dict[str, Any], scf_cache: Optional[Dict[int, Dict[str, Any]]] = None, ensure_fn: Optional[Callable] = None) -> Optional[Dict[str, Any]]:
    if scf_cache and fi in scf_cache:
        return scf_cache[fi]
    return ensure_fn(fi, megajson) if ensure_fn else None


def _atom_nums(atoms: List[Any]) -> List[int]:
    return [int(m.group()) for a in atoms if (m := re.search(r'\d+', str(a)))]


def _match_cand(
    m_src: Optional[Dict[str, Any]],
    p2n: Dict[int, int],
    pool: Dict[str, Dict[str, Any]],
    is_star: bool,
    fd_src: Optional[Dict[str, Any]] = None,
    k_src: Optional[str] = None,
    fd_cand: Optional[Dict[str, Any]] = None,
    R: Optional[np.ndarray] = None,
    is_backward: bool = False,
    src_labels: Optional[List[str]] = None,
    cand_labels: Optional[List[str]] = None
) -> Optional[str]:
    if not m_src or not m_src.get("atoms"):
        return None

    src_atoms = [str(a) for a in m_src.get("atoms", [])]
    if src_labels and cand_labels:
        mapped_atoms = [a for a in src_atoms if a in cand_labels]
        if not mapped_atoms:
            return None  # All atoms of this orbital have left the system
        s_mapped = set(mapped_atoms)
        at_src_indices = [src_labels.index(a) for a in src_atoms if a in src_labels]
        mapped_indices = [p2n.get(i) for i in at_src_indices]
    else:
        at_src_indices = _atom_nums(src_atoms)
        mapped_indices = [p2n.get(a) for a in at_src_indices]
        if any(x is None for x in mapped_indices):
            return None
        s_mapped = set(mapped_indices)
        if not s_mapped:
            return None

    # Check for analytical quantum-mechanical PySCF cross-overlap if wavefunctions are available
    qm_ovlps: Dict[str, float] = {}
    if fd_src and fd_cand and "mol" in fd_src and "mol" in fd_cand and "vecs" in fd_src and "vecs" in fd_cand and k_src and k_src in fd_src["vecs"]:
        try:
            from pyscf import gto
            v_src = fd_src["vecs"][k_src]
            s_cross = gto.intor_cross("int1e_ovlp", fd_src["mol"], fd_cand["mol"]) if not is_backward else gto.intor_cross("int1e_ovlp", fd_cand["mol"], fd_src["mol"])
            for k2, v2 in fd_cand["vecs"].items():
                if k2 in pool:
                    ov_qm = abs(float(v_src.T @ s_cross @ v2)) if not is_backward else abs(float(v2.T @ s_cross @ v_src))
                    qm_ovlps[k2] = ov_qm
        except Exception:
            pass

    best_k, best_score = None, -9999.0
    for k2, m2 in pool.items():
        if (m2.get("type") == "core") != (m_src.get("type") == "core"):
            continue
        if is_star != ("STAR" in k2 or "star" in m2.get("type", "")):
            continue

        if src_labels and cand_labels:
            s2 = set(str(a) for a in m2.get("atoms", []))
        else:
            s2 = set(_atom_nums(m2.get("atoms", [])))

        qm_ov = qm_ovlps.get(k2, 0.0)

        # Allow matching if atom sets overlap/subset, or if PySCF QM overlap is high
        atom_match = bool(s2 and (s_mapped == s2 or s2.issubset(s_mapped) or s_mapped.issubset(s2) or (s_mapped & s2)))
        if not atom_match and qm_ov < 0.65:
            continue

        if s_mapped == s2:
            sc = 120.0
        elif s2.issubset(s_mapped) or s_mapped.issubset(s2):
            sc = 60.0
        elif s_mapped & s2:
            sc = 30.0
        else:
            sc = 0.0

        if qm_ov > 0.0:
            sc += qm_ov * 250.0  # huge reward for verified physical overlap

        t1, t2 = m_src.get("type", ""), m2.get("type", "")
        if t1 == t2:
            sc += 40.0
        elif (t1 in ("pi", "lone_pair") and t2 in ("pi", "lone_pair")):
            sc += 20.0

        p1, p2 = m_src.get("pop"), m2.get("pop")
        if p1 and p2 and at_src_indices and mapped_indices:
            sc -= sum(abs(p1[a1] - p2[a2]) for a1, a2 in zip(at_src_indices, mapped_indices) if a1 is not None and a2 is not None and a1 < len(p1) and a2 < len(p2)) * 30.0

        e1, e2 = m_src.get("energy_ev"), m2.get("energy_ev")
        if e1 is not None and e2 is not None:
            sc -= abs(e1 - e2) * 1.0

        if fd_src and fd_cand and k_src and k2:
            if not is_backward:
                ov = compute_orbital_spatial_overlap(fd_src, k_src, fd_cand, k2, p2n, R)
            else:
                p2n_rev = {v: k for k, v in p2n.items()}
                ov = compute_orbital_spatial_overlap(fd_cand, k2, fd_src, k_src, p2n_rev, R)
            if ov > 0:
                if ov < 0.25:
                    sc -= 80.0
                else:
                    sc += ov * 120.0

        if sc > best_score:
            best_score, best_k = sc, k2

    return best_k if best_score > 20.0 else None


def build_track_frame_map(
    megajson: Dict[str, Any],
    ref_frame: int,
    orb_key: str,
    scf_cache: Optional[Dict[int, Dict[str, Any]]] = None,
    ensure_fn: Optional[Callable] = None,
    is_star: bool = False
) -> Tuple[List[Optional[str]], List[Optional[str]]]:
    """Map an orbital across all pathway frames via cross-frame continuity and junction stitching."""
    frames = megajson.get("frames", [])
    n_frames = len(frames)
    if n_frames == 0 or not orb_key:
        return [None] * n_frames, []

    groups = get_frame_groups(megajson)
    if orb_key.startswith("HOMO") or orb_key.startswith("LUMO") or "esp" in orb_key:
        k = "track_esp_charge" if "esp" in orb_key else orb_key
        return [k] * n_frames, [k] * len(groups)

    ref_g = next((i for i, (gs, ge) in enumerate(groups) if gs <= ref_frame <= ge), 0)
    g_keys = [None] * len(groups)
    g_keys[ref_g] = orb_key
    rpt = megajson.get("stitching_alignment_report", [])
    cands = [r["all_candidates"][0] for r in rpt if r.get("all_candidates")]

    def _meta_for(fi: int, k: str) -> Optional[Dict[str, Any]]:
        fd = _get_frame_dict(fi, megajson, scf_cache, ensure_fn)
        if not fd:
            return None
        pool = fd.get("antibonds_meta", {}) if ("STAR" in k or is_star) else fd.get("ibos_meta", {})
        return pool.get(k)

    # Forward tracking across junctions
    for g in range(ref_g, len(groups) - 1):
        if not g_keys[g]:
            break
        fi_src = groups[g][1]
        fi_next = groups[g + 1][0]
        m = _meta_for(fi_src, g_keys[g])
        if not m:
            break

        f_src = megajson["frames"][fi_src] if fi_src < len(megajson["frames"]) else {}
        f_next = megajson["frames"][fi_next] if fi_next < len(megajson["frames"]) else {}
        src_lbls = f_src.get("atom_labels") or []
        next_lbls = f_next.get("atom_labels") or []

        p2n = {}
        if src_lbls and next_lbls:
            for i, lbl in enumerate(src_lbls):
                if lbl in next_lbls:
                    p2n[i] = next_lbls.index(lbl)
        elif g < len(cands):
            cand_rpt = cands[g]
            p2n = dict(zip(cand_rpt.get("prev_indices", []), cand_rpt.get("next_indices", [])))

        R_mat = None
        if not (src_lbls and next_lbls) and g < len(cands) and "R" in cands[g]:
            R_mat = np.array(cands[g]["R"])

        fd_src = _get_frame_dict(fi_src, megajson, scf_cache, ensure_fn)
        fd_next = _get_frame_dict(fi_next, megajson, scf_cache, ensure_fn)
        if not fd_next:
            break
        pool = fd_next.get("antibonds_meta", {}) if ("STAR" in g_keys[g] or is_star) else fd_next.get("ibos_meta", {})
        g_keys[g + 1] = _match_cand(
            m, p2n, pool, is_star or ("STAR" in g_keys[g]),
            fd_src=fd_src, k_src=g_keys[g], fd_cand=fd_next, R=R_mat, is_backward=False,
            src_labels=src_lbls, cand_labels=next_lbls
        )

    # Backward tracking across junctions
    for g in range(ref_g, 0, -1):
        if not g_keys[g]:
            break
        fi_src = groups[g][0]
        fi_prev = groups[g - 1][1]
        m = _meta_for(fi_src, g_keys[g])
        if not m:
            break

        f_src = megajson["frames"][fi_src] if fi_src < len(megajson["frames"]) else {}
        f_prev = megajson["frames"][fi_prev] if fi_prev < len(megajson["frames"]) else {}
        src_lbls = f_src.get("atom_labels") or []
        prev_lbls = f_prev.get("atom_labels") or []

        n2p = {}
        if src_lbls and prev_lbls:
            for i, lbl in enumerate(src_lbls):
                if lbl in prev_lbls:
                    n2p[i] = prev_lbls.index(lbl)
        elif (g - 1) < len(cands):
            cand_rpt = cands[g - 1]
            n2p = dict(zip(cand_rpt.get("next_indices", []), cand_rpt.get("prev_indices", [])))

        R_mat = None
        if not (src_lbls and prev_lbls) and (g - 1) < len(cands) and "R" in cands[g - 1]:
            R_mat = np.array(cands[g - 1]["R"])

        fd_src = _get_frame_dict(fi_src, megajson, scf_cache, ensure_fn)
        fd_prev = _get_frame_dict(fi_prev, megajson, scf_cache, ensure_fn)
        if not fd_prev:
            break
        pool = fd_prev.get("antibonds_meta", {}) if ("STAR" in g_keys[g] or is_star) else fd_prev.get("ibos_meta", {})
        g_keys[g - 1] = _match_cand(
            m, n2p, pool, is_star or ("STAR" in g_keys[g]),
            fd_src=fd_src, k_src=g_keys[g], fd_cand=fd_prev, R=R_mat, is_backward=True,
            src_labels=src_lbls, cand_labels=prev_lbls
        )

    frame_map: List[Optional[str]] = [None] * n_frames
    for g, k_g in enumerate(g_keys):
        if k_g:
            for fi in range(groups[g][0], groups[g][1] + 1):
                frame_map[fi] = k_g
    return frame_map, g_keys


def resolve_track_roles(
    megajson: Dict[str, Any],
    frame_map: List[Optional[str]],
    default_title: str = "",
    scf_cache: Optional[Dict[int, Dict[str, Any]]] = None,
    ensure_fn: Optional[Callable] = None
) -> Tuple[str, str, Optional[str], Optional[str]]:
    """Determine start and end chemical roles to describe orbital transformations (e.g., σ → LP)."""
    if default_title and any(x in default_title.upper() for x in ("HOMO", "LUMO", "ESP")):
        return default_title, default_title, default_title, default_title

    valid_frames = [i for i, k in enumerate(frame_map) if k is not None]
    if not valid_frames:
        return default_title, default_title, None, None

    cached = [fi for fi in valid_frames if scf_cache and fi in scf_cache]
    f_first, f_last = (cached[0] if cached else valid_frames[0]), (cached[-1] if cached else valid_frames[-1])
    fd_first, fd_last = _get_frame_dict(f_first, megajson, scf_cache, None), _get_frame_dict(f_last, megajson, scf_cache, None)

    def _get_m(fd, k):
        if not fd or not k: return {}
        return {**fd.get("ibos_meta", {}), **fd.get("canonical_pi_meta", {}), **fd.get("antibonds_meta", {})}.get(k, {})

    s_role = clean_role_label(_get_m(fd_first, frame_map[f_first]))
    e_role = clean_role_label(_get_m(fd_last, frame_map[f_last]))

    groups = get_frame_groups(megajson)
    g_from = next((i for i, (gs, ge) in enumerate(groups) if gs <= f_last <= ge), 0)
    g_to = next((i for i, (gs, ge) in enumerate(groups) if gs <= f_first <= ge), 0)
    f_first_fr = megajson.get("frames", [{}])[f_first] if f_first < len(megajson.get("frames", [])) else {}
    f_last_fr = megajson.get("frames", [{}])[f_last] if f_last < len(megajson.get("frames", [])) else {}
    has_global_labels = bool(f_first_fr.get("atom_labels") and f_last_fr.get("atom_labels"))
    if not has_global_labels and e_role and g_from > g_to:
        cands = [r["all_candidates"][0] for r in megajson.get("stitching_alignment_report", []) if r.get("all_candidates")]
        n2p_c: Dict[int, int] = {}
        for g in range(g_from, g_to, -1):
            if (g - 1) < len(cands):
                n2p = dict(zip(cands[g - 1].get("next_indices", []), cands[g - 1].get("prev_indices", [])))
                n2p_c = {k: n2p.get(v, v) for k, v in n2p_c.items()} if n2p_c else n2p
        e_role = re.sub(r'([A-Za-z]+)(\d+)', lambda m: f"{m.group(1)}{n2p_c.get(int(m.group(2)), int(m.group(2)))}", e_role)

    s_type, e_type = s_role.split('(')[0].strip() if s_role else "", e_role.split('(')[0].strip() if e_role else ""
    s_atoms = set(re.findall(r'[A-Za-z]+\d+', s_role or ""))
    e_atoms = set(re.findall(r'[A-Za-z]+\d+', e_role or ""))
    if s_role and e_role and s_type != e_type:
        if s_atoms == e_atoms and {s_type, e_type} == {"π*", "σ*"}:
            # Avoid mislabeling persistent double-bond pi* antibonds as pi* -> sigma*
            label, trans = (s_role if s_type == "π*" else e_role), "π*"
        else:
            label, trans = f"{s_role} → {e_role}", f"{s_type} → {e_type}"
    else:
        label = s_role or e_role or default_title
        trans = s_type or e_type or (label.split("(")[0] if "(" in label else label)
    return label, trans, s_role, e_role


def build_global_track_catalog(
    megajson: Dict[str, Any],
    scf_cache: Optional[Dict[int, Dict[str, Any]]] = None,
    ensure_fn: Optional[Callable] = None,
    rebuild: bool = False
) -> Dict[str, Dict[str, Any]]:
    """Construct global electron tracks uniquely identifying each orbital across all reaction steps."""
    frames = megajson.get("frames", [])
    n_frames = len(frames)
    if n_frames == 0:
        return {}

    groups = get_frame_groups(megajson)
    n_g = len(groups)
    catalog: Dict[str, Dict[str, Any]] = {} if rebuild else (megajson.get("global_electron_tracks") or {})
    if not rebuild and catalog and len(catalog) > 3:
        return catalog

    if not catalog:
        catalog = {
            "trk_esp": {"track_id": "trk_esp", "name": "ESP Surface", "label": "⚡ Electrostatic Surface & Charges (ESP)", "category": "esp", "orbital_class": "esp", "color": "#c084fc", "atoms": [], "group_keys": ["track_esp_charge"] * n_g, "frame_map": ["track_esp_charge"] * n_frames},
            "trk_homo": {"track_id": "trk_homo", "name": "HOMO", "label": "HOMO", "category": "frontier", "orbital_class": "bonding", "color": "#38bdf8", "atoms": [], "group_keys": ["HOMO"] * n_g, "frame_map": ["HOMO"] * n_frames},
            "trk_lumo": {"track_id": "trk_lumo", "name": "LUMO", "label": "LUMO", "category": "frontier", "orbital_class": "antibonding", "color": "#fb923c", "atoms": [], "group_keys": ["LUMO"] * n_g, "frame_map": ["LUMO"] * n_frames}
        }

    def _can_merge_tracks(trk_gk: List[Optional[str]], new_gk: List[Optional[str]], trk_atoms: List[Any], new_atoms: List[Any]) -> bool:
        has_match = False
        for gi in range(min(len(trk_gk), len(new_gk))):
            k1, k2 = trk_gk[gi], new_gk[gi]
            if k1 is not None and k2 is not None:
                if k1 != k2:
                    return False  # Strict conflict: an orbital track cannot be two different IBOs in the same step
                has_match = True
        if not has_match:
            return False
        # Exact group key match always merges
        if all((k1 == k2 or k1 is None or k2 is None) for k1, k2 in zip(trk_gk, new_gk)):
            return True
        # Do not merge tracks if atoms are completely disjoint across the step
        s1 = set(str(a) for a in trk_atoms)
        s2 = set(str(a) for a in new_atoms)
        if s1 and s2 and not (s1 & s2) and (trk_gk[0] is not None and new_gk[0] is not None):
            return False
        return True

    # Ensure at least the start of each step is sampled if ensure_fn is available
    avail_frames = sorted(list((scf_cache or {}).keys()))
    if ensure_fn:
        step_starts = [g[0] for g in groups if g[0] not in avail_frames]
        for st in step_starts:
            ensure_fn(st, megajson)
        avail_frames = sorted(list((scf_cache or {}).keys()))

    for fi in avail_frames:
        fd = scf_cache[fi]
        ref_g = next((i for i, (gs, ge) in enumerate(groups) if gs <= fi <= ge), 0)
        pools = [("bonding", fd.get("ibos_meta", {})), ("pi", fd.get("canonical_pi_meta", {})), ("star", fd.get("antibonds_meta", {}))]
        for p_type, pool in pools:
            for k, m in pool.items():
                if any(trk.get("group_keys", [None] * n_g)[ref_g] == k for trk in catalog.values()):
                    continue
                is_star = ("STAR" in k) or (p_type == "star")
                f_map, g_keys = build_track_frame_map(megajson, fi, k, scf_cache, ensure_fn, is_star=is_star)

                matched_trk = None
                for tid, trk in list(catalog.items()):
                    if tid in ("trk_esp", "trk_homo", "trk_lumo"): continue
                    if _can_merge_tracks(trk.get("group_keys", []), g_keys, trk.get("atoms", []), m.get("atoms", [])):
                        matched_trk = trk
                        break

                if matched_trk:
                    mgk = matched_trk.setdefault("group_keys", [None] * n_g)
                    mfm = matched_trk.setdefault("frame_map", [None] * n_frames)
                    for gi, gk in enumerate(g_keys):
                        if gk: mgk[gi] = gk
                    for f_i, fk in enumerate(f_map):
                        if fk and not mfm[f_i]: mfm[f_i] = fk
                    cur_at = set(str(a) for a in matched_trk.get("atoms", []))
                    new_at = set(str(a) for a in m.get("atoms", []))
                    if new_at:
                        matched_trk["atoms"] = sorted(list(cur_at | new_at), key=lambda x: int(re.search(r'\d+', str(x)).group()) if re.search(r'\d+', str(x)) else str(x))
                    lbl, trans, _, _ = resolve_track_roles(megajson, mfm, matched_trk.get("label", m.get("label", k)), scf_cache, ensure_fn)
                    matched_trk["label"] = matched_trk["name"] = lbl
                    matched_trk["role_transition"] = trans
                else:
                    lbl, trans, _, _ = resolve_track_roles(megajson, f_map, m.get("label", k), scf_cache, ensure_fn)
                    t = m.get("type", "")
                    cat = "core" if t == "core" else ("antibonding" if (is_star or "star" in t) else ("lone_pairs" if t == "lone_pair" else ("delocalized" if "pi" in t else "bonding")))
                    col = "#94a3b8" if cat == "core" else ("#fb923c" if cat == "antibonding" else "#38bdf8")
                    track_count = len([k for k in catalog.keys() if k.startswith("trk_") and k[4:].isdigit()]) + 1
                    tid = f"trk_{track_count:03d}"
                    catalog[tid] = {
                        "id": tid, "track_id": tid, "orbital_key": tid, "default_on": True, "name": lbl, "label": lbl, "role_transition": trans,
                        "category": cat, "orbital_class": "core" if cat == "core" else ("antibonding" if cat == "antibonding" else "bonding"),
                        "atoms": m.get("atoms", []), "color": col, "palette": m.get("palette"), "group_keys": g_keys, "frame_map": f_map
                    }

    megajson["global_electron_tracks"] = catalog
    return catalog


def get_or_build_catalog(
    megajson: Dict[str, Any],
    scf_cache: Optional[Dict[int, Dict[str, Any]]] = None,
    ensure_fn: Optional[Callable] = None,
    rebuild: bool = False
) -> Dict[str, Dict[str, Any]]:
    return build_global_track_catalog(megajson, scf_cache, ensure_fn, rebuild=rebuild)
