"""
Orbital Indexer and Multi-Step Pathway Tracker
Maps localized orbitals continuously across multi-step IRC reaction coordinates
and stitched intermediate junctions using atom correspondence maps and global track IDs.
"""

import re
from typing import Any, Callable, Dict, List, Optional, Tuple
from .orbital_selection import clean_role_label


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


def _match_cand(m_src: Optional[Dict[str, Any]], p2n: Dict[int, int], pool: Dict[str, Dict[str, Any]], is_star: bool) -> Optional[str]:
    if not m_src or not m_src.get("atoms"):
        return None
    at_src = _atom_nums(m_src.get("atoms", []))
    mapped = [p2n.get(a) for a in at_src]
    if any(x is None for x in mapped):
        return None
    s_mapped = set(mapped)
    if not s_mapped:
        return None

    best_k, best_score = None, -9999.0
    for k2, m2 in pool.items():
        if (m2.get("type") == "core") != (m_src.get("type") == "core"):
            continue
        if is_star != ("STAR" in k2 or "star" in m2.get("type", "")):
            continue
        s2 = set(_atom_nums(m2.get("atoms", [])))
        if not s2 or (s_mapped != s2 and not (s2.issubset(s_mapped) or s_mapped.issubset(s2))):
            continue
        sc = 100.0 if s_mapped == s2 else 40.0
        t1, t2 = m_src.get("type", ""), m2.get("type", "")
        sc += 30.0 if t1 == t2 else (10.0 if (t1 in ("pi", "lone_pair") and t2 in ("pi", "lone_pair")) else 0.0)

        p1, p2 = m_src.get("pop"), m2.get("pop")
        if p1 and p2:
            sc -= sum(abs(p1[a1] - p2[a2]) for a1, a2 in zip(at_src, mapped) if a1 < len(p1) and a2 < len(p2)) * 40.0
        e1, e2 = m_src.get("energy_ev"), m2.get("energy_ev")
        if e1 is not None and e2 is not None:
            sc -= abs(e1 - e2) * 1.5

        if sc > best_score:
            best_score, best_k = sc, k2
    return best_k if best_score > 0 else None


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
        if not g_keys[g] or g >= len(cands):
            break
        m = _meta_for(groups[g][1], g_keys[g])
        p2n = dict(zip(cands[g].get("prev_indices", []), cands[g].get("next_indices", [])))
        fd_next = _get_frame_dict(groups[g + 1][0], megajson, scf_cache, ensure_fn)
        if not fd_next:
            break
        pool = fd_next.get("antibonds_meta", {}) if ("STAR" in g_keys[g] or is_star) else fd_next.get("ibos_meta", {})
        g_keys[g + 1] = _match_cand(m, p2n, pool, is_star or ("STAR" in g_keys[g]))

    # Backward tracking across junctions
    for g in range(ref_g, 0, -1):
        if not g_keys[g] or (g - 1) >= len(cands):
            break
        m = _meta_for(groups[g][0], g_keys[g])
        n2p = dict(zip(cands[g - 1].get("next_indices", []), cands[g - 1].get("prev_indices", [])))
        fd_prev = _get_frame_dict(groups[g - 1][1], megajson, scf_cache, ensure_fn)
        if not fd_prev:
            break
        pool = fd_prev.get("antibonds_meta", {}) if ("STAR" in g_keys[g] or is_star) else fd_prev.get("ibos_meta", {})
        g_keys[g - 1] = _match_cand(m, n2p, pool, is_star or ("STAR" in g_keys[g]))

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
    if e_role and g_from > g_to:
        cands = [r["all_candidates"][0] for r in megajson.get("stitching_alignment_report", []) if r.get("all_candidates")]
        n2p_c: Dict[int, int] = {}
        for g in range(g_from, g_to, -1):
            if (g - 1) < len(cands):
                n2p = dict(zip(cands[g - 1].get("next_indices", []), cands[g - 1].get("prev_indices", [])))
                n2p_c = {k: n2p.get(v, v) for k, v in n2p_c.items()} if n2p_c else n2p
        e_role = re.sub(r'([A-Za-z]+)(\d+)', lambda m: f"{m.group(1)}{n2p_c.get(int(m.group(2)), int(m.group(2)))}", e_role)

    s_type, e_type = s_role.split('(')[0].strip() if s_role else "", e_role.split('(')[0].strip() if e_role else ""
    if s_role and e_role and s_type != e_type:
        label, trans = f"{s_role} → {e_role}", f"{s_type} → {e_type}"
    else:
        label = s_role or e_role or default_title
        trans = s_type or e_type or (label.split("(")[0] if "(" in label else label)
    return label, trans, s_role, e_role


def build_global_track_catalog(
    megajson: Dict[str, Any],
    scf_cache: Optional[Dict[int, Dict[str, Any]]] = None,
    ensure_fn: Optional[Callable] = None
) -> Dict[str, Dict[str, Any]]:
    """Construct global electron tracks uniquely identifying each orbital across all reaction steps."""
    frames = megajson.get("frames", [])
    n_frames = len(frames)
    if n_frames == 0:
        return {}

    groups = get_frame_groups(megajson)
    n_g = len(groups)
    catalog: Dict[str, Dict[str, Any]] = megajson.get("global_electron_tracks") or {}
    if not catalog:
        catalog = {
            "trk_esp": {"track_id": "trk_esp", "name": "ESP Surface", "label": "⚡ Electrostatic Surface & Charges (ESP)", "category": "esp", "orbital_class": "esp", "color": "#c084fc", "atoms": [], "group_keys": ["track_esp_charge"] * n_g, "frame_map": ["track_esp_charge"] * n_frames},
            "trk_homo": {"track_id": "trk_homo", "name": "HOMO", "label": "HOMO", "category": "frontier", "orbital_class": "bonding", "color": "#38bdf8", "atoms": [], "group_keys": ["HOMO"] * n_g, "frame_map": ["HOMO"] * n_frames},
            "trk_lumo": {"track_id": "trk_lumo", "name": "LUMO", "label": "LUMO", "category": "frontier", "orbital_class": "antibonding", "color": "#fb923c", "atoms": [], "group_keys": ["LUMO"] * n_g, "frame_map": ["LUMO"] * n_frames}
        }

    avail_frames = sorted(list((scf_cache or {}).keys()))
    for fi in avail_frames:
        fd = scf_cache[fi]
        pools = [("bonding", fd.get("ibos_meta", {})), ("pi", fd.get("canonical_pi_meta", {})), ("star", fd.get("antibonds_meta", {}))]
        for p_type, pool in pools:
            for k, m in pool.items():
                is_star = ("STAR" in k) or (p_type == "star")
                f_map, g_keys = build_track_frame_map(megajson, fi, k, scf_cache, None, is_star=is_star)

                matched_trk = None
                for tid, trk in catalog.items():
                    if tid in ("trk_esp", "trk_homo", "trk_lumo"): continue
                    trk_gk = trk.get("group_keys", [])
                    if any(gk is not None and len(trk_gk) > gi and trk_gk[gi] == gk for gi, gk in enumerate(g_keys)):
                        matched_trk = trk
                        break

                if matched_trk:
                    mgk = matched_trk.setdefault("group_keys", [None] * n_g)
                    mfm = matched_trk.setdefault("frame_map", [None] * n_frames)
                    for gi, gk in enumerate(g_keys):
                        if gk: mgk[gi] = gk
                    for f_i, fk in enumerate(f_map):
                        if fk: mfm[f_i] = fk
                    lbl, trans, _, _ = resolve_track_roles(megajson, mfm, matched_trk.get("label", m.get("label", k)), scf_cache, None)
                    matched_trk["label"] = matched_trk["name"] = lbl
                    matched_trk["role_transition"] = trans
                else:
                    lbl, trans, _, _ = resolve_track_roles(megajson, f_map, m.get("label", k), scf_cache, None)
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


def get_or_build_catalog(megajson: Dict[str, Any], scf_cache: Optional[Dict[int, Dict[str, Any]]] = None, ensure_fn: Optional[Callable] = None) -> Dict[str, Dict[str, Any]]:
    return build_global_track_catalog(megajson, scf_cache, ensure_fn)
