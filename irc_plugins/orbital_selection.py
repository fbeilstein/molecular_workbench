"""
Orbital Selection and Nomenclature Module
Resolves orbital keys, enriches timeline tracks, and applies chemical taxonomy and color palettes.
"""

import re
from typing import Any, Dict, List, Optional
from rdkit import Chem


def get_system_charge(megajson: Dict[str, Any], fi: Optional[int] = None) -> int:
    """Infer molecular system net formal charge for a specific frame or overall."""
    if fi is not None and "frames" in megajson and 0 <= fi < len(megajson["frames"]):
        st_idx = megajson["frames"][fi].get("step_index")
        target_step = st_idx if isinstance(st_idx, int) else (int(st_idx.split("➔")[1]) if isinstance(st_idx, str) and "➔" in st_idx else None)
        if target_step is not None:
            for b in megajson.get("step_boundaries", []):
                if b.get("step_index") == target_step:
                    smi = (b.get("step_name") or b.get("reaction_smiles", "")).split(">>")[0]
                    try:
                        m = Chem.MolFromSmiles(smi)
                        if m:
                            return sum(a.GetFormalCharge() for a in m.GetAtoms())
                    except Exception:
                        pass

    for ck in ("charge", "net_charge", "formal_charge"):
        if ck in megajson:
            try:
                return int(megajson[ck])
            except (ValueError, TypeError):
                pass
    boundaries = megajson.get("step_boundaries", [])
    if boundaries and ">>" in (boundaries[0].get("reaction_smiles") or boundaries[0].get("step_name", "")):
        smi = (boundaries[0].get("reaction_smiles") or boundaries[0].get("step_name", "")).split(">>")[0]
        try:
            m = Chem.MolFromSmiles(smi)
            return sum(a.GetFormalCharge() for a in m.GetAtoms()) if m else 0
        except Exception:
            return 0
    return 0



def resolve_orbital_key(d: Dict[str, Any], orb_name: str, chem_key: Optional[str] = None) -> Optional[str]:
    """Resolve an orbital or chem_key against available molecular orbital metadata."""
    if not d or "vecs" not in d:
        return None
    all_metas = {
        **d.get("ibos_meta", {}),
        **d.get("canonical_pi_meta", {}),
        **d.get("antibonds_meta", {}),
        **d.get("meta", {})
    }
    if orb_name in d.get("meta", {}):
        return orb_name
    if chem_key and chem_key not in ("HOMO", "LUMO", "esp_charge", orb_name):
        for k, meta in all_metas.items():
            if (meta.get("chem_key") == chem_key or meta.get("label") == chem_key) and k in d.get("vecs", {}):
                return k
        return None
    if orb_name in d.get("vecs", {}):
        return orb_name
    return None


def clean_role_label(meta: Optional[Dict[str, Any]]) -> Optional[str]:
    """Produce clean chemical shorthand for orbital role (LP, σ, π, σ*, π*)."""
    if not meta:
        return None
    typ = meta.get("type", "")
    atoms = meta.get("atoms", [])
    atoms_str = "-".join(atoms) if atoms else ""
    if typ == "lone_pair":
        return f"LP({atoms[0]})" if atoms else "LP"
    elif "sigma_star" in typ or "sigmastar" in typ:
        return f"σ*({atoms_str})" if atoms_str else "σ*"
    elif "pistar" in typ or "pi_star" in typ or typ == "delocalized_pistar":
        return f"π*({atoms_str})" if atoms_str else "π*"
    elif "pi" in typ or typ == "delocalized_pi":
        return f"π({atoms_str})" if atoms_str else "π"
    elif "sigma" in typ or typ == "delocalized_sigma":
        return f"σ({atoms_str})" if atoms_str else "σ"
    return meta.get("label") or typ


def tag_orbital_classes(f_meta: Dict[str, Any], ibos_m: Dict[str, Any], star_m: Dict[str, Any], c_meta: Dict[str, Any]):
    """Tag every orbital with its chemical nature and distinct color palette."""
    for k, m in f_meta.items():
        name = m.get("name", "").upper()
        if "HOMO" in name:
            m["orbital_class"] = "bonding"
            m["palette"] = {"pos": "#38bdf8", "neg": "#1d4ed8"}
        elif "LUMO" in name:
            m["orbital_class"] = "antibonding"
            m["palette"] = {"pos": "#fb923c", "neg": "#db2777"}
        else:
            m["orbital_class"] = "bonding"
            m["palette"] = {"pos": "#38bdf8", "neg": "#1d4ed8"}

    for k, m in ibos_m.items():
        t = m.get("type", "")
        if t == "core":
            m["orbital_class"] = "core"
            m["palette"] = {"pos": "#94a3b8", "neg": "#475569"}
        else:
            m["orbital_class"] = "bonding"
            m["palette"] = {"pos": "#38bdf8", "neg": "#1d4ed8"}

    for k, m in star_m.items():
        m["orbital_class"] = "antibonding"
        m["palette"] = {"pos": "#fb923c", "neg": "#db2777"}

    for k, m in c_meta.items():
        m["orbital_class"] = "bonding"
        m["palette"] = {"pos": "#38bdf8", "neg": "#1d4ed8"}


def enrich_timeline_tracks(megajson: Dict[str, Any], tracks: list, scf_cache: Optional[Dict[int, Dict[str, Any]]] = None) -> list:
    """Enrich orbital timeline tracks with normal chemical nomenclature and role transitions."""
    if not tracks:
        return tracks

    cache = scf_cache or {}
    n_frames = len(megajson.get("frames", []))
    f0 = cache.get(0) if n_frames > 0 else None
    flast = cache.get(n_frames - 1) if n_frames > 0 else None

    for trk in tracks:
        if not isinstance(trk, dict):
            continue
        ok = trk.get("orbital_key", "")
        if ok == "track_esp_charge" or trk.get("id") == "trk_esp_charge":
            trk["label"] = "ESP Surface"
            trk["title"] = "ESP Surface"
            trk["short_label"] = "ESP"
            trk["role_transition"] = "ESP"
            trk["category"] = "ESP"
            continue

        start_role = trk.get("start_role")
        end_role = trk.get("end_role")

        if (not start_role or not end_role) and f0 is not None and flast is not None:
            r0 = resolve_orbital_key(f0, ok, trk.get("chem_key"))
            m0 = {**f0.get("ibos_meta", {}), **f0.get("canonical_pi_meta", {}), **f0.get("antibonds_meta", {})}.get(r0, {}) if r0 else {}
            rlast = resolve_orbital_key(flast, ok, trk.get("chem_key"))
            mlast = {**flast.get("ibos_meta", {}), **flast.get("canonical_pi_meta", {}), **flast.get("antibonds_meta", {})}.get(rlast, {}) if rlast else {}
            start_role = clean_role_label(m0)
            end_role = clean_role_label(mlast)

        ck = trk.get("chem_key", "")
        if not start_role and ck:
            if ck.startswith("lp_"):
                start_role = f"LP({ck.replace('lp_', '')})"
            elif ck.startswith("sigma_"):
                start_role = f"σ({ck.replace('sigma_', '')})"
            elif ck.startswith("pi_"):
                start_role = f"π({ck.replace('pi_', '')})"
            elif "pistar" in ck:
                start_role = f"π*({ck.split('_')[1] if '_' in ck else ''})"

        if start_role and end_role and start_role != end_role:
            label = f"{start_role} → {end_role}"
            role_trans = f"{start_role.split('(')[0]} → {end_role.split('(')[0]}"
        else:
            label = start_role or end_role or trk.get("title") or ok
            label = re.sub(r"^Non-Bonding\s+", "", label)
            label = re.sub(r"\[Antibond\]", "", label).strip()
            role_trans = label.split("(")[0] if "(" in label else label

        trk["label"] = label
        trk["title"] = label
        trk["role_transition"] = role_trans
        trk["short_label"] = role_trans
        if start_role:
            trk["start_role"] = start_role
        if end_role:
            trk["end_role"] = end_role
        trk["category"] = role_trans

    return tracks
