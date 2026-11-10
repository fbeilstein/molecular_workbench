"""
Orbital Selection and Nomenclature Module
Resolves orbital keys, enriches timeline tracks, and applies chemical taxonomy and color palettes.
"""

import re
from typing import Any, Dict, List, Optional
from rdkit import Chem

SOLVENT_DIELECTRICS: Dict[str, float] = {
    "water": 80.4,
    "h2o": 80.4,
    "acetonitrile": 37.5,
    "mecn": 37.5,
    "ch3cn": 37.5,
    "dmso": 46.7,
    "dmf": 37.2,
    "methanol": 32.7,
    "meoh": 32.7,
    "ethanol": 24.5,
    "etoh": 24.5,
    "acetone": 20.7,
    "ch2cl2": 8.93,
    "dcm": 8.93,
    "dichloromethane": 8.93,
    "thf": 7.58,
    "chloroform": 4.81,
    "chcl3": 4.81,
    "ether": 4.33,
    "diethylether": 4.33,
    "toluene": 2.38,
    "benzene": 2.27,
    "hexane": 1.88,
}


def get_system_solvent(megajson: Dict[str, Any]) -> Optional[str]:
    """Return normalized implicit solvent key from megajson/project, or None for gas/vacuum."""
    raw = megajson.get("solvent")
    if raw is None:
        try:
            from .project_manager import init_project
            raw = init_project().get("solvent", "water")
        except Exception:
            raw = "water"
    s = str(raw).strip().lower()
    if not s or s in ("gas", "none", "vacuum", "gas_phase"):
        return None
    return s


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
    if orb_name in d.get("meta", {}):
        return orb_name
    all_metas = {
        **d.get("ibos_meta", {}),
        **d.get("canonical_pi_meta", {}),
        **d.get("antibonds_meta", {}),
        **d.get("meta", {})
    }
    syms = set(d.get("symbols", []))
    if chem_key:
        target_syms = set(re.findall(r'[A-Z][a-z]?', chem_key))
        if target_syms and not target_syms.issubset(syms):
            return None

    if orb_name in d.get("vecs", {}):
        return orb_name

    if chem_key and chem_key not in ("HOMO", "LUMO", "esp_charge", orb_name):
        for k, meta in all_metas.items():
            if (meta.get("chem_key") == chem_key or meta.get("label") == chem_key) and k in d.get("vecs", {}):
                return k
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
        t = str(m.get("type", "")).lower()
        lbl = str(m.get("label", ""))
        if "star" in t or "π*" in lbl:
            m["orbital_class"] = "antibonding"
            m["palette"] = {"pos": "#fb923c", "neg": "#db2777"}
        else:
            m["orbital_class"] = "bonding"
            m["palette"] = {"pos": "#38bdf8", "neg": "#1d4ed8"}


def enrich_timeline_tracks(
    megajson: Dict[str, Any],
    tracks: list,
    _scf_cache: Optional[Dict[int, Dict[str, Any]]] = None,
    _ensure_fn: Optional[Any] = None,
) -> list:
    """Enrich orbital timeline tracks with exact track metadata, role transitions, and frame maps."""
    if not tracks:
        return tracks
    catalog = megajson.get("global_electron_tracks")
    if not catalog or len(catalog) <= 3:
        try:
            from .orbital_tracker import get_existing_orbital_store
            store = get_existing_orbital_store(megajson)
            catalog = store.get("registry", {}).get("tracks", {}) if store else {}
        except Exception:
            catalog = {}

    n_frames = len(megajson.get("frames", []))
    for trk in tracks:
        if not isinstance(trk, dict):
            continue
        tid = trk.get("track_id") or trk.get("id") or trk.get("orbital_key", "")
        if tid in ("track_esp_charge", "trk_esp") or trk.get("id") in ("trk_esp_charge", "trk_esp"):
            trk["id"] = "trk_esp"
            trk["track_id"] = "trk_esp"
            trk["orbital_key"] = "track_esp_charge"
            trk["label"] = trk["title"] = "ESP Surface"
            trk["short_label"] = trk["role_transition"] = trk["category"] = "ESP"
            trk["orbital_class"] = "esp"
            trk.setdefault("default_on", True)
            continue

        trk["id"] = tid
        trk["track_id"] = tid
        trk["orbital_key"] = tid
        trk.setdefault("default_on", True)

        if tid in catalog:
            c = catalog[tid]
            trk["frame_map"] = c["frame_map"]
            trk["label"] = trk["title"] = c["label"]
            trk["role_transition"] = trk["short_label"] = c.get("role_transition", "")
            trk["category"] = c["category"]
            trk["orbital_class"] = c["orbital_class"]
            trk["type"] = c.get("type", trk.get("type", ""))
            trk["atoms"] = c.get("atoms", [])
            continue

        # Migration lookup for legacy tracks (e.g. old IBO_xx or old track id by chemical match)
        migrated = None
        trk_atoms = set(trk.get("atoms", []))
        trk_type = trk.get("type", "")
        for c_tid, c in catalog.items():
            if trk_atoms and set(c.get("atoms", [])) == trk_atoms:
                if not trk_type or c.get("type") == trk_type or (("pi" in trk_type) == ("pi" in c.get("type", ""))):
                    migrated = c
                    break

        if migrated:
            m_tid = migrated["track_id"]
            trk["id"] = trk["track_id"] = trk["orbital_key"] = m_tid
            trk["frame_map"] = migrated["frame_map"]
            trk["label"] = trk["title"] = migrated["label"]
            trk["role_transition"] = trk["short_label"] = migrated.get("role_transition", "")
            trk["category"] = migrated["category"]
            trk["orbital_class"] = migrated["orbital_class"]
            trk["atoms"] = migrated.get("atoms", [])
        else:
            trk["frame_map"] = [tid if trk.get("ref_frame") == f else None for f in range(n_frames)]

    return tracks
