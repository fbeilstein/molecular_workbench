"""
Orbital Store Access & Persistence Interface — In-Memory Singleton & Fast Querying.
"""

import hashlib
import json
import os
import threading
from typing import Any, Dict, List, Optional, Tuple
import numpy as np
from pyscf import gto

from ..orbital_selection import get_system_charge, get_system_solvent
from .common import (
    FRONTIER_TRACK_MAP,
    NPZ_PATH,
    REGISTRY_PATH,
    SUBSCRIPTS,
)
from .matcher import _add_subscript
from .builder import _build_orbital_store

_MEMORY_STORE: Optional[Dict[str, Any]] = None
_STORE_LOCK = threading.Lock()


def compute_pathway_fingerprint(megajson: Dict[str, Any]) -> str:
    """SHA-256 hex digest (first 16 chars) capturing pathway geometry, charges, and solvent."""
    frames = megajson.get("frames", [])
    solvent = get_system_solvent(megajson) or "gas"
    data = []
    for i, f in enumerate(frames):
        syms = f.get("symbols") or megajson.get("symbols", [])
        lbls = f.get("atom_labels", [])
        coords = np.round(np.array(f.get("coords", f.get("positions", []))), 3).tolist()
        chg = get_system_charge(megajson, i)
        data.append((syms, lbls, coords, chg))
    raw = json.dumps({"solvent": solvent, "frames": data}, sort_keys=True)
    return hashlib.sha256(raw.encode("utf-8")).hexdigest()[:16]


def clear_memory_cache():
    """Clear only the in-memory singleton cache."""
    global _MEMORY_STORE
    with _STORE_LOCK:
        _MEMORY_STORE = None


def clear_orbital_store():
    """Clear in-memory and disk caches."""
    global _MEMORY_STORE
    with _STORE_LOCK:
        _MEMORY_STORE = None
        for p in (NPZ_PATH, REGISTRY_PATH):
            if os.path.exists(p):
                try:
                    os.remove(p)
                except Exception:
                    pass


def get_or_build_orbital_store(megajson: Dict[str, Any], force_rebuild: bool = False) -> Dict[str, Any]:
    """Retrieve in-memory or on-disk orbital store, building if missing or stale."""
    global _MEMORY_STORE
    p_hash = compute_pathway_fingerprint(megajson)
    with _STORE_LOCK:
        if not force_rebuild and _MEMORY_STORE is not None:
            if _MEMORY_STORE.get("pathway_hash") == p_hash:
                megajson["global_electron_tracks"] = _MEMORY_STORE["registry"]["tracks"]
                return _MEMORY_STORE

        if not force_rebuild and os.path.exists(REGISTRY_PATH) and os.path.exists(NPZ_PATH):
            try:
                with open(REGISTRY_PATH, "r") as f:
                    reg = json.load(f)
                if reg.get("pathway_hash") == p_hash and reg.get("total_frames") == len(megajson.get("frames", [])):
                    npz = np.load(NPZ_PATH)
                    store_dict = {
                        "pathway_hash": p_hash,
                        "total_frames": reg["total_frames"],
                        "registry": reg,
                        "npz": npz,
                    }
                    _MEMORY_STORE = store_dict
                    megajson["global_electron_tracks"] = reg["tracks"]
                    return store_dict
            except Exception:
                pass

        _MEMORY_STORE = _build_orbital_store(megajson, p_hash)
        return _MEMORY_STORE


def get_frame_pool_payload(megajson: Dict[str, Any], frame_idx: int) -> Dict[str, Any]:
    """Return instant get_frame_orbitals dictionary for the given frame index."""
    store = get_or_build_orbital_store(megajson)
    frames = megajson.get("frames", [])
    if frame_idx < 0 or frame_idx >= len(frames):
        return {"success": False, "error": f"Frame {frame_idx} out of bounds"}

    fr = frames[frame_idx]
    npz = store["npz"]
    registry = store["registry"]
    tracks = registry["tracks"]

    tids_f = npz[f"f{frame_idx}_tids"]
    pops_f = npz[f"f{frame_idx}_pops"]
    energies_f = npz[f"f{frame_idx}_energies"]
    charges_f = npz[f"f{frame_idx}_charges"]

    frontier_orbitals = {}
    all_ibos = []
    antibonds = []
    canonical_pi = []

    raw_entries = []
    frame_lbl_groups: Dict[str, List[str]] = {}
    for k, tid_raw in enumerate(tids_f):
        tid_s = str(tid_raw)
        trk = tracks.get(tid_s)
        if not trk:
            continue
        f_labels = trk.get("frame_labels", [])
        f_types = trk.get("frame_types", [])
        f_cats = trk.get("frame_categories", [])
        f_atoms_arr = trk.get("frame_atoms", [])

        raw_f_lbl = f_labels[frame_idx] if (frame_idx < len(f_labels) and f_labels[frame_idx]) else trk["label"]
        f_type = f_types[frame_idx] if (frame_idx < len(f_types) and f_types[frame_idx]) else trk["type"]
        f_cat = f_cats[frame_idx] if (frame_idx < len(f_cats) and f_cats[frame_idx]) else trk["category"]
        f_atoms = f_atoms_arr[frame_idx] if (frame_idx < len(f_atoms_arr) and f_atoms_arr[frame_idx] is not None) else trk["atoms"]

        grp_key = f"{trk['orbital_class']}::{f_type}::{raw_f_lbl}"
        frame_lbl_groups.setdefault(grp_key, []).append(tid_s)
        raw_entries.append((k, tid_s, trk, raw_f_lbl, f_type, f_cat, f_atoms, grp_key))

    disambig_map: Dict[str, str] = {}
    for grp_key, grp_tids in frame_lbl_groups.items():
        if len(grp_tids) >= 2:
            grp_tids_sorted = sorted(grp_tids, key=lambda x: (tracks[x].get("start_frame", 0), x))
            base_lbl = grp_key.split("::", 2)[2]
            for idx, tid_s in enumerate(grp_tids_sorted):
                disambig_map[tid_s] = _add_subscript(base_lbl, SUBSCRIPTS[idx % len(SUBSCRIPTS)])

    for k, tid_s, trk, raw_f_lbl, f_type, f_cat, f_atoms, _ in raw_entries:
        p_list = pops_f[k].tolist() if k < len(pops_f) else []
        e_val = float(energies_f[k]) if k < len(energies_f) else 0.0
        disp_f_lbl = disambig_map.get(tid_s, raw_f_lbl)

        item = {
            "track_id": tid_s,
            "id": tid_s,
            "key": tid_s,
            "local_key": tid_s,
            "orbital_key": tid_s,
            "name": tid_s,
            "index": k,
            "type": f_type,
            "category": f_cat,
            "orbital_class": trk["orbital_class"],
            "label": disp_f_lbl,
            "frame_label": disp_f_lbl,
            "track_label": trk["label"],
            "decoupled_label": trk["label"],
            "atoms": f_atoms,
            "track_atoms": trk["atoms"],
            "energy_ev": round(e_val, 2),
            "pop": p_list,
            "palette": (
                {"pos": "#38bdf8", "neg": "#1d4ed8"} if trk["orbital_class"] == "bonding"
                else ({"pos": "#fb923c", "neg": "#db2777"} if trk["orbital_class"] == "antibonding"
                      else {"pos": "#94a3b8", "neg": "#475569"})
            ),
        }

        if f_type == "frontier":
            f_key = next((fk for fk, mapped_tid in FRONTIER_TRACK_MAP.items() if mapped_tid == tid_s), tid_s)
            frontier_orbitals[f_key] = item
        elif f_type in ("delocalized_pi", "delocalized_pistar"):
            canonical_pi.append(item)
        elif f_type in ("sigma_star", "pi_star") or trk["orbital_class"] == "antibonding":
            antibonds.append(item)
        elif f_type != "core":
            all_ibos.append(item)

    return {
        "success": True,
        "frame_index": frame_idx,
        "charges": charges_f.tolist(),
        "frontier_orbitals": frontier_orbitals,
        "all_ibos": all_ibos,
        "antibonds": antibonds,
        "canonical_pi": canonical_pi,
        "coords": fr.get("coords", fr.get("positions", [])),
        "symbols": fr.get("symbols") or megajson.get("symbols", []),
        "atom_labels": fr.get("atom_labels"),
        "bonds": fr.get("bonds", []),
        "step_index": fr.get("step_index"),
    }


def get_track_vector_and_mol(
    megajson: Dict[str, Any], frame_idx: int, track_id: str
) -> Tuple[Optional[Any], Optional[np.ndarray], Optional[Dict[str, Any]]]:
    """Retrieve lightweight PySCF Mol, orbital vector, and metadata without recomputing SCF."""
    store = get_or_build_orbital_store(megajson)
    npz = store["npz"]
    key_tids = f"f{frame_idx}_tids"
    if key_tids not in npz:
        return None, None, None

    tids_f = npz[key_tids]
    matches = np.where(tids_f == track_id)[0]
    if len(matches) == 0:
        return None, None, None

    col_idx = matches[0]
    vec = npz[f"f{frame_idx}_vecs"][:, col_idx].astype(np.float64)

    fr = megajson["frames"][frame_idx]
    syms = fr.get("symbols") or megajson.get("symbols", [])
    coords = fr.get("coords", fr.get("positions", []))
    charge = get_system_charge(megajson, frame_idx)

    atom_str = "; ".join(f"{s} {p[0]} {p[1]} {p[2]}" for s, p in zip(syms, coords))
    mol = gto.M(atom=atom_str, basis="sto-3g", charge=charge, spin=0, verbose=0)
    meta = store.get("registry", {}).get("tracks", {}).get(track_id, {})

    return mol, vec, meta

