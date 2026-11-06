"""
Orbital Tracker and Persistent Electron Trajectory Store
Provides part-agnostic continuous wavefunction tracking, parallel frame SCF+IBO,
subspace-specific Hungarian assignment, and single-source .npz/.json persistence.
Knows zero about reaction steps, boundaries, or group keys.
"""

from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
import os
import sys
import threading
import time
from typing import Any, Dict, List, Optional, Set, Tuple

import numpy as np
from pyscf import gto, scf, lo
from scipy.optimize import linear_sum_assignment

from .ts_engine.orbital_engine import compute_frame_scf
from .ts_engine.orbital_ibo import compute_frame_ibos
from .ts_engine.orbital_frontier import extract_frontier_orbitals
from .orbital_selection import get_system_charge, tag_orbital_classes

ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
TRAJECTORIES_DIR = os.path.join(ROOT_DIR, "trajectories")
NPZ_PATH = os.path.join(TRAJECTORIES_DIR, "orbital_tracks.npz")
REGISTRY_PATH = os.path.join(TRAJECTORIES_DIR, "orbital_registry.json")

_MEMORY_STORE: Optional[Dict[str, Any]] = None
_STORE_LOCK = threading.Lock()

FRONTIER_KEYS = ["HOMO-3", "HOMO-2", "HOMO-1", "HOMO", "LUMO", "LUMO+1", "LUMO+2"]
FRONTIER_TRACK_MAP = {
    "HOMO-3": "trk_homo_3",
    "HOMO-2": "trk_homo_2",
    "HOMO-1": "trk_homo_1",
    "HOMO": "trk_homo",
    "LUMO": "trk_lumo",
    "LUMO+1": "trk_lumo_1",
    "LUMO+2": "trk_lumo_2",
}

SUBSCRIPTS = ["₁", "₂", "₃", "₄", "₅"]


def compute_pathway_fingerprint(megajson: Dict[str, Any]) -> str:
    """SHA-256 hex digest (first 16 chars) capturing pathway geometry and identity."""
    frames = megajson.get("frames", [])
    data = []
    for i, f in enumerate(frames):
        syms = f.get("symbols") or megajson.get("symbols", [])
        lbls = f.get("atom_labels", [])
        coords = np.round(np.array(f.get("coords", f.get("positions", []))), 3).tolist()
        chg = get_system_charge(megajson, i)
        data.append((syms, lbls, coords, chg))
    raw = json.dumps(data, sort_keys=True)
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


def _add_subscript(label: str, sub: str) -> str:
    for prefix in ("π*(", "π(", "σ*(", "σ(", "LP("):
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


def _worker_compute_frame(args: Tuple[int, List[str], List[List[float]], Optional[List[str]], int]) -> Tuple[int, List[float], Dict[str, Any], Dict[str, Any], Dict[str, np.ndarray], np.ndarray, Dict[str, Any], Dict[str, Any], Dict[str, np.ndarray]]:
    """Worker function for parallel SCF + IBO evaluation."""
    t, symbols, coords, atom_labels, charge = args
    mol, mf, mull = compute_frame_scf(symbols, coords, charge=charge, spin=0, basis="sto-3g")
    ibos_dict, star_dict, vecs_dict, pops, ibo_coeff = compute_frame_ibos(
        mol, mf, symbols, atom_labels=atom_labels
    )
    f_meta, c_meta, f_vecs = extract_frontier_orbitals(
        mol, mf, None, None, ibos_dict, star_dict, vecs_dict, atom_labels=atom_labels
    )
    tag_orbital_classes(f_meta, ibos_dict, star_dict, c_meta)
    return (t, mull, ibos_dict, star_dict, vecs_dict, pops, f_meta, c_meta, f_vecs)


def _build_orbital_store(megajson: Dict[str, Any], p_hash: str) -> Dict[str, Any]:
    """Run parallel frame SCF and sequential tracking to build persistent store."""
    frames = megajson.get("frames", [])
    N = len(frames)
    if N == 0:
        return {"pathway_hash": p_hash, "total_frames": 0, "registry": {"tracks": {}}}

    # Step 1: Prepare parallel evaluation tasks
    tasks = []
    for t in range(N):
        fr = frames[t]
        symbols = fr.get("symbols") or megajson.get("symbols", [])
        coords = fr.get("coords", fr.get("positions", []))
        atom_labels = fr.get("atom_labels")
        if not atom_labels or len(atom_labels) != len(symbols):
            atom_labels = [f"{symbols[i]}{i}" for i in range(len(symbols))]
        charge = get_system_charge(megajson, t)
        tasks.append((t, symbols, coords, atom_labels, charge))

    workers = min(8, max(1, os.cpu_count() or 4))
    with ProcessPoolExecutor(max_workers=workers) as ex:
        precomputed_list = list(ex.map(_worker_compute_frame, tasks))

    # Sort results by frame index
    precomputed_list.sort(key=lambda x: x[0])

    track_records: Dict[str, Dict[str, Any]] = {}
    frame_outputs: List[Dict[str, Any]] = []

    prev_mol: Optional[Any] = None
    active_core: List[Dict[str, Any]] = []
    active_valence: List[Dict[str, Any]] = []
    active_antibond: List[Dict[str, Any]] = []
    active_deloc_pi: List[Dict[str, Any]] = []
    prev_frontier_vecs: Dict[str, np.ndarray] = {}

    next_seq = 1
    next_core_seq = 0

    # Step 2: Sequential cross-frame propagation across subspaces
    for t in range(N):
        fr = frames[t]
        symbols = fr.get("symbols") or megajson.get("symbols", [])
        coords = fr.get("coords", fr.get("positions", []))
        charge = get_system_charge(megajson, t)

        atom_str = "; ".join(f"{s} {p[0]} {p[1]} {p[2]}" for s, p in zip(symbols, coords))
        mol = gto.M(atom=atom_str, basis="sto-3g", charge=charge, spin=0, verbose=0)
        ovlp = mol.intor_symmetric("int1e_ovlp")
        ao_atom = [int(l.split()[0]) for l in mol.ao_labels()]

        _, mull, ibos_dict, star_dict, vecs_dict, pops, f_meta, c_meta, f_vecs = precomputed_list[t]

        def _calc_pop(v: np.ndarray, is_virt: bool = False) -> np.ndarray:
            p_mat = np.outer(v, v) @ ovlp
            ap = np.zeros(mol.natm, dtype=np.float32)
            for mu in range(mol.nao):
                ap[ao_atom[mu]] += p_mat[mu, mu]
            ap *= (1.0 if is_virt else 2.0)
            return ap

        s_cross = None
        if t > 0 and prev_mol is not None:
            s_cross = gto.intor_cross("int1e_ovlp", prev_mol, mol)

        frame_tids: List[str] = []
        frame_vecs: List[np.ndarray] = []
        frame_pops: List[np.ndarray] = []
        frame_energies: List[float] = []

        def _record_track_instance(
            tid: str,
            vec: np.ndarray,
            pop_arr: np.ndarray,
            energy_ev: float,
            role_lbl: str,
            typ: str,
            cat: str,
            orb_class: str,
            atoms: List[str],
            subspace: str,
        ):
            frame_tids.append(tid)
            frame_vecs.append(vec)
            frame_pops.append(pop_arr)
            frame_energies.append(energy_ev)

            if tid not in track_records:
                track_records[tid] = {
                    "track_id": tid,
                    "subspace": subspace,
                    "category": cat,
                    "orbital_class": orb_class,
                    "type": typ,
                    "atoms_set": set(atoms),
                    "frame_present": [False] * N,
                    "frame_roles": [None] * N,
                    "frame_energies": [None] * N,
                    "frame_pops": [None] * N,
                }
            rec = track_records[tid]
            rec["frame_present"][t] = True
            rec["frame_roles"][t] = role_lbl
            rec["frame_energies"][t] = round(float(energy_ev), 2)
            rec["frame_pops"][t] = pop_arr.round(2).tolist()
            rec["atoms_set"].update(atoms)

        # ----------------------------------------------------
        # Subspace A: Frontier Orbitals
        # ----------------------------------------------------
        curr_frontier_vecs: Dict[str, np.ndarray] = {}
        for fk in FRONTIER_KEYS:
            if fk in f_meta and fk in f_vecs:
                tid = FRONTIER_TRACK_MAP[fk]
                v = f_vecs[fk].copy()
                if s_cross is not None and tid in prev_frontier_vecs:
                    if float(prev_frontier_vecs[tid] @ s_cross @ v) < 0.0:
                        v = -v
                curr_frontier_vecs[tid] = v
                meta = f_meta[fk]
                is_virt = "LUMO" in fk
                p_arr = _calc_pop(v, is_virt=is_virt)
                e_ev = float(meta.get("energy_ev", 0.0))
                _record_track_instance(
                    tid, v, p_arr, e_ev, fk, "frontier", meta.get("category", "frontier"),
                    meta.get("orbital_class", "bonding"), meta.get("atoms", []), "frontier"
                )
        prev_frontier_vecs = curr_frontier_vecs

        # ----------------------------------------------------
        # Subspace B: Core Orbitals
        # ----------------------------------------------------
        core_keys = [k for k, v in ibos_dict.items() if v.get("type") == "core"]
        core_keys.sort(key=lambda k: ibos_dict[k]["energy_ev"])
        curr_core_list = []
        for k in core_keys:
            v = vecs_dict[k].copy()
            meta = ibos_dict[k]
            p_arr = pops[meta["index"]] if meta["index"] < len(pops) else _calc_pop(v)
            curr_core_list.append({"vec": v, "pop": p_arr, "meta": meta})

        next_active_core: List[Dict[str, Any]] = []
        if t == 0 or s_cross is None or len(active_core) == 0:
            for item in curr_core_list:
                tid = f"trk_core_{next_core_seq:02d}"
                next_core_seq += 1
                meta = item["meta"]
                _record_track_instance(
                    tid, item["vec"], item["pop"], meta["energy_ev"], meta["label"],
                    "core", "core", "core", meta.get("atoms", []), "core"
                )
                next_active_core.append({"tid": tid, "vec": item["vec"]})
        else:
            prev_mat = np.column_stack([ac["vec"] for ac in active_core])
            curr_mat = np.column_stack([item["vec"] for item in curr_core_list])
            Omega = prev_mat.T @ s_cross @ curr_mat
            r_ind, c_ind = linear_sum_assignment(-np.abs(Omega))
            matched_curr = set()
            for r, c in zip(r_ind, c_ind):
                if abs(Omega[r, c]) >= 0.60:
                    matched_curr.add(c)
                    tid = active_core[r]["tid"]
                    v = curr_mat[:, c].copy()
                    if Omega[r, c] < 0.0:
                        v = -v
                    meta = curr_core_list[c]["meta"]
                    _record_track_instance(
                        tid, v, curr_core_list[c]["pop"], meta["energy_ev"], meta["label"],
                        "core", "core", "core", meta.get("atoms", []), "core"
                    )
                    next_active_core.append({"tid": tid, "vec": v})
            for c in range(len(curr_core_list)):
                if c not in matched_curr:
                    tid = f"trk_core_{next_core_seq:02d}"
                    next_core_seq += 1
                    meta = curr_core_list[c]["meta"]
                    _record_track_instance(
                        tid, curr_core_list[c]["vec"], curr_core_list[c]["pop"], meta["energy_ev"],
                        meta["label"], "core", "core", "core", meta.get("atoms", []), "core"
                    )
                    next_active_core.append({"tid": tid, "vec": curr_core_list[c]["vec"]})
        active_core = next_active_core

        # ----------------------------------------------------
        # Subspace C: Localized Valence IBOs
        # ----------------------------------------------------
        val_keys = [k for k, v in ibos_dict.items() if v.get("type") != "core"]
        curr_val_list = []
        for k in val_keys:
            v = vecs_dict[k].copy()
            meta = ibos_dict[k]
            p_arr = pops[meta["index"]] if meta["index"] < len(pops) else _calc_pop(v)
            curr_val_list.append({"vec": v, "pop": p_arr, "meta": meta})

        next_active_valence: List[Dict[str, Any]] = []
        if t == 0 or s_cross is None or len(active_valence) == 0:
            for item in curr_val_list:
                tid = f"trk_{next_seq:04d}"
                next_seq += 1
                meta = item["meta"]
                _record_track_instance(
                    tid, item["vec"], item["pop"], meta["energy_ev"], meta["label"],
                    meta["type"], meta["category"], meta["orbital_class"], meta.get("atoms", []), "valence_ibo"
                )
                next_active_valence.append({"tid": tid, "vec": item["vec"], "atoms": meta.get("atoms", []), "type": meta["type"]})
        else:
            prev_mat = np.column_stack([av["vec"] for av in active_valence])
            curr_mat = np.column_stack([item["vec"] for item in curr_val_list])

            # Degenerate triple-bond pi gauge pre-alignment
            Q = len(curr_val_list)
            for i1 in range(Q):
                for i2 in range(i1 + 1, Q):
                    m1, m2 = curr_val_list[i1]["meta"], curr_val_list[i2]["meta"]
                    if "pi" in m1.get("type", "") and "pi" in m2.get("type", ""):
                        ats1 = tuple(sorted(m1.get("atoms", [])))
                        ats2 = tuple(sorted(m2.get("atoms", [])))
                        if len(ats1) == 2 and ats1 == ats2:
                            prev_matches = [
                                pi_idx for pi_idx, av in enumerate(active_valence)
                                if "pi" in av.get("type", "") and tuple(sorted(av.get("atoms", []))) == ats1
                            ]
                            if len(prev_matches) == 2:
                                p_pair = prev_mat[:, prev_matches]
                                c_pair = curr_mat[:, [i1, i2]]
                                c_pair_rot = _align_degenerate_2d_pair(p_pair, c_pair, s_cross)
                                curr_mat[:, [i1, i2]] = c_pair_rot
                                curr_val_list[i1]["vec"] = c_pair_rot[:, 0]
                                curr_val_list[i2]["vec"] = c_pair_rot[:, 1]

            Omega = prev_mat.T @ s_cross @ curr_mat
            r_ind, c_ind = linear_sum_assignment(-np.abs(Omega))
            matched_curr = set()
            for r, c in zip(r_ind, c_ind):
                if abs(Omega[r, c]) >= 0.45:
                    matched_curr.add(c)
                    tid = active_valence[r]["tid"]
                    v = curr_mat[:, c].copy()
                    if Omega[r, c] < 0.0:
                        v = -v
                    meta = curr_val_list[c]["meta"]
                    _record_track_instance(
                        tid, v, curr_val_list[c]["pop"], meta["energy_ev"], meta["label"],
                        meta["type"], meta["category"], meta["orbital_class"], meta.get("atoms", []), "valence_ibo"
                    )
                    next_active_valence.append({"tid": tid, "vec": v, "atoms": meta.get("atoms", []), "type": meta["type"]})
            for c in range(len(curr_val_list)):
                if c not in matched_curr:
                    tid = f"trk_{next_seq:04d}"
                    next_seq += 1
                    meta = curr_val_list[c]["meta"]
                    _record_track_instance(
                        tid, curr_val_list[c]["vec"], curr_val_list[c]["pop"], meta["energy_ev"],
                        meta["label"], meta["type"], meta["category"], meta["orbital_class"],
                        meta.get("atoms", []), "valence_ibo"
                    )
                    next_active_valence.append({"tid": tid, "vec": curr_val_list[c]["vec"], "atoms": meta.get("atoms", []), "type": meta["type"]})
        active_valence = next_active_valence

        # ----------------------------------------------------
        # Subspace D: Localized Antibonds & Delocalized Pi
        # ----------------------------------------------------
        star_keys = list(star_dict.keys())
        curr_star_list = []
        for k in star_keys:
            v = vecs_dict[k].copy()
            meta = star_dict[k]
            p_arr = _calc_pop(v, is_virt=True)
            curr_star_list.append({"vec": v, "pop": p_arr, "meta": meta})

        next_active_antibond: List[Dict[str, Any]] = []
        if len(curr_star_list) > 0:
            if t == 0 or s_cross is None or len(active_antibond) == 0:
                for item in curr_star_list:
                    tid = f"trk_{next_seq:04d}"
                    next_seq += 1
                    meta = item["meta"]
                    _record_track_instance(
                        tid, item["vec"], item["pop"], meta["energy_ev"], meta["label"],
                        meta["type"], meta["category"], meta["orbital_class"], meta.get("atoms", []), "antibond"
                    )
                    next_active_antibond.append({"tid": tid, "vec": item["vec"], "atoms": meta.get("atoms", []), "type": meta["type"]})
            else:
                prev_smat = np.column_stack([ab["vec"] for ab in active_antibond])
                curr_smat = np.column_stack([item["vec"] for item in curr_star_list])

                # Degenerate pi* alignment
                Q_s = len(curr_star_list)
                for i1 in range(Q_s):
                    for i2 in range(i1 + 1, Q_s):
                        m1, m2 = curr_star_list[i1]["meta"], curr_star_list[i2]["meta"]
                        if "pi_star" in m1.get("type", "") and "pi_star" in m2.get("type", ""):
                            ats1 = tuple(sorted(m1.get("atoms", [])))
                            ats2 = tuple(sorted(m2.get("atoms", [])))
                            if len(ats1) == 2 and ats1 == ats2:
                                prev_smatches = [
                                    pi_idx for pi_idx, ab in enumerate(active_antibond)
                                    if "pi_star" in ab.get("type", "") and tuple(sorted(ab.get("atoms", []))) == ats1
                                ]
                                if len(prev_smatches) == 2:
                                    p_pair = prev_smat[:, prev_smatches]
                                    c_pair = curr_smat[:, [i1, i2]]
                                    c_pair_rot = _align_degenerate_2d_pair(p_pair, c_pair, s_cross)
                                    curr_smat[:, [i1, i2]] = c_pair_rot
                                    curr_star_list[i1]["vec"] = c_pair_rot[:, 0]
                                    curr_star_list[i2]["vec"] = c_pair_rot[:, 1]

                Omega_s = prev_smat.T @ s_cross @ curr_smat
                r_ind, c_ind = linear_sum_assignment(-np.abs(Omega_s))
                matched_scurr = set()
                for r, c in zip(r_ind, c_ind):
                    if abs(Omega_s[r, c]) >= 0.45:
                        matched_scurr.add(c)
                        tid = active_antibond[r]["tid"]
                        v = curr_smat[:, c].copy()
                        if Omega_s[r, c] < 0.0:
                            v = -v
                        meta = curr_star_list[c]["meta"]
                        _record_track_instance(
                            tid, v, curr_star_list[c]["pop"], meta["energy_ev"], meta["label"],
                            meta["type"], meta["category"], meta["orbital_class"], meta.get("atoms", []), "antibond"
                        )
                        next_active_antibond.append({"tid": tid, "vec": v, "atoms": meta.get("atoms", []), "type": meta["type"]})
                for c in range(len(curr_star_list)):
                    if c not in matched_scurr:
                        tid = f"trk_{next_seq:04d}"
                        next_seq += 1
                        meta = curr_star_list[c]["meta"]
                        _record_track_instance(
                            tid, curr_star_list[c]["vec"], curr_star_list[c]["pop"], meta["energy_ev"],
                            meta["label"], meta["type"], meta["category"], meta["orbital_class"],
                            meta.get("atoms", []), "antibond"
                        )
                        next_active_antibond.append({"tid": tid, "vec": curr_star_list[c]["vec"], "atoms": meta.get("atoms", []), "type": meta["type"]})
        active_antibond = next_active_antibond

        # Delocalized pi systems from c_meta (if any)
        curr_deloc_list = []
        for k, meta in c_meta.items():
            if k in f_vecs:
                v = f_vecs[k].copy()
                is_virt = "pistar" in meta.get("type", "")
                p_arr = _calc_pop(v, is_virt=is_virt)
                curr_deloc_list.append({"vec": v, "pop": p_arr, "meta": meta})

        next_active_deloc: List[Dict[str, Any]] = []
        if len(curr_deloc_list) > 0:
            if t == 0 or s_cross is None or len(active_deloc_pi) == 0:
                for item in curr_deloc_list:
                    tid = f"trk_{next_seq:04d}"
                    next_seq += 1
                    meta = item["meta"]
                    _record_track_instance(
                        tid, item["vec"], item["pop"], meta["energy_ev"], meta["label"],
                        meta["type"], meta["category"], meta["orbital_class"], meta.get("atoms", []), "deloc_pi"
                    )
                    next_active_deloc.append({"tid": tid, "vec": item["vec"]})
            else:
                prev_dmat = np.column_stack([ad["vec"] for ad in active_deloc_pi])
                curr_dmat = np.column_stack([item["vec"] for item in curr_deloc_list])
                Omega_d = prev_dmat.T @ s_cross @ curr_dmat
                r_ind, c_ind = linear_sum_assignment(-np.abs(Omega_d))
                matched_dcurr = set()
                for r, c in zip(r_ind, c_ind):
                    if abs(Omega_d[r, c]) >= 0.45:
                        matched_dcurr.add(c)
                        tid = active_deloc_pi[r]["tid"]
                        v = curr_dmat[:, c].copy()
                        if Omega_d[r, c] < 0.0:
                            v = -v
                        meta = curr_deloc_list[c]["meta"]
                        _record_track_instance(
                            tid, v, curr_deloc_list[c]["pop"], meta["energy_ev"], meta["label"],
                            meta["type"], meta["category"], meta["orbital_class"], meta.get("atoms", []), "deloc_pi"
                        )
                        next_active_deloc.append({"tid": tid, "vec": v})
                for c in range(len(curr_deloc_list)):
                    if c not in matched_dcurr:
                        tid = f"trk_{next_seq:04d}"
                        next_seq += 1
                        meta = curr_deloc_list[c]["meta"]
                        _record_track_instance(
                            tid, curr_deloc_list[c]["vec"], curr_deloc_list[c]["pop"], meta["energy_ev"],
                            meta["label"], meta["type"], meta["category"], meta["orbital_class"],
                            meta.get("atoms", []), "deloc_pi"
                        )
                        next_active_deloc.append({"tid": tid, "vec": curr_deloc_list[c]["vec"]})
        active_deloc_pi = next_active_deloc

        frame_outputs.append({
            "tids": np.array(frame_tids, dtype="U24"),
            "vecs": np.column_stack(frame_vecs).astype(np.float32) if frame_vecs else np.zeros((mol.nao, 0), dtype=np.float32),
            "pops": np.array(frame_pops, dtype=np.float32) if frame_pops else np.zeros((0, mol.natm), dtype=np.float32),
            "energies": np.array(frame_energies, dtype=np.float32),
            "charges": np.array(mull, dtype=np.float32),
        })

        prev_mol = mol

    # ----------------------------------------------------
    # Post-Pass: Smoothing and Summary Label Construction
    # ----------------------------------------------------
    registry_tracks: Dict[str, Dict[str, Any]] = {}

    for tid, rec in track_records.items():
        present = rec["frame_present"]
        roles = rec["frame_roles"]
        if not any(present):
            continue
        s_f = min(i for i, p in enumerate(present) if p)
        e_f = max(i for i, p in enumerate(present) if p)

        smoothed = list(roles)
        for i in range(s_f + 1, e_f):
            if smoothed[i - 1] is not None and smoothed[i + 1] is not None and smoothed[i] is not None:
                if smoothed[i - 1] == smoothed[i + 1] and smoothed[i] != smoothed[i - 1]:
                    smoothed[i] = smoothed[i - 1]

        # Collect distinct roles persisting for >= 3 frames
        filtered_roles = []
        cur_role = None
        cur_len = 0
        for i in range(s_f, e_f + 1):
            r = smoothed[i]
            if r == cur_role:
                cur_len += 1
            else:
                if cur_role is not None and cur_len >= 3:
                    if not filtered_roles or filtered_roles[-1] != cur_role:
                        filtered_roles.append(cur_role)
                cur_role = r
                cur_len = 1
        if cur_role is not None and cur_len >= 3:
            if not filtered_roles or filtered_roles[-1] != cur_role:
                filtered_roles.append(cur_role)

        if not filtered_roles:
            valid = [r for r in smoothed[s_f:e_f + 1] if r is not None]
            if valid:
                from collections import Counter
                filtered_roles = [Counter(valid).most_common(1)[0][0]]
            else:
                filtered_roles = [rec["type"]]

        r_first, r_last = filtered_roles[0], filtered_roles[-1]
        if r_first == r_last or len(filtered_roles) == 1:
            summary_label = r_first
            role_trans = r_first
        else:
            summary_label = f"{r_first} → {r_last}"
            role_trans = f"{r_first} → {r_last}"

        atoms_list = sorted(list(rec["atoms_set"]))
        registry_tracks[tid] = {
            "track_id": tid,
            "id": tid,
            "orbital_key": tid,
            "label": summary_label,
            "role_transition": role_trans,
            "short_label": role_trans,
            "category": rec["category"],
            "orbital_class": rec["orbital_class"],
            "type": rec["type"],
            "atoms": atoms_list,
            "start_frame": s_f,
            "end_frame": e_f,
            "frame_map": [tid if present[f] else None for f in range(N)],
            "frame_labels": smoothed,
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

    # ----------------------------------------------------
    # Serialization to .npz and registry.json
    # ----------------------------------------------------
    os.makedirs(TRAJECTORIES_DIR, exist_ok=True)
    npz_data: Dict[str, Any] = {"pathway_hash": np.array(p_hash)}
    for f in range(N):
        fo = frame_outputs[f]
        npz_data[f"f{f}_tids"] = fo["tids"]
        npz_data[f"f{f}_vecs"] = fo["vecs"]
        npz_data[f"f{f}_pops"] = fo["pops"]
        npz_data[f"f{f}_energies"] = fo["energies"]
        npz_data[f"f{f}_charges"] = fo["charges"]

    np.savez_compressed(NPZ_PATH, **npz_data)

    registry_json = {
        "pathway_hash": p_hash,
        "total_frames": N,
        "tracks": registry_tracks,
    }
    with open(REGISTRY_PATH, "w") as f:
        json.dump(registry_json, f, indent=2)

    megajson["global_electron_tracks"] = registry_tracks

    store_dict = {
        "pathway_hash": p_hash,
        "total_frames": N,
        "registry": registry_json,
        "npz": npz_data,
    }
    return store_dict


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

    for k, tid in enumerate(tids_f):
        trk = tracks.get(str(tid))
        if not trk:
            continue
        p_list = pops_f[k].tolist() if k < len(pops_f) else []
        e_val = float(energies_f[k]) if k < len(energies_f) else 0.0

        item = {
            "track_id": str(tid),
            "id": str(tid),
            "key": str(tid),
            "local_key": str(tid),
            "orbital_key": str(tid),
            "name": str(tid),
            "index": k,
            "type": trk["type"],
            "category": trk["category"],
            "orbital_class": trk["orbital_class"],
            "label": trk["label"],
            "frame_label": trk["frame_labels"][frame_idx] if frame_idx < len(trk["frame_labels"]) else trk["label"],
            "decoupled_label": trk["label"],
            "atoms": trk["atoms"],
            "energy_ev": round(e_val, 2),
            "pop": p_list,
            "palette": (
                {"pos": "#38bdf8", "neg": "#1d4ed8"} if trk["orbital_class"] == "bonding"
                else ({"pos": "#fb923c", "neg": "#db2777"} if trk["orbital_class"] == "antibonding"
                      else {"pos": "#94a3b8", "neg": "#475569"})
            ),
        }

        if trk["type"] == "frontier":
            f_key = next((fk for fk, mapped_tid in FRONTIER_TRACK_MAP.items() if mapped_tid == tid), tid)
            frontier_orbitals[f_key] = item
        elif trk["type"] in ("sigma_star", "pi_star") or trk["orbital_class"] == "antibonding":
            antibonds.append(item)
        elif trk["type"] in ("delocalized_pi", "delocalized_pistar"):
            canonical_pi.append(item)
        elif trk["type"] != "core":
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
