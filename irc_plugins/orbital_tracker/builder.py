"""
Orbital Store Builder — Parallel Frame SCF, Sequential Cross-Frame Assignment & Serialization.
"""

from concurrent.futures import ProcessPoolExecutor
import json
import os
from typing import Any, Dict, List, Optional, Tuple
import numpy as np
from pyscf import gto, lib as pyscf_lib

from ..ts_engine.orbital_engine import compute_frame_scf
from ..ts_engine.orbital_ibo import compute_frame_ibos
from ..ts_engine.orbital_frontier import extract_step0_planar_pi, extract_frontier_orbitals
from ..orbital_selection import get_system_charge, get_system_solvent, tag_orbital_classes
from .common import (
    FRONTIER_KEYS,
    FRONTIER_TRACK_MAP,
    NPZ_PATH,
    REGISTRY_PATH,
    TRAJECTORIES_DIR,
)
from .matcher import match_subspace_orbitals, smooth_and_label_tracks


def _worker_compute_frame(
    args: Tuple[int, List[str], List[List[float]], Optional[List[str]], int, Optional[str]]
) -> Tuple[int, List[float], Dict[str, Any], Dict[str, Any], Dict[str, np.ndarray], np.ndarray, Dict[str, Any], Dict[str, Any], Dict[str, np.ndarray]]:
    """Worker function for parallel solvated SCF + Step-0 conjugated pi separation + IBO evaluation."""
    t, symbols, coords, atom_labels, charge, solvent = args
    pyscf_lib.num_threads(1)
    mol, mf, mull = compute_frame_scf(
        symbols, coords, charge=charge, spin=0, basis="sto-3g", solvent=solvent
    )
    non_pi_occ, non_pi_virt, c_meta, c_vecs = extract_step0_planar_pi(
        mol, mf, atom_labels=atom_labels
    )
    ibos_dict, star_dict, vecs_dict, pops, ibo_coeff = compute_frame_ibos(
        mol, mf, symbols, atom_labels=atom_labels,
        non_pi_occ_coeff=non_pi_occ, non_pi_virt_coeff=non_pi_virt,
    )
    f_meta, c_meta, f_vecs = extract_frontier_orbitals(
        mol, mf, None, None, ibos_dict, star_dict, vecs_dict,
        atom_labels=atom_labels, step0_pi=(c_meta, c_vecs),
    )
    tag_orbital_classes(f_meta, ibos_dict, star_dict, c_meta)
    return (t, mull, ibos_dict, star_dict, vecs_dict, pops, f_meta, c_meta, f_vecs)


def _build_orbital_store(megajson: Dict[str, Any], p_hash: str) -> Dict[str, Any]:
    """Run parallel frame SCF and sequential tracking to build persistent store."""
    frames = megajson.get("frames", [])
    N = len(frames)
    if N == 0:
        return {"pathway_hash": p_hash, "total_frames": 0, "registry": {"tracks": {}}}

    solvent = get_system_solvent(megajson)

    # Step 1: Parallel evaluation
    tasks = []
    for t in range(N):
        fr = frames[t]
        symbols = fr.get("symbols") or megajson.get("symbols", [])
        coords = fr.get("coords", fr.get("positions", []))
        atom_labels = fr.get("atom_labels")
        if not atom_labels or len(atom_labels) != len(symbols):
            atom_labels = [f"{symbols[i]}{i}" for i in range(len(symbols))]
        charge = get_system_charge(megajson, t)
        tasks.append((t, symbols, coords, atom_labels, charge, solvent))

    workers = min(8, max(1, os.cpu_count() or 4))
    with ProcessPoolExecutor(max_workers=workers) as ex:
        precomputed_list = list(ex.map(_worker_compute_frame, tasks))

    precomputed_list.sort(key=lambda x: x[0])

    track_records: Dict[str, Dict[str, Any]] = {}
    frame_outputs: List[Dict[str, Any]] = []

    prev_mol: Optional[Any] = None
    prev_atom_labels: Optional[List[str]] = None
    prev_ao_atom: Optional[List[int]] = None
    active_core: List[Dict[str, Any]] = []
    active_valence: List[Dict[str, Any]] = []
    active_antibond: List[Dict[str, Any]] = []
    prev_frontier_vecs: Dict[str, np.ndarray] = {}

    next_seq = 1
    next_core_seq = 0

    # Step 2: Sequential cross-frame propagation across subspaces
    for t in range(N):
        fr = frames[t]
        symbols = fr.get("symbols") or megajson.get("symbols", [])
        coords = fr.get("coords", fr.get("positions", []))
        curr_atom_labels = fr.get("atom_labels")
        if not curr_atom_labels or len(curr_atom_labels) != len(symbols):
            curr_atom_labels = [f"{symbols[i]}{i}" for i in range(len(symbols))]
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
            if (
                prev_atom_labels is not None
                and prev_ao_atom is not None
                and set(prev_atom_labels) != set(curr_atom_labels)
            ):
                curr_set = set(curr_atom_labels)
                prev_set = set(prev_atom_labels)
                for mu, at_i in enumerate(prev_ao_atom):
                    if prev_atom_labels[at_i] not in curr_set:
                        s_cross[mu, :] = 0.0
                for nu, at_j in enumerate(ao_atom):
                    if curr_atom_labels[at_j] not in prev_set:
                        s_cross[:, nu] = 0.0

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
                    "frame_types": [None] * N,
                    "frame_categories": [None] * N,
                    "frame_atoms": [None] * N,
                    "frame_energies": [None] * N,
                    "frame_pops": [None] * N,
                }
            rec = track_records[tid]
            rec["frame_present"][t] = True
            rec["frame_roles"][t] = role_lbl
            rec["frame_types"][t] = typ
            rec["frame_categories"][t] = cat
            rec["frame_atoms"][t] = list(atoms)
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

        active_core, next_core_seq = match_subspace_orbitals(
            active_core, curr_core_list, s_cross, t,
            subspace="core", threshold=0.60, next_seq=next_core_seq,
            record_fn=_record_track_instance, prefix="trk_core_", pad_zeros=2
        )

        # ----------------------------------------------------
        # Subspace C: Occupied Valence Space (Localized IBOs + Step-0 Conjugated Pi)
        # ----------------------------------------------------
        val_keys = [k for k, v in ibos_dict.items() if v.get("type") != "core"]
        curr_val_list = []
        for k in val_keys:
            v = vecs_dict[k].copy()
            meta = ibos_dict[k]
            p_arr = pops[meta["index"]] if meta["index"] < len(pops) else _calc_pop(v)
            curr_val_list.append({"vec": v, "pop": p_arr, "meta": meta})
        for k, meta in c_meta.items():
            if k in f_vecs and "pistar" not in meta.get("type", ""):
                v = f_vecs[k].copy()
                p_arr = _calc_pop(v, is_virt=False)
                curr_val_list.append({"vec": v, "pop": p_arr, "meta": meta})

        active_valence, next_seq = match_subspace_orbitals(
            active_valence, curr_val_list, s_cross, t,
            subspace="valence_ibo", threshold=0.45, next_seq=next_seq,
            record_fn=_record_track_instance, prefix="trk_", pad_zeros=4,
            pi_type_keyword="pi"
        )

        # ----------------------------------------------------
        # Subspace D: Virtual Space (Localized Antibonds + Step-0 Conjugated Pi*)
        # ----------------------------------------------------
        star_keys = list(star_dict.keys())
        curr_virt_list = []
        for k in star_keys:
            v = vecs_dict[k].copy()
            meta = star_dict[k]
            p_arr = _calc_pop(v, is_virt=True)
            curr_virt_list.append({"vec": v, "pop": p_arr, "meta": meta})
        for k, meta in c_meta.items():
            if k in f_vecs and "pistar" in meta.get("type", ""):
                v = f_vecs[k].copy()
                p_arr = _calc_pop(v, is_virt=True)
                curr_virt_list.append({"vec": v, "pop": p_arr, "meta": meta})

        active_antibond, next_seq = match_subspace_orbitals(
            active_antibond, curr_virt_list, s_cross, t,
            subspace="antibond", threshold=0.45, next_seq=next_seq,
            record_fn=_record_track_instance, prefix="trk_", pad_zeros=4,
            pi_type_keyword="pi_star"
        )

        frame_outputs.append({
            "tids": np.array(frame_tids, dtype="U24"),
            "vecs": np.column_stack(frame_vecs).astype(np.float32) if frame_vecs else np.zeros((mol.nao, 0), dtype=np.float32),
            "pops": np.array(frame_pops, dtype=np.float32) if frame_pops else np.zeros((0, mol.natm), dtype=np.float32),
            "energies": np.array(frame_energies, dtype=np.float32),
            "charges": np.array(mull, dtype=np.float32),
        })

        prev_mol = mol
        prev_atom_labels = curr_atom_labels
        prev_ao_atom = ao_atom

    # Step 3: Track smoothing and disambiguated role labeling
    registry_tracks = smooth_and_label_tracks(track_records, N)

    # Step 4: Serialization to .npz and registry.json
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

    return {
        "pathway_hash": p_hash,
        "total_frames": N,
        "registry": registry_json,
        "npz": npz_data,
    }

