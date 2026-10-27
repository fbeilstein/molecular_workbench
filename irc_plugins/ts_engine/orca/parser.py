"""
ORCA Parser Module
Parses multi-frame trajectory files and builds trajectory JSON objects.
"""

import os
import re
from typing import Any, Dict, List, Optional
import numpy as np


def parse_orca_trj_file(filepath: str) -> List[Dict[str, Any]]:
    """Parse multi-frame .xyz trajectory from ORCA (_IRC_Full_trj.xyz or _MEP_trj.xyz)."""
    if not os.path.exists(filepath):
        return []
    with open(filepath, "r") as f:
        lines = [l.strip() for l in f if l.strip()]
    frames = []
    idx = 0
    while idx < len(lines):
        try:
            n_atoms = int(lines[idx])
        except ValueError:
            idx += 1
            continue
        if idx + 1 >= len(lines):
            break
        comment = lines[idx + 1]
        m = re.search(r"E\s+([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)", comment)
        energy_eh = float(m.group(1)) if m else 0.0
        pos = []
        syms = []
        for a_idx in range(idx + 2, min(len(lines), idx + 2 + n_atoms)):
            parts = lines[a_idx].split()
            if len(parts) >= 4:
                syms.append(parts[0])
                pos.append([float(parts[1]), float(parts[2]), float(parts[3])])
        if len(pos) == n_atoms:
            frames.append({
                "energy_eh": energy_eh,
                "energy_ev": energy_eh * 27.211386,
                "symbols": syms,
                "positions": pos
            })
        idx += 2 + n_atoms
    return frames


def build_trajectory_result(
    frames: List[Dict[str, Any]],
    ts_index: int,
    model: Any,
    atoms: Any,
    log: Any,
    energy_type: str = "electronic",
    thermo_info: Optional[Dict[str, Any]] = None
) -> Dict[str, Any]:
    """Package raw frames into full trajectory JSON format with arc length s and barriers."""
    EV_TO_KCAL = 23.060541945329334

    try:
        atom_masses = atoms.get_masses()
    except Exception:
        atom_masses = np.ones(model.n_atoms) * 12.0
    masses = np.array([atom_masses[k] if k < len(atom_masses) else 12.0 for k in range(model.n_atoms)])
    s_coords = [0.0] * len(frames)

    s_acc = 0.0
    for k in range(ts_index - 1, -1, -1):
        dr = np.array(frames[k]["positions"]) - np.array(frames[k + 1]["positions"])
        ds = float(np.sqrt(np.sum(masses[:, None] * (dr ** 2))))
        s_acc -= ds
        s_coords[k] = s_acc

    s_acc = 0.0
    for k in range(ts_index + 1, len(frames)):
        dr = np.array(frames[k]["positions"]) - np.array(frames[k - 1]["positions"])
        ds = float(np.sqrt(np.sum(masses[:, None] * (dr ** 2))))
        s_acc += ds
        s_coords[k] = s_acc

    e_r = frames[0]["energy_ev"]
    e_p = frames[-1]["energy_ev"]
    e_ts = frames[ts_index]["energy_ev"]
    e_complex = min(frames[k]["energy_ev"] for k in range(ts_index + 1))
    e_baseline = e_r

    out_frames = []
    for k, item in enumerate(frames):
        e_cur = item["energy_ev"]
        rel_kcal = (e_cur - e_baseline) * EV_TO_KCAL
        out_frames.append({
            "step": k - ts_index,
            "s_arc": round(s_coords[k], 4),
            "energy_ev": round(e_cur, 5),
            "rel_energy_kcal": round(rel_kcal, 3),
            "positions": item["positions"]
        })

    barrier_fwd = max(0.0, (e_ts - e_complex) * EV_TO_KCAL)
    barrier_rev = max(0.0, (e_ts - e_p) * EV_TO_KCAL)
    reaction_e = (e_p - e_r) * EV_TO_KCAL

    p_start = np.array(out_frames[0]["positions"])
    p_end = np.array(out_frames[-1]["positions"])
    sep_dist = float(np.linalg.norm(p_start - p_end))

    bonds = []
    ch_pairs = [(b[0], b[1]) for b in getattr(model, "order_changed", [])]
    for pair in sorted(getattr(model, "r_bonds", set()) | getattr(model, "p_bonds", set())):
        r_ord = int(round(getattr(model, "r_bond_orders", {}).get(pair, 0)))
        p_ord = int(round(getattr(model, "p_bond_orders", {}).get(pair, 0)))
        b_type = "forming" if pair in getattr(model, "formed", set()) else \
                 "breaking" if pair in getattr(model, "broken", set()) else \
                 "order_changed" if pair in ch_pairs else "spectator"
        bonds.append({"from": pair[0], "to": pair[1], "r_order": max(1, r_ord), "p_order": max(1, p_ord), "type": b_type})

    from ..fragment_relax import build_fragments_meta

    barriers = {
        "e_ts_ev": round(e_ts, 5),
        "e_reactant_ev": round(e_r, 5),
        "e_product_ev": round(e_p, 5),
        "forward_barrier_kcal": round(barrier_fwd, 2),
        "reverse_barrier_kcal": round(barrier_rev, 2),
        "reaction_energy_kcal": round(reaction_e, 2),
        "endpoint_separation_A": round(sep_dist, 2),
        "energy_type": energy_type
    }
    if thermo_info:
        barriers.update(thermo_info)

    result_json = {
        "reaction_smiles": f"{'.'.join(getattr(model, 'r_smiles', []))}>>{'.'.join(getattr(model, 'p_smiles', []))}",
        "symbols": model.symbols,
        "n_atoms": model.n_atoms,
        "bonds": bonds,
        "ts_frame_index": ts_index,
        "asymp_r_count": 0,
        "asymp_p_count": 0,
        "barriers": barriers,
        "fragments": build_fragments_meta(model),
        "frames": out_frames
    }
    if thermo_info:
        result_json["thermo_data"] = thermo_info

    log(f"\n✓ Reaction Path Generated: {len(out_frames)} frames.")
    log(f"  Forward Barrier  = {barriers['forward_barrier_kcal']:.2f} kcal/mol")
    log(f"  Reverse Barrier  = {barriers['reverse_barrier_kcal']:.2f} kcal/mol")
    log(f"  Reaction Energy  = {barriers['reaction_energy_kcal']:+.2f} kcal/mol")
    return result_json
