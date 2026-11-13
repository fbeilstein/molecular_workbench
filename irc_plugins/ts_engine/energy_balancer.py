import numpy as np
from typing import Any, Dict, List, Optional
from rdkit import Chem
from .common import (
    EV_TO_KCAL,
    clean_smiles as _clean_smiles,
    compute_isolated_fragment_energies,
    get_formal_charge,
    rcov,
)


def analyze_pathway_stoichiometry(
    trajectories: List[Dict[str, Any]], fragment_roles: Optional[Dict[str, str]] = None
) -> Dict[str, Any]:
    """Map global reaction stoichiometry with automatic & user-configurable roles."""
    roles = fragment_roles or {}
    n_steps = len(trajectories)
    step_info, all_external_reactants, all_external_byproducts = [], [], []

    for k, traj in enumerate(trajectories):
        frags = traj.get("fragments", {})
        symbols = traj.get("symbols", [])
        f0_pos = traj.get("frames", [{}])[0].get("positions", [])
        flast_pos = traj.get("frames", [{}])[-1].get("positions", [])

        r_list = [
            {
                "id": f"s{k+1}_r{r_idx}", "step_index": k + 1, "formula": r.get("formula", ""),
                "smiles": _clean_smiles(r.get("smiles", "")) or f"R_{k}_{r_idx}",
                "charge": get_formal_charge(r.get("smiles", "")),
                "num_atoms": r.get("num_atoms", len(r.get("atom_indices", []))),
                "atom_indices": r.get("atom_indices", []),
                "symbols": [symbols[i] for i in r.get("atom_indices", [])] if symbols else [],
                "positions": [f0_pos[i] for i in r.get("atom_indices", [])] if f0_pos else [],
                "role": roles.get(f"s{k+1}_r{r_idx}", roles.get(_clean_smiles(r.get("smiles", "")), "auto"))
            }
            for r_idx, r in enumerate(frags.get("reactants", []))
        ]

        p_list = [
            {
                "id": f"s{k+1}_p{p_idx}", "step_index": k + 1, "formula": p.get("formula", ""),
                "smiles": _clean_smiles(p.get("smiles", "")) or f"P_{k}_{p_idx}",
                "charge": get_formal_charge(p.get("smiles", "")),
                "num_atoms": p.get("num_atoms", len(p.get("atom_indices", []))),
                "atom_indices": p.get("atom_indices", []),
                "symbols": [symbols[i] for i in p.get("atom_indices", [])] if symbols else [],
                "positions": [flast_pos[i] for i in p.get("atom_indices", [])] if flast_pos else [],
                "role": roles.get(f"s{k+1}_p{p_idx}", roles.get(_clean_smiles(p.get("smiles", "")), "auto"))
            }
            for p_idx, p in enumerate(frags.get("products", []))
        ]

        step_info.append({"step_index": k + 1, "reactants": r_list, "products": p_list})

    # Forward resolution: resolve core vs incoming
    produced = {}
    for k, s in enumerate(step_info):
        s["incoming_external"] = []
        max_r_atoms = max((r["num_atoms"] for r in s["reactants"]), default=0)
        for r in s["reactants"]:
            if r["role"] == "auto":
                if r["smiles"] in produced and produced[r["smiles"]] > 0:
                    produced[r["smiles"]] -= 1
                    r["role"] = "core"
                elif k == 0:
                    r["role"] = "core" if (r["num_atoms"] == max_r_atoms and len(s["reactants"]) > 1) else "incoming"
                else:
                    r["role"] = "incoming"
            if r["role"] == "incoming":
                s["incoming_external"].append(r)
                all_external_reactants.append(r)
        for p in s["products"]:
            produced[p["smiles"]] = produced.get(p["smiles"], 0) + 1

    # Backward resolution: resolve core vs leaving
    consumed = {}
    for k in reversed(range(n_steps)):
        s = step_info[k]
        s["departing_byproducts"] = []
        max_p_atoms = max((p["num_atoms"] for p in s["products"]), default=0)
        for p in s["products"]:
            if p["role"] == "auto":
                if p["smiles"] in consumed and consumed[p["smiles"]] > 0:
                    consumed[p["smiles"]] -= 1
                    p["role"] = "core"
                elif k == n_steps - 1:
                    p["role"] = "core" if (p["num_atoms"] == max_p_atoms and len(s["products"]) > 1) else "leaving"
                else:
                    p["role"] = "leaving"
            if p["role"] == "leaving":
                s["departing_byproducts"].append(p)
                all_external_byproducts.append(p)
        for r in s["reactants"]:
            consumed[r["smiles"]] = consumed.get(r["smiles"], 0) + 1

    for k in range(n_steps):
        step_info[k]["waiting_reagents"] = [r for m in range(k + 1, n_steps) for r in step_info[m]["incoming_external"]]
        step_info[k]["departed_byproducts"] = [p for m in range(0, k) for p in step_info[m]["departing_byproducts"]]

    return {"steps": step_info, "all_external_reactants": all_external_reactants, "all_external_byproducts": all_external_byproducts}


def _smoothstep(x: float) -> float:
    x = float(np.clip(x, 0.0, 1.0))
    return x * x * (3.0 - 2.0 * x)


def _hermite_spline(t: np.ndarray, y0: float, y1: float, g0: float, g1: float) -> np.ndarray:
    h00 = 2.0 * t**3 - 3.0 * t**2 + 1.0
    h10 = t**3 - 2.0 * t**2 + t
    h01 = -2.0 * t**3 + 3.0 * t**2
    h11 = t**3 - t**2
    return h00 * y0 + h10 * g0 + h01 * y1 + h11 * g1


def balance_pathway_energies(
    megajson: Dict[str, Any],
    trajectories: List[Dict[str, Any]],
    recalculate_with_xtb: bool = True,
    solvent: str = "water",
    fragment_roles: Optional[Dict[str, str]] = None,
    truncate_at_minimum: bool = False
) -> Dict[str, Any]:
    """
    Assemble the multi-step cascade energy profile using rigid Hess's-Law vertical stage offsets.
    Strictly preserves authentic quantum free-energy barrier heights and thermodynamic profile shapes
    without synthetic arm-warping.
    """
    boundaries = megajson.get("step_boundaries", [])
    master_frames = megajson.get("frames", [])
    n_frames = len(master_frames)
    if n_frames == 0:
        return megajson

    stoich = analyze_pathway_stoichiometry(trajectories, fragment_roles)
    solvent = (solvent or megajson.get("solvent", "water") or "water").strip()
    eV_to_kcal = 23.0605

    rel_kcal = np.zeros(n_frames)
    cum_baseline = 0.0

    for k, b in enumerate(boundaries):
        start, end = b["start_frame"], b["end_frame"]
        traj = trajectories[k] if k < len(trajectories) else {}
        is_barrierless = bool(b.get("is_barrierless", False) or traj.get("is_barrierless", False))

        # Extract authentic step profile
        raw = np.array([float(master_frames[i].get("free_energy_kcal", master_frames[i].get("rel_energy_kcal", 0.0))) for i in range(start, end + 1)])
        if len(raw) > 0:
            raw = raw - raw[0]
        else:
            raw = np.zeros(max(1, end - start + 1))

        # Apply rigid vertical offset to preserve exact intrinsic barrier and curve shape
        for local_i, global_i in enumerate(range(start, end + 1)):
            rel_kcal[global_i] = cum_baseline + raw[local_i]

        step_rxn_kcal = raw[-1] if len(raw) > 0 else 0.0
        cum_baseline += step_rxn_kcal

        # Bridge frames to next step
        if k < len(boundaries) - 1:
            next_start = boundaries[k + 1]["start_frame"]
            if end + 1 < next_start:
                for br_i in range(end + 1, next_start):
                    rel_kcal[br_i] = cum_baseline

        # Handle TS marker and barrier metrics
        if is_barrierless:
            b["is_barrierless"] = True
            b["ts_frame"] = None
            for fr_i in range(start, end + 1):
                master_frames[fr_i]["is_ts"] = False
        else:
            raw_peak_idx = int(np.argmax(raw))
            b["is_barrierless"] = False
            b["ts_frame"] = start + raw_peak_idx
            for fr_i in range(start, end + 1):
                master_frames[fr_i]["is_ts"] = bool(fr_i == b["ts_frame"])

    balanced_ev = rel_kcal / eV_to_kcal

    step_metrics = []
    for k, b in enumerate(boundaries):
        start, end, ts = b["start_frame"], b["end_frame"], b["ts_frame"]
        is_bl = b.get("is_barrierless", False)
        step_raw = rel_kcal[start : end + 1]
        barrier = 0.0 if is_bl else float(max(0.0, np.max(step_raw) - rel_kcal[start]))
        rxn_e = float(rel_kcal[end] - rel_kcal[start])
        ts_e = None if is_bl or ts is None else round(float(rel_kcal[ts]), 2)

        m = {
            "step_index": k + 1,
            "step_name": b.get("step_name", f"Step {k+1}"),
            "start_frame": start,
            "ts_frame": ts,
            "end_frame": end,
            "is_barrierless": is_bl,
            "barrier_kcal": round(barrier, 2),
            "reaction_energy_kcal": round(rxn_e, 2),
            "ts_energy_kcal": ts_e,
            "reactant_energy_kcal": round(float(rel_kcal[start]), 2),
            "product_energy_kcal": round(float(rel_kcal[end]), 2)
        }
        b.update(m)
        step_metrics.append(m)

    for i, fr in enumerate(master_frames):
        fr["energy_kcal"] = round(float(rel_kcal[i]), 3)
        fr["rel_energy_kcal"] = round(float(rel_kcal[i]), 3)
        fr["free_energy_kcal"] = round(float(rel_kcal[i]), 3)
        fr["balanced_energy_ev"] = round(float(balanced_ev[i]), 5)

    megajson["overall_barrier_kcal"] = round(float(max(np.max(rel_kcal), 0.0)), 2)
    megajson["overall_reaction_energy_kcal"] = round(float(rel_kcal[-1]), 2)
    megajson["step_boundaries"] = boundaries
    megajson["energy_profile"] = {
        "frame_energies_kcal": rel_kcal.round(3).tolist(),
        "step_metrics": step_metrics,
        "overall_barrier_kcal": megajson["overall_barrier_kcal"],
        "overall_reaction_energy_kcal": megajson["overall_reaction_energy_kcal"],
        "initial_reactants": [r["formula"] for r in stoich.get("all_external_reactants", [])],
        "final_products": [p["formula"] for p in stoich.get("all_external_byproducts", [])]
    }
    megajson["stoichiometry_summary"] = stoich
    megajson["solvent"] = solvent
    return megajson
