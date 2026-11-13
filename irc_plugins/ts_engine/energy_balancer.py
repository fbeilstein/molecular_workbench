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
    """Calculate Hess's Law balanced Free Energy (G) profile with standard-state isolated fragment thermodynamics in chosen solvent."""
    boundaries = megajson.get("step_boundaries", [])
    master_frames = megajson.get("frames", [])
    n_frames = len(master_frames)
    if n_frames == 0:
        return megajson

    stoich = analyze_pathway_stoichiometry(trajectories, fragment_roles)
    solvent = (solvent or megajson.get("solvent", "water") or "water").strip()

    # Collect all chemical species involved across elementary steps
    all_species = []
    for s in stoich.get("steps", []):
        all_species.extend(s.get("reactants", []))
        all_species.extend(s.get("products", []))
    all_species.extend(stoich.get("all_external_reactants", []))
    all_species.extend(stoich.get("all_external_byproducts", []))

    frag_energies = {}
    if recalculate_with_xtb and all_species:
        try:
            frag_energies = compute_isolated_fragment_energies(all_species, solvent=solvent, relax=True)
        except Exception as fe_err:
            import logging
            logging.getLogger("irc_workbench").warning(f"Could not evaluate isolated fragment energies: {fe_err}")

    eV_to_kcal = 23.0605

    # Check if we have valid fragment energies to compute thermodynamic Hess baselines
    has_thermo = False
    step_delta_kcal = []
    if frag_energies and stoich.get("steps"):
        valid_steps = 0
        from .free_energy_engine import DELTA_G_ASSOC_KCAL
        for k, s in enumerate(stoich["steps"]):
            traj = trajectories[k] if k < len(trajectories) else {}
            barriers = traj.get("barriers", {})
            e_r = sum(frag_energies.get(_clean_smiles(r["smiles"]), 0.0) for r in s.get("reactants", []))
            e_p = sum(frag_energies.get(_clean_smiles(p["smiles"]), 0.0) for p in s.get("products", []))
            e_r_traj = barriers.get("e_reactant_ev")
            e_p_traj = barriers.get("e_product_ev")
            cross_engine = (
                (e_r_traj is not None and e_r != 0.0 and abs(float(e_r_traj) - e_r) > 20.0)
                or (e_p_traj is not None and e_p != 0.0 and abs(float(e_p_traj) - e_p) > 20.0)
                or str(traj.get("engine", "")).lower() == "orca"
            )
            if not cross_engine and e_r != 0.0 and e_p != 0.0:
                valid_steps += 1
                n_r_k = max(1, len(s.get("reactants", [])))
                n_p_k = max(1, len(s.get("products", [])))
                delta_n_entropy = (n_r_k - n_p_k) * DELTA_G_ASSOC_KCAL
                step_delta_kcal.append((e_p - e_r) * eV_to_kcal + delta_n_entropy)
            else:
                # Preserve the step trajectory's native ΔG_rxn when it comes from a different level
                # of theory (e.g. ORCA DFT mixed with xTB steps) or lacks matching xTB fragments.
                frs_k = traj.get("frames", [])
                if barriers.get("energy_type") == "gibbs" and barriers.get("reaction_energy_kcal") is not None:
                    step_delta_kcal.append(float(barriers["reaction_energy_kcal"]))
                    valid_steps += 1
                elif len(frs_k) >= 2:
                    e0 = float(frs_k[0].get("free_energy_kcal", frs_k[0].get("rel_energy_kcal", 0.0)))
                    e1 = float(frs_k[-1].get("free_energy_kcal", frs_k[-1].get("rel_energy_kcal", 0.0)))
                    step_delta_kcal.append(e1 - e0)
                    valid_steps += 1
                else:
                    step_delta_kcal.append(0.0)
        if valid_steps == len(stoich["steps"]) and len(boundaries) == len(stoich["steps"]):
            has_thermo = True

    rel_kcal = np.zeros(n_frames)

    if has_thermo:
        # Full Hess's Law balancing with standard-state thermodynamics in solvent
        cum_baseline = 0.0
        for k, b in enumerate(boundaries):
            start, end, ts = b["start_frame"], b["end_frame"], b["ts_frame"]
            traj = trajectories[k] if k < len(trajectories) else {}
            de_step = step_delta_kcal[k]

            e_start_target = cum_baseline
            e_end_target = cum_baseline + de_step

            # Extract raw IRC profile for this step
            s_len = end - start + 1
            raw = np.array([float(master_frames[i].get("rel_energy_kcal", 0.0)) for i in range(start, end + 1)])
            raw = raw - raw[0]

            # Forward barrier height from quantum calculations
            raw_barrier = float(np.max(raw)) if len(raw) > 0 else 0.0
            rep_barrier = float(traj.get("barriers", {}).get("forward_barrier_kcal", 0.0) or 0.0)
            barrier = max(raw_barrier, rep_barrier, 0.0)

            # Determine whether this step has an authentic barrier or is monotonic downhill / barrierless
            is_downhill = (raw_barrier <= 1e-4 and rep_barrier <= 1e-4) or (raw[-1] < -0.1 and raw_barrier <= 0.05)

            if is_downhill:
                # Strictly preserve authentic downhill IRC shape without any synthetic bump or barrier
                raw_delta = raw[-1]
                for local_i, global_i in enumerate(range(start, end + 1)):
                    if abs(raw_delta) > 1e-6:
                        progress = raw[local_i] / raw_delta
                        rel_kcal[global_i] = e_start_target + progress * de_step
                    else:
                        rel_kcal[global_i] = e_start_target + raw[local_i]
            else:
                # Step possesses an authentic TS barrier
                ts_rel_idx = ts - start if start <= ts <= end else int(np.argmax(raw))
                peak_rel_idx = int(np.argmax(raw))
                raw_at_ts = raw[ts_rel_idx] if 0 <= ts_rel_idx < len(raw) else 0.0
                if raw_barrier > 1.0 and (raw_at_ts < 0.75 * raw_barrier or raw_at_ts <= 1e-4):
                    ts_rel_idx = peak_rel_idx
                    b["ts_frame"] = start + ts_rel_idx
                    for fr_i in range(start, end + 1):
                        master_frames[fr_i]["is_ts"] = (fr_i == b["ts_frame"])

                scale_idx = peak_rel_idx if raw_barrier > 1e-4 else ts_rel_idx
                e_ts_target = e_start_target + barrier
                raw_ts = raw[scale_idx] if scale_idx < len(raw) else raw_barrier
                raw_denom_prod = raw[-1] - raw_ts

                for local_i, global_i in enumerate(range(start, end + 1)):
                    if local_i <= scale_idx:
                        # Follow authentic quantum IRC barrier curve up to apex
                        if raw_ts > 1e-4:
                            f = raw[local_i] / raw_ts
                        else:
                            f = local_i / max(1, scale_idx)
                        rel_kcal[global_i] = e_start_target + f * (e_ts_target - e_start_target)
                    else:
                        # Follow authentic quantum IRC descent curve from apex down to product
                        if abs(raw_denom_prod) > 1e-4:
                            f = (raw[local_i] - raw_ts) / raw_denom_prod
                        else:
                            denom_frames = max(1, s_len - 1 - scale_idx)
                            f = (local_i - scale_idx) / denom_frames
                        rel_kcal[global_i] = e_ts_target + f * (e_end_target - e_ts_target)

            cum_baseline = e_end_target

            # Fill bridge frames to next step smoothly at the intermediate baseline
            if k < len(boundaries) - 1:
                next_start = boundaries[k + 1]["start_frame"]
                if end + 1 < next_start:
                    for br_i in range(end + 1, next_start):
                        rel_kcal[br_i] = e_end_target
        balanced_ev = rel_kcal / eV_to_kcal
    else:
        # Fallback if fragment energies are unavailable: relative IRC curve
        cum_baseline = 0.0
        for k, b in enumerate(boundaries):
            start, end = b["start_frame"], b["end_frame"]
            raw = np.array([float(master_frames[i].get("rel_energy_kcal", 0.0)) for i in range(start, end + 1)])
            raw = raw - raw[0]
            for local_i, global_i in enumerate(range(start, end + 1)):
                rel_kcal[global_i] = cum_baseline + raw[local_i]
            cum_baseline += raw[-1]
            if k < len(boundaries) - 1:
                next_start = boundaries[k + 1]["start_frame"]
                for br_i in range(end + 1, next_start):
                    rel_kcal[br_i] = cum_baseline
        balanced_ev = rel_kcal / eV_to_kcal

    step_metrics = [
        {
            "step_index": k + 1, "step_name": b.get("step_name", f"Step {k+1}"),
            "start_frame": b["start_frame"], "ts_frame": b["ts_frame"], "end_frame": b["end_frame"],
            "barrier_kcal": round(float(max(0.0, np.max(rel_kcal[b["start_frame"]:b["end_frame"]+1]) - rel_kcal[b["start_frame"]])), 2),
            "reaction_energy_kcal": round(float(rel_kcal[b["end_frame"]] - rel_kcal[b["start_frame"]]), 2),
            "ts_energy_kcal": round(float(rel_kcal[b["ts_frame"]]), 2),
            "reactant_energy_kcal": round(float(rel_kcal[b["start_frame"]]), 2),
            "product_energy_kcal": round(float(rel_kcal[b["end_frame"]]), 2)
        }
        for k, b in enumerate(boundaries)
    ]
    for b, m in zip(boundaries, step_metrics):
        b.update(m)

    for i, fr in enumerate(master_frames):
        fr["energy_kcal"] = round(float(rel_kcal[i]), 3)
        fr["rel_energy_kcal"] = round(float(rel_kcal[i]), 3)
        fr["balanced_energy_ev"] = round(float(balanced_ev[i]), 5)

    megajson["overall_barrier_kcal"] = round(float(max(np.max(rel_kcal), 0.0)), 2)
    megajson["overall_reaction_energy_kcal"] = round(float(rel_kcal[-1]), 2)
    megajson["step_boundaries"] = boundaries
    megajson["energy_profile"] = {
        "frame_energies_kcal": rel_kcal.round(3).tolist(),
        "step_metrics": step_metrics,
        "overall_barrier_kcal": megajson["overall_barrier_kcal"],
        "overall_reaction_energy_kcal": megajson["overall_reaction_energy_kcal"],
        "initial_reactants": [r["formula"] for r in stoich["all_external_reactants"]],
        "final_products": [p["formula"] for p in stoich["all_external_byproducts"]]
    }
    megajson["stoichiometry_summary"], megajson["solvent"] = stoich, solvent
    return megajson
