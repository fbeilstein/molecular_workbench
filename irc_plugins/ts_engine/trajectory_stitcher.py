"""
Master Trajectory Stitcher
Orchestrates multistep IRC cascade alignment, topological reordering, and geodesic bridging.
Decomposed into modular subpackages under irc_plugins.ts_engine.alignment.
"""

import json
from typing import Any, Dict, List, Optional, Set, Tuple
import numpy as np
from rdkit import Chem
from .alignment import (
    kabsch_fit,
    compute_rmsd,
    align_coordinates,
    find_step_chemical_overlaps,
    generate_bridge_frames,
    build_multistep_bridge,
    find_optimal_overlap_window,
    extract_fragment_molblock,
)

# Backward-compatibility alias
find_all_hero_overlaps = find_step_chemical_overlaps


def get_canonical_smiles_set(smiles_str: str) -> Set[str]:
    """Extract canonical SMILES for individual molecules in a dot-separated SMILES string."""
    if not smiles_str:
        return set()
    res = set()
    for part in smiles_str.split("."):
        part = part.strip()
        if not part:
            continue
        try:
            m = Chem.MolFromSmiles(part)
            if m:
                res.add(Chem.MolToSmiles(Chem.RemoveHs(m)))
        except Exception:
            pass
    return res


def sort_trajectories_chemically(trajectories: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Sort trajectories into a contiguous cascade based on product -> reactant coincidence."""
    if len(trajectories) <= 1:
        return trajectories
    n = len(trajectories)
    step_data = []
    for traj in trajectories:
        rxn = traj.get("reaction_smiles", "")
        r_smi = get_canonical_smiles_set(rxn.split(">>")[0] if ">>" in rxn else "")
        p_smi = get_canonical_smiles_set(rxn.split(">>")[1] if ">>" in rxn else "")
        step_data.append({"traj": traj, "r": r_smi, "p": p_smi})

    chain = [0]
    used = {0}
    for _ in range(n - 1):
        curr_p = step_data[chain[-1]]["p"]
        best_next, best_score = -1, 0
        for cand_idx in range(n):
            if cand_idx in used:
                continue
            common = len(curr_p & step_data[cand_idx]["r"])
            if common > best_score:
                best_score = common
                best_next = cand_idx
        if best_next != -1 and best_score > 0:
            chain.append(best_next)
            used.add(best_next)
        else:
            for rem in range(n):
                if rem not in used:
                    chain.append(rem)
                    used.add(rem)
                    break
    return [trajectories[i] for i in chain]


def stitch_multistep_trajectories(
    trajectories: List[Dict[str, Any]],
    junction_choices: Optional[Dict[str, str]] = None,
    n_bridge_frames: int = 8,
    auto_sort: bool = False,
    auto_trim: bool = True,
    balance_energies: bool = True,
    solvent: str = "water"
) -> Dict[str, Any]:
    """Stitch multiple single-step IRC trajectories into a unified cascade pathway."""
    if not trajectories:
        return {"title": "Empty Pathway", "n_steps": 0, "total_frames": 0, "frames": []}

    steps = sort_trajectories_chemically(trajectories) if auto_sort else [json.loads(json.dumps(t)) for t in trajectories]
    choices = junction_choices or {}
    aligned_steps = [json.loads(json.dumps(s)) for s in steps]
    step_reports, all_bridges = [], []

    for k in range(1, len(aligned_steps)):
        prev_step, next_step = aligned_steps[k - 1], aligned_steps[k]
        orig_p_len, orig_n_len = len(prev_step["frames"]), len(next_step["frames"])

        win_report = find_optimal_overlap_window(prev_step, next_step, max_window=8)
        p_cut = win_report["best_prev_idx"] if auto_trim else orig_p_len - 1
        n_cut = win_report["best_next_idx"] if auto_trim else 0

        all_cands = find_step_chemical_overlaps(prev_step, next_step, prev_frame_idx=p_cut, next_frame_idx=n_cut)
        cand = next((c for c in all_cands if c["id"] == choices.get(str(k))), None) if all_cands and choices.get(str(k)) else (all_cands[0] if all_cands else None)

        out_next = aligned_steps[k]
        if cand:
            R, t = np.array(cand["R"]), np.array(cand["t"])
            idxs_p, idxs_n = cand["prev_indices"], cand["next_indices"]
            prev_syms, next_syms = prev_step.get("symbols", []), next_step.get("symbols", [])
            order = list(idxs_n) + [c for c in range(len(next_syms)) if c not in set(idxs_n)]

            out_next["symbols"] = [next_syms[c] for c in order]
            inv_order = {c: r for r, c in enumerate(order)}
            remapped_bonds = []
            for b in next_step.get("bonds", []):
                u = b.get("from", b[0] if isinstance(b, (list, tuple)) else 0)
                v = b.get("to", b[1] if isinstance(b, (list, tuple)) else 0)
                if u in inv_order and v in inv_order:
                    nb = dict(b)
                    nb["from"], nb["to"] = inv_order[u], inv_order[v]
                    remapped_bonds.append(nb)
            out_next["bonds"] = remapped_bonds

            if "fragments" in out_next:
                for side in ("reactants", "products"):
                    for fg in out_next["fragments"].get(side, []):
                        if "atom_indices" in fg:
                            fg["atom_indices"] = [inv_order.get(idx, idx) for idx in fg["atom_indices"]]

            for fr in out_next["frames"]:
                raw_pos = np.array(fr.get("coords", fr.get("positions", [])))
                fr["coords"] = ((raw_pos @ R.T + t)[order]).round(4).tolist()
                fr.pop("positions", None)
            b_rmsd = cand["rmsd"]
        else:
            b_rmsd, idxs_p, idxs_n = 0.0, [], []

        if p_cut < len(prev_step["frames"]) - 1:
            prev_step["frames"] = prev_step["frames"][:p_cut + 1]
        if n_cut > 0:
            out_next["frames"] = out_next["frames"][n_cut:]
            out_next["ts_frame_index"] = max(0, out_next.get("ts_frame_index", 0) - n_cut)

        prev_pos = np.array(prev_step["frames"][-1].get("coords", prev_step["frames"][-1].get("positions", [])))
        align_pos = np.array(out_next["frames"][0].get("coords", out_next["frames"][0].get("positions", [])))
        e_p = prev_step["frames"][-1].get("rel_energy_kcal", 0.0)
        e_n = out_next["frames"][0].get("rel_energy_kcal", 0.0)

        bridge = build_multistep_bridge(prev_pos, align_pos, e_p, e_n, prev_step, out_next, n_bridge_frames, idxs_p=idxs_p)

        raw_next_pos = np.array(aligned_steps[k]["frames"][n_cut].get("coords", aligned_steps[k]["frames"][n_cut].get("positions", [])))
        j_info = {
            "junction_index": k, "step_prev_index": k, "step_next_index": k + 1,
            "hero_formula": cand.get("formula", "Intermediate") if cand else "Intermediate",
            "prev_symbols": prev_step.get("symbols", []), "prev_coords": prev_pos.round(4).tolist(),
            "next_symbols": out_next.get("symbols", []), "next_coords": align_pos.round(4).tolist(),
            "next_coords_aligned": align_pos.round(4).tolist(), "next_coords_unaligned": raw_next_pos.round(4).tolist(),
            "matched_prev_indices": idxs_p, "matched_next_indices": list(range(len(idxs_p))) if cand else [],
            "all_candidates": all_cands, "selected_candidate_id": cand["id"] if cand else None,
            "n_matched_atoms": len(idxs_p), "joint_rmsd_angstrom": b_rmsd, "rmsd": b_rmsd,
            "bridge_frames": bridge, "n_bridge_frames": len(bridge),
            "prev_cut": p_cut, "next_cut": n_cut, "is_trimmed": (p_cut < orig_p_len - 1 or n_cut > 0)
        }
        step_reports.append(j_info)
        all_bridges.append(bridge)

    master_frames, step_boundaries, curr_idx, accum_s, curr_energy_offset = [], [], 0, 0.0, 0.0
    symbols = max((s.get("symbols", []) for s in aligned_steps), key=len, default=[])
    max_energy, final_energy = 0.0, 0.0
    master_bonds_dict = {}
    for step in aligned_steps:
        for b in step.get("bonds", []):
            p = (min(b.get("from", 0), b.get("to", 0)), max(b.get("from", 0), b.get("to", 0)))
            master_bonds_dict[p] = dict(b)

    for k, step in enumerate(aligned_steps):
        s_frames, ts_idx = step.get("frames", []), step.get("ts_frame_index", 0)
        step_boundaries.append({
            "step_index": k + 1, "step_name": step.get("reaction_smiles", f"Step {k + 1}"),
            "start_frame": curr_idx, "end_frame": curr_idx + len(s_frames) - 1, "ts_frame": curr_idx + ts_idx
        })
        k_bonds = step.get("bonds", [])
        raw_s = [float(fr.get("s_arc", idx * 0.1)) for idx, fr in enumerate(s_frames)]
        s_base, step_start_s = raw_s[0], (0.0 if k == 0 else accum_s + 0.15)
        raw_e = [float(fr.get("rel_energy_kcal", fr.get("energy_kcal", 0.0))) for fr in s_frames]
        e_base = raw_e[0] if raw_e else 0.0

        for f_i, fr in enumerate(s_frames):
            c = fr.get("coords", fr.get("positions", []))
            fr["coords"], fr["positions"], fr["step_index"], fr["is_ts"] = c, c, k + 1, (f_i == ts_idx)
            fr["symbols"] = step.get("symbols", [])
            e_k = curr_energy_offset + (raw_e[f_i] - e_base)
            fr["energy_kcal"], fr["rel_energy_kcal"], fr["s_arc"] = round(float(e_k), 2), round(float(e_k), 2), round(step_start_s + (raw_s[f_i] - s_base), 3)
            accum_s, max_energy, final_energy = fr["s_arc"], max(max_energy, e_k), e_k
            step_bonds = []
            for b in k_bonds:
                b_type = b.get("type")
                r_ord = 0 if b_type == "forming" else b.get("r_order", 1)
                p_ord = 0 if b_type == "breaking" else b.get("p_order", 1)
                b_ord = r_ord if f_i <= ts_idx else p_ord
                if b_ord > 0:
                    step_bonds.append({"from": b.get("from", 0), "to": b.get("to", 0), "order": b_ord})
            fr["bonds"] = step_bonds
            master_frames.append(fr)
            curr_idx += 1

        step_end_e = curr_energy_offset + (raw_e[-1] - e_base) if raw_e else curr_energy_offset

        if k < len(all_bridges) and all_bridges[k]:
            for b_fr in all_bridges[k]:
                c = b_fr.get("coords", b_fr.get("positions", []))
                b_fr["coords"], b_fr["positions"], b_fr["step_index"], b_fr["is_ts"] = c, c, f"Bridge {k + 1}➔{k + 2}", False
                accum_s += 0.15
                b_fr["energy_kcal"] = round(float(step_end_e), 2)
                b_fr["rel_energy_kcal"], b_fr["s_arc"] = b_fr["energy_kcal"], round(float(accum_s), 3)
                max_energy = max(max_energy, b_fr["energy_kcal"])
                master_frames.append(b_fr)
                curr_idx += 1

        curr_energy_offset = step_end_e

    res = {
        "title": "Master Multistep IRC Trajectory", "n_steps": len(aligned_steps), "total_frames": len(master_frames),
        "symbols": symbols, "step_symbols": [step.get("symbols", []) for step in aligned_steps],
        "bonds": list(master_bonds_dict.values()), "step_boundaries": step_boundaries,
        "stitching_alignment_report": step_reports, "frames": master_frames,
        "overall_barrier_kcal": round(float(max_energy), 2), "overall_reaction_energy_kcal": round(float(final_energy), 2)
    }

    if balance_energies and len(trajectories) >= 1:
        try:
            from .energy_balancer import balance_pathway_energies
            res = balance_pathway_energies(res, trajectories, solvent=solvent, recalculate_with_xtb=False)
        except Exception as eb_err:
            import logging
            logging.getLogger("irc_workbench").warning(f"Failed to balance pathway energies: {eb_err}")

    return res


# Re-export from alignment package for backward compatibility
extract_fragment_molblock = extract_fragment_molblock
