"""
Window Optimizer Module
Scans post-TS and pre-TS frames at step boundaries to find the minimal RMSD cut points.
"""

from typing import Any, Dict, List
from .graph_mapper import find_step_chemical_overlaps


def find_optimal_overlap_window(
    step_prev: Dict[str, Any],
    step_next: Dict[str, Any],
    max_window: int = 8
) -> Dict[str, Any]:
    """Scan candidate boundary frames across the junction window to find the frame pair with minimal Kabsch RMSD."""
    frames_p = step_prev.get("frames", [])
    frames_n = step_next.get("frames", [])
    if not frames_p or not frames_n:
        return {
            "best_prev_idx": -1, "best_next_idx": 0, "min_rmsd": 0.0, "default_rmsd": 0.0,
            "improvement": 0.0, "candidate_p_range": [0, 0], "candidate_n_range": [0, 0], "scan_results": []
        }

    n_p = len(frames_p)
    n_n = len(frames_n)
    ts_p = step_prev.get("ts_frame_index", 0)
    ts_n = step_next.get("ts_frame_index", n_n - 1)

    p_min_idx = max(ts_p + 1, n_p - max_window) if n_p > ts_p + 1 else max(0, n_p - max_window)
    e_p_last = float(frames_p[-1].get("free_energy_kcal", frames_p[-1].get("rel_energy_kcal", 0.0)))
    candidate_p = [
        idx for idx in range(p_min_idx, n_p)
        if abs(float(frames_p[idx].get("free_energy_kcal", frames_p[idx].get("rel_energy_kcal", 0.0))) - e_p_last) <= 0.10
    ] or [n_p - 1]

    n_max_idx = min(ts_n, max_window - 1) if ts_n > 0 else min(n_n - 1, max_window - 1)
    e_n_first = float(frames_n[0].get("free_energy_kcal", frames_n[0].get("rel_energy_kcal", 0.0)))
    candidate_n = [
        idx for idx in range(0, n_max_idx + 1)
        if abs(float(frames_n[idx].get("free_energy_kcal", frames_n[idx].get("rel_energy_kcal", 0.0))) - e_n_first) <= 0.10
    ] or [0]

    best_p, best_n = n_p - 1, 0
    best_rmsd, default_rmsd = 999.0, 999.0
    scan_results = []

    for fp in candidate_p:
        for fn in candidate_n:
            cands = find_step_chemical_overlaps(step_prev, step_next, prev_frame_idx=fp, next_frame_idx=fn)
            rmsd_val = cands[0]["rmsd"] if cands else 999.0

            if fp == n_p - 1 and fn == 0:
                default_rmsd = rmsd_val

            scan_results.append({"prev_idx": fp, "next_idx": fn, "rmsd": round(rmsd_val, 4)})

            if rmsd_val < best_rmsd:
                best_rmsd = rmsd_val
                best_p, best_n = fp, fn

    improvement = max(0.0, default_rmsd - best_rmsd) if default_rmsd < 900.0 else 0.0

    return {
        "best_prev_idx": best_p,
        "best_next_idx": best_n,
        "min_rmsd": round(float(best_rmsd if best_rmsd < 900.0 else 0.0), 4),
        "default_rmsd": round(float(default_rmsd if default_rmsd < 900.0 else 0.0), 4),
        "improvement": round(float(improvement), 4),
        "candidate_p_range": [candidate_p[0], candidate_p[-1]],
        "candidate_n_range": [candidate_n[0], candidate_n[-1]],
        "scan_results": scan_results
    }
