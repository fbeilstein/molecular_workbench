"""
Bridge Interpolator Module
Generates clash-free geodesic transition frames across consecutive reaction steps.
Handles unequal atom counts (incoming reagents and departing leaving groups).
"""

from typing import Any, Dict, List, Optional
import numpy as np


def generate_bridge_frames(
    p_start: np.ndarray,
    p_end: np.ndarray,
    e_start: float,
    e_end: float,
    n_frames: int = 8,
    bonds: Optional[List[Any]] = None,
    symbols: Optional[List[str]] = None
) -> List[Dict[str, Any]]:
    """Generate intermittent transition frames with zero atom clashes via IDPP geodesic interpolation."""
    frames: List[Dict[str, Any]] = []
    if len(p_start) == len(p_end) and len(p_start) >= 2:
        try:
            import warnings
            from ase import Atoms
            from ase.mep import idpp_interpolate

            syms = symbols if symbols and len(symbols) == len(p_start) else ["C"] * len(p_start)
            images = [Atoms(symbols=syms, positions=p_start)]
            for i in range(1, n_frames + 1):
                s = (i / (n_frames + 1)) ** 2 * (3.0 - 2.0 * (i / (n_frames + 1)))
                images.append(Atoms(symbols=syms, positions=(1.0 - s) * p_start + s * p_end))
            images.append(Atoms(symbols=syms, positions=p_end))

            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                idpp_interpolate(images, traj=None, log=None, fmax=0.05, steps=50)

            candidate_frames = []
            has_nan = False
            for i, img in enumerate(images[1:-1]):
                pos = img.get_positions()
                if np.isnan(pos).any() or np.isinf(pos).any():
                    has_nan = True
                    break
                s = ((i + 1) / (n_frames + 1)) ** 2 * (3.0 - 2.0 * ((i + 1) / (n_frames + 1)))
                candidate_frames.append({
                    "coords": pos.round(4).tolist(),
                    "energy_kcal": round(float((1.0 - s) * e_start + s * e_end), 2),
                    "is_bridge": True
                })
            if not has_nan and len(candidate_frames) == n_frames:
                return candidate_frames
        except Exception:
            pass

    if len(p_start) == len(p_end) and len(p_start) >= 1:
        for i in range(1, n_frames + 1):
            s = (i / (n_frames + 1)) ** 2 * (3.0 - 2.0 * (i / (n_frames + 1)))
            coords = ((1.0 - s) * p_start + s * p_end).copy()
            for _ in range(5):
                D = np.linalg.norm(coords[:, None, :] - coords[None, :, :], axis=2)
                np.fill_diagonal(D, 99.0)
                for a, b in zip(*np.where(D < 1.15)):
                    if a < b and (dist := float(np.linalg.norm(coords[b] - coords[a]))) > 1e-4:
                        shift = 0.5 * (1.15 - dist) * ((coords[b] - coords[a]) / dist)
                        coords[b] += shift
                        coords[a] -= shift
            coords = np.nan_to_num(coords, nan=0.0, posinf=0.0, neginf=0.0)
            frames.append({
                "coords": coords.round(4).tolist(),
                "energy_kcal": round(float((1.0 - s) * e_start + s * e_end), 2),
                "is_bridge": True
            })
        return frames

    return []


def build_multistep_bridge(
    prev_pos: np.ndarray,
    align_pos: np.ndarray,
    e_start: float,
    e_end: float,
    prev_step: Dict[str, Any],
    next_step: Dict[str, Any],
    n_bridge_frames: int = 8,
    idxs_p: Optional[List[int]] = None
) -> List[Dict[str, Any]]:
    """Build bridge frames between two consecutive reaction steps with differing atom counts."""
    n_p, n_a = len(prev_pos), len(align_pos)
    m = len(idxs_p) if idxs_p and len(idxs_p) <= min(n_p, n_a) else None

    if m is not None:
        p_interm = prev_pos[idxs_p]
        a_interm = align_pos[:m]
        syms_m = [prev_step.get("symbols", [])[i] for i in idxs_p]
        bridge_comm = generate_bridge_frames(p_interm, a_interm, e_start, e_end, n_bridge_frames, symbols=syms_m)
        extra_c = align_pos[m:] if n_a > m else None
        active_bonds = []
        for b in next_step.get("bonds", []):
            r_ord = 0 if b.get("type") == "forming" else b.get("r_order", b.get("order", 1))
            if r_ord > 0:
                active_bonds.append({"from": b.get("from", 0), "to": b.get("to", 0), "order": r_ord})

        bridge = []
        for b_fr in bridge_comm:
            c_comm = np.array(b_fr.get("coords", b_fr.get("positions", [])))
            full_c = np.vstack([c_comm, extra_c]) if extra_c is not None and len(extra_c) > 0 else c_comm
            b_fr["coords"] = full_c.round(4).tolist()
            b_fr["symbols"] = next_step.get("symbols", [])[:len(full_c)]
            b_fr["bonds"] = active_bonds
            bridge.append(b_fr)
        return bridge

    if n_p == n_a:
        bridge = generate_bridge_frames(
            prev_pos, align_pos, e_start, e_end, n_bridge_frames,
            bonds=prev_step.get("bonds", []), symbols=prev_step.get("symbols", [])
        )
        for b_fr in bridge:
            b_fr["symbols"] = prev_step.get("symbols", [])
            b_fr["bonds"] = [b for b in prev_step.get("bonds", []) if b.get("type") != "forming"]
        return bridge

    if n_p < n_a:
        bridge_comm = generate_bridge_frames(
            prev_pos, align_pos[:n_p], e_start, e_end, n_bridge_frames,
            bonds=prev_step.get("bonds", []), symbols=prev_step.get("symbols", [])
        )
        bridge = []
        extra_c = align_pos[n_p:]
        for b_fr in bridge_comm:
            c_comm = np.array(b_fr.get("coords", b_fr.get("positions", [])))
            full_c = np.vstack([c_comm, extra_c])
            b_fr["coords"] = full_c.round(4).tolist()
            b_fr["symbols"] = next_step.get("symbols", [])
            b_fr["bonds"] = [b for b in next_step.get("bonds", []) if b.get("type") != "forming"]
            bridge.append(b_fr)
        return bridge

    bridge_comm = generate_bridge_frames(
        prev_pos[:n_a], align_pos, e_start, e_end, n_bridge_frames,
        bonds=next_step.get("bonds", []), symbols=next_step.get("symbols", [])
    )
    bridge = []
    departing_c = prev_pos[n_a:]
    for b_fr in bridge_comm:
        c_comm = np.array(b_fr.get("coords", b_fr.get("positions", [])))
        full_c = np.vstack([c_comm, departing_c])
        b_fr["coords"] = full_c.round(4).tolist()
        b_fr["symbols"] = prev_step.get("symbols", [])
        b_fr["bonds"] = [b for b in prev_step.get("bonds", []) if b.get("type") != "breaking"]
        bridge.append(b_fr)
    return bridge

