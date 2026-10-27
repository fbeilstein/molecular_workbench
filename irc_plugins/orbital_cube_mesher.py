"""
Orbital Cube Mesher and Caching Module
Computes volumetric orbital grids, serializes Gaussian .cube files, and manages disk/memory caches.
"""

import os
import threading
from typing import Any, Callable, Dict, List, Optional, Tuple
import numpy as np

from .ts_engine import orbital_engine
from .orbital_selection import resolve_orbital_key

ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
OUT_DIR = os.path.join(ROOT_DIR, 'out')
CUBES_DIR = os.path.join(OUT_DIR, 'cubes')

_GRID_CACHE: Dict[str, Any] = {}
_CUBE_CACHE: Dict[str, str] = {}
_CACHE_LOCK = threading.Lock()


def clear_orbital_caches():
    """Clear memory and on-disk Gaussian cube files."""
    with _CACHE_LOCK:
        _CUBE_CACHE.clear()
        _GRID_CACHE.clear()
        if os.path.exists(CUBES_DIR):
            for cf in os.listdir(CUBES_DIR):
                if cf.endswith(".cube"):
                    try:
                        os.remove(os.path.join(CUBES_DIR, cf))
                    except Exception:
                        pass


def get_grid(megajson: Dict[str, Any]) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Retrieve or build the global bounding grid for the reaction coordinate."""
    if "lo" not in _GRID_CACHE:
        all_pos = [f.get("coords", f.get("positions", [])) for f in megajson["frames"]]
        _GRID_CACHE["lo"], _GRID_CACHE["step"], _GRID_CACHE["npts"] = orbital_engine.build_global_grid(
            all_pos, grid_points=32, margin_angstrom=2.5
        )
    return _GRID_CACHE["lo"], _GRID_CACHE["step"], _GRID_CACHE["npts"]


def get_or_compute_cube(
    frame_idx: int,
    ok: str,
    ck: Optional[str],
    iso: float,
    d: Dict[str, Any],
    megajson: Dict[str, Any]
) -> Dict[str, Any]:
    """Fetch cached Gaussian cube text or evaluate molecular orbital on grid."""
    coords = megajson["frames"][frame_idx].get("coords", megajson["frames"][frame_idx].get("positions", []))
    rk = resolve_orbital_key(d, ok, ck)
    if not rk or rk not in d.get("vecs", {}):
        return {"success": False, "available": False, "error": f"Orbital '{ck or ok}' not formed at frame {frame_idx}"}

    cache_key = f"{frame_idx}_{ok}"
    disk_cube_path = os.path.join(CUBES_DIR, f"{frame_idx}_{rk}.cube")

    # 1. Check memory cache
    if cache_key in _CUBE_CACHE:
        cube_text = _CUBE_CACHE[cache_key]
        if cube_text.startswith(f"F{frame_idx}_{rk}"):
            return {
                "success": True, "cached": True, "frame_index": frame_idx,
                "orbital": ok, "resolved_orbital": rk, "chem_key": ck,
                "isovalue": iso, "cube": cube_text, "filename": f"{ok}_f{frame_idx}.cube", "coords": coords
            }

    # 2. Check disk cache
    if os.path.exists(disk_cube_path):
        try:
            with open(disk_cube_path, "r") as cf:
                first_line = cf.readline().strip()
                if first_line == f"F{frame_idx}_{rk}":
                    cf.seek(0)
                    cube_text = cf.read()
                    _CUBE_CACHE[cache_key] = cube_text
                    return {
                        "success": True, "cached": True, "frame_index": frame_idx,
                        "orbital": ok, "resolved_orbital": rk, "chem_key": ck,
                        "isovalue": iso, "cube": cube_text, "filename": f"{ok}_f{frame_idx}.cube", "coords": coords
                    }
        except Exception:
            pass

    # 3. Evaluate MO grid and serialize cube
    lo, step, npts = get_grid(megajson)
    v = orbital_engine.evaluate_mo_grid(d["mol"], d["vecs"][rk], lo, step, npts)
    fr = megajson["frames"][frame_idx]
    syms = fr.get("symbols") or megajson.get("symbols", [])
    cube_text = orbital_engine.generate_cube_text(
        f"F{frame_idx}_{rk}", syms, coords, v, lo, step, npts
    )
    _CUBE_CACHE[cache_key] = cube_text

    try:
        os.makedirs(CUBES_DIR, exist_ok=True)
        with open(disk_cube_path, "w") as cf:
            cf.write(cube_text)
    except Exception:
        pass

    return {
        "success": True, "cached": False, "frame_index": frame_idx,
        "orbital": ok, "resolved_orbital": rk, "chem_key": ck,
        "isovalue": iso, "cube": cube_text, "filename": f"{ok}_f{frame_idx}.cube", "coords": coords
    }


def batch_precompute_cubes(
    megajson: Dict[str, Any],
    sf: int,
    sz: int,
    sel: List[Any],
    force: bool,
    ensure_frame_computed_fn: Callable[[int, Dict[str, Any]], Dict[str, Any]]
) -> Dict[str, Any]:
    """Precompute batches of cubes for selected tracks along the reaction path."""
    n_frames = len(megajson.get("frames", []))
    ef = min(n_frames - 1, sf + sz - 1)
    lo, step, npts = get_grid(megajson)
    batch_cubes = {}

    os.makedirs(CUBES_DIR, exist_ok=True)
    for fi in range(sf, ef + 1):
        d = ensure_frame_computed_fn(fi, megajson)
        coords = megajson["frames"][fi].get("coords", megajson["frames"][fi].get("positions", []))
        for trk in sel:
            ok = trk.get("orbital_key") if isinstance(trk, dict) else str(trk)
            ck = trk.get("chem_key") if isinstance(trk, dict) else None
            if not ok or ok == "track_esp_charge":
                continue
            cache_key = f"{fi}_{ok}"
            rk = resolve_orbital_key(d, ok, ck)
            if not rk or rk not in d.get("vecs", {}):
                continue

            disk_cube_path = os.path.join(CUBES_DIR, f"{fi}_{rk}.cube")
            if not force:
                if cache_key in _CUBE_CACHE and _CUBE_CACHE[cache_key].startswith(f"F{fi}_{rk}"):
                    batch_cubes[cache_key] = _CUBE_CACHE[cache_key]
                    continue
                if os.path.exists(disk_cube_path):
                    try:
                        with open(disk_cube_path, "r") as cf:
                            if cf.readline().strip() == f"F{fi}_{rk}":
                                cf.seek(0)
                                cube_text = cf.read()
                                _CUBE_CACHE[cache_key] = cube_text
                                batch_cubes[cache_key] = cube_text
                                continue
                    except Exception:
                        pass

            v = orbital_engine.evaluate_mo_grid(d["mol"], d["vecs"][rk], lo, step, npts)
            fr_i = megajson["frames"][fi]
            syms_i = fr_i.get("symbols") or megajson.get("symbols", [])
            cube_text = orbital_engine.generate_cube_text(
                f"F{fi}_{rk}", syms_i, coords, v, lo, step, npts
            )
            _CUBE_CACHE[cache_key] = cube_text
            batch_cubes[cache_key] = cube_text
            try:
                with open(disk_cube_path, "w") as cf:
                    cf.write(cube_text)
            except Exception:
                pass

    return {
        "success": True, "start_frame": sf, "end_frame": ef, "total_frames": n_frames,
        "percent": int(round((ef + 1) / n_frames * 100)),
        "is_complete": bool(ef >= n_frames - 1),
        "cubes": batch_cubes
    }


def get_cached_cubes(frames: List[int], tracks: List[Any], scf_cache: Dict[int, Any]) -> Dict[str, str]:
    """Retrieve pre-existing cached cubes across requested frames and tracks."""
    found = {}
    for fi in frames:
        d = scf_cache.get(fi)
        for trk in tracks:
            ok = trk.get("orbital_key") if isinstance(trk, dict) else str(trk)
            ck = trk.get("chem_key") if isinstance(trk, dict) else None
            if not ok or ok == "track_esp_charge":
                continue
            cache_key = f"{fi}_{ok}"
            if cache_key in _CUBE_CACHE:
                found[cache_key] = _CUBE_CACHE[cache_key]
                continue
            rk = resolve_orbital_key(d, ok, ck) if d else ok
            if rk:
                disk_cube_path = os.path.join(CUBES_DIR, f"{fi}_{rk}.cube")
                if os.path.exists(disk_cube_path):
                    try:
                        with open(disk_cube_path, "r") as cf:
                            if cf.readline().strip() == f"F{fi}_{rk}":
                                cf.seek(0)
                                cube_text = cf.read()
                                _CUBE_CACHE[cache_key] = cube_text
                                found[cache_key] = cube_text
                    except Exception:
                        pass
    return found
