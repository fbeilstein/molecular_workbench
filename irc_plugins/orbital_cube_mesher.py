"""
Orbital Cube Mesher and Caching Module
Computes volumetric orbital grids, serializes Gaussian .cube files, and manages disk/memory caches.
Keys strictly by invariant track_id with zero legacy key fallback.
"""

import os
import threading
from typing import Any, Callable, Dict, List, Optional, Tuple
import numpy as np

from .ts_engine import orbital_engine
from .orbital_tracker import get_track_vector_and_mol

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


def _read_cached_disk_cube(disk_cube_path: str, header_expected: str) -> Optional[str]:
    """Read cube file from disk if present and matching expected header."""
    if not os.path.exists(disk_cube_path):
        return None
    try:
        with open(disk_cube_path, "r") as cf:
            if cf.readline().strip() == header_expected:
                cf.seek(0)
                return cf.read()
    except Exception:
        pass
    return None


def get_or_compute_cube(
    frame_idx: int,
    ok: str,
    ck: Optional[str],
    iso: float,
    _d: Optional[Dict[str, Any]],
    megajson: Dict[str, Any],
    track_id: Optional[str] = None,
    _ensure_frame_computed_fn: Optional[Callable[[int, Dict[str, Any]], Dict[str, Any]]] = None,
) -> Dict[str, Any]:
    """Fetch cached Gaussian cube text or evaluate molecular orbital on grid using exact track_id."""
    tid = track_id or ok
    if not tid or tid in ("track_esp_charge", "trk_esp"):
        return {"success": False, "available": False, "track_id": tid}

    cache_key = f"{frame_idx}_{tid}"
    coords = megajson["frames"][frame_idx].get("coords", megajson["frames"][frame_idx].get("positions", []))
    disk_cube_path = os.path.join(CUBES_DIR, f"{cache_key}.cube")
    header_expected = f"F{frame_idx}_{tid}"

    # 1. Check memory cache
    if cache_key in _CUBE_CACHE:
        cube_text = _CUBE_CACHE[cache_key]
        if cube_text.startswith(header_expected):
            return {
                "success": True, "cached": True, "frame_index": frame_idx,
                "track_id": tid, "orbital": tid, "resolved_orbital": tid, "chem_key": ck,
                "isovalue": iso, "cube": cube_text, "filename": f"{tid}_f{frame_idx}.cube", "coords": coords
            }

    # 2. Check disk cache
    cached_disk = _read_cached_disk_cube(disk_cube_path, header_expected)
    if cached_disk:
        _CUBE_CACHE[cache_key] = cached_disk
        return {
            "success": True, "cached": True, "frame_index": frame_idx,
            "track_id": tid, "orbital": tid, "resolved_orbital": tid, "chem_key": ck,
            "isovalue": iso, "cube": cached_disk, "filename": f"{tid}_f{frame_idx}.cube", "coords": coords
        }

    # 3. Retrieve orbital vector and PySCF Mol directly from orbital_tracker
    mol, vec, meta = get_track_vector_and_mol(megajson, frame_idx, tid)
    if vec is None or mol is None:
        return {
            "success": False, "available": False, "track_id": tid,
            "message": f"Orbital '{tid}' not present at frame {frame_idx}"
        }

    # 4. Evaluate volumetric grid and serialize cube text
    lo, step, npts = get_grid(megajson)
    v = orbital_engine.evaluate_mo_grid(mol, vec, lo, step, npts)
    fr = megajson["frames"][frame_idx]
    syms = fr.get("symbols") or megajson.get("symbols", [])
    cube_text = orbital_engine.generate_cube_text(header_expected, syms, coords, v, lo, step, npts)
    _CUBE_CACHE[cache_key] = cube_text

    try:
        os.makedirs(CUBES_DIR, exist_ok=True)
        with open(disk_cube_path, "w") as cf:
            cf.write(cube_text)
    except Exception:
        pass

    return {
        "success": True, "cached": False, "frame_index": frame_idx,
        "track_id": tid, "orbital": tid, "resolved_orbital": tid, "chem_key": ck,
        "isovalue": iso, "cube": cube_text, "filename": f"{tid}_f{frame_idx}.cube", "coords": coords
    }


def batch_precompute_cubes(
    megajson: Dict[str, Any],
    sf: int,
    sz: int,
    sel: List[Any],
    force: bool,
    _ensure_frame_computed_fn: Optional[Callable[[int, Dict[str, Any]], Dict[str, Any]]] = None,
) -> Dict[str, Any]:
    """Precompute batches of cubes for selected tracks along the reaction path."""
    n_frames = len(megajson.get("frames", []))
    ef = min(n_frames - 1, sf + sz - 1)
    lo, step, npts = get_grid(megajson)
    batch_cubes = {}

    os.makedirs(CUBES_DIR, exist_ok=True)
    for fi in range(sf, ef + 1):
        coords = megajson["frames"][fi].get("coords", megajson["frames"][fi].get("positions", []))
        fr_i = megajson["frames"][fi]
        syms_i = fr_i.get("symbols") or megajson.get("symbols", [])

        for trk in sel:
            tid = trk.get("track_id") if isinstance(trk, dict) else (trk.get("orbital_key") if isinstance(trk, dict) else str(trk))
            if not tid or tid in ("track_esp_charge", "trk_esp"):
                continue

            cache_key = f"{fi}_{tid}"
            disk_cube_path = os.path.join(CUBES_DIR, f"{cache_key}.cube")
            header_expected = f"F{fi}_{tid}"

            if not force:
                if cache_key in _CUBE_CACHE and _CUBE_CACHE[cache_key].startswith(header_expected):
                    batch_cubes[cache_key] = _CUBE_CACHE[cache_key]
                    continue
                cached_disk = _read_cached_disk_cube(disk_cube_path, header_expected)
                if cached_disk:
                    _CUBE_CACHE[cache_key] = cached_disk
                    batch_cubes[cache_key] = cached_disk
                    continue

            mol, vec, meta = get_track_vector_and_mol(megajson, fi, tid)
            if vec is None or mol is None:
                continue

            v = orbital_engine.evaluate_mo_grid(mol, vec, lo, step, npts)
            cube_text = orbital_engine.generate_cube_text(header_expected, syms_i, coords, v, lo, step, npts)
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


def get_cached_cubes(
    frames: List[int],
    tracks: List[Any],
    megajson_or_cache: Any = None,
    megajson: Optional[Dict[str, Any]] = None
) -> Dict[str, str]:
    """Retrieve pre-existing cached cubes across requested frames and tracks."""
    found = {}
    for fi in frames:
        for trk in tracks:
            tid = trk.get("track_id") if isinstance(trk, dict) else (trk.get("orbital_key") if isinstance(trk, dict) else str(trk))
            if not tid or tid in ("track_esp_charge", "trk_esp"):
                continue
            cache_key = f"{fi}_{tid}"
            header_expected = f"F{fi}_{tid}"
            if cache_key in _CUBE_CACHE and _CUBE_CACHE[cache_key].startswith(header_expected):
                found[cache_key] = _CUBE_CACHE[cache_key]
                continue
            disk_cube_path = os.path.join(CUBES_DIR, f"{cache_key}.cube")
            cached_disk = _read_cached_disk_cube(disk_cube_path, header_expected)
            if cached_disk:
                _CUBE_CACHE[cache_key] = cached_disk
                found[cache_key] = cached_disk
    return found
