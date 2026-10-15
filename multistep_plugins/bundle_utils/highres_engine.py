import os
from typing import Any, Dict, List, Optional
from irc_plugins.ts_engine import orbital_engine
from ..step_4_orbitals import Step4Orbitals, _resolve_orbital_key

_HIGHRES_GRID_CACHE: Dict[int, Any] = {}
_HIGHRES_CUBE_CACHE: Dict[str, str] = {}
ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
CUBE_DIR = os.path.join(ROOT_DIR, "out", "cubes")
os.makedirs(CUBE_DIR, exist_ok=True)


def get_highres_grid(megajson: Dict[str, Any], resolution: int):
    """Compute and cache invariant 3D bounding box grid at specified resolution."""
    if resolution not in _HIGHRES_GRID_CACHE:
        all_pos = [f["coords"] for f in megajson["frames"]]
        lo, step, npts = orbital_engine.build_global_grid(all_pos, grid_points=resolution, margin_angstrom=2.5)
        _HIGHRES_GRID_CACHE[resolution] = (lo, step, npts)
    return _HIGHRES_GRID_CACHE[resolution]


def get_cached_highres_cube(frame_idx: int, orb_key: str, megajson: Dict[str, Any]) -> Optional[str]:
    """Retrieve high-res cube from memory cache or disk if available."""
    cache_k = f"{frame_idx}_{orb_key}"
    if cache_k in _HIGHRES_CUBE_CACHE:
        return _HIGHRES_CUBE_CACHE[cache_k]

    # Check on-disk cubes
    cube_file = os.path.join(CUBE_DIR, f"{cache_k}.cube")
    if os.path.exists(cube_file):
        try:
            with open(cube_file, "r", encoding="utf-8") as f:
                content = f.read()
                _HIGHRES_CUBE_CACHE[cache_k] = content
                return content
        except Exception:
            pass

    # Fallback to string if present in legacy state
    val = megajson.get("highres_cubes", {}).get(cache_k)
    if isinstance(val, str) and len(val) > 200:
        _HIGHRES_CUBE_CACHE[cache_k] = val
        return val
    return None


def recompute_highres_batch(
    megajson: Dict[str, Any],
    tracks: List[Any],
    start_frame: int,
    end_frame: int,
    resolution: int = 50,
    isovalue: float = 0.03
) -> Dict[str, Any]:
    """Compute high-resolution grid and volumetric orbital cubes for a range of frames."""
    lo, step, npts = get_highres_grid(megajson, resolution)
    s4 = Step4Orbitals()
    total_frames = len(megajson["frames"])
    ef = min(end_frame, total_frames - 1)
    sf = max(0, start_frame)
    computed_cubes = {}

    highres_dict = megajson.setdefault("highres_cubes", {})
    megajson["orbital_resolution"] = resolution

    for fi in range(sf, ef + 1):
        d = s4._ensure_frame_computed(fi, megajson)
        symbols = megajson["symbols"]
        coords = megajson["frames"][fi]["coords"]

        for trk in tracks:
            if isinstance(trk, dict):
                ok = trk.get("orbital_key")
                ck = trk.get("chem_key")
            else:
                ok = str(trk)
                ck = None

            if not ok or ok == "track_esp_charge":
                continue

            rk = _resolve_orbital_key(d, ok, ck)
            if rk and rk in d.get("vecs", {}):
                kr = f"{resolution}_{fi}_{rk}"
                if kr not in _HIGHRES_CUBE_CACHE:
                    v = orbital_engine.evaluate_mo_grid(d["mol"], d["vecs"][rk], lo, step, npts)
                    cube_str = orbital_engine.generate_cube_text(
                        f"F{fi}_{rk}_res{resolution}", symbols, coords, v, lo, step, npts
                    )
                    _HIGHRES_CUBE_CACHE[kr] = cube_str

                cube_content = _HIGHRES_CUBE_CACHE[kr]
                _HIGHRES_CUBE_CACHE[f"{fi}_{ok}"] = cube_content
                cube_filename = f"{fi}_{ok}.cube"
                cube_disk_path = os.path.join(CUBE_DIR, cube_filename)
                try:
                    with open(cube_disk_path, "w", encoding="utf-8") as cf:
                        cf.write(cube_content)
                except Exception as ce:
                    print(f"[HighresEngine] Warning: could not write {cube_disk_path}: {ce}")
                highres_dict[f"{fi}_{ok}"] = f"trajectories/cubes/{cube_filename}"
                computed_cubes[f"{fi}_{ok}"] = True

    pct = int(round((ef + 1) / total_frames * 100))
    return {
        "success": True,
        "resolution": resolution,
        "start_frame": sf,
        "end_frame": ef,
        "total_frames": total_frames,
        "percent": pct,
        "is_complete": bool(ef >= total_frames - 1),
        "cubes_count": len(computed_cubes)
    }
