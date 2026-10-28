import json
import os
import sys
import zipfile
from typing import Any, Dict, List, Optional, Tuple
from .highres_engine import _HIGHRES_CUBE_CACHE, CUBE_DIR

TOOLS_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "tools"))
if TOOLS_DIR not in sys.path:
    sys.path.insert(0, TOOLS_DIR)
from cube_converter import CubeConverter

_FAILED_MESH_KEYS = set()


def convert_and_pack_trajectory_orbitals(
    megajson: Dict[str, Any],
    zf: zipfile.ZipFile,
    isovalue: float = 0.035,
    global_state: Optional[Dict[str, Any]] = None
) -> Tuple[List[str], int]:
    """Extract and pack orbital surface meshes for all frames across configured timeline tracks."""
    timeline = megajson.get("orbital_timeline") or (global_state.get("step_6_orbitals", {}).get("orbital_timeline") if global_state else None) or {}
    raw_tracks = timeline.get("tracks", [])
    frames = megajson.get("frames", [])
    total_frames = len(frames)

    active_tracks = []
    for trk in raw_tracks:
        if isinstance(trk, dict):
            ok = trk.get("track_id") or trk.get("orbital_key")
            ck = trk.get("chem_key")
            clips = trk.get("clips", [])
            f_map = trk.get("frame_map")
        else:
            ok, ck, clips, f_map = str(trk), None, [], None
        if not ok or str(ok).startswith("track_esp") or str(ok) == "trk_esp":
            continue
        active_tracks.append((ok, ck, clips, f_map))

    if not active_tracks and frames:
        for tree_cat in megajson.get("orbitals_tree", []):
            for it in tree_cat.get("items", []):
                k = it.get("track_id") or it.get("key") or it.get("orbital_key")
                if k and k not in [a[0] for a in active_tracks]:
                    active_tracks.append((k, it.get("chem_key"), [], it.get("frame_map")))

    needed_targets: List[Tuple[int, str, Optional[str]]] = []
    for ok, ck, clips, f_map in active_tracks:
        if clips:
            f_set = set()
            for c in clips:
                s, e = max(0, int(c.get("start", 0))), min(total_frames - 1, int(c.get("end", total_frames - 1)))
                f_set.update(range(s, e + 1))
        else:
            f_set = set(range(total_frames))
        for fi in sorted(f_set):
            if f_map and 0 <= fi < len(f_map) and f_map[fi] is None:
                continue
            needed_targets.append((fi, ok, ck))

    packed_keys: List[str] = []
    uncompressed_bytes = 0
    written_entries = set(zf.namelist())

    def _pack_mesh_data(key: str, p_bytes: Optional[bytes], n_bytes: Optional[bytes]):
        nonlocal uncompressed_bytes
        has_m = False
        p_name, n_name = f"orbitals/{key}_pos.json", f"orbitals/{key}_neg.json"
        if p_bytes and p_name not in written_entries:
            zf.writestr(p_name, p_bytes)
            written_entries.add(p_name)
            uncompressed_bytes += len(p_bytes)
            has_m = True
        if n_bytes and n_name not in written_entries:
            zf.writestr(n_name, n_bytes)
            written_entries.add(n_name)
            uncompressed_bytes += len(n_bytes)
            has_m = True
        if has_m and key not in packed_keys:
            packed_keys.append(key)

    missing_targets = list(needed_targets)
    sources: Dict[str, str] = {}
    if os.path.exists(CUBE_DIR):
        for f in os.listdir(CUBE_DIR):
            if f.endswith(".cube"):
                sources[f[:-5]] = os.path.join(CUBE_DIR, f)

    for k, text in _HIGHRES_CUBE_CACHE.items():
        if isinstance(text, str) and "\n" in text:
            parts = k.split("_")
            clean_k = f"{parts[1]}_{'_'.join(parts[2:])}" if (len(parts) > 2 and parts[0].isdigit() and parts[1].isdigit()) else k
            sources.setdefault(clean_k, text)

    try:
        from ..step_6_orbitals import _CUBE_CACHE
        for k, text in _CUBE_CACHE.items():
            if k not in sources:
                sources[k] = text
    except Exception:
        pass

    still_missing: List[Tuple[int, str, Optional[str]]] = []
    for fi, ok, ck in missing_targets:
        key = f"{fi}_{ok}"
        if key in sources:
            try:
                vol, origin, step = CubeConverter.read_cube(sources[key])
                pos = CubeConverter.extract_mesh(vol, origin, step, isovalue)
                neg = CubeConverter.extract_mesh(vol, origin, step, -isovalue)
                _pack_mesh_data(key, json.dumps(pos, separators=(',', ':')).encode('utf-8') if pos else None, json.dumps(neg, separators=(',', ':')).encode('utf-8') if neg else None)
            except Exception as e:
                still_missing.append((fi, ok, ck))
        else:
            still_missing.append((fi, ok, ck))

    if still_missing and total_frames > 0:
        try:
            from ..step_6_orbitals import Step6Orbitals
            s6 = Step6Orbitals()
            for fi, ok, ck in still_missing:
                key = f"{fi}_{ok}"
                try:
                    res = s6.custom_action({"action_type": "get_orbital_mesh", "frame_index": fi, "track_id": ok, "orbital": ok, "chem_key": ck, "isovalue": isovalue}, global_state or {"step_6_orbitals": {"megajson": megajson}})
                    cube_text = res.get("cube")
                    if cube_text:
                        vol, origin, step = CubeConverter.read_cube(cube_text)
                        pos = CubeConverter.extract_mesh(vol, origin, step, isovalue)
                        neg = CubeConverter.extract_mesh(vol, origin, step, -isovalue)
                        _pack_mesh_data(key, json.dumps(pos, separators=(',', ':')).encode('utf-8') if pos else None, json.dumps(neg, separators=(',', ':')).encode('utf-8') if neg else None)
                    elif res.get("pos") or res.get("neg"):
                        _pack_mesh_data(key, json.dumps(res["pos"], separators=(',', ':')).encode('utf-8') if res.get("pos") else None, json.dumps(res["neg"], separators=(',', ':')).encode('utf-8') if res.get("neg") else None)
                    else:
                        _FAILED_MESH_KEYS.add(key)
                except Exception:
                    _FAILED_MESH_KEYS.add(key)
        except Exception:
            pass

    print(f"[OrbitalPackager] Successfully packed {len(packed_keys)} orbital surface tracks into bundle.")
    return packed_keys, uncompressed_bytes


_SURFACE_CACHE: Dict[str, Tuple[Optional[Dict[str, Any]], Optional[Dict[str, Any]]]] = {}


def get_or_extract_surface_mesh(frame_idx: int, ok: str, iso: float, megajson: Dict[str, Any], global_state: Optional[Dict[str, Any]] = None):
    iso_val = round(float(iso), 4)
    cache_key = f"{frame_idx}_{ok}_{iso_val}"
    if cache_key in _SURFACE_CACHE:
        return _SURFACE_CACHE[cache_key]

    root_dir = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
    bundle_zip = os.path.join(root_dir, "out", "bundles", "reaction_cascade_bundle.zip")
    if os.path.exists(bundle_zip):
        try:
            with zipfile.ZipFile(bundle_zip, "r") as zf:
                zip_iso = 0.03
                if "bundle.json" in zf.namelist():
                    try:
                        zip_iso = float(json.loads(zf.read("bundle.json").decode("utf-8")).get("orbital_isovalue") or 0.03)
                    except Exception:
                        pass
                if abs(iso_val - zip_iso) <= 0.003:
                    p_fn, n_fn = f"orbitals/{frame_idx}_{ok}_pos.json", f"orbitals/{frame_idx}_{ok}_neg.json"
                    pos = json.loads(zf.read(p_fn).decode("utf-8")) if p_fn in zf.namelist() else None
                    neg = json.loads(zf.read(n_fn).decode("utf-8")) if n_fn in zf.namelist() else None
                    if pos or neg:
                        _SURFACE_CACHE[cache_key] = (pos, neg)
                        return pos, neg
        except Exception:
            pass

    from .highres_engine import get_cached_highres_cube
    cube_text = get_cached_highres_cube(frame_idx, ok, megajson)
    if not cube_text:
        try:
            from ..step_6_orbitals import _CUBE_CACHE
            cube_text = _CUBE_CACHE.get(f"{frame_idx}_{ok}")
        except Exception:
            pass

    if not cube_text:
        try:
            from ..step_6_orbitals import Step6Orbitals
            res = Step6Orbitals().custom_action({"action_type": "get_orbital_mesh", "frame_index": frame_idx, "track_id": ok, "orbital": ok, "isovalue": iso_val}, global_state or {"step_6_orbitals": {"megajson": megajson}})
            if isinstance(res, dict):
                if res.get("cube"): cube_text = res["cube"]
                elif res.get("pos") or res.get("neg"):
                    _SURFACE_CACHE[cache_key] = (res.get("pos"), res.get("neg"))
                    return res.get("pos"), res.get("neg")
        except Exception:
            pass

    if cube_text:
        try:
            vol, origin, step = CubeConverter.read_cube(cube_text)
            pos = CubeConverter.extract_mesh(vol, origin, step, iso_val)
            neg = CubeConverter.extract_mesh(vol, origin, step, -iso_val)
            _SURFACE_CACHE[cache_key] = (pos, neg)
            return pos, neg
        except Exception:
            pass

    return None, None


def save_state_to_file(global_state: Dict[str, Any]):
    """Persist pipeline state to project state."""
    from ..project_manager import save_step_state
    try:
        for sid, sdata in global_state.items():
            if isinstance(sdata, dict):
                save_step_state(sid, sdata)
    except Exception:
        pass
