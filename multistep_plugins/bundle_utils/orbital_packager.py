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
    """Extract and pack orbital surface meshes for all frames across all configured timeline tracks."""
    timeline = megajson.get("orbital_timeline") or (global_state.get("step_4_orbitals", {}).get("orbital_timeline") if global_state else None) or {}
    raw_tracks = timeline.get("tracks", [])
    frames = megajson.get("frames", [])
    total_frames = len(frames)

    active_tracks = []
    for trk in raw_tracks:
        if isinstance(trk, dict):
            ok = trk.get("orbital_key")
            ck = trk.get("chem_key")
            clips = trk.get("clips", [])
        else:
            ok, ck, clips = str(trk), None, []
        if not ok or str(ok).startswith("track_esp"):
            continue
        active_tracks.append((ok, ck, clips))

    if not active_tracks and frames:
        for tree_cat in megajson.get("orbitals_tree", []):
            for it in tree_cat.get("items", []):
                k = it.get("key") or it.get("orbital_key")
                if k and k not in [a[0] for a in active_tracks]:
                    active_tracks.append((k, it.get("chem_key"), []))

    needed_targets: List[Tuple[int, str, Optional[str]]] = []
    for ok, ck, clips in active_tracks:
        if clips:
            f_set = set()
            for c in clips:
                s = max(0, int(c.get("start", 0)))
                e = min(total_frames - 1, int(c.get("end", total_frames - 1)))
                f_set.update(range(s, e + 1))
        else:
            f_set = set(range(total_frames))
        for fi in sorted(f_set):
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

    ws_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
    search_zips = [
        os.path.join(ws_root, "bundles", "tetrahedral_intermediate_CO.zip"),
        os.path.join(ws_root, "out", "tetrahedral_intermediate_CO.zip"),
        os.path.join(ws_root, "out", "bundles", "multistep_reaction_bundle_maxcomp.zip"),
    ]
    reusable_meshes: Dict[str, bytes] = {}
    for zp in search_zips:
        if os.path.exists(zp):
            try:
                with zipfile.ZipFile(zp, "r") as ez:
                    for name in ez.namelist():
                        if name.startswith("orbitals/") and name.endswith(("_pos.json", "_neg.json")):
                            if name not in reusable_meshes:
                                reusable_meshes[name] = ez.read(name)
            except Exception:
                pass

    missing_targets: List[Tuple[int, str, Optional[str]]] = []
    for fi, ok, ck in needed_targets:
        key = f"{fi}_{ok}"
        if key in _FAILED_MESH_KEYS:
            continue
        p_name, n_name = f"orbitals/{key}_pos.json", f"orbitals/{key}_neg.json"
        if p_name in reusable_meshes or n_name in reusable_meshes:
            _pack_mesh_data(key, reusable_meshes.get(p_name), reusable_meshes.get(n_name))
        else:
            missing_targets.append((fi, ok, ck))

    sources: Dict[str, Any] = {}
    legacy_d = os.path.abspath(os.path.join(ws_root, "trajectories", "cubes"))
    for d in (CUBE_DIR, legacy_d):
        if os.path.exists(d):
            sources.update({f[:-5]: os.path.join(d, f) for f in os.listdir(d) if f.endswith(".cube")})

    for k, text in _HIGHRES_CUBE_CACHE.items():
        if "_" in k and not k.startswith("res"):
            parts = k.split("_")
            clean_k = f"{parts[1]}_{'_'.join(parts[2:])}" if (len(parts) > 2 and parts[0].isdigit() and parts[1].isdigit()) else k
            sources.setdefault(clean_k, text)

    for k, val in megajson.get("highres_cubes", {}).items():
        if k not in sources:
            if isinstance(val, str) and "\n" in val:
                sources[k] = val
            elif isinstance(val, str) and val.endswith(".cube"):
                cp = os.path.abspath(os.path.join(ws_root, val))
                if os.path.exists(cp):
                    sources[k] = cp

    try:
        from ..step_4_orbitals import _CUBE_CACHE
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
                p_b = json.dumps(pos, separators=(',', ':')).encode('utf-8') if pos else None
                n_b = json.dumps(neg, separators=(',', ':')).encode('utf-8') if neg else None
                _pack_mesh_data(key, p_b, n_b)
            except Exception as e:
                print(f"[OrbitalPackager] Error extracting cube for {key}: {e}")
                still_missing.append((fi, ok, ck))
        else:
            still_missing.append((fi, ok, ck))

    if still_missing and total_frames > 0:
        print(f"[OrbitalPackager] Computing {len(still_missing)} trajectory orbital meshes on demand...")
        try:
            from ..step_4_orbitals import Step4Orbitals
            s4 = Step4Orbitals()
            for fi, ok, ck in still_missing:
                key = f"{fi}_{ok}"
                try:
                    res = s4.custom_action({
                        "action_type": "get_orbital_mesh",
                        "frame_index": fi,
                        "orbital": ok,
                        "chem_key": ck,
                        "isovalue": isovalue
                    }, global_state or {"step_4_orbitals": {"megajson": megajson}})
                    cube_text = res.get("cube")
                    if cube_text:
                        vol, origin, step = CubeConverter.read_cube(cube_text)
                        pos = CubeConverter.extract_mesh(vol, origin, step, isovalue)
                        neg = CubeConverter.extract_mesh(vol, origin, step, -isovalue)
                        p_b = json.dumps(pos, separators=(',', ':')).encode('utf-8') if pos else None
                        n_b = json.dumps(neg, separators=(',', ':')).encode('utf-8') if neg else None
                        _pack_mesh_data(key, p_b, n_b)
                    else:
                        _FAILED_MESH_KEYS.add(key)
                except Exception as ex:
                    _FAILED_MESH_KEYS.add(key)
                    print(f"[OrbitalPackager] Warning frame {fi} {ok}: {ex}")
        except Exception as e:
            print(f"[OrbitalPackager] Error during MO computation: {e}")


    for key, src in sources.items():
        if key not in packed_keys:
            try:
                vol, origin, step = CubeConverter.read_cube(src)
                pos = CubeConverter.extract_mesh(vol, origin, step, isovalue)
                neg = CubeConverter.extract_mesh(vol, origin, step, -isovalue)
                p_b = json.dumps(pos, separators=(',', ':')).encode('utf-8') if pos else None
                n_b = json.dumps(neg, separators=(',', ':')).encode('utf-8') if neg else None
                _pack_mesh_data(key, p_b, n_b)
            except Exception:
                pass

    print(f"[OrbitalPackager] Successfully packed {len(packed_keys)} orbital surface tracks into bundle.")
    return packed_keys, uncompressed_bytes


_SURFACE_CACHE: Dict[str, Tuple[Optional[Dict[str, Any]], Optional[Dict[str, Any]]]] = {}


def get_or_extract_surface_mesh(frame_idx: int, ok: str, iso: float, megajson: Dict[str, Any]):
    cache_key = f"{frame_idx}_{ok}_{iso}"
    if cache_key in _SURFACE_CACHE:
        return _SURFACE_CACHE[cache_key]

    root_dir = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
    for sub in ("out/bundles", "static_multistep/bundles"):
        bundle_zip = os.path.join(root_dir, sub, "multistep_reaction_bundle_maxcomp.zip")
        if os.path.exists(bundle_zip):
            try:
                with zipfile.ZipFile(bundle_zip, "r") as zf:
                    p_fn, n_fn = f"orbitals/{frame_idx}_{ok}_pos.json", f"orbitals/{frame_idx}_{ok}_neg.json"
                    pos = json.loads(zf.read(p_fn).decode("utf-8")) if p_fn in zf.namelist() else None
                    neg = json.loads(zf.read(n_fn).decode("utf-8")) if n_fn in zf.namelist() else None
                    if pos or neg:
                        _SURFACE_CACHE[cache_key] = (pos, neg)
                        return pos, neg
            except Exception:
                pass

    cube_text = get_cached_highres_cube(frame_idx, ok, megajson)
    if not cube_text:
        try:
            from ..step_4_orbitals import _CUBE_CACHE
            cube_text = _CUBE_CACHE.get(f"{frame_idx}_{ok}")
        except Exception:
            pass

    if cube_text:
        try:
            vol, origin, step = CubeConverter.read_cube(cube_text)
            pos = CubeConverter.extract_mesh(vol, origin, step, iso)
            neg = CubeConverter.extract_mesh(vol, origin, step, -iso)
            _SURFACE_CACHE[cache_key] = (pos, neg)
            return pos, neg
        except Exception as e:
            print(f"[Step5] Warning: mesh extraction failed for {frame_idx}_{ok}: {e}")

    return None, None


def save_state_to_file(global_state: Dict[str, Any]):
    """Persist pipeline state to out/multistep_state.json."""
    root_dir = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
    out_dir = os.path.join(root_dir, "out")
    os.makedirs(out_dir, exist_ok=True)
    try:
        with open(os.path.join(out_dir, "multistep_state.json"), "w") as f:
            json.dump(global_state, f, indent=2)
    except Exception as e:
        print(f"[Step5] Warning: could not persist state to file: {e}")
