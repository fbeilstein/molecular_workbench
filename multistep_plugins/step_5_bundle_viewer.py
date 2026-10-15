import json
import os
import sys
import zipfile
from typing import Any, Dict, Optional, Tuple
from .base import PipelineStep, StepResult
from .step_4_orbitals import Step4Orbitals
from .bundle_utils import (
    recompute_highres_batch,
    get_cached_highres_cube,
    generate_xyz,
    create_maxcomp_bundle_zip,
    unpack_molecule_bundle,
    purge_cube_files,
)

_SURFACE_CACHE: Dict[str, Tuple[Optional[Dict[str, Any]], Optional[Dict[str, Any]]]] = {}


def _get_or_extract_surface_mesh(frame_idx: int, ok: str, iso: float, megajson: Dict[str, Any]):
    cache_key = f"{frame_idx}_{ok}_{iso}"
    if cache_key in _SURFACE_CACHE:
        return _SURFACE_CACHE[cache_key]

    # 1. Try reading precomputed surface from maxcomp bundle zip
    bundle_zip = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "static_multistep", "bundles", "multistep_reaction_bundle_maxcomp.zip")
    if os.path.exists(bundle_zip):
        try:
            with zipfile.ZipFile(bundle_zip, "r") as zf:
                p_fn = f"orbitals/{frame_idx}_{ok}_pos.json"
                n_fn = f"orbitals/{frame_idx}_{ok}_neg.json"
                pos = json.loads(zf.read(p_fn).decode("utf-8")) if p_fn in zf.namelist() else None
                neg = json.loads(zf.read(n_fn).decode("utf-8")) if n_fn in zf.namelist() else None
                if pos or neg:
                    _SURFACE_CACHE[cache_key] = (pos, neg)
                    return pos, neg
        except Exception:
            pass

    # 2. Try extracting from available cube
    cube_text = get_cached_highres_cube(frame_idx, ok, megajson)
    if not cube_text:
        try:
            from .step_4_orbitals import _CUBE_CACHE
            cube_text = _CUBE_CACHE.get(f"{frame_idx}_{ok}")
        except Exception:
            pass

    if cube_text:
        try:
            TOOLS_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "tools"))
            if TOOLS_DIR not in sys.path:
                sys.path.insert(0, TOOLS_DIR)
            from cube_converter import CubeConverter
            vol, origin, step = CubeConverter.read_cube(cube_text)
            pos = CubeConverter.extract_mesh(vol, origin, step, iso)
            neg = CubeConverter.extract_mesh(vol, origin, step, -iso)
            _SURFACE_CACHE[cache_key] = (pos, neg)
            return pos, neg
        except Exception as e:
            print(f"[Step5] Warning: mesh extraction failed for {frame_idx}_{ok}: {e}")

    return None, None


def _save_state_to_file(global_state: Dict[str, Any]):
    """Persist pipeline state to multistep_state.json."""
    state_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "multistep_state.json")
    try:
        with open(state_file, "w") as f:
            json.dump(global_state, f, indent=2)
    except Exception as e:
        print(f"[Step5] Warning: could not persist state to file: {e}")


class Step5BundleViewer(PipelineStep):
    """Step 5: Master Synchronized 3D Pathway + Orbitals + SVG Mechanism + Max-Compression Bundle Repacker."""

    @property
    def step_id(self) -> str:
        return "step_5_bundle_viewer"

    @property
    def step_name(self) -> str:
        return "5. Bundle & Viewer"

    def _get_megajson(self, current_state: Dict[str, Any]) -> Dict[str, Any]:
        s5 = current_state.get("step_5_bundle_viewer", {})
        if s5.get("megajson"):
            return s5["megajson"]
        s4 = current_state.get("step_4_orbitals", {})
        mj = s4.get("megajson")
        if not mj:
            s3 = current_state.get("step_3_energy_stitch", {})
            mj = s3.get("megajson")
        if not mj:
            mj = current_state.get("step_2_stitch_align", {}).get("megajson")
        return mj or {}

    def get_ui_schema(self, current_state: Dict[str, Any]) -> Dict[str, Any]:
        megajson = self._get_megajson(current_state)
        if not megajson or not megajson.get("frames"):
            return {
                "type": "bundle_viewer_ui",
                "has_megajson": False,
                "message": "No stitched pathway found. Complete earlier steps first."
            }


        frames = megajson.get("frames", [])
        timeline = megajson.get("orbital_timeline") or current_state.get("step_4_orbitals", {}).get("orbital_timeline", {})
        highres_cubes = megajson.get("highres_cubes", {})
        molecule_bundles = megajson.get("molecule_bundles", {})
        molecule_summary = {
            m_id: {
                "filename": data.get("filename", "bundle.zip"),
                "orbitals": data.get("orbitals", []),
                "count": len(data.get("orbitals", []))
            }
            for m_id, data in molecule_bundles.items()
        }

        ui_megajson = dict(megajson)
        if "highres_cubes" in ui_megajson:
            ui_megajson["highres_cubes"] = {
                k: (v if (isinstance(v, str) and "\n" not in v) else True)
                for k, v in ui_megajson["highres_cubes"].items()
            }

        return {
            "type": "bundle_viewer_ui",
            "has_megajson": True,
            "megajson": ui_megajson,
            "total_frames": len(frames),
            "step_boundaries": megajson.get("step_boundaries", []),
            "overall_barrier_kcal": megajson.get("overall_barrier_kcal", 0.0),
            "overall_reaction_energy_kcal": megajson.get("overall_reaction_energy_kcal", 0.0),
            "energy_profile": megajson.get("energy_profile", {}),
            "symbols": megajson.get("symbols", []),
            "orbital_timeline": timeline,
            "orbital_schedule": megajson.get("orbital_schedule", {}),
            "landmarks": current_state.get("step_4_orbitals", {}).get("landmarks", []),
            "mechanism_svg": megajson.get("mechanism_svg", ""),
            "svg_wiring": megajson.get("svg_wiring", {}),
            "orbital_resolution": megajson.get("orbital_resolution", 32),
            "highres_count": len(highres_cubes),
            "molecule_bundles": molecule_summary,
            "message": "Synchronized 3D pathway, orbitals, interactive reaction mechanism SVG, and PES profile ready."
        }

    def custom_action(self, user_input: Dict[str, Any], global_state: Dict[str, Any] = None) -> Dict[str, Any]:
        global_state = global_state or {}
        megajson = self._get_megajson(global_state)
        if not megajson:
            return {"error": "Pathway data is missing"}

        action = user_input.get("action_type")

        # 1. High-resolution orbital recomputation
        if action == "recompute_highres_cubes":
            resolution = int(user_input.get("resolution", 50))
            tracks = user_input.get("selected_tracks", [])
            sf = int(user_input.get("start_frame", 0))
            ef = int(user_input.get("end_frame", sf + 15))
            iso = float(user_input.get("isovalue", 0.03))

            res = recompute_highres_batch(megajson, tracks, sf, ef, resolution, iso)
            _save_state_to_file(global_state)
            return res

        # 2. Orbital Surface Mesh / Cube Retrieval (Prefers precomputed lightweight meshes)
        if action in ["get_cube", "get_orbital_mesh"]:
            frame_idx = int(user_input.get("frame_index", 0))
            ok = str(user_input.get("orbital", "HOMO"))
            iso = float(user_input.get("isovalue", 0.03))

            pos, neg = _get_or_extract_surface_mesh(frame_idx, ok, iso, megajson)
            if pos or neg:
                return {
                    "success": True,
                    "frame_index": frame_idx,
                    "orbital": ok,
                    "isovalue": iso,
                    "is_mesh": True,
                    "pos": pos,
                    "neg": neg,
                    "coords": megajson["frames"][frame_idx]["coords"]
                }

            highres_cube = get_cached_highres_cube(frame_idx, ok, megajson)
            if highres_cube:
                return {
                    "success": True,
                    "frame_index": frame_idx,
                    "orbital": ok,
                    "isovalue": iso,
                    "is_highres": True,
                    "cube": highres_cube,
                    "filename": f"{ok}_f{frame_idx}_highres.cube",
                    "coords": megajson["frames"][frame_idx]["coords"]
                }

            s4 = Step4Orbitals()
            return s4.custom_action(user_input, global_state)

        if action in ["get_cached_cubes", "get_cached_meshes"]:
            frames = [int(x) for x in user_input.get("frames", [])]
            tracks = user_input.get("tracks", [])
            iso = float(user_input.get("isovalue", 0.03))
            meshes = {}
            for fi in frames:
                for trk in tracks:
                    ok = trk.get("orbital_key") if isinstance(trk, dict) else str(trk)
                    if ok == "track_esp_charge":
                        continue
                    pos, neg = _get_or_extract_surface_mesh(fi, ok, iso, megajson)
                    if pos or neg:
                        meshes[f"{fi}_{ok}"] = {"pos": pos, "neg": neg}

            if meshes:
                return {"success": True, "meshes": meshes, "is_mesh": True}

            s4 = Step4Orbitals()
            res = s4.custom_action(user_input, global_state)
            if isinstance(res, dict) and res.get("success"):
                cubes = res.setdefault("cubes", {})
                for fi in frames:
                    for trk in tracks:
                        ok = trk.get("orbital_key") if isinstance(trk, dict) else str(trk)
                        if f"{fi}_{ok}" not in cubes:
                            hrc = get_cached_highres_cube(fi, ok, megajson)
                            if hrc:
                                cubes[f"{fi}_{ok}"] = hrc
            return res

        # 3. Mechanism SVG & Wiring persistence
        if action == "save_mechanism_svg":
            svg_text = user_input.get("svg_text", "")
            megajson["mechanism_svg"] = svg_text
            for s_key in ("step_5_bundle_viewer", "step_4_orbitals", "step_3_energy_stitch"):
                if s_key in global_state and isinstance(global_state[s_key], dict):
                    s_mj = global_state[s_key].get("megajson")
                    if isinstance(s_mj, dict):
                        s_mj["mechanism_svg"] = svg_text
            _save_state_to_file(global_state)
            return {"success": True, "svg_length": len(svg_text)}

        if action == "save_svg_wiring":
            wiring = user_input.get("wiring", {})
            megajson["svg_wiring"] = wiring
            for s_key in ("step_5_bundle_viewer", "step_4_orbitals", "step_3_energy_stitch"):
                if s_key in global_state and isinstance(global_state[s_key], dict):
                    s_mj = global_state[s_key].get("megajson")
                    if isinstance(s_mj, dict):
                        s_mj["svg_wiring"] = wiring
            _save_state_to_file(global_state)
            return {"success": True, "wired_count": len(wiring)}

        # 4. Molecular Workbench ZIP bundle upload & orbital catalog
        if action == "upload_molecule_bundle":
            m_id = str(user_input.get("molecule_id", "molecule"))
            filename = str(user_input.get("filename", "bundle.zip"))
            data_b64 = user_input.get("data_b64", "")

            try:
                orbitals_found, files_dict = unpack_molecule_bundle(data_b64)
            except Exception as e:
                return {"success": False, "error": f"Failed to parse ZIP bundle: {e}"}

            bundles = megajson.setdefault("molecule_bundles", {})
            bundles[m_id] = {
                "filename": filename,
                "orbitals": orbitals_found,
                "files": files_dict,
                "data_b64": data_b64
            }
            _save_state_to_file(global_state)
            return {
                "success": True,
                "molecule_id": m_id,
                "orbitals": orbitals_found,
                "count": len(orbitals_found)
            }

        if action == "get_molecule_bundle":
            m_id = str(user_input.get("molecule_id", ""))
            b = megajson.get("molecule_bundles", {}).get(m_id)
            if not b:
                return {"success": False, "error": f"No bundle found for molecule '{m_id}'"}
            return {
                "success": True,
                "molecule_id": m_id,
                "filename": b.get("filename"),
                "orbitals": b.get("orbitals", []),
                "files": b.get("files", {})
            }

        # 5. Maximal Compression Bundle Repacker
        if action in ["repack_bundle", "export_bundle"]:
            bundle_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "static_multistep", "bundles")
            os.makedirs(bundle_dir, exist_ok=True)
            bundle_filename = "multistep_reaction_bundle_maxcomp.zip"
            bundle_path = os.path.join(bundle_dir, bundle_filename)

            zip_bytes, uncompressed_bytes = create_maxcomp_bundle_zip(megajson, global_state)
            with open(bundle_path, "wb") as f:
                f.write(zip_bytes)

            ratio = (1.0 - (len(zip_bytes) / max(1, uncompressed_bytes))) * 100.0
            return {
                "success": True,
                "download_url": f"/static_multistep/bundles/{bundle_filename}",
                "filename": bundle_filename,
                "size_bytes": len(zip_bytes),
                "uncompressed_bytes": uncompressed_bytes,
                "compression_ratio_pct": round(ratio, 1)
            }

        if action == "purge_cubes":
            purged = purge_cube_files()
            megajson.pop("highres_cubes", None)
            _save_state_to_file(global_state)
            return {"success": True, "purged": purged}

        if action == "get_xyz":
            return {"success": True, "xyz": generate_xyz(megajson), "filename": "multistep_pathway.xyz"}

        if action == "get_megajson":
            return {"success": True, "megajson": megajson, "filename": "multistep_pathway.json"}

        return {"error": f"Unknown action: {action}"}

    def run(self, global_state: Dict[str, Any], user_input: Dict[str, Any]) -> StepResult:
        megajson = self._get_megajson(global_state)
        if not megajson:
            return StepResult(success=False, error_message="No stitched pathway available.")

        return StepResult(
            success=True,
            data={
                "completed": True,
                "megajson": megajson
            }
        )

    def validate(self, step_data: Dict[str, Any]) -> bool:
        return bool(step_data and step_data.get("completed"))
