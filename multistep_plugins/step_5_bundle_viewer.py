import json
import os
from typing import Any, Dict
from .base import PipelineStep, StepResult
from .step_4_orbitals import Step4Orbitals
from .bundle_utils import (
    recompute_highres_batch,
    get_cached_highres_cube,
    generate_xyz,
    create_maxcomp_bundle_zip,
    unpack_molecule_bundle,
)


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

        return {
            "type": "bundle_viewer_ui",
            "has_megajson": True,
            "megajson": megajson,
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

        # 2. Orbital Cube Retrieval (Checks high-res first, then Step 4 standard)
        if action in ["get_cube", "get_orbital_mesh"]:
            frame_idx = int(user_input.get("frame_index", 0))
            ok = str(user_input.get("orbital", "HOMO"))
            iso = float(user_input.get("isovalue", 0.03))

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

        # 3. Mechanism SVG & Wiring persistence
        if action == "save_mechanism_svg":
            svg_text = user_input.get("svg_text", "")
            megajson["mechanism_svg"] = svg_text
            _save_state_to_file(global_state)
            return {"success": True, "svg_length": len(svg_text)}

        if action == "save_svg_wiring":
            wiring = user_input.get("wiring", {})
            megajson["svg_wiring"] = wiring
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
