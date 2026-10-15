import io
import json
import os
import zipfile
from typing import Any, Dict, List
from .base import PipelineStep, StepResult
from irc_plugins.ts_engine import orbital_engine


class Step5BundleViewer(PipelineStep):
    """Step 5: Master Synchronized 3D Pathway + Orbitals + PES Viewer and Bundle Exporter."""

    @property
    def step_id(self) -> str:
        return "step_5_bundle_viewer"

    @property
    def step_name(self) -> str:
        return "5. Bundle & Viewer"

    def _get_megajson(self, current_state: Dict[str, Any]) -> Dict[str, Any]:
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
            "message": "Synchronized 3D pathway, orbitals, and Free Energy profile ready for analysis and export."
        }

    def _generate_xyz(self, megajson: Dict[str, Any]) -> str:
        symbols = megajson.get("symbols", [])
        n_atoms = len(symbols)
        frames = megajson.get("frames", [])
        lines = []
        for f_idx, frame in enumerate(frames):
            coords = frame.get("coords", frame.get("positions", []))
            energy = frame.get("energy_kcal", 0.0)
            step_idx = frame.get("step_index", 0)
            lines.append(f"{n_atoms}")
            lines.append(f"Frame {f_idx} | Step {step_idx} | Free Energy: {energy:.2f} kcal/mol")
            for sym, pos in zip(symbols, coords):
                lines.append(f"{sym:<2} {pos[0]:12.6f} {pos[1]:12.6f} {pos[2]:12.6f}")
        return "\n".join(lines) + "\n"

    def _build_bundle_manifest(self, megajson: Dict[str, Any], global_state: Dict[str, Any]) -> Dict[str, Any]:
        boundaries = megajson.get("step_boundaries", [])
        title = " → ".join(b.get("step_name", f"Step {i+1}") for i, b in enumerate(boundaries)) or "Multistep Reaction"
        rxn_smi = " -> ".join(b.get("reaction_smiles", "") for b in boundaries)
        return {
            "title": title,
            "reaction_smiles": rxn_smi,
            "total_frames": len(megajson.get("frames", [])),
            "overall_barrier_kcal": megajson.get("overall_barrier_kcal", 0.0),
            "overall_reaction_energy_kcal": megajson.get("overall_reaction_energy_kcal", 0.0),
            "solvent": megajson.get("solvent", "water"),
            "step_boundaries": boundaries,
            "symbols": megajson.get("symbols", []),
            "energy_profile": megajson.get("energy_profile", {}),
            "orbital_schedule": megajson.get("orbital_schedule", {}),
            "files": {
                "trajectory_xyz": "multistep_pathway.xyz",
                "trajectory_json": "multistep_pathway.json",
                "orbitals_directory": "orbitals/"
            }
        }

    def _create_bundle_zip(self, megajson: Dict[str, Any], global_state: Dict[str, Any]) -> bytes:
        zip_buf = io.BytesIO()
        manifest = self._build_bundle_manifest(megajson, global_state)
        xyz_content = self._generate_xyz(megajson)
        json_content = json.dumps(megajson, indent=2)

        with zipfile.ZipFile(zip_buf, "w", zipfile.ZIP_DEFLATED) as zf:
            zf.writestr("bundle.json", json.dumps(manifest, indent=2))
            zf.writestr("multistep_pathway.xyz", xyz_content)
            zf.writestr("multistep_pathway.json", json_content)

            # Store pre-generated landmark orbital meshes
            landmark_meshes = megajson.get("landmark_orbital_meshes", {})
            for f_str, orbs in landmark_meshes.items():
                for orb_name, mesh in orbs.items():
                    if mesh.get("pos"):
                        zf.writestr(f"orbitals/frame_{int(f_str):03d}_{orb_name}_pos.json", json.dumps(mesh["pos"]))
                    if mesh.get("neg"):
                        zf.writestr(f"orbitals/frame_{int(f_str):03d}_{orb_name}_neg.json", json.dumps(mesh["neg"]))

        zip_buf.seek(0)
        return zip_buf.getvalue()

    def custom_action(self, user_input: Dict[str, Any], global_state: Dict[str, Any] = None) -> Dict[str, Any]:
        global_state = global_state or {}
        megajson = self._get_megajson(global_state)
        if not megajson:
            return {"error": "Pathway data is missing"}

        action = user_input.get("action_type")
        if action == "get_xyz":
            return {"success": True, "xyz": self._generate_xyz(megajson), "filename": "multistep_pathway.xyz"}

        if action == "get_megajson":
            return {"success": True, "megajson": megajson, "filename": "multistep_pathway.json"}

        if action == "export_bundle":
            bundle_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "static_multistep", "bundles")
            os.makedirs(bundle_dir, exist_ok=True)
            bundle_filename = "multistep_reaction_bundle.zip"
            bundle_path = os.path.join(bundle_dir, bundle_filename)
            zip_bytes = self._create_bundle_zip(megajson, global_state)
            with open(bundle_path, "wb") as f:
                f.write(zip_bytes)
            return {
                "success": True,
                "download_url": f"/static_multistep/bundles/{bundle_filename}",
                "filename": bundle_filename,
                "size_bytes": len(zip_bytes)
            }

        if action in ["get_cube", "get_orbital_mesh", "get_frame_orbitals", "batch_precompute_cubes"]:
            from .step_4_orbitals import Step4Orbitals
            s4 = Step4Orbitals()
            return s4.custom_action(user_input, global_state)

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
