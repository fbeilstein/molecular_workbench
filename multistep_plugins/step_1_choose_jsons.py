import glob
import json
import os
from typing import Any, Dict, List
from .base import PipelineStep, StepResult

TRAJ_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "trajectories")


def _extract_summary(file_path: str) -> Dict[str, Any]:
    """Extract lightweight metadata from an IRC trajectory JSON file for quick preview."""
    try:
        with open(file_path, "r") as f:
            data = json.load(f)
        frames = data.get("frames", [])
        barriers = data.get("barriers", {})
        fragments = data.get("fragments", {})
        return {
            "filename": os.path.basename(file_path),
            "filepath": file_path,
            "reaction_smiles": data.get("reaction_smiles", "Unknown Reaction"),
            "n_atoms": data.get("n_atoms", len(data.get("symbols", []))),
            "symbols": data.get("symbols", []),
            "total_frames": len(frames),
            "ts_frame_index": data.get("ts_frame_index", 0),
            "forward_barrier_kcal": barriers.get("forward_barrier_kcal", 0.0),
            "reverse_barrier_kcal": barriers.get("reverse_barrier_kcal", 0.0),
            "reaction_energy_kcal": barriers.get("reaction_energy_kcal", 0.0),
            "has_fragments": bool(fragments),
            "reactants_meta": fragments.get("reactants", []),
            "products_meta": fragments.get("products", []),
            "valid": len(frames) > 0 and "reaction_smiles" in data
        }
    except Exception as err:
        return {
            "filename": os.path.basename(file_path),
            "filepath": file_path,
            "valid": False,
            "error": str(err)
        }


def scan_available_trajectories() -> List[Dict[str, Any]]:
    """Scan the trajectories directory and workspace root for valid IRC JSON files."""
    os.makedirs(TRAJ_DIR, exist_ok=True)
    paths = glob.glob(os.path.join(TRAJ_DIR, "*.json"))
    root_traj = os.path.join(os.path.dirname(TRAJ_DIR), "irc_trajectory.json")
    if os.path.exists(root_traj) and root_traj not in paths:
        paths.append(root_traj)

    summaries = []
    seen_names = set()
    for p in sorted(paths):
        s = _extract_summary(p)
        if s.get("valid") and s["filename"] not in seen_names:
            summaries.append(s)
            seen_names.add(s["filename"])
    return summaries


class Step1ChooseJSONs(PipelineStep):
    """Step 1: Choose, preview, and arrange IRC trajectories into a multistep sequence."""

    @property
    def step_id(self) -> str:
        return "step_1_choose_jsons"

    @property
    def step_name(self) -> str:
        return "1. Choose Trajectories"

    def get_ui_schema(self, current_state: Dict[str, Any]) -> Dict[str, Any]:
        available = scan_available_trajectories()
        step_data = current_state.get(self.step_id, {})
        selected_sequence = step_data.get("sequence", [])

        return {
            "type": "choose_jsons_ui",
            "available_trajectories": available,
            "selected_sequence": selected_sequence,
            "message": "Select at least 2 elementary IRC trajectories and arrange their reaction sequence."
        }


    def _handle_get_trajectory(self, user_input: Dict[str, Any]) -> Dict[str, Any]:
        """Load full trajectory details for interactive 3D preview."""
        filename = user_input.get("filename")
        if not filename:
            return {"error": "Filename is required"}

        target_path = os.path.join(TRAJ_DIR, filename)
        if not os.path.exists(target_path):
            target_path = os.path.join(os.path.dirname(TRAJ_DIR), filename)

        if not os.path.exists(target_path):
            return {"error": f"Trajectory file '{filename}' not found"}

        try:
            with open(target_path, "r") as f:
                traj_data = json.load(f)
            for fr in traj_data.get("frames", []):
                if "coords" not in fr and "positions" in fr:
                    fr["coords"] = fr["positions"]
                elif "positions" not in fr and "coords" in fr:
                    fr["positions"] = fr["coords"]
                if "energy_kcal" not in fr:
                    fr["energy_kcal"] = fr.get("rel_energy_kcal", 0.0)
            return {"success": True, "trajectory": traj_data, "filename": filename}
        except Exception as e:
            return {"error": f"Failed to load trajectory: {e}"}


    def _handle_upload(self, user_input: Dict[str, Any]) -> Dict[str, Any]:
        """Save an uploaded trajectory JSON directly to the trajectories directory."""
        filename = user_input.get("filename")
        content = user_input.get("content")
        if not filename or not content:
            return {"error": "Filename and JSON content are required"}

        if not filename.endswith(".json"):
            filename += ".json"

        target_path = os.path.join(TRAJ_DIR, os.path.basename(filename))
        try:
            if isinstance(content, str):
                data = json.loads(content)
            else:
                data = content

            if "frames" not in data or "reaction_smiles" not in data:
                return {"error": "Invalid format: JSON must contain 'frames' and 'reaction_smiles'"}

            with open(target_path, "w") as f:
                json.dump(data, f, indent=2)

            summary = _extract_summary(target_path)
            return {"success": True, "summary": summary, "available": scan_available_trajectories()}
        except Exception as e:
            return {"error": f"Upload failed: {e}"}

    def custom_action(self, user_input: Dict[str, Any], global_state: Dict[str, Any] = None) -> Dict[str, Any]:
        action_type = user_input.get("action_type")
        if action_type == "get_trajectory":
            return self._handle_get_trajectory(user_input)
        if action_type == "upload":
            return self._handle_upload(user_input)
        if action_type == "rescan":
            return {"success": True, "available": scan_available_trajectories()}
        return {"error": f"Unknown action: {action_type}"}

    def run(self, global_state: Dict[str, Any], user_input: Dict[str, Any]) -> StepResult:
        sequence_files = user_input.get("sequence", [])
        if not sequence_files or len(sequence_files) < 2:
            return StepResult(
                success=False,
                error_message="Please select at least 2 consecutive reaction trajectories to stitch."
            )

        # Validate that all files exist and load their full trajectory objects
        loaded_trajectories = []
        for fn in sequence_files:
            p = os.path.join(TRAJ_DIR, fn)
            if not os.path.exists(p):
                p = os.path.join(os.path.dirname(TRAJ_DIR), fn)
            if not os.path.exists(p):
                return StepResult(success=False, error_message=f"Missing trajectory file: {fn}")
            try:
                with open(p, "r") as f:
                    loaded_trajectories.append(json.load(f))
            except Exception as e:
                return StepResult(success=False, error_message=f"Error reading {fn}: {e}")

        return StepResult(
            success=True,
            data={
                "sequence": sequence_files,
                "trajectories": loaded_trajectories
            }
        )

    def validate(self, step_data: Dict[str, Any]) -> bool:
        seq = step_data.get("sequence", []) if step_data else []
        return len(seq) >= 2
