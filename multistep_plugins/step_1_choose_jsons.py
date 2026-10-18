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
        is_stitched = bool(data.get("n_steps", 1) > 1 or data.get("title") == "Master Multistep IRC Trajectory" or data.get("step_boundaries"))
        n_steps = data.get("n_steps", len(data.get("step_boundaries", [])) or 1)
        return {
            "filename": os.path.basename(file_path),
            "filepath": file_path,
            "reaction_smiles": data.get("reaction_smiles") or data.get("title") or os.path.splitext(os.path.basename(file_path))[0],
            "n_atoms": data.get("n_atoms", len(data.get("symbols", []))),
            "symbols": data.get("symbols", []),
            "total_frames": len(frames),
            "n_steps": n_steps,
            "is_stitched": is_stitched,
            "ts_frame_index": data.get("ts_frame_index", 0),
            "forward_barrier_kcal": barriers.get("forward_barrier_kcal", data.get("overall_barrier_kcal", 0.0)),
            "reverse_barrier_kcal": barriers.get("reverse_barrier_kcal", 0.0),
            "reaction_energy_kcal": barriers.get("reaction_energy_kcal", data.get("overall_reaction_energy_kcal", 0.0)),
            "has_fragments": bool(fragments),
            "reactants_meta": fragments.get("reactants", []),
            "products_meta": fragments.get("products", []),
            "valid": len(frames) > 0
        }
    except Exception as err:
        return {"filename": os.path.basename(file_path), "filepath": file_path, "valid": False, "error": str(err)}


def scan_available_trajectories() -> List[Dict[str, Any]]:
    """Scan trajectories/ and out/ directories for valid elementary IRC JSON files."""
    os.makedirs(TRAJ_DIR, exist_ok=True)
    paths = list(glob.glob(os.path.join(TRAJ_DIR, "*.json")))
    out_dir = os.path.join(os.path.dirname(TRAJ_DIR), "out")
    if os.path.exists(out_dir):
        for p in sorted(glob.glob(os.path.join(out_dir, "*.json"))):
            if not p.endswith(("_state.json", "_wiring.json", "optimization.log")):
                if p not in paths:
                    paths.append(p)

    summaries = []
    seen_names = set()
    for p in sorted(paths):
        s = _extract_summary(p)
        fn = s.get("filename")
        if s.get("valid") and fn not in seen_names:
            summaries.append(s)
            seen_names.add(fn)
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
            target_path = os.path.join(os.path.dirname(TRAJ_DIR), "out", filename)
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
                if "positions" not in fr and "coords" in fr:
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

            if "frames" not in data or not isinstance(data["frames"], list) or len(data["frames"]) == 0:
                return {"error": "Invalid format: JSON must contain non-empty 'frames' list"}

            if not data.get("reaction_smiles"):
                data["reaction_smiles"] = data.get("title") or os.path.splitext(os.path.basename(filename))[0]

            with open(target_path, "w") as f:
                json.dump(data, f, indent=2)

            summary = _extract_summary(target_path)
            return {"success": True, "summary": summary, "available": scan_available_trajectories()}
        except Exception as e:
            return {"error": f"Upload failed: {e}"}

    def _handle_delete(self, user_input: Dict[str, Any]) -> Dict[str, Any]:
        """Delete an unwanted trajectory from trajectories/ or out/."""
        filename = user_input.get("filename")
        if not filename:
            return {"error": "Filename is required"}
        p = os.path.join(TRAJ_DIR, os.path.basename(filename))
        if not os.path.exists(p):
            p = os.path.join(os.path.dirname(TRAJ_DIR), "out", os.path.basename(filename))
        if os.path.exists(p):
            try:
                os.remove(p)
                return {"success": True, "message": f"Deleted {filename}", "available": scan_available_trajectories()}
            except Exception as e:
                return {"error": f"Failed to delete {filename}: {e}"}
        return {"error": f"File '{filename}' not found"}

    def custom_action(self, user_input: Dict[str, Any], global_state: Dict[str, Any] = None) -> Dict[str, Any]:
        action_type = user_input.get("action_type")
        if action_type == "get_trajectory":
            return self._handle_get_trajectory(user_input)
        if action_type == "upload":
            return self._handle_upload(user_input)
        if action_type == "delete":
            return self._handle_delete(user_input)
        if action_type == "rescan":
            return {"success": True, "available": scan_available_trajectories()}
        return {"error": f"Unknown action: {action_type}"}

    def run(self, global_state: Dict[str, Any], user_input: Dict[str, Any]) -> StepResult:
        sequence_files = user_input.get("sequence", [])
        if not sequence_files:
            return StepResult(success=False, error_message="Please select trajectory file(s).")

        loaded_trajectories = []
        for fn in sequence_files:
            p = os.path.join(TRAJ_DIR, fn)
            if not os.path.exists(p):
                p = os.path.join(os.path.dirname(TRAJ_DIR), "out", fn)
            if not os.path.exists(p):
                p = os.path.join(os.path.dirname(TRAJ_DIR), fn)
            if not os.path.exists(p):
                return StepResult(success=False, error_message=f"Missing trajectory file: {fn}")
            try:
                with open(p, "r") as f:
                    t_data = json.load(f)
                for fr in t_data.get("frames", []):
                    if "coords" not in fr and "positions" in fr:
                        fr["coords"] = fr["positions"]
                    if "positions" not in fr and "coords" in fr:
                        fr["positions"] = fr["coords"]
                loaded_trajectories.append(t_data)
            except Exception as e:
                return StepResult(success=False, error_message=f"Error reading {fn}: {e}")

        is_single_cascade = False
        megajson = None
        if len(loaded_trajectories) == 1:
            t0 = loaded_trajectories[0]
            if t0.get("n_steps", 1) > 1 or t0.get("title") == "Master Multistep IRC Trajectory" or t0.get("step_boundaries"):
                is_single_cascade = True
                megajson = t0

        if len(loaded_trajectories) < 2 and not is_single_cascade:
            return StepResult(
                success=False,
                error_message="Please select at least 2 elementary trajectories, or 1 pre-stitched cascade trajectory."
            )

        return StepResult(
            success=True,
            data={
                "sequence": sequence_files,
                "trajectories": loaded_trajectories,
                "is_stitched_cascade": is_single_cascade,
                "megajson": megajson
            }
        )

    def validate(self, step_data: Dict[str, Any]) -> bool:
        seq = step_data.get("sequence", []) if step_data else []
        return len(seq) >= 2 or bool(step_data and step_data.get("is_stitched_cascade"))
