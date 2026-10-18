import glob
import json
import os
from typing import Any, Dict, List
from .base import PipelineStep, StepResult
from .ts_engine.trajectory_stitcher import stitch_multistep_trajectories

ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
OUT_DIR = os.path.join(ROOT_DIR, 'out')
TRAJECTORIES_DIR = os.path.join(ROOT_DIR, 'trajectories')
PROLONGED_PATH = os.path.join(TRAJECTORIES_DIR, 'prolonged_pathway.json')


class Step5Cascade(PipelineStep):
    """Step 5: Inspect and tune multi-step reaction cascades with rotamer alignment & bridge frames."""

    @property
    def step_id(self) -> str:
        return "step_5_cascade"

    @property
    def step_name(self) -> str:
        return "5. Multistep Cascade"

    def _discover_trajectories(self) -> List[Dict[str, Any]]:
        cand_files = []
        for p in [os.path.join(OUT_DIR, "irc_trajectory_prev.json"), os.path.join(OUT_DIR, "irc_trajectory.json")]:
            if os.path.exists(p) and p not in cand_files:
                cand_files.append(p)
        if os.path.exists(TRAJECTORIES_DIR):
            for p in sorted(glob.glob(os.path.join(TRAJECTORIES_DIR, "*.json"))):
                if p not in cand_files and not p.endswith("prolonged_pathway.json"):
                    cand_files.append(p)

        results = []
        seen_smiles = set()
        for p in cand_files:
            try:
                with open(p, "r") as f:
                    tj = json.load(f)
                if not tj.get("frames"):
                    continue
                if tj.get("n_steps", 1) > 1 or tj.get("title") == "Master Multistep IRC Trajectory":
                    continue
                smi = tj.get("reaction_smiles") or os.path.basename(p)
                if smi in seen_smiles:
                    continue
                seen_smiles.add(smi)
                results.append({
                    "path": p,
                    "filename": os.path.basename(p),
                    "reaction_smiles": tj.get("reaction_smiles", ""),
                    "n_frames": len(tj.get("frames", [])),
                    "n_atoms": len(tj.get("symbols", [])),
                    "n_steps": tj.get("n_steps", 1),
                    "barriers": tj.get("barriers", {}),
                    "data": tj
                })
            except Exception:
                continue
        return results

    def get_ui_schema(self, current_state: Dict[str, Any]) -> Dict[str, Any]:
        available = self._discover_trajectories()
        step_state = current_state.get(self.step_id, {})
        megajson = step_state.get("megajson")

        if not megajson and os.path.exists(PROLONGED_PATH):
            try:
                with open(PROLONGED_PATH, "r") as f:
                    megajson = json.load(f)
            except Exception:
                pass

        if not megajson and len(available) >= 2:
            try:
                megajson = stitch_multistep_trajectories([a["data"] for a in available[:2]], n_bridge_frames=8)
                with open(PROLONGED_PATH, "w") as f:
                    json.dump(megajson, f, indent=2)
            except Exception:
                megajson = None

        stitching_report = megajson.get("stitching_alignment_report", []) if megajson else []
        step_boundaries = megajson.get("step_boundaries", []) if megajson else []

        traj_meta = [{
            "filename": a["filename"],
            "path": a["path"],
            "reaction_smiles": a["reaction_smiles"],
            "n_frames": a["n_frames"],
            "n_atoms": a["n_atoms"]
        } for a in available]

        return {
            "type": "cascade_ui",
            "trajectories": traj_meta,
            "megajson": megajson,
            "has_cascade": bool(megajson),
            "stitching_report": stitching_report,
            "step_boundaries": step_boundaries,
            "n_bridge_frames": step_state.get("n_bridge_frames", 8),
            "junction_choices": step_state.get("junction_choices", {}),
            "message": "Inspect aligned reaction sequence, select rotamer/scaffold modes, and tune bridge interpolation."
        }

    def _handle_stitch_cascade(self, global_state: Dict[str, Any], user_input: Dict[str, Any]) -> Dict[str, Any]:
        user_input = user_input or {}
        junction_choices = user_input.get("junction_choices", {})
        n_bridge = int(user_input.get("n_bridge_frames", 8))
        ordered_files = user_input.get("ordered_files", [])

        available = {a["filename"]: a["data"] for a in self._discover_trajectories()}
        if not available:
            return {"error": "No trajectories found to stitch."}

        to_stitch = []
        if ordered_files:
            for fn in ordered_files:
                if fn in available:
                    to_stitch.append(available[fn])
        if len(to_stitch) < 2:
            to_stitch = [a["data"] for a in self._discover_trajectories()][:2]

        if len(to_stitch) < 2:
            return {"error": "Need at least 2 trajectories to stitch a cascade."}

        try:
            megajson = stitch_multistep_trajectories(to_stitch, junction_choices=junction_choices, n_bridge_frames=n_bridge)
            os.makedirs(TRAJECTORIES_DIR, exist_ok=True)
            with open(PROLONGED_PATH, "w") as f:
                json.dump(megajson, f, indent=2)

            # Also update out/irc_trajectory.json master
            try:
                out_traj = os.path.join(OUT_DIR, "irc_trajectory.json")
                with open(out_traj, "w") as f:
                    json.dump(megajson, f, indent=2)
            except Exception:
                pass

            return {
                "success": True,
                "megajson": megajson,
                "stitching_report": megajson.get("stitching_alignment_report", []),
                "step_boundaries": megajson.get("step_boundaries", []),
                "message": f"Successfully stitched {megajson.get('n_steps', len(to_stitch))} steps ({megajson.get('total_frames', len(megajson['frames']))} frames)!"
            }
        except Exception as e:
            return {"error": f"Stitching failed: {e}"}

    def custom_action(self, user_input: Dict[str, Any], global_state: Dict[str, Any] = None) -> Dict[str, Any]:
        global_state = global_state or {}
        action = user_input.get("action_type")
        if action == "stitch_cascade":
            return self._handle_stitch_cascade(global_state, user_input)
        return {"error": f"Unknown action: {action}"}

    def run(self, global_state: Dict[str, Any], user_input: Dict[str, Any]) -> StepResult:
        megajson = user_input.get("megajson")
        if not megajson and os.path.exists(PROLONGED_PATH):
            try:
                with open(PROLONGED_PATH, "r") as f:
                    megajson = json.load(f)
            except Exception:
                pass
        return StepResult(success=True, data={"megajson": megajson})

    def validate(self, step_data: Dict[str, Any]) -> bool:
        return bool(step_data and step_data.get("megajson"))
