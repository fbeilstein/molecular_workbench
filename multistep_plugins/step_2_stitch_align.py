from typing import Any, Dict
from .base import PipelineStep, StepResult
from irc_plugins.ts_engine.trajectory_stitcher import stitch_multistep_trajectories


class Step2StitchAlign(PipelineStep):
    """Step 2: Detect Hero intermediate overlap, align coordinates with Kabsch fit, and stitch."""

    @property
    def step_id(self) -> str:
        return "step_2_stitch_align"

    @property
    def step_name(self) -> str:
        return "2. Align & Stitch"

    def get_ui_schema(self, current_state: Dict[str, Any]) -> Dict[str, Any]:
        step_1_data = current_state.get("step_1_choose_jsons", {})
        trajectories = step_1_data.get("trajectories", [])
        existing_megajson = current_state.get(self.step_id, {}).get("megajson")

        report = []
        junction_previews = []
        if existing_megajson:
            report = existing_megajson.get("stitching_alignment_report", [])
            junction_previews = existing_megajson.get("junction_previews", [])

        return {
            "type": "stitch_align_ui",
            "n_steps": len(trajectories),
            "sequence": step_1_data.get("sequence", []),
            "stitching_report": report,
            "junction_previews": junction_previews,
            "has_megajson": bool(existing_megajson),
            "message": "Inspect boundary 3D alignment, choose Hero intermediate molecules, and verify zero RMSD continuity."
        }

    def _handle_preview_stitch(self, global_state: Dict[str, Any], user_input: Dict[str, Any] = None) -> Dict[str, Any]:
        step_1_data = global_state.get("step_1_choose_jsons", {})
        trajectories = step_1_data.get("trajectories", [])
        if len(trajectories) < 2:
            return {"error": "Need at least 2 trajectories from Step 1."}

        user_input = user_input or {}
        junction_choices = user_input.get("junction_choices", {})
        n_bridge = int(user_input.get("n_bridge_frames", 8))

        try:
            megajson = stitch_multistep_trajectories(trajectories, junction_choices=junction_choices, n_bridge_frames=n_bridge)
            return {
                "success": True,
                "megajson": megajson,
                "report": megajson.get("stitching_alignment_report", []),
                "junction_previews": megajson.get("junction_previews", [])
            }
        except Exception as e:
            return {"error": f"Stitching failed: {e}"}

    def custom_action(self, user_input: Dict[str, Any], global_state: Dict[str, Any] = None) -> Dict[str, Any]:
        global_state = global_state or {}
        action_type = user_input.get("action_type")
        if action_type == "preview_stitch":
            return self._handle_preview_stitch(global_state, user_input)
        return {"error": f"Unknown action: {action_type}"}

    def run(self, global_state: Dict[str, Any], user_input: Dict[str, Any]) -> StepResult:
        step_1_data = global_state.get("step_1_choose_jsons", {})
        trajectories = step_1_data.get("trajectories", [])
        if len(trajectories) < 2:
            return StepResult(success=False, error_message="At least 2 trajectories are required.")

        junction_choices = user_input.get("junction_choices", {})
        n_bridge = int(user_input.get("n_bridge_frames", 8))
        try:
            megajson = stitch_multistep_trajectories(trajectories, junction_choices=junction_choices, n_bridge_frames=n_bridge)
            return StepResult(
                success=True,
                data={
                    "megajson": megajson,
                    "report": megajson.get("stitching_alignment_report", []),
                    "junction_previews": megajson.get("junction_previews", []),
                    "junction_choices": junction_choices,
                    "n_bridge_frames": n_bridge
                }
            )
        except Exception as e:
            return StepResult(success=False, error_message=f"Stitching execution failed: {e}")


    def validate(self, step_data: Dict[str, Any]) -> bool:
        return bool(step_data and step_data.get("megajson"))
