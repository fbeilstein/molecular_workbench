import os
from typing import Any, Dict
from .base import PipelineStep, StepResult


class Step3MegaViewer(PipelineStep):
    """Step 3: Interactive 3D visualization and analysis of the stitched reaction pathway."""

    @property
    def step_id(self) -> str:
        return "step_3_mega_viewer"

    @property
    def step_name(self) -> str:
        return "3. Pathway Viewer"

    def get_ui_schema(self, current_state: Dict[str, Any]) -> Dict[str, Any]:
        step_2_data = current_state.get("step_2_stitch_align", {})
        megajson = step_2_data.get("megajson")

        if not megajson:
            return {
                "type": "mega_viewer_ui",
                "has_megajson": False,
                "message": "No stitched trajectory found. Complete Step 2 first."
            }

        frames = megajson.get("frames", [])
        step_boundaries = megajson.get("step_boundaries", [])
        overall_barrier = megajson.get("overall_barrier_kcal", 0.0)
        overall_reaction_energy = megajson.get("overall_reaction_energy_kcal", 0.0)

        return {
            "type": "mega_viewer_ui",
            "has_megajson": True,
            "megajson": megajson,
            "total_frames": len(frames),
            "step_boundaries": step_boundaries,
            "overall_barrier_kcal": overall_barrier,
            "overall_reaction_energy_kcal": overall_reaction_energy,
            "symbols": megajson.get("symbols", []),
            "message": "Continuous multistep reaction pathway ready for playback and export."
        }

    def _generate_xyz(self, megajson: Dict[str, Any]) -> str:
        """Convert MegaJSON frames to standard multi-frame XYZ format."""
        symbols = megajson.get("symbols", [])
        n_atoms = len(symbols)
        frames = megajson.get("frames", [])
        lines = []

        for f_idx, frame in enumerate(frames):
            coords = frame.get("coords", [])
            energy = frame.get("energy_kcal", 0.0)
            step_idx = frame.get("step_index", 0)
            lines.append(f"{n_atoms}")
            lines.append(f"Frame {f_idx} | Step {step_idx} | Energy: {energy:.2f} kcal/mol")
            for sym, pos in zip(symbols, coords):
                lines.append(f"{sym:<2} {pos[0]:12.6f} {pos[1]:12.6f} {pos[2]:12.6f}")

        return "\n".join(lines) + "\n"

    def custom_action(self, user_input: Dict[str, Any], global_state: Dict[str, Any] = None) -> Dict[str, Any]:
        global_state = global_state or {}
        action_type = user_input.get("action_type")
        megajson = global_state.get("step_2_stitch_align", {}).get("megajson")

        if not megajson:
            return {"error": "Stitched pathway data is missing"}

        if action_type == "get_xyz":
            xyz_str = self._generate_xyz(megajson)
            return {"success": True, "xyz": xyz_str, "filename": "multistep_pathway.xyz"}

        if action_type == "get_megajson":
            return {"success": True, "megajson": megajson, "filename": "multistep_pathway.json"}

        return {"error": f"Unknown action: {action_type}"}

    def run(self, global_state: Dict[str, Any], user_input: Dict[str, Any]) -> StepResult:
        megajson = global_state.get("step_2_stitch_align", {}).get("megajson")
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
