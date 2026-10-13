from typing import Any, Dict
from .base import PipelineStep, StepResult
from irc_plugins.ts_engine.energy_balancer import balance_pathway_energies


class Step3EnergyStitch(PipelineStep):
    """Step 3: Collect general stoichiometry, calculate isolated fragment energies, and balance PES via Hess's Law."""

    @property
    def step_id(self) -> str:
        return "step_3_energy_stitch"

    @property
    def step_name(self) -> str:
        return "3. Hess's Law PES"

    def get_ui_schema(self, current_state: Dict[str, Any]) -> Dict[str, Any]:
        step_1_data = current_state.get("step_1_choose_jsons", {})
        step_2_data = current_state.get("step_2_stitch_align", {})
        trajectories = step_1_data.get("trajectories", [])
        megajson = step_2_data.get("megajson")

        if not megajson:
            return {
                "type": "energy_stitch_ui",
                "has_megajson": False,
                "message": "Complete Step 2 alignment first."
            }

        existing_data = current_state.get(self.step_id, {})
        balanced_megajson = existing_data.get("megajson")
        saved_roles = existing_data.get("fragment_roles", {})
        saved_solvent = existing_data.get("solvent", "water")

        if not balanced_megajson or "energy_profile" not in balanced_megajson:
            balanced_megajson = balance_pathway_energies(
                megajson, trajectories, recalculate_with_xtb=True, solvent=saved_solvent, fragment_roles=saved_roles
            )

        return {
            "type": "energy_stitch_ui",
            "has_megajson": True,
            "megajson": balanced_megajson,
            "stoichiometry_summary": balanced_megajson.get("stoichiometry_summary", {}),
            "energy_profile": balanced_megajson.get("energy_profile", {}),
            "step_boundaries": balanced_megajson.get("step_boundaries", []),
            "solvent": balanced_megajson.get("solvent", saved_solvent),
            "fragment_roles": saved_roles,
            "message": "Specify fragment roles (Hero Core, Incoming Reagent, Leaving Group Excluded) and recalculate PES."
        }

    def _handle_recalculate(self, global_state: Dict[str, Any], user_input: Dict[str, Any] = None) -> Dict[str, Any]:
        user_input = user_input or {}
        use_xtb = bool(user_input.get("recalculate_with_xtb", True))
        solvent = user_input.get("solvent", "water")
        fragment_roles = user_input.get("fragment_roles", {})

        step_1_data = global_state.get("step_1_choose_jsons", {})
        step_2_data = global_state.get("step_2_stitch_align", {})
        trajectories = step_1_data.get("trajectories", [])
        megajson = step_2_data.get("megajson")

        if not megajson or not trajectories:
            return {"error": "Missing stitched trajectory or input steps."}

        try:
            balanced = balance_pathway_energies(
                megajson, trajectories, recalculate_with_xtb=use_xtb, solvent=solvent, fragment_roles=fragment_roles
            )
            return {
                "success": True,
                "megajson": balanced,
                "energy_profile": balanced.get("energy_profile", {}),
                "stoichiometry_summary": balanced.get("stoichiometry_summary", {}),
                "solvent": solvent,
                "fragment_roles": fragment_roles
            }
        except Exception as e:
            return {"error": f"Hess's Law energy calculation failed: {e}"}

    def custom_action(self, user_input: Dict[str, Any], global_state: Dict[str, Any] = None) -> Dict[str, Any]:
        global_state = global_state or {}
        action = user_input.get("action_type")
        if action == "recalculate_energies":
            return self._handle_recalculate(global_state, user_input)
        return {"error": f"Unknown action: {action}"}

    def run(self, global_state: Dict[str, Any], user_input: Dict[str, Any]) -> StepResult:
        step_1_data = global_state.get("step_1_choose_jsons", {})
        step_2_data = global_state.get("step_2_stitch_align", {})
        trajectories = step_1_data.get("trajectories", [])
        megajson = step_2_data.get("megajson")

        if not megajson:
            return StepResult(success=False, error_message="Step 2 trajectory missing.")

        use_xtb = bool(user_input.get("recalculate_with_xtb", True))
        solvent = user_input.get("solvent", "water")
        fragment_roles = user_input.get("fragment_roles", {})

        try:
            balanced = balance_pathway_energies(
                megajson, trajectories, recalculate_with_xtb=use_xtb, solvent=solvent, fragment_roles=fragment_roles
            )
            return StepResult(
                success=True,
                data={
                    "megajson": balanced,
                    "energy_profile": balanced.get("energy_profile", {}),
                    "stoichiometry_summary": balanced.get("stoichiometry_summary", {}),
                    "solvent": solvent,
                    "fragment_roles": fragment_roles
                }
            )
        except Exception as e:
            return StepResult(success=False, error_message=f"Energy balancing execution failed: {e}")

    def validate(self, step_data: Dict[str, Any]) -> bool:
        return bool(step_data and step_data.get("megajson") and step_data.get("energy_profile"))
