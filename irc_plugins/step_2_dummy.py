from typing import Any, Dict
from .base import PipelineStep, StepResult

class DummyStep2(PipelineStep):
    @property
    def step_id(self) -> str:
        return "dummy_step_2"
        
    @property
    def step_name(self) -> str:
        return "2. Review & Correction (Dummy)"

    def get_ui_schema(self, current_state: Dict[str, Any]) -> Dict[str, Any]:
        prev_result = current_state.get("dummy_step_1", {}).get("step_1_result", "None")
        return {
            "type": "dummy_review",
            "message": f"Previous step produced: {prev_result}. Approve or correct it.",
            "fields": [
                {"name": "correction", "label": "Correction (optional)", "type": "text"}
            ]
        }

    def run(self, global_state: Dict[str, Any], user_input: Dict[str, Any]) -> StepResult:
        prev_result = global_state.get("dummy_step_1", {}).get("step_1_result", "None")
        correction = user_input.get("correction", "")
        
        if correction:
            final_val = f"Corrected to: {correction}"
        else:
            final_val = f"Approved: {prev_result}"
            
        return StepResult(success=True, data={"step_2_final": final_val})
        
    def validate(self, step_data: Dict[str, Any]) -> bool:
        return "step_2_final" in step_data
