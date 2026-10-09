from typing import Any, Dict
from .base import PipelineStep, StepResult

class DummyStep1(PipelineStep):
    @property
    def step_id(self) -> str:
        return "dummy_step_1"
        
    @property
    def step_name(self) -> str:
        return "1. Initial Data Input (Dummy)"

    def get_ui_schema(self, current_state: Dict[str, Any]) -> Dict[str, Any]:
        return {
            "type": "dummy_input",
            "message": "Please provide an initial value.",
            "fields": [
                {"name": "test_value", "label": "Test Value", "type": "text"}
            ]
        }

    def run(self, global_state: Dict[str, Any], user_input: Dict[str, Any]) -> StepResult:
        # Simulate processing and needing user review
        val = user_input.get("test_value", "")
        if not val:
            return StepResult(success=False, error_message="test_value is required")
        
        # In a real tool, we might generate a guess here and ask for correction.
        # Here we just pass it along.
        return StepResult(success=True, data={"step_1_result": f"Processed: {val}"})
        
    def validate(self, step_data: Dict[str, Any]) -> bool:
        return "step_1_result" in step_data
