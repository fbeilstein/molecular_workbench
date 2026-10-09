from typing import Any, Dict
from .base import PipelineStep, StepResult

class Step1Input(PipelineStep):
    @property
    def step_id(self) -> str:
        return "step_1_input"
        
    @property
    def step_name(self) -> str:
        return "1. Reaction Input"

    def get_ui_schema(self, current_state: Dict[str, Any]) -> Dict[str, Any]:
        return {
            "type": "ketcher",
            "message": "Draw the reactants on the left and products on the right, separated by a reaction arrow. (No need to map atoms yet!)"
        }

    def run(self, global_state: Dict[str, Any], user_input: Dict[str, Any]) -> StepResult:
        smiles = user_input.get("smiles", "")
        rxnfile = user_input.get("rxnfile", "")
        
        if not smiles or ">>" not in smiles:
            return StepResult(success=False, error_message="Please draw a valid reaction (with an arrow).")
            
        return StepResult(success=True, data={
            "raw_smiles": smiles,
            "raw_rxnfile": rxnfile
        })
        
    def validate(self, step_data: Dict[str, Any]) -> bool:
        return "raw_smiles" in step_data
