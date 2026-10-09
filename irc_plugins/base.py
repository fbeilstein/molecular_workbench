from typing import Any, Dict, Optional
from dataclasses import dataclass

@dataclass
class StepResult:
    success: bool
    data: Optional[Dict[str, Any]] = None
    error_message: Optional[str] = None
    warnings: Optional[list[str]] = None

class PipelineStep:
    """Base interface that every IRC pipeline phase must implement."""
    
    @property
    def step_id(self) -> str:
        """Unique identifier (e.g., 'phase_1_mapping')."""
        raise NotImplementedError
        
    @property
    def step_name(self) -> str:
        """Human-readable name (e.g., '1. 2D Atom Mapping')."""
        raise NotImplementedError

    def get_ui_schema(self, current_state: Dict[str, Any]) -> Dict[str, Any]:
        """
        Returns the data required to render this step's UI on the frontend.
        """
        return {}

    def run(self, global_state: Dict[str, Any], user_input: Dict[str, Any]) -> StepResult:
        """
        Executes the computational payload for this step.
        """
        raise NotImplementedError
        
    def validate(self, step_data: Dict[str, Any]) -> bool:
        """
        Fast check to verify if the data produced by this step is valid.
        """
        return True
