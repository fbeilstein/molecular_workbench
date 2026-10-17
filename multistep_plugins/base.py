from dataclasses import dataclass
from typing import Any, Dict, List, Optional


@dataclass
class StepResult:
    success: bool
    data: Optional[Dict[str, Any]] = None
    error_message: Optional[str] = None
    warnings: Optional[List[str]] = None


class PipelineStep:
    """Base interface for Multistep Reaction Pathway Wizard steps."""

    @property
    def step_id(self) -> str:
        """Unique identifier for this pipeline step."""
        raise NotImplementedError

    @property
    def step_name(self) -> str:
        """Human-readable name displayed on the wizard tab."""
        raise NotImplementedError

    def get_ui_schema(self, current_state: Dict[str, Any]) -> Dict[str, Any]:
        """Returns the payload required to render this step's frontend UI."""
        return {}

    def run(self, global_state: Dict[str, Any], user_input: Dict[str, Any]) -> StepResult:
        """Executes the computation or validation for this step."""
        raise NotImplementedError

    def validate(self, step_data: Dict[str, Any]) -> bool:
        """Verifies if the step's output is complete to enable the Next button."""
        return True

    def __eq__(self, other: object) -> bool:
        if isinstance(other, PipelineStep):
            return self.step_id == other.step_id
        return False

    def __hash__(self) -> int:
        return hash(self.step_id)
