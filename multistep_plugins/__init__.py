from typing import Dict, List
from .base import PipelineStep
from .step_1_choose_jsons import Step1ChooseJSONs
from .step_2_stitch_align import Step2StitchAlign
from .step_3_energy_stitch import Step3EnergyStitch
from .step_4_mega_viewer import Step4MegaViewer

STEPS: List[PipelineStep] = [
    Step1ChooseJSONs(),
    Step2StitchAlign(),
    Step3EnergyStitch(),
    Step4MegaViewer()
]

STEP_MAP: Dict[str, PipelineStep] = {s.step_id: s for s in STEPS}


def get_all_steps() -> List[PipelineStep]:
    """Returns the ordered list of pipeline steps."""
    return STEPS


def get_step_by_id(step_id: str) -> PipelineStep:
    """Returns a specific pipeline step by ID."""
    return STEP_MAP.get(step_id)
