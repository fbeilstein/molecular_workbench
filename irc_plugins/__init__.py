from .base import PipelineStep, StepResult

_PLUGINS = []

def register_plugin(plugin: PipelineStep):
    _PLUGINS.append(plugin)

def get_plugins():
    return _PLUGINS

# Import and register your plugins here
from .step_1_input import Step1Input
from .step_2_mapping import Step2Mapping
from .step_3_ts import Step3TS
from .step_4_irc import Step4IRC
from .step_5_cascade import Step5Cascade
from .step_6_orbitals import Step6Orbitals
from .step_7_bundle import Step7Bundle

register_plugin(Step1Input())
register_plugin(Step2Mapping())
register_plugin(Step3TS())
register_plugin(Step4IRC())
register_plugin(Step5Cascade())
register_plugin(Step6Orbitals())
register_plugin(Step7Bundle())



