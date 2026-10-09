from .base import PipelineStep, StepResult

_PLUGINS = []

def register_plugin(plugin: PipelineStep):
    _PLUGINS.append(plugin)

def get_plugins():
    return _PLUGINS

# Import and register your plugins here
from .step_1_input import Step1Input
from .step_2_mapping import Step2Mapping

register_plugin(Step1Input())
register_plugin(Step2Mapping())
