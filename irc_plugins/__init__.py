from .base import PipelineStep, StepResult

_PLUGINS = []

def register_plugin(plugin: PipelineStep):
    _PLUGINS.append(plugin)

def get_plugins():
    return _PLUGINS

# Import and register your plugins here
from .step_1_dummy import DummyStep1
from .step_2_dummy import DummyStep2

register_plugin(DummyStep1())
register_plugin(DummyStep2())
