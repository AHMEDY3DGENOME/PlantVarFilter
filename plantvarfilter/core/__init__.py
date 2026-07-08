"""
PlantOmicsGWAS Core Engine
"""

from .checkpoint import CheckpointManager
from .context import PipelineContext
from .results import PipelineResult, StepResult
from .runner import PipelineRunner
from .step import PipelineStep

from .pipeline_factory import (
    create_step,
    create_steps,
    list_available_steps,
    step_is_installed,
)

__version__ = "1.0.0"

__all__ = [
    "CheckpointManager",
    "PipelineContext",
    "PipelineResult",
    "StepResult",
    "PipelineStep",
    "PipelineRunner",
    "create_step",
    "create_steps",
    "list_available_steps",
    "step_is_installed",
]