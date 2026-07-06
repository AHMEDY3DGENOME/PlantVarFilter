"""
PlantOmicsGWAS Core Engine
"""

from .checkpoint import CheckpointManager
from .context import PipelineContext
from .results import PipelineResult, StepResult
from .runner import PipelineRunner
from .step import PipelineStep

from .pipeline_factory import (
    PIPELINE_STEP_REGISTRY,
    create_step,
    create_steps,
    list_available_steps,
)

__version__ = "1.0.0"

__all__ = [
    "CheckpointManager",
    "PipelineContext",
    "PipelineResult",
    "StepResult",
    "PipelineStep",
    "PipelineRunner",
    "PIPELINE_STEP_REGISTRY",
    "create_step",
    "create_steps",
    "list_available_steps",
]