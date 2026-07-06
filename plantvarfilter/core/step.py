"""
PlantOmicsGWAS Core

Base Pipeline Step
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Optional

from .context import PipelineContext
from .results import StepResult


class PipelineStep(ABC):
    """
    Base class for every analysis step.

    Every pipeline (QC, Alignment, GWAS, Prediction, ...)
    should inherit from this class.
    """

    step_id: str = "unknown"

    name: str = "Unknown Step"

    description: str = ""

    version: str = "1.0"

    def __init__(self):

        self.result = StepResult(self.step_id)

    # --------------------------------------------------

    def before_run(self, context: PipelineContext) -> None:
        """
        Optional hook before execution.
        """
        pass

    # --------------------------------------------------

    @abstractmethod
    def execute(self, context: PipelineContext) -> None:
        """
        Main implementation.

        Every child class MUST implement this.
        """
        raise NotImplementedError

    # --------------------------------------------------

    def after_run(self, context: PipelineContext) -> None:
        """
        Optional hook after execution.
        """
        pass

    # --------------------------------------------------

    def run(self, context: PipelineContext) -> StepResult:

        try:

            self.before_run(context)

            self.execute(context)

            self.after_run(context)

            self.result.finish(
                success=True,
                message="Completed successfully."
            )

        except Exception as e:

            self.result.add_error(str(e))

            self.result.finish(
                success=False,
                return_code=1,
                message=str(e)
            )

        return self.result

    # --------------------------------------------------

    def skip(self, reason: str) -> StepResult:

        self.result.skipped = True

        self.result.finish(
            success=True,
            message=reason
        )

        return self.result

    # --------------------------------------------------

    def add_output(self, name: str, value: str):

        self.result.add_output(name, value)

    # --------------------------------------------------

    def warning(self, message: str):

        self.result.add_warning(message)

    # --------------------------------------------------

    def error(self, message: str):

        self.result.add_error(message)

    # --------------------------------------------------

    def __repr__(self):

        return f"<PipelineStep {self.step_id}>"