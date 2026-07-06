"""
PlantOmicsGWAS Core

Pipeline Result Objects
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional
import time


@dataclass
class StepResult:
    """
    Represents the result of one workflow step.
    """

    step_id: str
    success: bool = False
    skipped: bool = False
    return_code: int = 0
    start_time: float = field(default_factory=time.time)
    end_time: float = 0.0
    runtime: float = 0.0
    command: Optional[List[str]] = None
    outputs: Dict[str, str] = field(default_factory=dict)
    warnings: List[str] = field(default_factory=list)
    errors: List[str] = field(default_factory=list)
    log_file: Optional[Path] = None
    message: str = ""

    def finish(
        self,
        success: bool,
        return_code: int = 0,
        message: str = "",
    ) -> None:
        self.success = success
        self.return_code = return_code
        self.message = message
        self.end_time = time.time()
        self.runtime = self.end_time - self.start_time

    def skip(self, message: str = "Skipped.") -> None:
        self.success = True
        self.skipped = True
        self.return_code = 0
        self.message = message
        self.end_time = time.time()
        self.runtime = self.end_time - self.start_time

    def add_output(self, key: str, value: str) -> None:
        self.outputs[key] = str(value)

    def add_warning(self, message: str) -> None:
        self.warnings.append(message)

    def add_error(self, message: str) -> None:
        self.errors.append(message)

    @property
    def status(self) -> str:
        if self.skipped:
            return "SKIPPED"
        if self.success:
            return "SUCCESS"
        return "FAILED"


@dataclass
class PipelineResult:
    """
    Represents the result of the complete pipeline.
    """

    project_name: str
    start_time: float = field(default_factory=time.time)
    end_time: float = 0.0
    runtime: float = 0.0
    steps: List[StepResult] = field(default_factory=list)
    success: bool = True

    def add_step(self, result: StepResult) -> None:
        self.steps.append(result)

        if not result.success and not result.skipped:
            self.success = False

    def finish(self) -> None:
        self.end_time = time.time()
        self.runtime = self.end_time - self.start_time

    @property
    def successful_steps(self) -> int:
        return sum(1 for step in self.steps if step.success and not step.skipped)

    @property
    def failed_steps(self) -> int:
        return sum(1 for step in self.steps if not step.success and not step.skipped)

    @property
    def skipped_steps(self) -> int:
        return sum(1 for step in self.steps if step.skipped)

    def summary(self) -> Dict:
        return {
            "project": self.project_name,
            "success": self.success,
            "runtime_seconds": round(self.runtime, 2),
            "total_steps": len(self.steps),
            "successful_steps": self.successful_steps,
            "failed_steps": self.failed_steps,
            "skipped_steps": self.skipped_steps,
        }