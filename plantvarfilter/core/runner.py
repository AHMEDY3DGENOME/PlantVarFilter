"""
PlantOmicsGWAS Core

Pipeline Runner
"""

from __future__ import annotations

from typing import Iterable, List

from .checkpoint import CheckpointManager
from .context import PipelineContext
from .results import PipelineResult, StepResult
from .step import PipelineStep


class PipelineRunner:
    """
    Executes the workflow while supporting resume/checkpoints.
    """

    def __init__(
        self,
        context: PipelineContext,
        steps: Iterable[PipelineStep],
        stop_on_error: bool = True,
    ):
        self.context = context
        self.steps: List[PipelineStep] = list(steps)
        self.stop_on_error = stop_on_error
        self.result = PipelineResult(project_name=context.project_name)

        resume_enabled = bool(
            context.config.get("compute", {}).get("resume", True)
        )

        self.checkpoints = CheckpointManager(
            output_dir=context.output_dir,
            enabled=resume_enabled,
        )

    def _update_context_outputs(self, step_result: StepResult) -> None:
        if not step_result.outputs:
            return

        self.context.config.setdefault("input", {})

        for key, value in step_result.outputs.items():
            self.context.config["input"][key] = value

    def _print_step_result(self, step_result: StepResult) -> None:
        print(
            f"Status : {step_result.status} | "
            f"Runtime: {step_result.runtime:.2f}s"
        )

        if step_result.message:
            print(f"Message: {step_result.message}")

        if step_result.errors:
            for error in step_result.errors:
                print(f"Error  : {error}")

        if step_result.outputs:
            print("Outputs:")
            for key, value in step_result.outputs.items():
                print(f"  - {key}: {value}")

    def run(self) -> PipelineResult:
        self.context.prepare_directories()

        print("\nPlantOmicsGWAS Core Runner")
        print("=" * 50)
        print(f"Project: {self.context.project_name}")
        print(f"Output : {self.context.output_dir}")
        print(f"Steps  : {len(self.steps)}")
        print("=" * 50)

        for index, step in enumerate(self.steps, start=1):
            print(f"\n[{index}/{len(self.steps)}] {step.step_id} - {step.name}")

            if self.checkpoints.exists(step.step_id):
                outputs = self.checkpoints.get_outputs(step.step_id)

                step_result = StepResult(step_id=step.step_id)
                step_result.outputs.update(outputs)
                step_result.skip(
                    message=f"Skipped because checkpoint exists for '{step.step_id}'."
                )

                if outputs:
                    self.context.config.setdefault("input", {})
                    self.context.config["input"].update(outputs)

                self.result.add_step(step_result)
                self._print_step_result(step_result)
                continue

            step_result: StepResult = step.run(self.context)
            self.result.add_step(step_result)

            if step_result.success:
                self._update_context_outputs(step_result)

                self.checkpoints.save(
                    step_id=step.step_id,
                    status=step_result.status,
                    message=step_result.message,
                    outputs=step_result.outputs,
                    runtime=step_result.runtime,
                )

            self._print_step_result(step_result)

            if not step_result.success and self.stop_on_error:
                print("\nPipeline stopped because a step failed.")
                break

        self.result.finish()

        print("\nPipeline Summary")
        print("=" * 50)
        print(f"Success          : {self.result.success}")
        print(f"Total steps      : {len(self.result.steps)}")
        print(f"Successful steps : {self.result.successful_steps}")
        print(f"Failed steps     : {self.result.failed_steps}")
        print(f"Skipped steps    : {self.result.skipped_steps}")
        print(f"Runtime seconds  : {self.result.runtime:.2f}")
        print("=" * 50)

        return self.result