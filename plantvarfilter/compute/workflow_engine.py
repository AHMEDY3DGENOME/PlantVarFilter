from pathlib import Path
from typing import Any, Dict, List

from plantvarfilter.core import (
    PipelineContext,
    PipelineRunner,
)

from plantvarfilter.core.pipeline_factory import create_steps

from .config_loader import get_enabled_steps, load_config
from .workflow_registry import WORKFLOW_STEPS, WorkflowStep


class WorkflowEngine:
    def __init__(self, config_path: str):
        self.config_path = config_path
        self.config: Dict[str, Any] = load_config(config_path)

        self.enabled_step_ids: List[str] = get_enabled_steps(self.config)

        self.workflow_steps: List[WorkflowStep] = [
            WORKFLOW_STEPS[step_id]
            for step_id in self.enabled_step_ids
        ]

        self.output_dir = Path(self.config["output"]["dir"])

        self.context = PipelineContext(
            config=self.config,
            config_path=self.config_path,
        )

    def prepare_directories(self) -> None:
        self.context.prepare_directories()

    def print_plan(self) -> None:
        print("\nPlantOmicsGWAS Compute Engine")
        print("=" * 60)
        print(f"Project : {self.context.project_name}")
        print(f"Config  : {self.config_path}")
        print(f"Output  : {self.context.output_dir}")

        print("\nWorkflow:")

        for idx, step in enumerate(self.workflow_steps, start=1):
            print(
                f"{idx:02d}. "
                f"[{step.category}] "
                f"{step.step_id} "
                f"- {step.name}"
            )

        print("=" * 60)

    def dry_run(self) -> None:
        self.prepare_directories()
        self.print_plan()

        print("\nExecution Plan")
        print("=" * 60)

        for idx, step_id in enumerate(self.enabled_step_ids, start=1):
            step = WORKFLOW_STEPS[step_id]

            print(
                f"{idx:02d}. "
                f"{step.step_id}"
            )

            print(f"     Module : {step.module}")

            if step.dependencies:
                print(
                    f"     Depends: "
                    f"{', '.join(step.dependencies)}"
                )

            print()

        print("=" * 60)
        print("Dry-run completed successfully.")

    def run(self) -> None:
        self.prepare_directories()
        self.print_plan()

        pipeline_steps = create_steps(self.enabled_step_ids)

        runner = PipelineRunner(
            context=self.context,
            steps=pipeline_steps,
            stop_on_error=True,
        )

        result = runner.run()

        if not result.success:
            raise RuntimeError(
                "Workflow failed. "
                "See the log above."
            )

        print("\nWorkflow completed successfully.")


def run_workflow(
    config_path: str,
    dry_run: bool = False,
) -> None:

    engine = WorkflowEngine(config_path)

    if dry_run:
        engine.dry_run()
    else:
        engine.run()