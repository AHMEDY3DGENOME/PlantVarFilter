"""
PlantOmicsGWAS Compute Engine - Step Executor

Runs workflow commands and stores logs for each step.
"""

import subprocess
from pathlib import Path
from typing import List, Optional


class StepExecutionError(RuntimeError):
    pass


class StepExecutor:
    def __init__(self, logs_dir: str | Path):
        self.logs_dir = Path(logs_dir)
        self.logs_dir.mkdir(parents=True, exist_ok=True)

    def run(
        self,
        step_id: str,
        command: Optional[List[str]],
        dry_run: bool = False,
    ) -> bool:
        if command is None:
            print(f"[SKIP] {step_id}: no command implemented yet.")
            return True

        print(f"\n[RUN] {step_id}")
        print(" ".join(command))

        if dry_run:
            return True

        log_file = self.logs_dir / f"{step_id}.log"

        with log_file.open("w", encoding="utf-8") as log:
            log.write(f"Step: {step_id}\n")
            log.write(f"Command: {' '.join(command)}\n")
            log.write("=" * 80 + "\n\n")

            process = subprocess.run(
                command,
                stdout=log,
                stderr=subprocess.STDOUT,
                text=True,
            )

        if process.returncode != 0:
            raise StepExecutionError(
                f"Step '{step_id}' failed. Check log file: {log_file}"
            )

        print(f"[DONE] {step_id} | log: {log_file}")
        return True