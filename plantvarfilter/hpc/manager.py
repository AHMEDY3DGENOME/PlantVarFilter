"""
PlantOmicsGWAS HPC Manager

High-level HPC operations:

- write job script
- submit job
- submit dependency workflow
- show cluster/scheduler info
- check job status
- cancel jobs

This module is scheduler-aware and supports:
- SLURM
- PBS/Torque
- LSF
- local fallback
"""

from __future__ import annotations

import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional

from plantvarfilter.hpc.job_writer import write_job_script
from plantvarfilter.hpc.scheduler import (
    SchedulerInfo,
    SchedulerType,
    scheduler_from_config,
)


@dataclass
class HPCSubmissionResult:
    scheduler: str
    job_script: str
    submitted: bool
    job_id: Optional[str] = None
    command: Optional[List[str]] = None
    stdout: str = ""
    return_code: int = 0


class HPCManager:
    def __init__(self, config: Dict, config_path: str):
        self.config = config
        self.config_path = str(Path(config_path).expanduser().resolve())
        self.scheduler_info: SchedulerInfo = scheduler_from_config(config)

    # ------------------------------------------------------------------
    # Core helpers
    # ------------------------------------------------------------------

    def write_job(self, output_path: Optional[str] = None) -> str:
        return write_job_script(
            config=self.config,
            config_path=self.config_path,
            output_path=output_path,
        )

    def submit(
        self,
        output_path: Optional[str] = None,
        no_submit: bool = False,
    ) -> HPCSubmissionResult:
        job_script = self.write_job(output_path=output_path)

        if no_submit:
            return HPCSubmissionResult(
                scheduler=self.scheduler_info.scheduler.value,
                job_script=job_script,
                submitted=False,
                command=None,
                stdout="Submission skipped because no_submit=True.",
                return_code=0,
            )

        if self.scheduler_info.scheduler == SchedulerType.LOCAL:
            return HPCSubmissionResult(
                scheduler="local",
                job_script=job_script,
                submitted=False,
                command=None,
                stdout="Local scheduler selected. Run the script manually.",
                return_code=0,
            )

        if not self.scheduler_info.submit_command:
            raise RuntimeError("No submit command available for selected scheduler.")

        cmd = [self.scheduler_info.submit_command, job_script]

        proc = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )

        stdout = proc.stdout or ""
        job_id = self._parse_job_id(stdout)

        return HPCSubmissionResult(
            scheduler=self.scheduler_info.scheduler.value,
            job_script=job_script,
            submitted=proc.returncode == 0,
            job_id=job_id,
            command=cmd,
            stdout=stdout,
            return_code=proc.returncode,
        )

    # ------------------------------------------------------------------
    # Dependency submission
    # ------------------------------------------------------------------

    def submit_dependency_chain(
        self,
        job_scripts: List[str],
        no_submit: bool = False,
    ) -> List[HPCSubmissionResult]:
        """
        Submit multiple job scripts as a dependency chain.

        SLURM:
            sbatch job1.sh
            sbatch --dependency=afterok:<jobid1> job2.sh

        PBS:
            qsub job1.sh
            qsub -W depend=afterok:<jobid1> job2.sh

        LSF:
            bsub < job1.sh
            bsub -w "done(<jobid1>)" < job2.sh

        Note:
        This function expects already-written job script paths.
        """

        results: List[HPCSubmissionResult] = []

        if not job_scripts:
            return results

        previous_job_id: Optional[str] = None

        for script in job_scripts:
            script_path = str(Path(script).expanduser().resolve())

            if no_submit:
                results.append(
                    HPCSubmissionResult(
                        scheduler=self.scheduler_info.scheduler.value,
                        job_script=script_path,
                        submitted=False,
                        job_id=None,
                        command=None,
                        stdout="Dependency submission skipped because no_submit=True.",
                        return_code=0,
                    )
                )
                continue

            cmd = self._dependency_submit_command(
                script_path=script_path,
                previous_job_id=previous_job_id,
            )

            proc = subprocess.run(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
            )

            stdout = proc.stdout or ""
            job_id = self._parse_job_id(stdout)

            result = HPCSubmissionResult(
                scheduler=self.scheduler_info.scheduler.value,
                job_script=script_path,
                submitted=proc.returncode == 0,
                job_id=job_id,
                command=cmd,
                stdout=stdout,
                return_code=proc.returncode,
            )

            results.append(result)

            if proc.returncode != 0:
                break

            previous_job_id = job_id

        return results

    def _dependency_submit_command(
        self,
        script_path: str,
        previous_job_id: Optional[str],
    ) -> List[str]:
        scheduler = self.scheduler_info.scheduler

        if scheduler == SchedulerType.SLURM:
            if previous_job_id:
                return [
                    "sbatch",
                    f"--dependency=afterok:{previous_job_id}",
                    script_path,
                ]
            return ["sbatch", script_path]

        if scheduler == SchedulerType.PBS:
            if previous_job_id:
                return [
                    "qsub",
                    "-W",
                    f"depend=afterok:{previous_job_id}",
                    script_path,
                ]
            return ["qsub", script_path]

        if scheduler == SchedulerType.LSF:
            if previous_job_id:
                return [
                    "bsub",
                    "-w",
                    f"done({previous_job_id})",
                    "<",
                    script_path,
                ]
            return ["bsub", "<", script_path]

        raise RuntimeError("Dependency jobs are only supported on SLURM, PBS, and LSF.")

    # ------------------------------------------------------------------
    # Scheduler utilities
    # ------------------------------------------------------------------

    def cluster_info(self) -> str:
        info = self.scheduler_info

        lines = [
            "PlantOmicsGWAS HPC Cluster Information",
            "=" * 45,
            f"Scheduler      : {info.scheduler.value}",
            f"Detected       : {info.detected}",
            f"Submit command : {info.submit_command}",
            f"Status command : {info.status_command}",
            f"Cancel command : {info.cancel_command}",
            f"Array support  : {info.array_supported}",
            "=" * 45,
        ]

        return "\n".join(lines)

    def status(self, job_id: Optional[str] = None) -> str:
        if not self.scheduler_info.status_command:
            return "No status command available for local scheduler."

        cmd = [self.scheduler_info.status_command]

        if job_id:
            cmd.append(str(job_id))

        proc = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )

        return proc.stdout or ""

    def cancel(self, job_id: str) -> str:
        if not self.scheduler_info.cancel_command:
            return "No cancel command available for local scheduler."

        cmd = [self.scheduler_info.cancel_command, str(job_id)]

        proc = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )

        return proc.stdout or ""

    # ------------------------------------------------------------------
    # Job ID parsing
    # ------------------------------------------------------------------

    def _parse_job_id(self, stdout: str) -> Optional[str]:
        text = stdout.strip()

        if not text:
            return None

        scheduler = self.scheduler_info.scheduler

        # SLURM: Submitted batch job 123456
        if scheduler == SchedulerType.SLURM:
            for token in text.split():
                if token.isdigit():
                    return token

        # PBS: usually prints job id directly, e.g. 12345.server
        if scheduler == SchedulerType.PBS:
            return text.split()[0]

        # LSF: Job <12345> is submitted...
        if scheduler == SchedulerType.LSF:
            if "<" in text and ">" in text:
                return text.split("<", 1)[1].split(">", 1)[0]

        return None