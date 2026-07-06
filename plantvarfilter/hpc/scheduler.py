"""
PlantOmicsGWAS HPC Scheduler

Detects and describes supported execution backends:

- local
- SLURM
- PBS/Torque
- LSF

This module does not submit jobs directly.
It provides scheduler detection and command definitions for the HPC layer.
"""

from __future__ import annotations

import shutil
from dataclasses import dataclass
from enum import Enum
from typing import Dict, Optional


class SchedulerType(str, Enum):
    LOCAL = "local"
    SLURM = "slurm"
    PBS = "pbs"
    LSF = "lsf"


@dataclass(frozen=True)
class SchedulerInfo:
    scheduler: SchedulerType
    submit_command: Optional[str]
    status_command: Optional[str]
    cancel_command: Optional[str]
    array_supported: bool
    detected: bool


SCHEDULER_COMMANDS: Dict[SchedulerType, Dict[str, Optional[str]]] = {
    SchedulerType.LOCAL: {
        "submit": None,
        "status": None,
        "cancel": None,
    },
    SchedulerType.SLURM: {
        "submit": "sbatch",
        "status": "squeue",
        "cancel": "scancel",
    },
    SchedulerType.PBS: {
        "submit": "qsub",
        "status": "qstat",
        "cancel": "qdel",
    },
    SchedulerType.LSF: {
        "submit": "bsub",
        "status": "bjobs",
        "cancel": "bkill",
    },
}


def _command_exists(command: Optional[str]) -> bool:
    if not command:
        return False
    return shutil.which(command) is not None


def detect_scheduler(preferred: Optional[str] = None) -> SchedulerInfo:
    if preferred:
        preferred_type = SchedulerType(preferred.lower())
        return get_scheduler_info(preferred_type)

    for scheduler in (
        SchedulerType.SLURM,
        SchedulerType.PBS,
        SchedulerType.LSF,
    ):
        info = get_scheduler_info(scheduler)
        if info.detected:
            return info

    return get_scheduler_info(SchedulerType.LOCAL)


def get_scheduler_info(scheduler: SchedulerType) -> SchedulerInfo:
    commands = SCHEDULER_COMMANDS[scheduler]

    submit_command = commands["submit"]
    status_command = commands["status"]
    cancel_command = commands["cancel"]

    detected = (
        scheduler == SchedulerType.LOCAL
        or _command_exists(submit_command)
    )

    return SchedulerInfo(
        scheduler=scheduler,
        submit_command=submit_command,
        status_command=status_command,
        cancel_command=cancel_command,
        array_supported=scheduler in {
            SchedulerType.SLURM,
            SchedulerType.PBS,
            SchedulerType.LSF,
        },
        detected=detected,
    )


def scheduler_from_config(config: dict) -> SchedulerInfo:
    hpc_cfg = config.get("hpc", {}) or {}
    scheduler_name = hpc_cfg.get("scheduler")

    if not scheduler_name:
        return detect_scheduler()

    return detect_scheduler(str(scheduler_name))


def is_hpc_scheduler(scheduler: SchedulerType) -> bool:
    return scheduler in {
        SchedulerType.SLURM,
        SchedulerType.PBS,
        SchedulerType.LSF,
    }


def require_scheduler(scheduler: SchedulerType) -> SchedulerInfo:
    info = get_scheduler_info(scheduler)

    if not info.detected:
        cmd = SCHEDULER_COMMANDS[scheduler]["submit"]
        raise RuntimeError(
            f"Requested scheduler '{scheduler.value}' was not detected. "
            f"Missing submit command: {cmd}"
        )

    return info