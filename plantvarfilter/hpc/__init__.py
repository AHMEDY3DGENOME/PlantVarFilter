"""
PlantOmicsGWAS HPC Layer

Utilities for running PlantOmicsGWAS workflows on High Performance Clusters.

Supported scheduler targets:
- Local
- SLURM
- PBS/Torque
- LSF
"""

from .scheduler import (
    SchedulerInfo,
    SchedulerType,
    detect_scheduler,
    get_scheduler_info,
    is_hpc_scheduler,
    require_scheduler,
    scheduler_from_config,
)

from .job_writer import (
    build_job_script,
    preview_job_script,
    write_job_script,
)

from .array import (
    ArrayJob,
    build_array_job,
    slurm_array_directive,
    pbs_array_directive,
    lsf_array_directive,
    array_summary,
)

from .array_resolver import (
    ArrayTask,
    get_array_task_id,
    resolve_array_task,
    apply_array_task_to_config,
    array_task_summary,
)

from .manager import (
    HPCManager,
    HPCSubmissionResult,
)

__all__ = [
    "SchedulerInfo",
    "SchedulerType",
    "detect_scheduler",
    "get_scheduler_info",
    "is_hpc_scheduler",
    "require_scheduler",
    "scheduler_from_config",
    "build_job_script",
    "preview_job_script",
    "write_job_script",
    "ArrayJob",
    "build_array_job",
    "slurm_array_directive",
    "pbs_array_directive",
    "lsf_array_directive",
    "array_summary",
    "ArrayTask",
    "get_array_task_id",
    "resolve_array_task",
    "apply_array_task_to_config",
    "array_task_summary",
    "HPCManager",
    "HPCSubmissionResult",
]