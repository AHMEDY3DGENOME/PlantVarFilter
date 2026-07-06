"""
PlantOmicsGWAS HPC Job Writer
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Dict

from plantvarfilter.hpc.scheduler import SchedulerType, scheduler_from_config
from plantvarfilter.hpc.templates import (
    build_local_script,
    build_lsf_script,
    build_pbs_script,
    build_slurm_script,
)


def _default_job_filename(scheduler: SchedulerType) -> str:
    if scheduler == SchedulerType.SLURM:
        return "plantomicsgwas.slurm"
    if scheduler == SchedulerType.PBS:
        return "plantomicsgwas.pbs"
    if scheduler == SchedulerType.LSF:
        return "plantomicsgwas.lsf"
    return "plantomicsgwas.sh"


def build_job_script(config: Dict, config_path: str) -> str:
    scheduler_info = scheduler_from_config(config)
    scheduler = scheduler_info.scheduler

    if scheduler == SchedulerType.SLURM:
        return build_slurm_script(config, config_path)

    if scheduler == SchedulerType.PBS:
        return build_pbs_script(config, config_path)

    if scheduler == SchedulerType.LSF:
        return build_lsf_script(config, config_path)

    return build_local_script(config_path)


def write_job_script(
    config: Dict,
    config_path: str,
    output_path: str | None = None,
) -> str:
    scheduler = scheduler_from_config(config).scheduler

    output_dir = Path(config.get("output", {}).get("dir", "."))
    output_dir.mkdir(parents=True, exist_ok=True)

    if output_path is None:
        output_path = str(output_dir / _default_job_filename(scheduler))

    script = build_job_script(config, config_path)

    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(script, encoding="utf-8")
    os.chmod(path, 0o755)

    return str(path)


def preview_job_script(config: Dict, config_path: str) -> str:
    return build_job_script(config, config_path)