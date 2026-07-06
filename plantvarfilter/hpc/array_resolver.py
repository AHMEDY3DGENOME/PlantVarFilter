"""
PlantOmicsGWAS HPC Array Resolver

Resolves the current scheduler array task into a concrete work item.

Supported environment variables:
- SLURM_ARRAY_TASK_ID
- PBS_ARRAY_INDEX
- PBS_ARRAYID
- LSB_JOBINDEX
- PLANTOMICS_ARRAY_TASK_ID
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Dict, List, Optional

from plantvarfilter.hpc.array import ArrayJob, build_array_job


@dataclass(frozen=True)
class ArrayTask:
    enabled: bool
    mode: Optional[str]
    task_id: Optional[int]
    items: List[str]
    chunk_size: int

    @property
    def is_active(self) -> bool:
        return self.enabled and self.task_id is not None

    @property
    def first_item(self) -> Optional[str]:
        return self.items[0] if self.items else None


def get_array_task_id() -> Optional[int]:
    env_names = [
        "PLANTOMICS_ARRAY_TASK_ID",
        "SLURM_ARRAY_TASK_ID",
        "PBS_ARRAY_INDEX",
        "PBS_ARRAYID",
        "LSB_JOBINDEX",
    ]

    for name in env_names:
        value = os.environ.get(name)

        if value is None or str(value).strip() == "":
            continue

        try:
            task_id = int(value)
        except ValueError:
            continue

        # LSF arrays are usually 1-based.
        if name == "LSB_JOBINDEX":
            task_id -= 1

        return task_id

    return None


def resolve_array_task(config: Dict) -> ArrayTask:
    array: Optional[ArrayJob] = build_array_job(config)

    if array is None:
        return ArrayTask(
            enabled=False,
            mode=None,
            task_id=None,
            items=[],
            chunk_size=1,
        )

    task_id = get_array_task_id()

    if task_id is None:
        return ArrayTask(
            enabled=True,
            mode=array.mode,
            task_id=None,
            items=[],
            chunk_size=array.chunk_size,
        )

    start = task_id * array.chunk_size
    end = start + array.chunk_size

    selected_items = array.items[start:end]

    return ArrayTask(
        enabled=True,
        mode=array.mode,
        task_id=task_id,
        items=selected_items,
        chunk_size=array.chunk_size,
    )


def apply_array_task_to_config(config: Dict) -> Dict:
    task = resolve_array_task(config)

    if not task.is_active:
        return config

    config.setdefault("hpc", {})
    config["hpc"]["resolved_array_task_id"] = task.task_id
    config["hpc"]["resolved_array_mode"] = task.mode
    config["hpc"]["resolved_array_items"] = task.items

    if task.mode == "traits":
        config.setdefault("gwas", {})
        config["gwas"]["trait_columns"] = task.items

    elif task.mode == "chromosomes":
        config.setdefault("input", {})
        config["input"]["chromosomes"] = task.items

    elif task.mode == "samples":
        config.setdefault("input", {})
        config["input"]["samples"] = task.items

    elif task.mode == "custom":
        config.setdefault("hpc", {})
        config["hpc"]["custom_tasks"] = task.items

    return config


def array_task_summary(config: Dict) -> str:
    task = resolve_array_task(config)

    if not task.enabled:
        return "HPC array mode is disabled."

    if task.task_id is None:
        return "HPC array mode is enabled, but no array task ID was detected."

    return (
        f"HPC array task resolved | "
        f"mode={task.mode} | "
        f"task_id={task.task_id} | "
        f"items={task.items}"
    )