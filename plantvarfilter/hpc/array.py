"""
PlantOmicsGWAS HPC Job Arrays

Utilities for generating scheduler-independent array jobs.

Supported schedulers:
- SLURM
- PBS/Torque
- LSF

The Compute Engine can use this information to split analyses
across traits, chromosomes, samples, or user-defined tasks.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional


@dataclass
class ArrayJob:
    """
    Description of one HPC array job.
    """

    enabled: bool
    mode: str
    items: List[str]
    chunk_size: int = 1

    @property
    def size(self) -> int:
        return len(self.items)

    @property
    def task_count(self) -> int:
        if self.chunk_size <= 1:
            return self.size

        return (self.size + self.chunk_size - 1) // self.chunk_size


def _to_list(value) -> List[str]:
    if value is None:
        return []

    if isinstance(value, list):
        return [str(v) for v in value]

    return [str(value)]


def build_array_job(config: Dict) -> Optional[ArrayJob]:
    """
    Build an ArrayJob from the YAML configuration.

    Supported modes:

    traits
        Run one GWAS job per phenotype.

    chromosomes
        Run one GWAS job per chromosome.

    samples
        Run preprocessing/alignment per sample.

    custom
        User-provided task list.
    """

    hpc = config.get("hpc", {})

    if not hpc.get("array", False):
        return None

    mode = hpc.get("array_mode", "traits")

    chunk_size = int(hpc.get("chunk_size", 1))

    items: List[str] = []

    if mode == "traits":
        gwas = config.get("gwas", {})
        items = _to_list(gwas.get("trait_columns"))

    elif mode == "chromosomes":
        items = _to_list(
            hpc.get(
                "chromosomes",
                [str(i) for i in range(1, 13)],
            )
        )

    elif mode == "samples":
        items = _to_list(
            hpc.get("samples")
        )

    elif mode == "custom":
        items = _to_list(
            hpc.get("tasks")
        )

    else:
        raise ValueError(
            f"Unsupported HPC array mode: {mode}"
        )

    return ArrayJob(
        enabled=True,
        mode=mode,
        items=items,
        chunk_size=chunk_size,
    )


def slurm_array_directive(array: ArrayJob) -> str:
    """
    Returns:

    #SBATCH --array=0-49
    """

    if array.task_count == 0:
        return ""

    return f"#SBATCH --array=0-{array.task_count-1}"


def pbs_array_directive(array: ArrayJob) -> str:
    """
    Returns:

    #PBS -J 0-49
    """

    if array.task_count == 0:
        return ""

    return f"#PBS -J 0-{array.task_count-1}"


def lsf_array_directive(array: ArrayJob) -> str:
    """
    Returns:

    #BSUB -J myjob[1-50]
    """

    if array.task_count == 0:
        return ""

    return f"[1-{array.task_count}]"


def array_summary(array: Optional[ArrayJob]) -> str:
    if array is None:
        return "No HPC array configured."

    return (
        f"Mode: {array.mode} | "
        f"Items: {array.size} | "
        f"Tasks: {array.task_count} | "
        f"Chunk Size: {array.chunk_size}"
    )