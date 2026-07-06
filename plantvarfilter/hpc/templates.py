"""
PlantOmicsGWAS HPC Job Templates

Generates scheduler-specific job scripts.

Supported:
- Local
- SLURM
- PBS/Torque
- LSF
"""

from __future__ import annotations

from pathlib import Path
from textwrap import dedent
from typing import Dict

from plantvarfilter.hpc.array import (
    build_array_job,
    lsf_array_directive,
    pbs_array_directive,
    slurm_array_directive,
)


def _abs(path: str) -> str:
    return str(Path(path).expanduser().resolve())


def _conda_activate(hpc: Dict) -> str:
    if hpc.get("conda_env"):
        return f"conda activate {hpc['conda_env']}"
    return ""


def _array_env_block(config: Dict) -> str:
    array = build_array_job(config)

    if array is None:
        return ""

    return dedent(
        f"""\
        export PLANTOMICS_ARRAY_MODE="{array.mode}"
        export PLANTOMICS_ARRAY_TASK_COUNT="{array.task_count}"
        export PLANTOMICS_ARRAY_CHUNK_SIZE="{array.chunk_size}"
        """
    )


def build_slurm_script(config: Dict, config_path: str) -> str:
    hpc = config.get("hpc", {})

    cfg = _abs(config_path)
    workdir = str(Path(cfg).parent)
    array = build_array_job(config)

    lines = [
        "#!/bin/bash",
        f"#SBATCH --job-name={hpc.get('job_name', 'plantomicsgwas')}",
        f"#SBATCH --partition={hpc.get('partition', 'compute')}",
        f"#SBATCH --nodes={hpc.get('nodes', 1)}",
        f"#SBATCH --ntasks={hpc.get('tasks_per_node', 1)}",
        f"#SBATCH --cpus-per-task={hpc.get('cpus_per_task', 8)}",
        f"#SBATCH --mem={hpc.get('memory', '32G')}",
        f"#SBATCH --time={hpc.get('time', '24:00:00')}",
    ]

    if array:
        directive = slurm_array_directive(array)
        if directive:
            lines.append(directive)

    if hpc.get("account"):
        lines.append(f"#SBATCH --account={hpc['account']}")

    if hpc.get("email"):
        lines.append(f"#SBATCH --mail-user={hpc['email']}")
        if hpc.get("mail_type"):
            lines.append(f"#SBATCH --mail-type={hpc['mail_type']}")

    body = f"""
set -e

echo "======================================"
echo " PlantOmicsGWAS HPC Job"
echo " Scheduler: SLURM"
echo "======================================"

cd "{workdir}"

source ~/.bashrc

{_conda_activate(hpc)}

{_array_env_block(config)}

plantomicsgwas-compute run --config "{cfg}"

echo "Finished."
"""

    return "\n".join(lines) + "\n\n" + dedent(body)


def build_pbs_script(config: Dict, config_path: str) -> str:
    hpc = config.get("hpc", {})

    cfg = _abs(config_path)
    workdir = str(Path(cfg).parent)
    array = build_array_job(config)

    lines = [
        "#!/bin/bash",
        f"#PBS -N {hpc.get('job_name', 'plantomicsgwas')}",
        f"#PBS -l nodes={hpc.get('nodes', 1)}:ppn={hpc.get('cpus_per_task', 8)}",
        f"#PBS -l mem={hpc.get('memory', '32G')}",
        f"#PBS -l walltime={hpc.get('time', '24:00:00')}",
    ]

    if array:
        directive = pbs_array_directive(array)
        if directive:
            lines.append(directive)

    body = f"""
set -e

echo "======================================"
echo " PlantOmicsGWAS HPC Job"
echo " Scheduler: PBS"
echo "======================================"

cd "{workdir}"

source ~/.bashrc

{_conda_activate(hpc)}

{_array_env_block(config)}

plantomicsgwas-compute run --config "{cfg}"

echo "Finished."
"""

    return "\n".join(lines) + "\n\n" + dedent(body)


def build_lsf_script(config: Dict, config_path: str) -> str:
    hpc = config.get("hpc", {})
    job_name = hpc.get("job_name", "plantomicsgwas")

    cfg = _abs(config_path)
    workdir = str(Path(cfg).parent)
    array = build_array_job(config)

    if array:
        job_name = f"{job_name}{lsf_array_directive(array)}"

    lines = [
        "#!/bin/bash",
        f"#BSUB -J {job_name}",
        f"#BSUB -n {hpc.get('cpus_per_task', 8)}",
        f"#BSUB -M {hpc.get('memory', '32000')}",
        f"#BSUB -W {hpc.get('time', '24:00')}",
    ]

    body = f"""
set -e

echo "======================================"
echo " PlantOmicsGWAS HPC Job"
echo " Scheduler: LSF"
echo "======================================"

cd "{workdir}"

source ~/.bashrc

{_conda_activate(hpc)}

{_array_env_block(config)}

plantomicsgwas-compute run --config "{cfg}"

echo "Finished."
"""

    return "\n".join(lines) + "\n\n" + dedent(body)


def build_local_script(config_path: str) -> str:
    cfg = _abs(config_path)
    workdir = str(Path(cfg).parent)

    return dedent(
        f"""\
        #!/bin/bash

        set -e

        cd "{workdir}"

        plantomicsgwas-compute run --config "{cfg}"
        """
    )