"""
PlantOmicsGWAS - HPC environment setup command.

This provides the `plantomicsgwas-setup-env` console command, which runs
the bundled install_plantomicsgwas.sh script shipped inside this package.

Why this exists:
Researchers on a fresh HPC account often have no conda, no bcftools/
samtools/plink, and no Python environment prepared. This command
bootstraps all of that in one step, entirely from the package installed
via PyPI -- no external download, no GitHub link required.

Usage:
    pip install plantomicsgwas
    plantomicsgwas-setup-env [environment_name]
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


def main() -> None:
    script_path = Path(__file__).resolve().parent / "install_plantomicsgwas.sh"

    if not script_path.exists():
        print(
            f"[ERROR] Bundled installer script not found at: {script_path}\n"
            f"This is a packaging issue -- please report it.",
            file=sys.stderr,
        )
        sys.exit(1)

    cmd = ["bash", str(script_path), *sys.argv[1:]]

    result = subprocess.run(cmd)
    sys.exit(result.returncode)


if __name__ == "__main__":
    main()