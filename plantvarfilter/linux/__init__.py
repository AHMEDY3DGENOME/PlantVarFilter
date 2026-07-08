# plantvarfilter/linux/__init__.py
from pathlib import Path
import os
import stat
import subprocess
from shutil import which

BIN_DIR = Path(__file__).resolve().parent
LIB_DIR = BIN_DIR / "lib"  # optional vendored .so files (libhts.so.3, etc.)


def _ensure_exec(p: Path) -> None:
    try:
        mode = p.stat().st_mode
        p.chmod(mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    except Exception:
        pass


def _bundled_binary_runs(candidate: Path) -> bool:
    """
    Confirm the bundled binary actually starts (e.g. its shared
    libraries resolve), instead of assuming presence == working.
    """
    env = os.environ.copy()

    if LIB_DIR.exists():
        env["LD_LIBRARY_PATH"] = (
            f"{LIB_DIR}:{env.get('LD_LIBRARY_PATH', '')}"
        )

    try:
        proc = subprocess.run(
            [str(candidate), "--version"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            env=env,
            timeout=10,
        )
        return proc.returncode == 0
    except Exception:
        return False


def resolve_tool(name: str) -> str | None:
    """
    Resolution order (most reliable first):

    1. A working copy already on the system PATH
       (e.g. installed via the conda/mamba environment
       described in the README). This is what most users
       will have, and it carries its own correct dependencies.

    2. The bundled binary in plantvarfilter/linux/<name>,
       but ONLY if it actually runs on this machine
       (some systems are missing shared libraries the
       bundled binary was built against).

    Returns None, with a clear error surfaced by the caller,
    if neither option works.
    """
    sys_path = which(name)
    if sys_path:
        return sys_path

    candidate = BIN_DIR / name

    if candidate.exists():
        _ensure_exec(candidate)

        if _bundled_binary_runs(candidate):
            return str(candidate)

        raise RuntimeError(
            f"\nThe bundled '{name}' binary could not run on this system "
            f"(likely a missing shared library, e.g. libhts.so.3), and no "
            f"working system installation of '{name}' was found on PATH.\n\n"
            f"Fix with either:\n"
            f"  1) Install it via conda/mamba (recommended):\n"
            f"       mamba install -c bioconda {name}\n"
            f"  2) Or install the missing system library, e.g. (Rocky/RHEL):\n"
            f"       dnf install -y htslib htslib-libs\n"
        )

    return None