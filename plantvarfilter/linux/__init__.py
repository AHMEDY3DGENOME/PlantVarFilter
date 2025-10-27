# plantvarfilter/linux/__init__.py
from pathlib import Path
import os, stat
from shutil import which

BIN_DIR = Path(__file__).resolve().parent

def _ensure_exec(p: Path):
    try:
        mode = p.stat().st_mode
        p.chmod(mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    except Exception:
        pass

def resolve_tool(name: str) -> str | None:
    """Prefer bundled binary in plantvarfilter/linux/<name>, else fall back to PATH."""
    candidate = BIN_DIR / name
    if candidate.exists():
        _ensure_exec(candidate)
        return str(candidate)

    sys_path = which(name)
    return sys_path
