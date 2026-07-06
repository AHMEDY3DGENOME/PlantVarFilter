"""
PlantOmicsGWAS Core Checkpoint Manager

Provides resume/checkpoint support for long-running workflows.

Each completed step writes a small JSON checkpoint file:

results/
  checkpoints/
    alignment.done.json
    variant_calling.done.json

On resume, PipelineRunner can skip completed steps safely.
"""

from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Any, Dict, List, Optional


class CheckpointManager:
    def __init__(self, output_dir: str | Path, enabled: bool = True):
        self.output_dir = Path(output_dir)
        self.enabled = bool(enabled)
        self.checkpoint_dir = self.output_dir / "checkpoints"

        if self.enabled:
            self.checkpoint_dir.mkdir(parents=True, exist_ok=True)

    def checkpoint_path(self, step_id: str) -> Path:
        safe_step_id = step_id.replace("/", "_").replace("\\", "_")
        return self.checkpoint_dir / f"{safe_step_id}.done.json"

    def exists(self, step_id: str) -> bool:
        if not self.enabled:
            return False
        return self.checkpoint_path(step_id).exists()

    def save(
        self,
        step_id: str,
        status: str = "success",
        message: str = "",
        outputs: Optional[Dict[str, Any]] = None,
        runtime: Optional[float] = None,
    ) -> str:
        if not self.enabled:
            return ""

        payload = {
            "step_id": step_id,
            "status": status,
            "message": message,
            "outputs": outputs or {},
            "runtime": runtime,
            "timestamp": time.time(),
        }

        path = self.checkpoint_path(step_id)
        path.write_text(
            json.dumps(payload, indent=2),
            encoding="utf-8",
        )

        return str(path)

    def load(self, step_id: str) -> Dict[str, Any]:
        path = self.checkpoint_path(step_id)

        if not path.exists():
            raise FileNotFoundError(f"Checkpoint not found: {path}")

        return json.loads(path.read_text(encoding="utf-8"))

    def get_outputs(self, step_id: str) -> Dict[str, Any]:
        if not self.exists(step_id):
            return {}

        data = self.load(step_id)
        outputs = data.get("outputs", {})

        if isinstance(outputs, dict):
            return outputs

        return {}

    def delete(self, step_id: str) -> bool:
        path = self.checkpoint_path(step_id)

        if path.exists():
            path.unlink()
            return True

        return False

    def clear(self) -> int:
        if not self.checkpoint_dir.exists():
            return 0

        count = 0

        for path in self.checkpoint_dir.glob("*.done.json"):
            path.unlink()
            count += 1

        return count

    def list_checkpoints(self) -> List[Dict[str, Any]]:
        if not self.checkpoint_dir.exists():
            return []

        checkpoints: List[Dict[str, Any]] = []

        for path in sorted(self.checkpoint_dir.glob("*.done.json")):
            try:
                data = json.loads(path.read_text(encoding="utf-8"))
                data["_path"] = str(path)
                checkpoints.append(data)
            except Exception:
                checkpoints.append(
                    {
                        "step_id": path.stem.replace(".done", ""),
                        "status": "unknown",
                        "_path": str(path),
                    }
                )

        return checkpoints