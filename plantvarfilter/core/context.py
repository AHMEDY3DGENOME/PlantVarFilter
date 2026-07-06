from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Optional


@dataclass
class PipelineContext:
    config: Dict[str, Any]
    config_path: Optional[str] = None

    @property
    def project_name(self) -> str:
        return self.config.get("project", {}).get("name", "plantomicsgwas_project")

    @property
    def input(self) -> Dict[str, Any]:
        return self.config.get("input", {})

    @property
    def output(self) -> Dict[str, Any]:
        return self.config.get("output", {})

    @property
    def compute(self) -> Dict[str, Any]:
        return self.config.get("compute", {})

    @property
    def tools(self) -> Dict[str, Any]:
        return self.config.get("tools", {})

    @property
    def output_dir(self) -> Path:
        return Path(self.output.get("dir", "results"))

    @property
    def logs_dir(self) -> Path:
        return Path(self.output.get("logs_dir", self.output_dir / "logs"))

    @property
    def reports_dir(self) -> Path:
        return Path(self.output.get("reports_dir", self.output_dir / "reports"))

    @property
    def plots_dir(self) -> Path:
        return Path(self.output.get("plots_dir", self.output_dir / "plots"))

    @property
    def threads(self) -> int:
        return int(self.compute.get("threads", 1))

    @property
    def resume(self) -> bool:
        return bool(self.compute.get("resume", False))

    @property
    def overwrite(self) -> bool:
        return bool(self.compute.get("overwrite", False))

    def prepare_directories(self) -> None:
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.logs_dir.mkdir(parents=True, exist_ok=True)
        self.reports_dir.mkdir(parents=True, exist_ok=True)
        self.plots_dir.mkdir(parents=True, exist_ok=True)

    def get_input(self, key: str, default: Any = None) -> Any:
        return self.input.get(key, default)

    def get_output(self, key: str, default: Any = None) -> Any:
        return self.output.get(key, default)

    def get_tool(self, key: str, default: Any = None) -> Any:
        return self.tools.get(key, default)