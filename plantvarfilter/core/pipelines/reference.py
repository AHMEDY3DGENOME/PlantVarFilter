"""
PlantOmicsGWAS Core Pipeline - Reference Indexing

Wrapper for:
plantvarfilter.preanalysis.reference_manager.ReferenceManager
"""

from __future__ import annotations

from pathlib import Path

from plantvarfilter.core.context import PipelineContext
from plantvarfilter.core.step import PipelineStep
from plantvarfilter.preanalysis.reference_manager import ReferenceManager


class ReferenceIndexingStep(PipelineStep):
    step_id = "reference_indexing"
    name = "Reference Genome Indexing"
    description = "Build and validate reference genome indexes."

    def execute(self, context: PipelineContext) -> None:
        reference = (
            context.get_input("reference_fasta")
            or context.get_input("reference")
            or context.get_input("genome")
        )

        if not reference:
            raise ValueError("Missing required input: reference_fasta")

        reference_path = Path(reference).expanduser()

        if not reference_path.exists():
            raise FileNotFoundError(f"Reference FASTA not found: {reference_path}")

        out_dir = context.output_dir / "reference_index"
        out_dir.mkdir(parents=True, exist_ok=True)

        ref_cfg = context.config.get("reference", {})

        build_mmi = bool(ref_cfg.get("build_mmi", True))
        build_bt2 = bool(ref_cfg.get("build_bt2", True))
        build_dict = bool(ref_cfg.get("build_dict", True))

        manager = ReferenceManager(
            logger=print,
            workspace=str(context.output_dir),
        )

        status = manager.build_indices(
            fasta=str(reference_path),
            out_dir=str(out_dir),
            build_mmi=build_mmi,
            build_bt2=build_bt2,
            build_dict=build_dict,
        )

        self.add_output("reference_dir", status.reference_dir)
        self.add_output("reference_fasta", status.fasta)

        if status.faidx:
            self.add_output("faidx", status.faidx)

        if status.dict:
            self.add_output("dict", status.dict)

        if status.mmi:
            self.add_output("mmi", status.mmi)

        if status.bt2_prefix:
            self.add_output("bt2_prefix", status.bt2_prefix)

        self.result.message = (
            "Reference indexing completed successfully."
            if status.ok
            else "Reference indexing completed with missing optional indexes."
        )


def create_step() -> ReferenceIndexingStep:
    return ReferenceIndexingStep()