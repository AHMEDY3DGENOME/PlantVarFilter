"""
PlantOmicsGWAS Core Pipeline - BAM Processing

Wrapper for:
plantvarfilter.samtools_utils.Samtools
"""

from __future__ import annotations

from pathlib import Path

from plantvarfilter.core.context import PipelineContext
from plantvarfilter.core.step import PipelineStep
from plantvarfilter.samtools_utils import Samtools


class BamProcessingStep(PipelineStep):
    step_id = "bam_processing"
    name = "BAM Processing"
    description = "Sort, fixmate, mark duplicates, index, and generate BAM QC statistics."

    def _find_bam(self, bam_dir: Path) -> str:
        bam_files = sorted(bam_dir.glob("*.bam"))

        if not bam_files:
            raise FileNotFoundError(f"No BAM files found in: {bam_dir}")

        # Prefer alignment output if present
        preferred = [
            p for p in bam_files
            if "alignment" in p.name.lower() or "sorted" in p.name.lower()
        ]

        return str(preferred[0] if preferred else bam_files[0])

    def execute(self, context: PipelineContext) -> None:
        bam_file = (
            context.get_input("final_bam")
            or context.get_input("bam_file")
            or context.get_input("alignment_bam")
            or context.get_input("bam")
        )

        bam_dir = (
            context.get_input("bam_dir")
            or str(context.output_dir / "alignment")
        )

        if not bam_file:
            bam_dir_path = Path(bam_dir)

            if not bam_dir_path.exists():
                raise FileNotFoundError(f"BAM directory not found: {bam_dir_path}")

            bam_file = self._find_bam(bam_dir_path)

        bam_path = Path(bam_file)

        if not bam_path.exists():
            raise FileNotFoundError(f"BAM file not found: {bam_path}")

        out_dir = context.output_dir / "bam_processing"
        out_dir.mkdir(parents=True, exist_ok=True)

        out_prefix = str(out_dir / bam_path.stem)

        remove_dups = bool(context.compute.get("remove_duplicates", False))
        compute_stats = bool(context.compute.get("compute_bam_stats", True))
        keep_temps = bool(context.compute.get("keep_temps", False))

        samtools = Samtools()

        result = samtools.preprocess(
            input_path=str(bam_path),
            out_prefix=out_prefix,
            threads=context.threads,
            remove_dups=remove_dups,
            compute_stats=compute_stats,
            log=print,
            keep_temps=keep_temps,
        )

        self.add_output("bam_processing_dir", str(out_dir))
        self.add_output("final_bam", result.final_bam)
        self.add_output("sorted_bam", result.final_bam)
        self.add_output("bam_file", result.final_bam)

        if result.bai:
            self.add_output("bai", result.bai)
            self.add_output("bam_index", result.bai)

        for key, value in result.stats_files.items():
            self.add_output(key, value)

        self.result.message = "BAM processing completed successfully."


def create_step() -> BamProcessingStep:
    return BamProcessingStep()