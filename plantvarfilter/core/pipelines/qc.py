"""
PlantOmicsGWAS Core Pipeline - FASTQ Quality Control

This wrapper connects the Core Engine with the existing FASTQ QC module:
plantvarfilter.preanalysis.fastq_qc.run_fastq_qc
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

from plantvarfilter.core.context import PipelineContext
from plantvarfilter.core.step import PipelineStep
from plantvarfilter.preanalysis.fastq_qc import run_fastq_qc


class FastqQCStep(PipelineStep):
    step_id = "fastq_qc"
    name = "FASTQ Quality Control"
    description = "Run quality control on raw FASTQ reads."

    def _find_fastq_pair(self, fastq_dir: Path) -> tuple[str, Optional[str]]:
        fastq_files = sorted(
            list(fastq_dir.glob("*.fastq"))
            + list(fastq_dir.glob("*.fq"))
            + list(fastq_dir.glob("*.fastq.gz"))
            + list(fastq_dir.glob("*.fq.gz"))
        )

        if not fastq_files:
            raise FileNotFoundError(f"No FASTQ files found in: {fastq_dir}")

        reads1 = str(fastq_files[0])
        reads2 = str(fastq_files[1]) if len(fastq_files) > 1 else None

        return reads1, reads2

    def execute(self, context: PipelineContext) -> None:
        fastq_dir = context.get_input("fastq_dir")
        reads1 = context.get_input("reads1")
        reads2 = context.get_input("reads2")

        out_dir = context.output_dir / "fastq_qc"

        if not reads1:
            if not fastq_dir:
                raise ValueError("Missing required input: fastq_dir or reads1")

            fastq_path = Path(fastq_dir)

            if not fastq_path.exists():
                raise FileNotFoundError(f"FASTQ directory not found: {fastq_path}")

            reads1, reads2 = self._find_fastq_pair(fastq_path)

        platform = context.get_input("platform", "illumina")
        sample_max = int(context.compute.get("sample_max", 1_000_000))
        use_fastqc = bool(context.compute.get("use_fastqc_if_available", True))

        out_dir.mkdir(parents=True, exist_ok=True)

        report = run_fastq_qc(
            reads1=reads1,
            reads2=reads2,
            platform=platform,
            out_dir=str(out_dir),
            sample_max=sample_max,
            use_fastqc_if_available=use_fastqc,
            logger=print,
        )

        self.add_output("fastq_qc_dir", str(out_dir))
        self.add_output("report_txt", report.report_txt)
        self.add_output("details_json", report.details_json)

        if report.length_hist_png:
            self.add_output("length_hist_png", report.length_hist_png)

        if report.gc_hist_png:
            self.add_output("gc_hist_png", report.gc_hist_png)

        if report.per_cycle_q_mean_png:
            self.add_output("per_cycle_q_mean_png", report.per_cycle_q_mean_png)

        if report.fastqc_html:
            self.add_output("fastqc_html", report.fastqc_html)

        self.result.message = f"FASTQ QC completed with verdict: {report.verdict}"


def create_step() -> FastqQCStep:
    return FastqQCStep()