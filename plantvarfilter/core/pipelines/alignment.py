"""
PlantOmicsGWAS Core Pipeline - Read Alignment

Wrapper for:
plantvarfilter.preanalysis.aligner.Aligner
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

from plantvarfilter.core.context import PipelineContext
from plantvarfilter.core.step import PipelineStep
from plantvarfilter.preanalysis.aligner import Aligner


class AlignmentStep(PipelineStep):
    step_id = "alignment"
    name = "Read Alignment"
    description = "Align sequencing reads to the reference genome."

    def _find_fastq_pair(self, fastq_dir: Path) -> tuple[str, Optional[str]]:
        fastq_files = sorted(
            list(fastq_dir.glob("*.fastq"))
            + list(fastq_dir.glob("*.fq"))
            + list(fastq_dir.glob("*.fastq.gz"))
            + list(fastq_dir.glob("*.fq.gz"))
        )

        if not fastq_files:
            raise FileNotFoundError(f"No FASTQ files found in: {fastq_dir}")

        reads1_candidates = [
            p for p in fastq_files
            if "_R1" in p.name or "_1" in p.name or "R1" in p.name
        ]
        reads2_candidates = [
            p for p in fastq_files
            if "_R2" in p.name or "_2" in p.name or "R2" in p.name
        ]

        if reads1_candidates:
            reads1 = reads1_candidates[0]
            reads2 = reads2_candidates[0] if reads2_candidates else None
        else:
            reads1 = fastq_files[0]
            reads2 = fastq_files[1] if len(fastq_files) > 1 else None

        return str(reads1), str(reads2) if reads2 else None

    def _resolve_reference(self, context: PipelineContext, platform: str) -> str:
        if platform.lower() in {"ont", "nanopore", "pb", "pacbio", "hifi", "long"}:
            reference = (
                context.get_input("mmi")
                or context.get_input("reference_mmi")
                or context.get_input("reference_fasta")
                or context.get_input("reference")
            )
        else:
            reference = (
                context.get_input("bt2_prefix")
                or context.get_input("reference_index")
                or context.get_input("reference_fasta")
                or context.get_input("reference")
            )

        if not reference:
            raise ValueError("Missing reference input.")

        ref_path = Path(reference)

        if platform.lower() in {"ont", "nanopore", "pb", "pacbio", "hifi", "long"}:
            if not ref_path.exists():
                raise FileNotFoundError(f"Reference FASTA/MMI not found: {ref_path}")
        else:
            if str(reference).endswith((".fa", ".fasta", ".fna")):
                if not ref_path.exists():
                    raise FileNotFoundError(f"Reference FASTA not found: {ref_path}")
            else:
                # bowtie2 uses prefix, so prefix itself may not exist as a file.
                bt2_files = list(ref_path.parent.glob(ref_path.name + "*.bt2")) + list(
                    ref_path.parent.glob(ref_path.name + "*.bt2l")
                )
                if not bt2_files and not ref_path.exists():
                    raise FileNotFoundError(
                        f"Bowtie2 index prefix not found: {reference}"
                    )

        return str(reference)

    def execute(self, context: PipelineContext) -> None:
        reads1 = context.get_input("reads1")
        reads2 = context.get_input("reads2")
        fastq_dir = context.get_input("fastq_dir")

        if not reads1:
            if not fastq_dir:
                raise ValueError("Missing required input: fastq_dir or reads1")

            fastq_path = Path(fastq_dir)

            if not fastq_path.exists():
                raise FileNotFoundError(f"FASTQ directory not found: {fastq_path}")

            reads1, reads2 = self._find_fastq_pair(fastq_path)

        if not Path(reads1).exists():
            raise FileNotFoundError(f"reads1 not found: {reads1}")

        if reads2 and not Path(reads2).exists():
            raise FileNotFoundError(f"reads2 not found: {reads2}")

        platform = (
            context.get_input("platform")
            or context.get_tool("platform")
            or "illumina"
        )

        reference = self._resolve_reference(context, platform)

        out_dir = context.output_dir / "alignment"
        out_dir.mkdir(parents=True, exist_ok=True)

        save_sam = bool(context.compute.get("save_sam", False))

        aligner = Aligner(
            logger=print,
            workspace=str(context.output_dir),
        )

        result = aligner.align(
            platform=platform,
            reference=reference,
            reads1=str(reads1),
            reads2=str(reads2) if reads2 else None,
            threads=context.threads,
            save_sam=save_sam,
            out_dir=str(out_dir),
            out_prefix="plantomicsgwas_alignment",
        )

        self.add_output("alignment_dir", str(out_dir))
        self.add_output("bam", result.bam)
        self.add_output("alignment_bam", result.bam)
        self.add_output("bam_file", result.bam)
        self.add_output("bai", result.bai)
        self.add_output("flagstat", result.flagstat)

        if result.sam:
            self.add_output("sam", result.sam)

        self.result.command = result.cmdline.split()
        self.result.message = (
            f"Alignment completed using {result.tool} "
            f"in {result.elapsed_sec:.2f} seconds."
        )


def create_step() -> AlignmentStep:
    return AlignmentStep()