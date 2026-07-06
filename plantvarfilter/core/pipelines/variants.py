"""
PlantOmicsGWAS Core Pipeline - Variant Calling and VCF Processing
"""

from __future__ import annotations

from pathlib import Path
from typing import List

from plantvarfilter.bcftools_utils import BCFtools
from plantvarfilter.core.context import PipelineContext
from plantvarfilter.core.step import PipelineStep
from plantvarfilter.variant_caller_utils import VariantCaller


class VariantCallingStep(PipelineStep):
    step_id = "variant_calling"
    name = "Variant Calling"
    description = "Call genomic variants from processed BAM files."

    def _find_bams(self, bam_dir: Path) -> List[str]:
        bam_files = sorted(bam_dir.glob("*.bam"))

        if not bam_files:
            raise FileNotFoundError(f"No BAM files found in: {bam_dir}")

        preferred = [
            p for p in bam_files
            if "dedup" in p.name.lower()
            or "final" in p.name.lower()
            or "sorted" in p.name.lower()
        ]

        return [str(p) for p in (preferred if preferred else bam_files)]

    def execute(self, context: PipelineContext) -> None:
        reference = (
            context.get_input("reference_fasta")
            or context.get_input("reference")
            or context.get_input("genome")
        )

        if not reference:
            raise ValueError("Missing required input: reference_fasta")

        reference_path = Path(reference)
        if not reference_path.exists():
            raise FileNotFoundError(f"Reference FASTA not found: {reference_path}")

        bam_file = (
            context.get_input("final_bam")
            or context.get_input("sorted_bam")
            or context.get_input("bam_file")
            or context.get_input("bam")
        )

        bam_dir = (
            context.get_input("sorted_bam_dir")
            or context.get_input("bam_dir")
            or str(context.output_dir / "bam_processing")
        )

        if bam_file:
            bams = [str(bam_file)]
        else:
            bam_dir_path = Path(bam_dir)

            if not bam_dir_path.exists():
                raise FileNotFoundError(f"BAM directory not found: {bam_dir_path}")

            bams = self._find_bams(bam_dir_path)

        for bam in bams:
            if not Path(bam).exists():
                raise FileNotFoundError(f"BAM file not found: {bam}")

        out_dir = context.output_dir / "variants"
        out_dir.mkdir(parents=True, exist_ok=True)

        out_prefix = str(out_dir / "raw_variants")

        regions_bed = context.get_input("regions_bed")
        min_baseq = int(context.compute.get("min_baseq", 13))
        min_mapq = int(context.compute.get("min_mapq", 20))
        ploidy = int(context.compute.get("ploidy", 2))
        split_after_calling = bool(context.compute.get("split_after_calling", False))

        caller = VariantCaller()

        vcf_gz, tbi = caller.call_bcftools(
            bams=bams,
            ref_fasta=str(reference_path),
            out_prefix=out_prefix,
            regions_bed=regions_bed,
            threads=context.threads,
            min_baseq=min_baseq,
            min_mapq=min_mapq,
            ploidy=ploidy,
            log=print,
            split_after_calling=split_after_calling,
        )

        self.add_output("variants_dir", str(out_dir))
        self.add_output("raw_vcf", vcf_gz)
        self.add_output("vcf", vcf_gz)
        self.add_output("raw_vcf_index", tbi)

        split_paths = caller.get_last_split()
        if split_paths:
            for key, value in split_paths.items():
                self.add_output(f"{key}_vcf", value)

        self.result.message = "Variant calling completed successfully."


class BCFtoolsProcessingStep(PipelineStep):
    step_id = "bcftools_processing"
    name = "BCFtools Variant Processing"
    description = "Normalize, sort, annotate IDs, filter, compress, and index VCF files."

    def execute(self, context: PipelineContext) -> None:
        input_vcf = (
            context.get_input("raw_vcf")
            or context.get_input("vcf")
            or str(context.output_dir / "variants" / "raw_variants.vcf.gz")
        )

        input_vcf_path = Path(input_vcf)

        if not input_vcf_path.exists():
            raise FileNotFoundError(f"Input VCF not found: {input_vcf_path}")

        out_dir = context.output_dir / "variants"
        out_dir.mkdir(parents=True, exist_ok=True)

        out_prefix = str(out_dir / "filtered_variants")

        ref_fasta = (
            context.get_input("reference_fasta")
            or context.get_input("reference")
            or context.get_input("genome")
        )

        regions_bed = context.get_input("regions_bed")
        filter_expr = (
            context.compute.get("filter_expr")
            or context.config.get("vcf", {}).get("filter_expr")
        )

        split_multiallelic = bool(context.compute.get("split_multiallelic", True))
        left_align = bool(context.compute.get("left_align", True))
        do_sort = bool(context.compute.get("sort_vcf", True))
        set_id_from_fields = bool(context.compute.get("set_id_from_fields", True))
        remove_filtered = bool(context.compute.get("remove_filtered", False))
        compress_output = bool(context.compute.get("compress_output", True))
        index_output = bool(context.compute.get("index_output", True))
        keep_temps = bool(context.compute.get("keep_temps", False))
        fill_tags = bool(context.compute.get("fill_tags", False))

        bcftools = BCFtools()

        final_vcf, stats_path = bcftools.preprocess(
            input_vcf=str(input_vcf_path),
            out_prefix=out_prefix,
            log=print,
            ref_fasta=ref_fasta,
            regions_bed=regions_bed,
            split_multiallelic=split_multiallelic,
            left_align=left_align,
            do_sort=do_sort,
            set_id_from_fields=set_id_from_fields,
            filter_expr=filter_expr,
            remove_filtered=remove_filtered,
            compress_output=compress_output,
            index_output=index_output,
            keep_temps=keep_temps,
            fill_tags=fill_tags,
        )

        self.add_output("variants_dir", str(out_dir))
        self.add_output("filtered_vcf", final_vcf)
        self.add_output("vcf", final_vcf)

        if stats_path:
            self.add_output("bcftools_stats", stats_path)

        self.result.message = "BCFtools VCF processing completed successfully."


def create_variant_calling_step() -> VariantCallingStep:
    return VariantCallingStep()


def create_bcftools_processing_step() -> BCFtoolsProcessingStep:
    return BCFtoolsProcessingStep()