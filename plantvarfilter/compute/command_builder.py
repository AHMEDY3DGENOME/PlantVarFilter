"""
PlantOmicsGWAS Compute Engine - Command Builder

Converts Compute Engine workflow steps into existing `plantomicsgwas` CLI commands.
This allows the HPC/server layer to reuse the current stable CLI instead of
duplicating pipeline logic.
"""

from pathlib import Path
from typing import Any, Dict, List, Optional


class CommandBuilder:
    def __init__(self, config: Dict[str, Any]):
        self.config = config
        self.input = config.get("input", {})
        self.output = config.get("output", {})
        self.compute = config.get("compute", {})
        self.tools = config.get("tools", {})

        self.output_dir = Path(self.output.get("dir", "results"))
        self.logs_dir = Path(self.output.get("logs_dir", self.output_dir / "logs"))
        self.plots_dir = Path(self.output.get("plots_dir", self.output_dir / "plots"))

        self.threads = str(self.compute.get("threads", 1))

    def build(self, step_id: str) -> Optional[List[str]]:
        builders = {
            "reference_indexing": self._reference_indexing,
            "fastq_qc": self._fastq_qc,
            "alignment": self._alignment,
            "bam_processing": self._bam_processing,
            "variant_calling": self._variant_calling,
            "bcftools_processing": self._bcftools_processing,
            "vcf_quality": self._vcf_quality,
            "annotation": self._annotation,
            "gwas": self._gwas,
            "batch_gwas": self._batch_gwas,
            "genomic_prediction": self._genomic_prediction,
            "ld_analysis": self._ld_analysis,
            "pav_matrix": self._pav_matrix,
            "vcf_pav": self._vcf_pav,
            "pangwas": self._pangwas,
            "pangwas_plots": self._pangwas_plots,
            "plots": self._plots,
        }

        builder = builders.get(step_id)

        if builder is None:
            return None

        return builder()

    def _get(self, key: str, default: Optional[str] = None) -> Optional[str]:
        return self.input.get(key) or self.output.get(key) or default

    def _reference_indexing(self) -> List[str]:
        reference = self._get("reference_fasta")
        out_dir = self.output_dir / "reference_index"

        return [
            "plantomicsgwas",
            "ref-index",
            "--fasta",
            str(reference),
            "--out",
            str(out_dir),
        ]

    def _fastq_qc(self) -> List[str]:
        fastq_dir = self._get("fastq_dir")
        out_dir = self.output_dir / "fastq_qc"

        return [
            "plantomicsgwas",
            "fastq-qc",
            "--input",
            str(fastq_dir),
            "--out",
            str(out_dir),
        ]

    def _alignment(self) -> List[str]:
        fastq_dir = self._get("fastq_dir")
        reference = self._get("reference_fasta")
        out_dir = self.output_dir / "alignment"
        aligner = self.tools.get("aligner", "bowtie2")

        return [
            "plantomicsgwas",
            "align",
            "--reads",
            str(fastq_dir),
            "--reference",
            str(reference),
            "--out",
            str(out_dir),
            "--aligner",
            str(aligner),
            "--threads",
            self.threads,
        ]

    def _bam_processing(self) -> List[str]:
        bam_dir = self._get("bam_dir") or str(self.output_dir / "alignment")
        out_dir = self.output_dir / "bam_processing"

        return [
            "plantomicsgwas",
            "preprocess-bam",
            "--bam-dir",
            str(bam_dir),
            "--out",
            str(out_dir),
            "--threads",
            self.threads,
        ]

    def _variant_calling(self) -> List[str]:
        bam_dir = self._get("sorted_bam_dir") or str(self.output_dir / "bam_processing")
        reference = self._get("reference_fasta")
        out_vcf = self.output_dir / "variants" / "raw_variants.vcf.gz"

        return [
            "plantomicsgwas",
            "call",
            "--bam-dir",
            str(bam_dir),
            "--reference",
            str(reference),
            "--out",
            str(out_vcf),
            "--threads",
            self.threads,
        ]

    def _bcftools_processing(self) -> List[str]:
        raw_vcf = self._get("raw_vcf") or str(self.output_dir / "variants" / "raw_variants.vcf.gz")
        out_vcf = self.output_dir / "variants" / "filtered_variants.vcf.gz"

        return [
            "plantomicsgwas",
            "vcf-prep",
            "--vcf",
            str(raw_vcf),
            "--out",
            str(out_vcf),
        ]

    def _vcf_quality(self) -> List[str]:
        filtered_vcf = self._get("filtered_vcf") or str(
            self.output_dir / "variants" / "filtered_variants.vcf.gz"
        )
        out_json = self.output_dir / "quality" / "vcf_qc.json"

        return [
            "plantomicsgwas",
            "vcf-qc",
            "--vcf",
            str(filtered_vcf),
            "--out-json",
            str(out_json),
        ]

    def _annotation(self) -> List[str]:
        gwas_csv = self.output_dir / "gwas" / "gwas_results.csv"
        gff_file = self._get("gff_file")
        out_csv = self.output_dir / "annotation" / "annotated_gwas.csv"

        command = [
            "plantomicsgwas",
            "annotate",
            "--gwas",
            str(gwas_csv),
            "--out",
            str(out_csv),
        ]

        if gff_file:
            command.extend(["--gff", str(gff_file)])

        return command

    def _gwas(self) -> List[str]:
        filtered_vcf = self._get("filtered_vcf") or str(
            self.output_dir / "variants" / "filtered_variants.vcf.gz"
        )
        phenotype = self._get("phenotype_file")
        out_dir = self.output_dir / "gwas"

        method = self.tools.get("gwas_method", "fastlmm")

        return [
            "plantomicsgwas",
            "gwas",
            "--vcf",
            str(filtered_vcf),
            "--pheno",
            str(phenotype),
            "--algo",
            str(method),
            "--out",
            str(out_dir),
        ]

    def _batch_gwas(self) -> List[str]:
        filtered_vcf = self._get("filtered_vcf") or str(
            self.output_dir / "variants" / "filtered_variants.vcf.gz"
        )
        phenotype = self._get("phenotype_file")
        out_dir = self.output_dir / "batch_gwas"

        return [
            "plantomicsgwas",
            "batch-gwas",
            "--vcf",
            str(filtered_vcf),
            "--pheno",
            str(phenotype),
            "--out",
            str(out_dir),
        ]

    def _genomic_prediction(self) -> List[str]:
        filtered_vcf = self._get("filtered_vcf") or str(
            self.output_dir / "variants" / "filtered_variants.vcf.gz"
        )
        phenotype = self._get("phenotype_file")
        out_dir = self.output_dir / "prediction"

        model = self.config.get("genomic_prediction", {}).get(
            "model",
            self.tools.get("prediction_model", "random_forest"),
        )

        return [
            "plantomicsgwas",
            "predict",
            "--vcf",
            str(filtered_vcf),
            "--pheno",
            str(phenotype),
            "--model",
            str(model),
            "--out",
            str(out_dir),
        ]

    def _ld_analysis(self) -> List[str]:
        filtered_vcf = self._get("filtered_vcf") or str(
            self.output_dir / "variants" / "filtered_variants.vcf.gz"
        )
        out_dir = self.output_dir / "ld"

        return [
            "plantomicsgwas",
            "ld-decay",
            "--vcf",
            str(filtered_vcf),
            "--out",
            str(out_dir),
        ]

    def _pav_matrix(self) -> Optional[List[str]]:
        return None

    def _vcf_pav(self) -> Optional[List[str]]:
        return None

    def _pangwas(self) -> Optional[List[str]]:
        return None

    def _pangwas_plots(self) -> Optional[List[str]]:
        return None

    def _plots(self) -> List[str]:
        gwas_csv = self.output_dir / "gwas" / "gwas_results.csv"
        out_dir = self.plots_dir

        return [
            "plantomicsgwas",
            "gwas-plot",
            "--gwas",
            str(gwas_csv),
            "--out",
            str(out_dir),
        ]


def build_command(step_id: str, config: Dict[str, Any]) -> Optional[List[str]]:
    builder = CommandBuilder(config)
    return builder.build(step_id)