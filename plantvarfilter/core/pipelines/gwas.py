"""
PlantOmicsGWAS Core Pipeline - GWAS

Stable wrapper for the existing GWAS CLI.
Later, this can be refactored to call plantvarfilter.gwas_pipeline.GWAS directly.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

from plantvarfilter.core.context import PipelineContext
from plantvarfilter.core.step import PipelineStep


class GWASStep(PipelineStep):
    step_id = "gwas"
    name = "GWAS Analysis"
    description = "Run genome-wide association analysis."

    def execute(self, context: PipelineContext) -> None:
        vcf = (
            context.get_input("filtered_vcf")
            or context.get_input("vcf")
            or str(context.output_dir / "variants" / "filtered_variants.vcf.gz")
        )

        phenotype = (
            context.get_input("phenotype_file")
            or context.get_input("phenotype")
            or context.get_input("pheno")
        )

        if not phenotype:
            raise ValueError("Missing required input: phenotype_file")

        if not Path(vcf).exists():
            raise FileNotFoundError(f"VCF file not found: {vcf}")

        if not Path(phenotype).exists():
            raise FileNotFoundError(f"Phenotype file not found: {phenotype}")

        out_dir = context.output_dir / "gwas"
        out_dir.mkdir(parents=True, exist_ok=True)

        method = (
            context.config.get("gwas", {}).get("method")
            or context.get_tool("gwas_method")
            or "fastlmm"
        )

        cmd = [
            "plantomicsgwas",
            "gwas",
            "--vcf",
            str(vcf),
            "--pheno",
            str(phenotype),
            "--algo",
            str(method),
            "--out",
            str(out_dir),
        ]

        covariates = (
            context.get_input("covariates_file")
            or context.get_input("covar")
        )

        if covariates:
            cmd.extend(["--covar", str(covariates)])

        self.result.command = cmd

        proc = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )

        if proc.returncode != 0:
            raise RuntimeError(proc.stdout)

        self.add_output("gwas_dir", str(out_dir))
        self.add_output("gwas_results", str(out_dir / "gwas_results.csv"))

        self.result.message = "GWAS analysis completed successfully."


class BatchGWASStep(PipelineStep):
    step_id = "batch_gwas"
    name = "Batch GWAS"
    description = "Run GWAS across multiple traits."

    def execute(self, context: PipelineContext) -> None:
        vcf = (
            context.get_input("filtered_vcf")
            or context.get_input("vcf")
            or str(context.output_dir / "variants" / "filtered_variants.vcf.gz")
        )

        phenotype = (
            context.get_input("phenotype_file")
            or context.get_input("phenotype")
            or context.get_input("pheno")
        )

        if not phenotype:
            raise ValueError("Missing required input: phenotype_file")

        if not Path(vcf).exists():
            raise FileNotFoundError(f"VCF file not found: {vcf}")

        if not Path(phenotype).exists():
            raise FileNotFoundError(f"Phenotype file not found: {phenotype}")

        out_dir = context.output_dir / "batch_gwas"
        out_dir.mkdir(parents=True, exist_ok=True)

        cmd = [
            "plantomicsgwas",
            "batch-gwas",
            "--vcf",
            str(vcf),
            "--pheno",
            str(phenotype),
            "--out",
            str(out_dir),
        ]

        covariates = (
            context.get_input("covariates_file")
            or context.get_input("covar")
        )

        if covariates:
            cmd.extend(["--covar", str(covariates)])

        self.result.command = cmd

        proc = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )

        if proc.returncode != 0:
            raise RuntimeError(proc.stdout)

        self.add_output("batch_gwas_dir", str(out_dir))
        self.result.message = "Batch GWAS completed successfully."


def create_gwas_step() -> GWASStep:
    return GWASStep()


def create_batch_gwas_step() -> BatchGWASStep:
    return BatchGWASStep()