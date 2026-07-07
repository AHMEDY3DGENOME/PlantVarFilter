"""
PlantOmicsGWAS Core Pipeline - GWAS

Stable wrapper for the existing GWAS CLI.

Supports two entry points:
- PLINK BED prefix: --bed
- VCF / VCF.GZ: --vcf
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

    def _get_phenotype(self, context: PipelineContext) -> str:
        phenotype = (
            context.get_input("phenotype_file")
            or context.get_input("phenotype")
            or context.get_input("pheno")
        )

        if not phenotype:
            raise ValueError("Missing required input: phenotype_file")

        if not Path(phenotype).exists():
            raise FileNotFoundError(f"Phenotype file not found: {phenotype}")

        return str(phenotype)

    def _build_input_args(self, context: PipelineContext) -> list[str]:
        plink_prefix = (
            context.get_input("plink_prefix")
            or context.get_input("bed_prefix")
        )

        if plink_prefix:
            bed_path = Path(str(plink_prefix) + ".bed")
            bim_path = Path(str(plink_prefix) + ".bim")
            fam_path = Path(str(plink_prefix) + ".fam")

            missing = [
                str(path)
                for path in (bed_path, bim_path, fam_path)
                if not path.exists()
            ]

            if missing:
                raise FileNotFoundError(
                    "Missing PLINK files: " + ", ".join(missing)
                )

            return ["--bed", str(plink_prefix)]

        vcf = (
            context.get_input("filtered_vcf")
            or context.get_input("vcf")
            or str(context.output_dir / "variants" / "filtered_variants.vcf.gz")
        )

        if not Path(vcf).exists():
            raise FileNotFoundError(f"VCF file not found: {vcf}")

        return ["--vcf", str(vcf)]

    def execute(self, context: PipelineContext) -> None:
        phenotype = self._get_phenotype(context)

        out_dir = context.output_dir / "gwas"
        out_dir.mkdir(parents=True, exist_ok=True)

        method = (
            context.config.get("gwas", {}).get("method")
            or context.get_tool("gwas_method")
            or "FaST-LMM"
        )

        cmd = [
            "plantomicsgwas",
            "gwas",
            *self._build_input_args(context),
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
            if not Path(covariates).exists():
                raise FileNotFoundError(f"Covariates file not found: {covariates}")
            cmd.extend(["--covar", str(covariates)])

        self.result.command = cmd

        proc = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )

        if proc.stdout:
            print(proc.stdout)

        if proc.returncode != 0:
            raise RuntimeError(proc.stdout)

        self.add_output("gwas_dir", str(out_dir))
        self.add_output("gwas_results", str(out_dir / "gwas_results.csv"))
        self.add_output("manhattan_plot", str(out_dir / "manhattan_plot.png"))
        self.add_output("qq_plot", str(out_dir / "qq_plot.png"))

        self.result.message = "GWAS analysis completed successfully."


class BatchGWASStep(PipelineStep):
    step_id = "batch_gwas"
    name = "Batch GWAS"
    description = "Run GWAS across multiple traits."

    def _build_input_args(self, context: PipelineContext) -> list[str]:
        plink_prefix = (
            context.get_input("plink_prefix")
            or context.get_input("bed_prefix")
        )

        if plink_prefix:
            bed_path = Path(str(plink_prefix) + ".bed")
            bim_path = Path(str(plink_prefix) + ".bim")
            fam_path = Path(str(plink_prefix) + ".fam")

            missing = [
                str(path)
                for path in (bed_path, bim_path, fam_path)
                if not path.exists()
            ]

            if missing:
                raise FileNotFoundError(
                    "Missing PLINK files: " + ", ".join(missing)
                )

            return ["--bed", str(plink_prefix)]

        vcf = (
            context.get_input("filtered_vcf")
            or context.get_input("vcf")
            or str(context.output_dir / "variants" / "filtered_variants.vcf.gz")
        )

        if not Path(vcf).exists():
            raise FileNotFoundError(f"VCF file not found: {vcf}")

        return ["--vcf", str(vcf)]

    def execute(self, context: PipelineContext) -> None:
        phenotype = (
            context.get_input("phenotype_file")
            or context.get_input("phenotype")
            or context.get_input("pheno")
        )

        if not phenotype:
            raise ValueError("Missing required input: phenotype_file")

        if not Path(phenotype).exists():
            raise FileNotFoundError(f"Phenotype file not found: {phenotype}")

        out_dir = context.output_dir / "batch_gwas"
        out_dir.mkdir(parents=True, exist_ok=True)

        cmd = [
            "plantomicsgwas",
            "batch-gwas",
            *self._build_input_args(context),
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
            if not Path(covariates).exists():
                raise FileNotFoundError(f"Covariates file not found: {covariates}")
            cmd.extend(["--covar", str(covariates)])

        self.result.command = cmd

        proc = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )

        if proc.stdout:
            print(proc.stdout)

        if proc.returncode != 0:
            raise RuntimeError(proc.stdout)

        self.add_output("batch_gwas_dir", str(out_dir))
        self.result.message = "Batch GWAS completed successfully."


def create_gwas_step() -> GWASStep:
    return GWASStep()


def create_batch_gwas_step() -> BatchGWASStep:
    return BatchGWASStep()