"""
PlantOmicsGWAS Core Pipeline - Pangenome / panGWAS

Wrappers for:
- plantvarfilter.pangenome_module.pav_matrix_builder.GenePresenceAbsenceBuilder
- plantvarfilter.pangenome_module.vcf_pav_builder.VariantPAVBuilder
- plantvarfilter.pangenome_module.pangwas_pipeline
- plantvarfilter.pangenome_module.pangwas_plots
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from plantvarfilter.core.context import PipelineContext
from plantvarfilter.core.step import PipelineStep
from plantvarfilter.pangenome_module.pav_matrix_builder import GenePresenceAbsenceBuilder
from plantvarfilter.pangenome_module.vcf_pav_builder import VariantPAVBuilder
from plantvarfilter.pangenome_module.pangwas_pipeline import (
    load_covariates,
    load_phenotype,
    run_pangwas,
)
from plantvarfilter.pangenome_module.pangwas_plots import plot_pangwas_manhattan


class PAVMatrixStep(PipelineStep):
    step_id = "pav_matrix"
    name = "PAV Matrix Construction"
    description = "Build gene presence/absence matrix from GFF/GTF files."

    def execute(self, context: PipelineContext) -> None:
        gff_dir = (
            context.get_input("gff_dir")
            or context.get_input("pav_input")
        )

        if not gff_dir:
            raise ValueError("Missing required input: gff_dir or pav_input")

        gff_dir_path = Path(gff_dir)
        if not gff_dir_path.exists():
            raise FileNotFoundError(f"GFF/GTF directory not found: {gff_dir_path}")

        out_dir = context.output_dir / "pangenome"
        out_dir.mkdir(parents=True, exist_ok=True)

        pangenome_cfg = context.config.get("pangenome", {})

        sample_list = pangenome_cfg.get("sample_list")
        gene_id_field = pangenome_cfg.get("gene_id_field", "ID")
        max_genes = pangenome_cfg.get("max_genes")
        sort_genes = bool(pangenome_cfg.get("sort_genes", True))

        out_matrix = out_dir / "gene_pav_matrix.csv"

        builder = GenePresenceAbsenceBuilder(
            gff_dir=str(gff_dir_path),
            sample_list=sample_list,
            gene_id_field=gene_id_field,
            max_genes=max_genes,
            sort_genes=sort_genes,
        )

        matrix_path = builder.save_matrix(str(out_matrix))

        self.add_output("pangenome_dir", str(out_dir))
        self.add_output("pav_matrix", matrix_path)
        self.result.message = "Gene PAV matrix generated successfully."


class VCFPAVStep(PipelineStep):
    step_id = "vcf_pav"
    name = "VCF to PAV Matrix"
    description = "Build variant presence/absence matrix from multi-sample VCF."

    def execute(self, context: PipelineContext) -> None:
        vcf = (
            context.get_input("filtered_vcf")
            or context.get_input("vcf")
        )

        if not vcf:
            raise ValueError("Missing required input: filtered_vcf or vcf")

        vcf_path = Path(vcf)
        if not vcf_path.exists():
            raise FileNotFoundError(f"VCF not found: {vcf_path}")

        out_dir = context.output_dir / "pangenome"
        out_dir.mkdir(parents=True, exist_ok=True)

        pangenome_cfg = context.config.get("pangenome", {})
        max_variants = int(pangenome_cfg.get("max_variants", 0))

        out_matrix = out_dir / "variant_pav_matrix.csv"

        builder = VariantPAVBuilder(
            vcf_path=str(vcf_path),
            max_variants=max_variants,
        )

        df = builder.build_matrix()
        df.to_csv(out_matrix)

        self.add_output("pangenome_dir", str(out_dir))
        self.add_output("pav_matrix", str(out_matrix))
        self.result.message = "Variant PAV matrix generated successfully."


class PanGWASStep(PipelineStep):
    step_id = "pangwas"
    name = "Pan-GWAS Analysis"
    description = "Run pan-genome-wide association analysis using a PAV matrix."

    def execute(self, context: PipelineContext) -> None:
        pav_matrix = (
            context.get_input("pav_matrix")
            or str(context.output_dir / "pangenome" / "gene_pav_matrix.csv")
        )

        phenotype = (
            context.get_input("phenotype_file")
            or context.get_input("phenotype")
        )

        if not phenotype:
            raise ValueError("Missing required input: phenotype_file")

        pav_path = Path(pav_matrix)
        pheno_path = Path(phenotype)

        if not pav_path.exists():
            raise FileNotFoundError(f"PAV matrix not found: {pav_path}")

        if not pheno_path.exists():
            raise FileNotFoundError(f"Phenotype file not found: {pheno_path}")

        out_dir = context.output_dir / "pangenome"
        out_dir.mkdir(parents=True, exist_ok=True)

        pangenome_cfg = context.config.get("pangenome", {})

        trait = (
            pangenome_cfg.get("trait")
            or context.config.get("gwas", {}).get("trait_column")
        )

        if not trait:
            raise ValueError("Missing pangenome.trait or gwas.trait_column")

        category = pangenome_cfg.get("category", "baseline")
        method = pangenome_cfg.get("method", "ttest")

        covariates_file = (
            context.get_input("covariates_file")
            or pangenome_cfg.get("covariates_file")
        )

        pav = pd.read_csv(pav_path, index_col=0)
        y = load_phenotype(str(pheno_path), trait=trait)

        covariates = None
        if covariates_file:
            covariates = load_covariates(str(covariates_file))

        results = run_pangwas(
            pav=pav,
            y=y,
            category=category,
            method=method,
            covariates=covariates,
        )

        out_csv = out_dir / "pangwas_results.csv"
        results.to_csv(out_csv, index=False)

        self.add_output("pangenome_dir", str(out_dir))
        self.add_output("pangwas_results", str(out_csv))
        self.result.message = f"Pan-GWAS completed using {category}/{method}."


class PanGWASPlotsStep(PipelineStep):
    step_id = "pangwas_plots"
    name = "Pan-GWAS Plots"
    description = "Generate Manhattan plot for panGWAS results."

    def execute(self, context: PipelineContext) -> None:
        results_csv = (
            context.get_input("pangwas_results")
            or str(context.output_dir / "pangenome" / "pangwas_results.csv")
        )

        results_path = Path(results_csv)
        if not results_path.exists():
            raise FileNotFoundError(f"panGWAS results not found: {results_path}")

        out_dir = context.output_dir / "pangenome"
        out_dir.mkdir(parents=True, exist_ok=True)

        pangenome_cfg = context.config.get("pangenome", {})

        p_col = pangenome_cfg.get("p_col", "p_value")
        marker_col = pangenome_cfg.get("marker_col", "marker")
        alpha = float(pangenome_cfg.get("alpha", 0.05))

        out_png = out_dir / "pangwas_manhattan.png"

        plot_path = plot_pangwas_manhattan(
            results_csv=str(results_path),
            out_png=str(out_png),
            p_col=p_col,
            marker_col=marker_col,
            alpha=alpha,
        )

        self.add_output("pangenome_dir", str(out_dir))

        if plot_path:
            self.add_output("pangwas_manhattan", plot_path)

        self.result.message = "Pan-GWAS plot generated successfully."


def create_pav_matrix_step() -> PAVMatrixStep:
    return PAVMatrixStep()


def create_vcf_pav_step() -> VCFPAVStep:
    return VCFPAVStep()


def create_pangwas_step() -> PanGWASStep:
    return PanGWASStep()


def create_pangwas_plots_step() -> PanGWASPlotsStep:
    return PanGWASPlotsStep()