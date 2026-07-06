"""
PlantOmicsGWAS Core Pipeline - Plots and Statistics
"""

from __future__ import annotations

from pathlib import Path

from pysnptools.snpreader import Bed, Pheno

from plantvarfilter.core.context import PipelineContext
from plantvarfilter.core.step import PipelineStep
from plantvarfilter.pipeline_plots import Plot


class PlottingStep(PipelineStep):
    step_id = "plots"
    name = "Pipeline Visualization"
    description = "Generate phenotype and genotype statistics plots."

    def execute(self, context: PipelineContext) -> None:
        out_dir = context.plots_dir
        out_dir.mkdir(parents=True, exist_ok=True)

        phenotype = (
            context.get_input("phenotype_file")
            or context.get_input("phenotype")
            or context.get_input("pheno")
        )

        bed_prefix = (
            context.get_input("plink_prefix")
            or context.get_input("bed_prefix")
        )

        if bed_prefix:
            bed_prefix = str(bed_prefix).replace(".bed", "")

        plotter = Plot()
        generated = 0

        if phenotype:
            phenotype_path = Path(phenotype)

            if phenotype_path.exists():
                pheno_pdf = out_dir / "phenotype_statistics.pdf"

                plotter.plot_pheno_statistics(
                    pheno_file=str(phenotype_path),
                    pheno_stats_name=str(pheno_pdf),
                )

                self.add_output("phenotype_statistics_pdf", str(pheno_pdf))
                generated += 1

        if bed_prefix and phenotype:
            phenotype_path = Path(phenotype)

            bed_path = Path(bed_prefix + ".bed")
            bim_path = Path(bed_prefix + ".bim")
            fam_path = Path(bed_prefix + ".fam")

            if bed_path.exists() and bim_path.exists() and fam_path.exists() and phenotype_path.exists():
                geno_pdf = out_dir / "genotype_statistics.pdf"

                bed = Bed(bed_prefix, count_A1=False)
                pheno = Pheno(str(phenotype_path))

                plotter.plot_geno_statistics(
                    bed=bed,
                    pheno=pheno,
                    geno_stats_name=str(geno_pdf),
                )

                self.add_output("genotype_statistics_pdf", str(geno_pdf))
                generated += 1

        self.add_output("plots_dir", str(out_dir))
        self.result.message = f"Pipeline plots generated successfully ({generated} reports)."


def create_step() -> PlottingStep:
    return PlottingStep()