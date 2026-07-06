"""
PlantOmicsGWAS Core Pipeline - Genomic Prediction
"""

from __future__ import annotations

from pathlib import Path

from pysnptools.snpreader import Bed, Pheno

from plantvarfilter.core.context import PipelineContext
from plantvarfilter.core.step import PipelineStep
from plantvarfilter.genomic_prediction_pipeline import GenomicPrediction
from plantvarfilter.helpers import HELPERS


class GenomicPredictionStep(PipelineStep):
    step_id = "genomic_prediction"
    name = "Genomic Prediction"
    description = "Run genomic prediction models."

    def execute(self, context: PipelineContext) -> None:
        bed_prefix = context.get_input("plink_prefix") or context.get_input("bed_prefix")
        phenotype = context.get_input("phenotype_file") or context.get_input("phenotype")

        if not bed_prefix:
            raise ValueError("Missing PLINK BED prefix")

        if not phenotype:
            raise ValueError("Missing phenotype file")

        bed_file = str(bed_prefix) if str(bed_prefix).endswith(".bed") else str(bed_prefix) + ".bed"

        for ext in [".bed", ".bim", ".fam"]:
            path = str(bed_prefix).replace(".bed", "") + ext
            if not Path(path).exists():
                raise FileNotFoundError(path)

        if not Path(phenotype).exists():
            raise FileNotFoundError(phenotype)

        out_dir = context.output_dir / "prediction"
        out_dir.mkdir(parents=True, exist_ok=True)

        cfg = context.config.get("prediction", {})

        algorithm = cfg.get("algorithm", "XGBoost (AI)")
        test_size = float(cfg.get("test_size", 0.3))
        estimators = int(cfg.get("estimators", 100))
        model_nr = int(cfg.get("model_nr", 5))
        max_depth = int(cfg.get("max_depth", 3))
        nr_jobs = int(cfg.get("nr_jobs", context.threads))
        alpha = float(cfg.get("alpha", 1.0))

        out_csv = str(out_dir / "genomic_prediction_results.csv")
        bland_png = str(out_dir / "Bland_Altman_plot.png")
        scatter_png = str(out_dir / "GP_scatter_plot.png")

        helper = HELPERS()
        chrom_mapping = helper.replace_with_integers(str(bed_prefix).replace(".bed", ".bim"))

        bed = Bed(str(bed_prefix).replace(".bed", ""), count_A1=False, chrom_map=chrom_mapping)
        pheno = Pheno(str(phenotype))
        bed_fixed = bed

        predictor = GenomicPrediction()

        if algorithm == "GP_LMM":
            gp_df = predictor.run_lmm_gp(
                bed_fixed, pheno, out_csv, model_nr, print, bed_file, chrom_mapping
            )

        elif algorithm == "Random Forest (AI)":
            gp_df = predictor.run_gp_rf(
                bed_fixed, pheno, bed_file, test_size, estimators,
                out_csv, chrom_mapping, print, model_nr, nr_jobs
            )

        elif algorithm == "Ridge Regression":
            gp_df = predictor.run_gp_ridge(
                bed_fixed, pheno, bed_file, test_size, alpha,
                out_csv, chrom_mapping, print, model_nr
            )

        else:
            gp_df = predictor.run_gp_xg(
                bed_fixed, pheno, bed_file, test_size, estimators,
                out_csv, chrom_mapping, print, model_nr, max_depth, nr_jobs
            )

        predictor.plot_gp(gp_df, bland_png, algorithm)
        predictor.plot_gp_scatter(gp_df, scatter_png, algorithm)

        self.add_output("prediction_dir", str(out_dir))
        self.add_output("genomic_prediction_results", out_csv)
        self.add_output("bland_altman_plot", bland_png)
        self.add_output("prediction_scatter_plot", scatter_png)
        self.add_output("prediction_algorithm", algorithm)

        self.result.message = f"Genomic prediction completed using {algorithm}."


def create_step() -> GenomicPredictionStep:
    return GenomicPredictionStep()