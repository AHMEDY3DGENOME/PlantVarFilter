"""
PlantOmicsGWAS Core Pipeline - Variant/GWAS Annotation
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from plantvarfilter.annotation_utils import Annotator
from plantvarfilter.core.context import PipelineContext
from plantvarfilter.core.step import PipelineStep


class AnnotationStep(PipelineStep):
    step_id = "annotation"
    name = "Variant Annotation"
    description = "Annotate GWAS results with nearest genes from GTF/GFF annotation."

    def execute(self, context: PipelineContext) -> None:
        gff_file = (
            context.get_input("gff_file")
            or context.get_input("gtf_file")
            or context.get_input("annotation_file")
        )

        gwas_file = (
            context.get_input("gwas_results")
            or context.get_input("gwas_file")
            or str(context.output_dir / "gwas" / "gwas_results.csv")
        )

        if not gff_file:
            raise ValueError("Missing required input: gff_file / gtf_file / annotation_file")

        gff_path = Path(gff_file)
        gwas_path = Path(gwas_file)

        if not gff_path.exists():
            raise FileNotFoundError(f"Annotation file not found: {gff_path}")

        if not gwas_path.exists():
            raise FileNotFoundError(f"GWAS results file not found: {gwas_path}")

        out_dir = context.output_dir / "annotation"
        out_dir.mkdir(parents=True, exist_ok=True)

        annotation_cfg = context.config.get("annotation", {})

        window_kb = int(annotation_cfg.get("window_kb", 50))
        chr_col = annotation_cfg.get("chr_col", "Chr")
        pos_col = annotation_cfg.get("pos_col", "ChrPos")
        region = annotation_cfg.get("region")

        gwas_df = pd.read_csv(gwas_path)

        if region:
            gwas_df = Annotator.filter_region(
                gwas_df,
                region=region,
                chr_col=chr_col,
                pos_col=pos_col,
            )

        annotator = Annotator()
        annotator.load_gtf_or_gff(str(gff_path))
        annotator.build_index()

        out_csv = out_dir / "gwas_annotated.csv"

        annotated_path = annotator.annotate_and_save(
            gwas_df=gwas_df,
            out_csv=str(out_csv),
            window_kb=window_kb,
            chr_col=chr_col,
            pos_col=pos_col,
        )

        self.add_output("annotation_dir", str(out_dir))
        self.add_output("annotated_gwas", annotated_path)

        self.result.message = "GWAS annotation completed successfully."


def create_step() -> AnnotationStep:
    return AnnotationStep()