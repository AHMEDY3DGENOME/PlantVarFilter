"""
PlantOmicsGWAS Core Pipeline - VCF Quality Assessment
"""

from __future__ import annotations

import json
from pathlib import Path

from plantvarfilter.core.context import PipelineContext
from plantvarfilter.core.step import PipelineStep
from plantvarfilter.vcf_quality import VCFQualityChecker


class VCFQualityStep(PipelineStep):
    step_id = "vcf_quality"
    name = "VCF Quality Assessment"
    description = "Evaluate VCF quality and generate a QC report."

    def execute(self, context: PipelineContext) -> None:
        input_vcf = (
            context.get_input("filtered_vcf")
            or context.get_input("vcf")
            or str(context.output_dir / "variants" / "filtered_variants.vcf.gz")
        )

        input_vcf_path = Path(input_vcf)

        if not input_vcf_path.exists():
            raise FileNotFoundError(f"VCF not found: {input_vcf_path}")

        out_dir = context.output_dir / "quality"
        out_dir.mkdir(parents=True, exist_ok=True)

        qc_cfg = context.config.get("vcf_quality", {})

        max_sites_scan = int(qc_cfg.get("max_sites_scan", 200_000))
        min_sites_required = int(qc_cfg.get("min_sites_required", 5_000))

        checker = VCFQualityChecker(
            max_sites_scan=max_sites_scan,
            min_sites_required=min_sites_required,
        )

        report = checker.evaluate(
            str(input_vcf_path),
            log_fn=print,
        )

        report_txt = out_dir / "vcf_quality_report.txt"
        report_json = out_dir / "vcf_quality_report.json"

        with open(report_txt, "w", encoding="utf-8") as fw:
            fw.write(
                checker.to_text(
                    report,
                    str(input_vcf_path),
                )
            )

        with open(report_json, "w", encoding="utf-8") as fw:
            json.dump(
                {
                    "score": report.score,
                    "verdict": report.verdict,
                    "metrics": report.metrics,
                    "recommendations": report.recommendations,
                    "hard_fail_reasons": report.hard_fail_reasons,
                    "data_type": report.data_type,
                    "distributions_available": bool(report.dists),
                },
                fw,
                indent=2,
            )

        self.add_output("quality_dir", str(out_dir))
        self.add_output("vcf_quality_report", str(report_txt))
        self.add_output("vcf_quality_json", str(report_json))
        self.add_output("vcf_quality_score", str(report.score))
        self.add_output("vcf_quality_verdict", report.verdict)

        self.result.message = (
            f"VCF Quality Score = {report.score:.1f} ({report.verdict})"
        )


def create_step() -> VCFQualityStep:
    return VCFQualityStep()