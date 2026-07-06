"""
PlantOmicsGWAS Core Pipeline - Linkage Disequilibrium
"""

from __future__ import annotations

from pathlib import Path

from plantvarfilter.core.context import PipelineContext
from plantvarfilter.core.step import PipelineStep
from plantvarfilter.ld_utils import LDAnalyzer


class LDAnalysisStep(PipelineStep):
    step_id = "ld_analysis"
    name = "Linkage Disequilibrium Analysis"
    description = "Run LD decay, LD heatmap and diversity analysis."

    def execute(self, context: PipelineContext) -> None:
        plink_prefix = (
            context.get_input("plink_prefix")
            or context.get_input("bed_prefix")
        )

        vcf = (
            context.get_input("filtered_vcf")
            or context.get_input("vcf")
        )

        if plink_prefix:
            plink_prefix = str(plink_prefix).replace(".bed", "")
            for ext in [".bed", ".bim", ".fam"]:
                if not Path(plink_prefix + ext).exists():
                    raise FileNotFoundError(plink_prefix + ext)

        if vcf and not Path(vcf).exists():
            raise FileNotFoundError(f"VCF not found: {vcf}")

        if not plink_prefix and not vcf:
            raise ValueError("Either plink_prefix/bed_prefix or filtered_vcf/vcf must be provided.")

        out_dir = context.output_dir / "ld"
        out_dir.mkdir(parents=True, exist_ok=True)

        cfg = context.config.get("ld", {})

        region = cfg.get("region")
        keep_samples = cfg.get("keep_samples")
        window_kb = int(cfg.get("window_kb", 500))
        window_snps = int(cfg.get("window_snps", 5000))
        max_dist_kb = int(cfg.get("max_dist_kb", 1000))
        min_r2 = float(cfg.get("min_r2", 0.1))
        max_snps_for_plot = int(cfg.get("max_snps_for_plot", 500))

        do_decay = bool(cfg.get("do_decay", True))
        do_heatmap = bool(cfg.get("do_heatmap", True))
        do_diversity = bool(cfg.get("do_diversity", True))

        analyzer = LDAnalyzer()
        all_outputs = {}

        if do_decay:
            decay = analyzer.ld_decay(
                out_prefix=str(out_dir / "ld_decay"),
                bfile_prefix=plink_prefix,
                vcf_path=vcf,
                window_kb=window_kb,
                window_snps=window_snps,
                max_dist_kb=max_dist_kb,
                min_r2=min_r2,
                region=region,
                keep_samples=keep_samples,
            )
            all_outputs.update({f"decay_{k}": v for k, v in decay.items()})

        if do_heatmap:
            heat = analyzer.ld_heatmap(
                out_prefix=str(out_dir / "ld_heatmap"),
                bfile_prefix=plink_prefix,
                vcf_path=vcf,
                region=region,
                window_snps=window_snps,
                min_r2=min_r2,
                keep_samples=keep_samples,
                max_snps_for_plot=max_snps_for_plot,
            )
            all_outputs.update({f"heatmap_{k}": v for k, v in heat.items()})

        if do_diversity:
            diversity = analyzer.diversity(
                out_prefix=str(out_dir / "diversity"),
                bfile_prefix=plink_prefix,
                vcf_path=vcf,
                keep_samples=keep_samples,
                region=region,
            )
            all_outputs.update({f"diversity_{k}": v for k, v in diversity.items()})

        self.add_output("ld_dir", str(out_dir))

        for key, value in all_outputs.items():
            self.add_output(key, value)

        self.result.message = "LD analysis completed successfully."


def create_step() -> LDAnalysisStep:
    return LDAnalysisStep()