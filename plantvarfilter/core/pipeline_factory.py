"""
PlantOmicsGWAS Core Pipeline Factory

Converts workflow step IDs from config into executable PipelineStep objects.

IMPORTANT (fixed):
Each step is imported lazily, only when actually instantiated.
This means importing this module (and therefore `plantvarfilter.core`,
and therefore `plantvarfilter.compute.compute_cli`) never requires
optional heavy dependencies like fastlmm, pysnptools, xgboost, or cyvcf2
unless a step that needs them is actually enabled in the user's config.
"""

from __future__ import annotations

from importlib import import_module
from typing import Dict, List, Tuple

from plantvarfilter.core.step import PipelineStep


# step_id -> (module_path, class_name)
# NOTE: nothing here is imported yet. This is just a lookup table.
_STEP_IMPORT_MAP: Dict[str, Tuple[str, str]] = {
    "reference_indexing": ("plantvarfilter.core.pipelines.reference", "ReferenceIndexingStep"),
    "fastq_qc": ("plantvarfilter.core.pipelines.qc", "FastqQCStep"),
    "alignment": ("plantvarfilter.core.pipelines.alignment", "AlignmentStep"),
    "bam_processing": ("plantvarfilter.core.pipelines.bam", "BamProcessingStep"),
    "variant_calling": ("plantvarfilter.core.pipelines.variants", "VariantCallingStep"),
    "bcftools_processing": ("plantvarfilter.core.pipelines.variants", "BCFtoolsProcessingStep"),
    "vcf_quality": ("plantvarfilter.core.pipelines.vcf_quality", "VCFQualityStep"),
    "annotation": ("plantvarfilter.core.pipelines.annotation", "AnnotationStep"),
    "gwas": ("plantvarfilter.core.pipelines.gwas", "GWASStep"),
    "batch_gwas": ("plantvarfilter.core.pipelines.gwas", "BatchGWASStep"),
    "genomic_prediction": ("plantvarfilter.core.pipelines.prediction", "GenomicPredictionStep"),
    "ld_analysis": ("plantvarfilter.core.pipelines.ld", "LDAnalysisStep"),
    "plots": ("plantvarfilter.core.pipelines.plots", "PlottingStep"),
    "pav_matrix": ("plantvarfilter.core.pipelines.pangenome", "PAVMatrixStep"),
    "vcf_pav": ("plantvarfilter.core.pipelines.pangenome", "VCFPAVStep"),
    "pangwas": ("plantvarfilter.core.pipelines.pangenome", "PanGWASStep"),
    "pangwas_plots": ("plantvarfilter.core.pipelines.pangenome", "PanGWASPlotsStep"),
}

# Friendly hint: which pip extra fixes a missing dependency for each step.
_STEP_EXTRA_HINT: Dict[str, str] = {
    "gwas": "gwas",
    "batch_gwas": "gwas",
    "genomic_prediction": "all",
    "plots": "gwas",
    "vcf_pav": "all",
    "pav_matrix": "all",
    "pangwas": "all",
    "pangwas_plots": "all",
}


def create_step(step_id: str) -> PipelineStep:
    if step_id not in _STEP_IMPORT_MAP:
        available = ", ".join(sorted(_STEP_IMPORT_MAP.keys()))
        raise KeyError(
            f"Unknown pipeline step: {step_id}. "
            f"Available steps: {available}"
        )

    module_path, class_name = _STEP_IMPORT_MAP[step_id]

    try:
        module = import_module(module_path)
        step_cls = getattr(module, class_name)
    except ImportError as exc:
        extra = _STEP_EXTRA_HINT.get(step_id, "all")
        raise ImportError(
            f"\nStep '{step_id}' needs an optional dependency that is not installed.\n"
            f"Fix with:\n"
            f"    pip install \"plantomicsgwas[{extra}]\"\n"
            f"(original error: {exc})\n"
        ) from exc

    return step_cls()


def create_steps(step_ids: List[str]) -> List[PipelineStep]:
    return [create_step(step_id) for step_id in step_ids]


def list_available_steps() -> List[str]:
    return sorted(_STEP_IMPORT_MAP.keys())


def step_is_installed(step_id: str) -> bool:
    """
    Check whether a step's dependencies are importable, WITHOUT
    instantiating the step or raising an exception.
    """
    if step_id not in _STEP_IMPORT_MAP:
        return False

    module_path, _ = _STEP_IMPORT_MAP[step_id]

    try:
        import_module(module_path)
        return True
    except ImportError:
        return False