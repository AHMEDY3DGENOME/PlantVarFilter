"""
PlantOmicsGWAS Core Pipeline Factory

Converts workflow step IDs from config into executable PipelineStep objects.
"""

from __future__ import annotations

from typing import Dict, List, Type

from plantvarfilter.core.step import PipelineStep
from plantvarfilter.core.pipelines.reference import ReferenceIndexingStep
from plantvarfilter.core.pipelines.qc import FastqQCStep
from plantvarfilter.core.pipelines.alignment import AlignmentStep
from plantvarfilter.core.pipelines.bam import BamProcessingStep
from plantvarfilter.core.pipelines.variants import (
    VariantCallingStep,
    BCFtoolsProcessingStep,
)
from plantvarfilter.core.pipelines.vcf_quality import VCFQualityStep
from plantvarfilter.core.pipelines.annotation import AnnotationStep
from plantvarfilter.core.pipelines.gwas import GWASStep, BatchGWASStep
from plantvarfilter.core.pipelines.prediction import GenomicPredictionStep
from plantvarfilter.core.pipelines.ld import LDAnalysisStep
from plantvarfilter.core.pipelines.plots import PlottingStep
from plantvarfilter.core.pipelines.pangenome import (
    PAVMatrixStep,
    VCFPAVStep,
    PanGWASStep,
    PanGWASPlotsStep,
)


PIPELINE_STEP_REGISTRY: Dict[str, Type[PipelineStep]] = {
    "reference_indexing": ReferenceIndexingStep,
    "fastq_qc": FastqQCStep,
    "alignment": AlignmentStep,
    "bam_processing": BamProcessingStep,
    "variant_calling": VariantCallingStep,
    "bcftools_processing": BCFtoolsProcessingStep,
    "vcf_quality": VCFQualityStep,
    "annotation": AnnotationStep,
    "gwas": GWASStep,
    "batch_gwas": BatchGWASStep,
    "genomic_prediction": GenomicPredictionStep,
    "ld_analysis": LDAnalysisStep,
    "plots": PlottingStep,
    "pav_matrix": PAVMatrixStep,
    "vcf_pav": VCFPAVStep,
    "pangwas": PanGWASStep,
    "pangwas_plots": PanGWASPlotsStep,
}


def create_step(step_id: str) -> PipelineStep:
    step_cls = PIPELINE_STEP_REGISTRY.get(step_id)

    if step_cls is None:
        available = ", ".join(sorted(PIPELINE_STEP_REGISTRY.keys()))
        raise KeyError(
            f"Unknown pipeline step: {step_id}. "
            f"Available steps: {available}"
        )

    return step_cls()


def create_steps(step_ids: List[str]) -> List[PipelineStep]:
    return [create_step(step_id) for step_id in step_ids]


def list_available_steps() -> List[str]:
    return sorted(PIPELINE_STEP_REGISTRY.keys())