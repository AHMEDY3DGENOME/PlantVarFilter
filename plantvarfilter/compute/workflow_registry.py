"""
PlantOmicsGWAS Compute Engine - Workflow Registry
"""

from dataclasses import dataclass
from typing import Dict, List, Optional


@dataclass(frozen=True)
class WorkflowStep:
    step_id: str
    name: str
    module: str
    category: str
    description: str
    dependencies: List[str]
    required_inputs: List[str]
    optional_inputs: List[str]
    outputs: List[str]
    hpc_ready: bool = True
    enabled_by_default: bool = True


WORKFLOW_STEPS: Dict[str, WorkflowStep] = {
    "reference_indexing": WorkflowStep(
        "reference_indexing", "Reference Genome Indexing",
        "plantvarfilter.preanalysis.reference_manager", "preprocessing",
        "Prepare and index the reference genome.",
        [], ["reference_fasta"], ["index_prefix"], ["reference_index"],
    ),

    "fastq_qc": WorkflowStep(
        "fastq_qc", "FASTQ Quality Control",
        "plantvarfilter.preanalysis.fastq_qc", "preprocessing",
        "Perform quality control on raw FASTQ reads.",
        [], ["fastq_dir"], ["reads1", "reads2", "threads"], ["qc_reports"],
    ),

    "alignment": WorkflowStep(
        "alignment", "Read Alignment",
        "plantvarfilter.preanalysis.aligner", "preprocessing",
        "Align sequencing reads to the reference genome.",
        ["reference_indexing", "fastq_qc"],
        ["fastq_dir", "reference_fasta"],
        ["reads1", "reads2", "platform", "threads"],
        ["alignment_bam", "bam_file"],
    ),

    "bam_processing": WorkflowStep(
        "bam_processing", "BAM Processing",
        "plantvarfilter.samtools_utils", "variant_processing",
        "Sort, fixmate, mark duplicates, index, and generate BAM statistics.",
        ["alignment"],
        ["bam_file"],
        ["bam_dir", "threads", "remove_duplicates"],
        ["final_bam", "bam_index"],
    ),

    "variant_calling": WorkflowStep(
        "variant_calling", "Variant Calling",
        "plantvarfilter.variant_caller_utils", "variant_processing",
        "Call genomic variants from BAM files.",
        ["bam_processing"],
        ["final_bam", "reference_fasta"],
        ["regions_bed", "min_baseq", "min_mapq", "ploidy"],
        ["raw_vcf"],
    ),

    "bcftools_processing": WorkflowStep(
        "bcftools_processing", "BCFtools Variant Processing",
        "plantvarfilter.bcftools_utils", "variant_processing",
        "Normalize, sort, filter, compress, and index VCF files.",
        ["variant_calling"],
        ["raw_vcf"],
        ["reference_fasta", "regions_bed", "filter_expr"],
        ["filtered_vcf"],
    ),

    "vcf_quality": WorkflowStep(
        "vcf_quality", "VCF Quality Assessment",
        "plantvarfilter.vcf_quality", "quality_control",
        "Evaluate VCF quality metrics.",
        ["bcftools_processing"],
        ["filtered_vcf"],
        ["max_sites_scan", "min_sites_required"],
        ["vcf_quality_report"],
    ),

    "annotation": WorkflowStep(
        "annotation", "GWAS Annotation",
        "plantvarfilter.annotation_utils", "annotation",
        "Annotate GWAS results with nearest genes from GTF/GFF.",
        ["gwas"],
        ["gff_file", "gwas_results"],
        ["window_kb", "chr_col", "pos_col", "region"],
        ["annotated_gwas"],
    ),

    "gwas": WorkflowStep(
        "gwas", "GWAS Analysis",
        "plantvarfilter.gwas_pipeline", "association_analysis",
        "Run genome-wide association analysis.",
        [],
        ["plink_prefix", "phenotype_file"],
        ["filtered_vcf", "covariates_file", "method", "threads"],
        ["gwas_results", "manhattan_plot", "qq_plot"],
    ),

    "batch_gwas": WorkflowStep(
        "batch_gwas", "Batch GWAS",
        "plantvarfilter.batch_gwas", "association_analysis",
        "Run GWAS across multiple traits.",
        ["bcftools_processing"],
        ["filtered_vcf", "phenotype_file"],
        ["trait_columns", "covariates_file", "threads"],
        ["batch_gwas_results"],
    ),

    "genomic_prediction": WorkflowStep(
        "genomic_prediction", "Genomic Prediction",
        "plantvarfilter.genomic_prediction_pipeline", "machine_learning",
        "Run genomic prediction models.",
        ["gwas"],
        ["plink_prefix", "phenotype_file"],
        ["algorithm", "test_size", "estimators", "model_nr", "max_depth"],
        ["prediction_results", "model_metrics"],
    ),

    "ld_analysis": WorkflowStep(
        "ld_analysis", "LD Analysis",
        "plantvarfilter.ld_utils", "association_analysis",
        "Calculate LD decay, LD heatmap, and diversity metrics.",
        ["bcftools_processing"],
        ["filtered_vcf"],
        ["plink_prefix", "region", "window_kb", "window_snps", "min_r2"],
        ["ld_results"],
    ),

    "plots": WorkflowStep(
        "plots", "Pipeline Visualization",
        "plantvarfilter.pipeline_plots", "visualization",
        "Generate phenotype and genotype statistics plots.",
        [],
        [],
        ["plink_prefix", "phenotype_file"],
        ["plots"],
    ),

    "pav_matrix": WorkflowStep(
        "pav_matrix", "PAV Matrix Construction",
        "plantvarfilter.pangenome_module.pav_matrix_builder", "pangenome",
        "Build gene presence/absence matrix from GFF/GTF files.",
        [],
        ["gff_dir"],
        ["sample_list", "gene_id_field", "max_genes"],
        ["pav_matrix"],
    ),

    "vcf_pav": WorkflowStep(
        "vcf_pav", "VCF to PAV Matrix",
        "plantvarfilter.pangenome_module.vcf_pav_builder", "pangenome",
        "Build variant PAV matrix from multi-sample VCF.",
        ["bcftools_processing"],
        ["filtered_vcf"],
        ["max_variants"],
        ["pav_matrix"],
    ),

    "pangwas": WorkflowStep(
        "pangwas", "Pan-GWAS Analysis",
        "plantvarfilter.pangenome_module.pangwas_pipeline", "pangenome",
        "Run pan-genome-wide association analysis.",
        ["pav_matrix"],
        ["pav_matrix", "phenotype_file"],
        ["trait", "category", "method", "covariates_file"],
        ["pangwas_results"],
    ),

    "pangwas_plots": WorkflowStep(
        "pangwas_plots", "Pan-GWAS Plots",
        "plantvarfilter.pangenome_module.pangwas_plots", "visualization",
        "Generate Manhattan plot for panGWAS.",
        ["pangwas"],
        ["pangwas_results"],
        ["p_col", "marker_col", "alpha"],
        ["pangwas_manhattan"],
    ),
}

DEFAULT_HPC_WORKFLOW: List[str] = [
    "reference_indexing",
    "fastq_qc",
    "alignment",
    "bam_processing",
    "variant_calling",
    "bcftools_processing",
    "vcf_quality",
    "gwas",
    "annotation",
    "batch_gwas",
    "genomic_prediction",
    "ld_analysis",
    "plots",
]

PANGENOME_WORKFLOW: List[str] = [
    "pav_matrix",
    "vcf_pav",
    "pangwas",
    "pangwas_plots",
]

FULL_WORKFLOW: List[str] = DEFAULT_HPC_WORKFLOW + PANGENOME_WORKFLOW


def get_step(step_id: str) -> WorkflowStep:
    if step_id not in WORKFLOW_STEPS:
        available = ", ".join(WORKFLOW_STEPS.keys())
        raise KeyError(f"Unknown workflow step: {step_id}. Available steps: {available}")
    return WORKFLOW_STEPS[step_id]


def list_steps(category: Optional[str] = None) -> List[WorkflowStep]:
    steps = list(WORKFLOW_STEPS.values())
    if category:
        steps = [step for step in steps if step.category == category]
    return steps


def get_workflow(workflow_name: str = "full") -> List[WorkflowStep]:
    workflow_name = workflow_name.lower()

    if workflow_name == "default":
        step_ids = DEFAULT_HPC_WORKFLOW
    elif workflow_name == "pangenome":
        step_ids = PANGENOME_WORKFLOW
    elif workflow_name == "full":
        step_ids = FULL_WORKFLOW
    else:
        raise ValueError("workflow_name must be one of: default, pangenome, full")

    return [get_step(step_id) for step_id in step_ids]


def validate_workflow(step_ids: List[str]) -> bool:
    for step_id in step_ids:
        get_step(step_id)
    return True
