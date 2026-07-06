
from pathlib import Path
from typing import Any, Dict, List

import yaml

from plantvarfilter.hpc.array_resolver import apply_array_task_to_config

from .workflow_registry import WORKFLOW_STEPS, validate_workflow


REQUIRED_TOP_LEVEL_KEYS = ["project", "input", "output", "compute", "steps"]


def load_config(config_path: str) -> Dict[str, Any]:
    path = Path(config_path)

    if not path.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")

    with path.open("r", encoding="utf-8") as f:
        config = yaml.safe_load(f)

    if not isinstance(config, dict):
        raise ValueError("Config file must contain a valid YAML dictionary.")

    config = apply_array_task_to_config(config)

    validate_config(config)
    return config


def validate_config(config: Dict[str, Any]) -> bool:
    missing = [key for key in REQUIRED_TOP_LEVEL_KEYS if key not in config]
    if missing:
        raise ValueError(f"Missing required config sections: {', '.join(missing)}")

    steps = config.get("steps", {})
    if not isinstance(steps, dict):
        raise ValueError("'steps' section must be a dictionary.")

    enabled_steps = get_enabled_steps(config)
    validate_workflow(enabled_steps)

    produced_outputs = set()

    for step_id in enabled_steps:
        step = WORKFLOW_STEPS[step_id]

        for dependency_id in step.dependencies:
            if dependency_id not in enabled_steps:
                raise ValueError(
                    f"Step '{step_id}' depends on '{dependency_id}', "
                    f"but '{dependency_id}' is not enabled."
                )

        available_inputs = produced_outputs.copy()

        for required_input in step.required_inputs:
            if has_config_input(config, required_input):
                continue

            if required_input in available_inputs:
                continue

            raise ValueError(
                f"Step '{step_id}' requires input '{required_input}', "
                f"but it was not found in config and is not produced by previous steps."
            )

        produced_outputs.update(step.outputs)

    return True


def get_enabled_steps(config: Dict[str, Any]) -> List[str]:
    steps = config.get("steps", {})
    return [step_id for step_id, enabled in steps.items() if enabled]


def has_config_input(config: Dict[str, Any], input_name: str) -> bool:
    input_section = config.get("input", {})
    output_section = config.get("output", {})

    aliases = {
        "reference_fasta": ["reference_fasta", "reference", "genome"],
        "fastq_dir": ["fastq_dir", "reads_dir", "raw_reads"],
        "reads1": ["reads1", "r1", "fastq_r1"],
        "reads2": ["reads2", "r2", "fastq_r2"],

        "bam_file": ["bam_file", "bam", "alignment_bam", "final_bam"],
        "final_bam": ["final_bam", "sorted_bam", "bam_file", "bam"],
        "raw_vcf": ["raw_vcf", "vcf"],
        "filtered_vcf": ["filtered_vcf", "vcf"],

        "phenotype_file": ["phenotype_file", "phenotype", "pheno"],
        "covariates_file": ["covariates_file", "covar"],

        "plink_prefix": ["plink_prefix", "bed_prefix"],
        "bed_prefix": ["bed_prefix", "plink_prefix"],

        "gff_file": ["gff_file", "gtf_file", "annotation_file"],
        "gwas_results": ["gwas_results", "gwas_file"],

        "gff_dir": ["gff_dir", "pav_input"],
        "pav_input": ["pav_input", "gff_dir"],
        "pav_matrix": ["pav_matrix"],
        "pangwas_results": ["pangwas_results"],

        "chromosomes": ["chromosomes"],
        "samples": ["samples"],
    }

    possible_names = aliases.get(input_name, [input_name])

    for name in possible_names:
        if input_section.get(name):
            return True
        if output_section.get(name):
            return True

    return False