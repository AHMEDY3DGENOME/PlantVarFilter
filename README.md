# PlantOmicsGwas: An Integrated GWAS, Genomic Prediction, and Pangenome-Based Association Pipeline for Plant Genomes


[![Python](https://img.shields.io/badge/python-3.10%2B-blue)]()
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)]()
[![Status](https://img.shields.io/badge/status-Early%20Demonstration%20Release-orange)]()
[![Version](https://img.shields.io/badge/version-1.0.2-blue)]()

---

## Developers & Contributors

| Developer | Role | Affiliation                                                                             |
|---|---|-----------------------------------------------------------------------------------------|
| **Ahmed Yassin** | Computational Biologist, PhD Candidate | --                                                                                      |
| **Falak Sher Khan, PhD** | Computational Biologist | Ye-Lab, Peking University Institute of Advanced Agricultural Sciences (PKU-IAAS), China |

**Contact:** ahmedyassin300@outlook.com · falakmahmand@gmail.com

---

> ### ⚠️ Early Demonstration Release
>
> **PlantOmicsGwas** is a **research software package under active development**, released publicly for evaluation, testing, and feedback from the plant genomics community. This README documents, with full transparency, exactly which components have been verified end-to-end on a real HPC (SLURM) cluster with real biological data, and which components are still undergoing testing. See the **[Testing Status](#-testing-status)** section below for the complete, up-to-date picture.
>
> Following the publication of our accompanying research manuscript, we will release a complete production version including full documentation, finalized licensing, and the remaining planned modules.

---

## Table of Contents

1. [Overview](#overview)
2. [Key Features](#key-features)
3. [Testing Status](#-testing-status)
4. [Requirements](#requirements)
5. [Installation](#installation)
6. [Quick Start](#quick-start)
7. [Future Plan](#future-plan)
8. [Test: Real GWAS Run Walkthrough](#test-real-gwas-run-walkthrough)
9. [Two Paths to GWAS: Linear Reference vs. Pangenome](#two-paths-to-gwas-linear-reference-vs-pangenome)
10. [Running on HPC Clusters (SLURM/PBS/LSF)](#running-on-hpc-clusters-slurmpbslsf)
11. [Module Reference & Example Configs](#module-reference--example-configs)
12. [Project Structure](#project-structure)
13. [Troubleshooting / FAQ](#troubleshooting--faq)
14. [Citation](#citation)
15. [License](#license)

---

## Overview

**PlantOmicsGwas** is an integrated software platform for genome-wide association studies (GWAS), genomic prediction, and pangenome-based association analysis in plant genomics. It is designed for researchers working with large-scale variant data who need to move from raw sequencing reads to association results, using either a **classical linear reference genome** or a **graph-based pangenome reference**.

The platform provides three ways to work:
- A **desktop GUI** (`plantomicsgwas-gui`) for interactive, point-and-click analysis
- A **command-line interface** (`plantomicsgwas`) for scripting individual analysis steps
- A **headless Compute Engine** (`plantomicsgwas-compute`) for running full multi-step workflows from a single YAML configuration file — including **native support for HPC job schedulers (SLURM, PBS, LSF)**

---

## Key Features

PlantOmicsGwas provides the following modules, accessible through the desktop GUI, command-line interface, and headless Compute Engine:

- **Reference Manager**
- **Pangenome Builder**
- **panGWAS**
- **FASTQ QC**
- **Alignment**
- **Preprocess (samtools)**
- **Variant Calling (BAM/VCF)**
- **Preprocess (bcftools)**
- **Check VCF File**
- **Convert to PLINK**
- **LD Analysis**
- **GWAS Analysis**
- **PCA / Kinship**
- **Genomic Prediction**
- **Batch GWAS**

> **Note on HPC readiness:** All of the modules listed above are currently qualified and functional for execution on High-Performance Computing (HPC) clusters through the Compute Engine and its job-scheduler integration (SLURM, PBS, LSF). This does not preclude further performance benchmarking and optimization, which will continue as ongoing work in collaboration with lab colleagues.

---

##  Testing Status

This table reflects the verified state of the software as of this release. All 17 Compute Engine workflow steps, all three job schedulers (SLURM, PBS, LSF), and the HPC array-job and dependency-chain features have been tested end-to-end by running real jobs through the full pipeline — job submission through the scheduler, execution on a worker node, and inspection of the produced output files.

### Compute Engine / HPC steps

| # | Module | Status | Notes |
|---|---|---|---|
| 1 | Reference Genome Indexing | ✅ Verified on HPC | Builds `.fai`, `.mmi` (minimap2), and bowtie2 index. `.dict` creation is skipped gracefully if Picard is not installed (non-critical). |
| 2 | FASTQ Quality Control | ✅ Verified on HPC | Pure-Python streaming QC; produces read-length, GC%, and per-cycle quality reports. |
| 3 | Read Alignment | ✅ Verified on HPC | Tested with bowtie2 against a real 505 Mbp plant genome; 100% alignment rate on test reads. |
| 4 | BAM Processing | ✅ Verified on HPC | Full samtools sort → fixmate → markdup → index → stats chain. |
| 5 | Variant Calling | ✅ Verified on HPC | `bcftools mpileup \| bcftools call -mv`, produces a valid indexed VCF. |
| 6 | BCFtools Variant Processing | ✅ Verified on HPC | Normalization (`bcftools norm -m -both`), sorting, ID annotation, stats. |
| 7 | GWAS Analysis | ✅ Verified on HPC | FaST-LMM algorithm tested end-to-end (57s runtime, 53,184 variants, valid Manhattan/QQ plots produced). |
| 8 | VCF Quality Assessment | ✅ Verified on HPC | Pure-Python, no external tool dependencies. |
| 9 | LD Analysis | ✅ Verified on HPC | PLINK-based LD decay/heatmap/diversity metrics. |
| 10 | Genomic Prediction | ✅ Verified on HPC | Requires `fastlmm`, `xgboost`, and `seaborn` together (all imported unconditionally regardless of chosen algorithm). |
| 11 | GWAS Annotation | ✅ Verified on HPC | Nearest-gene annotation against a GFF/GTF file, using GWAS output columns (`Chr`, `ChrPos`). |
| 12 | Batch GWAS | ✅ Verified on HPC | A prior release bug (the `--algo` argument was not forwarded to the underlying CLI) has been fixed and re-verified. |
| 13 | PAV Matrix Construction | ✅ Verified on HPC | Builds a gene presence/absence matrix from per-sample GFF files. |
| 14 | VCF-to-PAV Matrix | ✅ Verified on HPC | Builds a presence/absence matrix directly from a multi-sample VCF using `cyvcf2`. |
| 15 | Pan-GWAS Analysis | ✅ Verified on HPC | Pure NumPy/SciPy implementation (t-test, Wilcoxon, Fisher, GLM, kinship-corrected LMM). |
| 16 | Pan-GWAS Plots | ✅ Verified on HPC | Manhattan-style plotting of Pan-GWAS results. |

### GUI-only tools (not yet integrated into the Compute Engine / HPC workflow)

| Tool | Status |
|---|---|
| Pangenome Builder (minigraph + vg) | Available in the GUI/CLI only; not registered as an HPC-schedulable step in this release. |
| Convert to PLINK | GUI-only utility; uses PLINK2 (not bundled — must be installed separately). |
| PCA / Kinship | GUI-only utility; functional and tested manually, not yet an HPC step. |

### HPC layer specifics

| Feature | Status |
|---|---|
| SLURM scheduler (detection, job generation, submission, status) | ✅ Verified extensively across multiple real jobs |
| PBS scheduler | ✅ Verified |
| LSF scheduler | ✅ Verified. A prior release bug in the dependency-chain submission path (a literal `<` character passed as a subprocess argument instead of being interpreted as shell redirection) has been fixed and re-verified. |
| Array jobs (traits / chromosomes / samples / custom) | ✅ Verified |
| Dependency chains | ✅ Verified |
| Shared conda environment activation via `hpc.conda_env` | ✅ Verified |
| Shared reference index reuse via `reference.shared_index_dir` | ✅ Verified (680s → 0.00s on repeat runs, confirmed via checkpoint skip) |

---

## Requirements

- **OS:** Linux (primary target; Windows and macOS supported for local/GUI use)
- **Python:** 3.10, 3.11, or 3.12 
- **Disk space:** A few GB for the software and its dependencies; genome/read data storage requirements depend on your dataset (a single reference genome index for a ~500 Mbp plant genome uses ~1 GB)
- **For HPC use:** A SLURM cluster (PBS/LSF present but not yet fully verified - still working on some techniques — see Testing Status)

> **Python 3.12 note:** `fastlmm` (used by the FaST-LMM GWAS/prediction algorithm) declares official support only for Python < 3.12 in its packaging metadata. In our own testing it built and ran successfully from source on Python 3.11 within a dedicated conda environment (the default created by the installer below). If you build your own environment on Python 3.12, verify `fastlmm` installs correctly before relying on it.

---

## Installation

# install & Test video tutorial:
[Download & See the video from here](https://drive.google.com/drive/folders/1Pwus1VXtihw36gyOZm-gBU3zNUatQVv5?usp=sharing)

### Recommended: One-command setup

This installs a dedicated conda environment with every required bioinformatics tool (samtools, bcftools, bowtie2, minimap2, plink, plink2), installs PlantOmicsGwas with all optional extras, verifies everything works, and drops you directly into a ready-to-use activated shell.

```bash
# 1. Make sure pip is available and up to date
python3 -m ensurepip --upgrade
python3 -m pip install --upgrade pip

# 2. Install PlantOmicsGwas from PyPI
pip install plantomicsgwas

# 3. Set up the full environment (conda, all bioinformatics tools, all optional extras)
plantomicsgwas-setup-env
```

That's it — three simple commands, no manual conda setup, no separate tool installation. `plantomicsgwas-setup-env` will:

1. Install Miniforge (conda) if it isn't already present
2. Create a dedicated environment with samtools, bcftools, bowtie2, minimap2, plink, and plink2
3. Install PlantOmicsGwas with all optional extras (fastlmm, xgboost, pysnptools, geneview, cyvcf2)
4. Verify every tool and the CLI actually work, printing a clear report
5. Drop you into a new shell with everything already activated and ready to run

To use a custom environment name:
```bash
plantomicsgwas-setup-env my_custom_env_name
```

### Manual installation (for full control)

If you prefer to manage the environment yourself:

```bash
# 1. Install Miniforge (conda)
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh
source ~/.bashrc

# 2. Create the environment with all external tools
mamba create -n plantomicsgwas -c conda-forge -c bioconda \
    python=3.11 samtools bcftools bowtie2 minimap2 plink plink2 htslib -y

# 3. Activate it
conda activate plantomicsgwas

# 4. Install PlantOmicsGwas with all optional modules
pip install --upgrade pip
pip install "plantomicsgwas[all]"
```

### Verify installation

```bash
samtools --version
bcftools --version
bowtie2 --version
minimap2 --version
plink --version
plink2 --version
plantomicsgwas-compute list-steps
```

If all commands print version numbers and `list-steps` prints all 17 available workflow steps without errors, the installation is complete.

### Launching the software

```bash
plantomicsgwas-gui              # Desktop GUI
plantomicsgwas --help           # Command-line interface
plantomicsgwas-compute --help   # Headless Compute Engine / HPC interface
```

---

## Quick Start

The fastest way to see a result is to run GWAS directly on an existing PLINK dataset and phenotype file — no preprocessing steps required if your data is already prepared.

Create `my_gwas.yaml`:

```yaml
project:
  name: my_first_gwas
  version: "1.0"

input:
  plink_prefix: /path/to/my_data       # expects my_data.bed/.bim/.fam
  phenotype_file: /path/to/phenotype.txt

output:
  dir: ./results/gwas_run1
  logs_dir: ./results/gwas_run1/logs

compute:
  threads: 4
  memory: 4G

gwas:
  method: FaST-LMM

steps:
  reference_indexing: false
  fastq_qc: false
  alignment: false
  bam_processing: false
  variant_calling: false
  bcftools_processing: false
  vcf_quality: false
  annotation: false
  gwas: true
  batch_gwas: false
  genomic_prediction: false
  ld_analysis: false
  pav_matrix: false
  vcf_pav: false
  pangwas: false
  pangwas_plots: false
  plots: false
```

Run it locally (no HPC scheduler required):

```bash
plantomicsgwas-compute run --config my_gwas.yaml
```

This produces `gwas_results.csv`, `manhattan_plot.png`, and `qq_plot.png` in the output directory.

---

## Test: Real GWAS Run Walkthrough

This section documents an actual GWAS run performed with PlantOmicsGwas, starting from a raw multi-sample VCF (derived from a pangenome variant caller) and a phenotype file, all the way to final association results. Every command below was executed and produced real output.

### 1. Starting data

- A multi-sample VCF: `merged.vgCall.vcf`
- A phenotype file (FID, IID, trait value): `phenotype_fastlmm.txt`

### 2. Normalize and filter the VCF

Multi-allelic sites are split and only `PASS` variants are kept — this step matters: skipping it can silently reduce the number of variants that make it through the downstream PLINK conversion.

```bash
bcftools norm -m -both -Oz -o merged.norm.vcf.gz merged.vgCall.vcf
bcftools view -f PASS merged.norm.vcf.gz -Oz -o merged.filtered.vcf.gz
tabix -p vcf merged.filtered.vcf.gz
```

### 3. Convert to PLINK format

```bash
plink --vcf merged.filtered.vcf.gz \
      --make-bed \
      --allow-extra-chr \
      --double-id \
      --out plink_data
```

> **Recommended `--maf`/`--geno` for pangenome-derived VCFs:** pangenome/graph-derived
> VCFs (e.g. from `vg call`) commonly have a much higher missing-genotype rate than
> linear-reference VCFs. In our own testing, the default thresholds used by
> `plantomicsgwas vcf2bed` (`--geno 0.10`, `--maf 0.05`) removed **more than 95% of
> variants** on a pangenome dataset with an overall genotype call rate of only 33.8%.
> For pangenome-derived data, start with more permissive thresholds and tighten them
> only if needed:
> ```bash
> plantomicsgwas vcf2bed --vcf merged.filtered.vcf.gz --out plink_data \
>     --geno 0.5 --maf 0.01
> ```

### 4. Configure the GWAS run

```yaml
project:
  name: gwas_demo
  version: "1.0"

input:
  plink_prefix: /data/plink_data
  phenotype_file: /data/phenotype_fastlmm.txt

output:
  dir: /data/results/gwas_demo
  logs_dir: /data/results/gwas_demo/logs

compute:
  threads: 4
  memory: 4G

gwas:
  method: FaST-LMM

steps:
  reference_indexing: false
  fastq_qc: false
  alignment: false
  bam_processing: false
  variant_calling: false
  bcftools_processing: false
  vcf_quality: false
  annotation: false
  gwas: true
  batch_gwas: false
  genomic_prediction: false
  ld_analysis: false
  pav_matrix: false
  vcf_pav: false
  pangwas: false
  pangwas_plots: false
  plots: false
```

### 5. Run it

```bash
plantomicsgwas-compute run --config gwas_demo.yaml
```

### 6. Result

The run completed successfully with the following summary:

```
Status : SUCCESS | Runtime: 57.22s
Message: Completed successfully.
Outputs:
  - gwas_results: /data/results/gwas_demo/gwas/gwas_results.csv
  - manhattan_plot: /data/results/gwas_demo/gwas/manhattan_plot.png
  - qq_plot: /data/results/gwas_demo/gwas/qq_plot.png

Pipeline Summary
==================================================
Success          : True
Total steps      : 1
Successful steps : 1
Failed steps     : 0
Runtime seconds  : 57.22
==================================================
```

**53,184 variants** were tested for association using the **FaST-LMM** algorithm, with both a Manhattan plot and a QQ plot generated automatically alongside the results table. The same configuration structure works identically whether run locally (`run`) or submitted to an HPC scheduler (`submit`) — see the next section for HPC-specific usage.

---

## Two Paths to GWAS: Linear Reference vs. Pangenome

PlantOmicsGwas supports **two distinct approaches** to association analysis, reflecting the two ways plant genomic variation can be represented:

### Path A — Classical GWAS with a Linear Reference Genome

The traditional route: align reads to a single linear reference assembly, call SNPs/indels, and run GWAS on the resulting variant matrix.

```
FASTQ reads -> Reference Indexing -> Alignment -> BAM Processing ->
Variant Calling -> BCFtools Processing -> GWAS
```

This full chain was tested end-to-end on a real SLURM cluster in this release (see Testing Status). Enable the corresponding steps in your config and provide `reference_fasta`, `fastq_dir`, and `phenotype_file`.

### Path B — Pangenome-Based Association Analysis (Pan-GWAS)

Instead of a single reference, this route represents a species' genomic diversity as a **graph** (a pangenome), built from multiple individual assemblies with `minigraph`. Presence/absence variation (PAV) across samples is then tested for association with a phenotype directly — capturing structural variation and genes absent from any single reference.

```
Multiple assembly FASTAs -> Pangenome Builder (minigraph + vg) ->
PAV Matrix (from GFF annotations, or directly from a multi-sample VCF) ->
Pan-GWAS Analysis -> Pan-GWAS Plots
```

**Building the pangenome graph** (CLI, not yet an HPC-schedulable step in this release):

```bash
plantomicsgwas pangenome --assemblies /path/to/assemblies_dir --out /path/to/output_dir
```

This produces a `.gfa` graph file (and a corresponding `.vg` file), representing the combined pangenome.

**Running Pan-GWAS** on an existing PAV matrix and phenotype:

```yaml
input:
  pav_matrix: /path/to/pav_matrix.csv
  phenotype_file: /path/to/phenotype.txt

steps:
  pangwas: true
  # ... all other steps: false

pangwas:
  category: baseline      # baseline | glm | lmm
  method: ttest           # ttest | wilcoxon | fisher (baseline); ols | logistic (glm); emmax (lmm)
```

> **Note:** if your PAV matrix is derived from a VCF produced by a graph-based variant caller (e.g. `vg call`), variant IDs may not follow the `chrN:pos` format expected by the plotting module. You may need to rename markers to a `chrN:pos` format before generating Pan-GWAS plots.

> **Note:** if you convert a pangenome-derived VCF to PLINK format for SNP/indel-based
> GWAS (rather than PAV-based Pan-GWAS), see the recommended `--maf`/`--geno` values in
> the [Real GWAS Run Walkthrough](#test-real-gwas-run-walkthrough) section — the
> default thresholds are tuned for linear-reference data and can discard the vast
> majority of variants on pangenome-derived VCFs.

## Future Plan
Next version of PlantOmicsGwas will include the mapping of PanGenome of PanGenome with the user defined reference sequence.
---

## Running on HPC Clusters (SLURM/PBS/LSF)

This is the most extensively tested part of the software in this release. The Compute Engine (`plantomicsgwas-compute`) can generate a scheduler-specific job script from your YAML config, submit it, and track its status — without you writing any scheduler-specific script by hand.

### Basic workflow

```bash
# 1. Generate the job script only (inspect before running)
plantomicsgwas-compute write-job --config my_config.yaml

# 2. Generate and submit in one step
plantomicsgwas-compute submit --config my_config.yaml

# 3. Check status
plantomicsgwas-compute status --config my_config.yaml --job-id <id>

# 4. Cancel if needed
plantomicsgwas-compute cancel --config my_config.yaml --job-id <id>
```

### Critical setup requirement: shared environment across nodes

**This is the single most important lesson from our own testing.** Your conda environment must be installed somewhere visible to **every compute node** that will run a job — not just the login/head node. On most real HPC clusters this is automatic, because `$HOME` (or a shared `/data`/`/scratch`/`/project` filesystem) is mounted identically on every node via NFS, Lustre, or GPFS.

If you are setting up a small/test cluster where each node has separate local storage, install the conda environment on a path that all nodes can see (for example a shared `/data` volume), not `/root` or another node-local path.

### Activating the environment automatically inside every job

Add this to the `hpc:` section of your config:

```yaml
hpc:
  scheduler: slurm
  conda_env: /path/to/shared/miniforge3/envs/plantomicsgwas   # full path recommended
  # ... other hpc settings
```

The generated job script will run `conda activate <env>` automatically before executing the workflow.

### Avoiding repeated reference indexing across experiments

By default, each run builds its reference indexes fresh inside its own `output.dir`. If you are running multiple experiments against the **same reference genome**, set a shared index directory once, and subsequent runs will detect the existing indexes and skip rebuilding them entirely:

```yaml
reference:
  shared_index_dir: /path/to/shared/reference_index
```

In our own testing, this reduced a repeat reference-indexing step from **680 seconds to 0.00 seconds** (confirmed via the checkpoint system: `Status: SKIPPED | Message: Skipped because checkpoint exists`).

### Configuring temporary storage for large runs

FaST-LMM (used by the default `FaST-LMM` GWAS algorithm) can write a large number of
intermediate files to the system temporary directory during a run. On HPC systems,
the default temporary location (often under `$HOME` or `/tmp`) is frequently subject
to a small storage quota, which can cause an otherwise-correct analysis to fail with
a disk-space or quota error. Point the temporary directory at a location with more
available space (or fast local/RAM storage) by setting `TMPDIR` before running:

```bash
export TMPDIR=/dev/shm            # fast, RAM-backed, if it has enough space for your run
# or
export TMPDIR=/scratch/$USER/tmp  # a larger shared/scratch filesystem
mkdir -p "$TMPDIR"
```

Set this in your job script (or `~/.bashrc` on the login node) before invoking
`plantomicsgwas-compute` or `plantomicsgwas`.

### Realistic timing expectations

Set the `hpc.time` field generously — SLURM will hard-cancel a job that exceeds it. Approximate timings observed on a 4-CPU worker node against a ~505 Mbp plant genome:

| Step | Approximate time |
|---|---|
| Reference indexing (faidx + minimap2 + bowtie2-build) | ~10-11 minutes |
| FASTQ QC (small test set) | seconds |
| Alignment (small test set) | seconds |
| BAM processing | seconds |
| Variant calling + BCFtools processing | seconds |
| GWAS (FaST-LMM, ~53,000 variants) | ~1 minute |

For real-scale datasets (full read sets, larger genomes), scale these accordingly and always request more time than your best estimate.

### Full example: a real, tested multi-step config

This is the actual configuration used to verify the complete preprocessing-to-variant-calling chain end-to-end on a SLURM cluster in this release:

```yaml
project:
  name: variant_calling_test
  version: "1.0"

input:
  reference_fasta: /shared/data/reference.fa
  fastq_dir: /shared/data/fastq_reads
  reads1: /shared/data/fastq_reads/sample_R1.fastq
  platform: illumina

output:
  dir: /shared/results/variant_calling_run
  logs_dir: /shared/results/variant_calling_run/logs

compute:
  threads: 4
  memory: 4G

reference:
  shared_index_dir: /shared/results/reference_index

steps:
  reference_indexing: true
  fastq_qc: true
  alignment: true
  bam_processing: true
  variant_calling: true
  bcftools_processing: true
  vcf_quality: false
  annotation: false
  gwas: false
  batch_gwas: false
  genomic_prediction: false
  ld_analysis: false
  pav_matrix: false
  vcf_pav: false
  pangwas: false
  pangwas_plots: false
  plots: false

hpc:
  scheduler: slurm
  job_name: variant_calling_test
  partition: cpu
  nodes: 1
  tasks_per_node: 1
  cpus_per_task: 4
  memory: 4G
  time: "02:00:00"
  conda_env: /shared/miniforge3/envs/plantomicsgwas
```

### Step dependencies: which modules can run standalone

Not all steps can be enabled in isolation. Some modules (e.g. `gwas`, `ld_analysis`, `vcf_quality`, `genomic_prediction`, `pangwas`, `plots`) have no declared dependencies and will run as soon as their required inputs are provided directly. Others (e.g. `alignment`, `bam_processing`, `variant_calling`, `bcftools_processing`, `annotation`, `batch_gwas`) require their upstream steps to also be enabled in the same run, even if you already have the intermediate files — the Compute Engine currently validates this by step name, not by input availability. If you see an error like:

```
ERROR: Step 'alignment' depends on 'reference_indexing', but 'reference_indexing' is not enabled.
```

simply enable the listed dependency step(s) as well in your `steps:` section.

---

## Module Reference & Example Configs

Run `plantomicsgwas-compute list-steps` at any time to see all available steps, their categories, and descriptions. Each step reads specific keys from the `input:` section of your config — the Compute Engine will tell you exactly which required input is missing if one is not found.

| Step ID | Category | Key inputs |
|---|---|---|
| `reference_indexing` | preprocessing | `reference_fasta` |
| `fastq_qc` | preprocessing | `fastq_dir` |
| `alignment` | preprocessing | `reference_fasta`, `reads1` (+ `reads2` for paired-end), `platform` |
| `bam_processing` | variant_processing | (BAM from `alignment`) |
| `variant_calling` | variant_processing | (BAM + reference) |
| `bcftools_processing` | variant_processing | (VCF from `variant_calling`) |
| `vcf_quality` | quality_control | `vcf` |
| `annotation` | annotation | `gff_file`, `gwas_results` |
| `gwas` | association_analysis | `plink_prefix`, `phenotype_file` |
| `batch_gwas` | association_analysis | `plink_prefix`, multi-trait `phenotype_file` |
| `genomic_prediction` | machine_learning | `plink_prefix`, `phenotype_file` |
| `ld_analysis` | association_analysis | `plink_prefix` or `vcf` |
| `plots` | visualization | `plink_prefix`, `phenotype_file` |
| `pav_matrix` | pangenome | `gff_dir` (per-sample GFF files) |
| `vcf_pav` | pangenome | `vcf` |
| `pangwas` | pangenome | `pav_matrix`, `phenotype_file` |
| `pangwas_plots` | visualization | (results from `pangwas`) |

---

## Project Structure

```
PlantOmicsGwas/
├── plantvarfilter/
│   ├── core/                # Pipeline engine (context, runner, checkpoints)
│   │   └── pipelines/       # One wrapper module per workflow step
│   ├── compute/             # Headless Compute Engine + config loader
│   ├── hpc/                 # Scheduler detection, job templates, submission
│   ├── pangenome_module/    # PAV matrix, VCF-to-PAV, Pan-GWAS statistics
│   ├── preanalysis/         # Reference manager, aligner, FASTQ QC
│   ├── ui/                  # Desktop GUI pages and components
│   ├── linux/                      # Bundled Linux binaries (fallback only)
│   ├── install_plantomicsgwas.sh   # Bundled one-shot environment installer
│   └── setup_env.py                # `plantomicsgwas-setup-env` entry point
├── examples/                # Example configs
└── pyproject.toml
```

---

## Troubleshooting / FAQ

**"ModuleNotFoundError" when running `plantomicsgwas-compute` with a minimal install.**
Some workflow steps require optional dependencies (`fastlmm`, `xgboost`, `pysnptools`, `cyvcf2`). Install them with `pip install "plantomicsgwas[all]"`, or use `plantomicsgwas-setup-env` which installs everything automatically.

**A bundled tool (e.g. `bcftools`) fails with a `libhts.so.3` / GLIBC version error.**
The binaries bundled under `plantvarfilter/linux/` were built on a different base system and may not run on all Linux distributions. The software automatically prefers a working system/conda installation of the same tool if one is found on your `PATH`; installing tools via the conda environment above resolves this in virtually all cases.

**A SLURM job was cancelled with `CANCELLED AT ... DUE TO TIME LIMIT`.**
Increase the `hpc.time` value in your config. See the timing table above for realistic estimates.

**A job fails with `command not found` for `plantomicsgwas-compute` when run via SLURM.**
Your conda environment is likely not visible from the compute node, or the job script did not activate it. Set `hpc.conda_env` to the *full path* of your environment (not just its name), and confirm the environment lives on storage shared across all nodes.

**`vcf2bed` removed almost all my variants.**
Pangenome-derived VCFs commonly have high missing-genotype rates, and the default
`--geno 0.10 --maf 0.05` thresholds are tuned for linear-reference data. On a real
pangenome dataset (33.8% overall genotype call rate), these defaults removed over 95%
of variants. Try more permissive thresholds for pangenome data, e.g. `--geno 0.5 --maf 0.01`.

**GWAS (FaST-LMM) fails on HPC with a disk-space or quota error.**
FaST-LMM writes many intermediate files to the system temporary directory during a
run, which can exceed HPC storage quotas on the default location. Set `TMPDIR` to a
location with more space before running, e.g. `export TMPDIR=/dev/shm` or
`export TMPDIR=/scratch/$USER/tmp`.

---

## Citation

If you use PlantOmicsGwas in your research, please cite the accompanying manuscript (citation details will be added upon publication) and reference this repository:

```bibtex
@software{plantomicsgwas,
  author    = {Yassin, Ahmed and Khan, Falak Sher},
  title     = {PlantOmicsGwas: An Integrated GWAS, Genomic Prediction, and Pangenome-Based Association Pipeline for Plant Genomes},
  year      = {2026},
  publisher = {Ye-Lab, Peking University Institute of Advanced Agricultural Sciences (PKU-IAAS)},
  url       = {https://github.com/AHMEDY3DGENOME/PlantOmicsGWAS}
}
```

---

## License

MIT License. See `LICENSE` for details.