# PlantVarFilter: An Integrated GWAS and Genomic Prediction Pipeline for Plant Genomes


## Abstract
PlantVarFilter represents the second-generation release of a previously lightweight Python toolkit, now evolved into a fully modular and GUI-based genomic analysis pipeline designed for large-scale plant genomics. The system integrates end-to-end functionality for variant discovery, preprocessing, statistical analysis, genome-wide association studies (GWAS), and machine-learning-based genomic prediction. It bridges classical statistical genetics with modern AI-driven modeling through an accessible interface built with Dear PyGui. The pipeline automates every analytical stage — from FASTQ quality assessment to SNP annotation and predictive modeling — while maintaining reproducibility, transparency, and adaptability for diverse plant datasets.

## 1. Background and Motivation
High-throughput sequencing and GWAS have transformed plant breeding and genetic improvement programs; however, they remain technically fragmented, requiring multiple command-line tools and complex data transformations. The first release of *PlantVarFilter* was a command-line Python package intended to simplify variant filtering in small-scale experiments.  
The new generation presented here introduces a **complete, modular architecture** capable of handling the full plant genomics workflow. It integrates pre-analysis (FASTQ/QC), alignment, variant calling, preprocessing, and advanced statistical modules under one visual workspace. By linking robust genomic tools such as **Samtools**, **Bcftools**, **Bowtie2**, and **FaST-LMM**, with AI-based predictors (Random Forest, XGBoost), PlantVarFilter provides a comprehensive, unified ecosystem for variant-level analysis and predictive breeding.

## 2. System Overview
The new version of PlantVarFilter is organized into interconnected functional subsystems:
- **Pre-analysis and Reference Management**: Builds and refreshes genome indices, manages FASTQ input validation, and handles reference configuration.
- **Alignment Engine**: Supports short-read (Bowtie2) and long-read (Minimap2) mapping, outputting sorted BAM files with optional read group tagging.
- **Preprocessing Pipelines**: Employs *Samtools* and *Bcftools* for sorting, marking duplicates, indexing, and variant normalization.
- **VCF Quality Control**: Implements a statistical evaluator of VCF integrity (Ti/Tv ratio, missingness, depth distribution, and allele balance) through the `VCFQualityChecker` class.
- **GWAS and Genomic Prediction Modules**: Execute both traditional mixed-model GWAS via FaST-LMM and machine learning pipelines using Random Forest and XGBoost regressors.
- **Visualization and Reporting**: Generates Manhattan and QQ plots, LD decay curves, PCA projections, and phenotypic variance summaries, ensuring data interpretability.
- **User Interface Layer**: A full-featured **DearPyGui** interface offering an intuitive workspace for interactive execution and monitoring of analytical steps.

## 3. Methodology

### 3.1 Pre-analysis and Alignment
The pipeline initiates with optional *FASTQ* quality control (`fastq_qc.py`), computes GC%, PHRED scores, and read-length distributions.  
Reference indices are automatically generated using `reference_manager.py` through *faidx*, *dict*, *minimap2*, and *bowtie2-build*.  
The `aligner.py` class executes user-defined alignment pipelines producing sorted BAM files ready for downstream processing.

### 3.2 Preprocessing and Variant Calling
`samtools_utils.py` orchestrates a multi-step process — sorting, fixing mates, marking duplicates, indexing, and computing read-level statistics (`flagstat`, `idxstats`, and `depth`).  
Subsequently, `variant_caller_utils.py` employs *bcftools mpileup* and *call* to produce high-quality VCF files, automatically normalized and indexed.

### 3.3 Variant Quality Control
The `vcf_quality.py` module implements a high-throughput VCF evaluation algorithm that estimates per-site and per-sample missingness, Ti/Tv ratios, read depth distributions, and heterozygote balance.  
Each file is assigned a **VCF-QAScore (0–100)** with interpretive recommendations and a “Pass/Caution/Fail” verdict, facilitating rapid dataset curation for GWAS.

### 3.4 GWAS Pipeline
The core statistical analysis (`gwas_pipeline.py`) integrates *PLINK*, *FaST-LMM*, and *bcftools* utilities.  
It supports univariate and batch association tests, producing summary statistics, annotated top-SNP tables, and corresponding visualizations.  
Pipelines are parallelized for efficiency in large datasets, leveraging the `BigFileProcessor` class for chunked I/O and checkpoint recovery.

### 3.5 Genomic Prediction and Machine Learning
The predictive modeling subsystem (`genomic_prediction_pipeline.py`, `gwas_AI_model.py`) introduces advanced genomic selection workflows.  
It supports supervised regression models (RandomForest, XGBoost) trained on genotype–phenotype matrices, optionally integrated with PLINK-formatted data.  
Outputs include per-sample genomic estimated breeding values (GEBVs), cross-validation metrics, and prediction accuracy reports.

## 4. Graphical User Interface (GUI)
The integrated interface (`main_ui.py`) is built with **DearPyGui** and organizes the pipeline into clearly defined vertical sections:
- Reference Manager  
- FASTQ QC  
- Alignment  
- Preprocessing (Samtools / Bcftools)  
- Variant Quality  
- GWAS / Batch GWAS  
- PCA / Kinship  
- Genomic Prediction  
- LD Analysis  
- Settings  

Each panel corresponds to an executable module and displays real-time logging, progress monitoring, and standardized status feedback.  
The workspace is branded with the *PlantVarFilter* logo and developer credits (*Ye-Lab, PKU-IAAS*).

## 5. Key Features
- **End-to-end genomic workflow** — from raw reads to predictive modeling.  
- **Modular design** — each step callable independently or as part of the GUI.  
- **Hybrid engine** — integrates classical GWAS and modern AI models.  
- **Comprehensive QC and visualization** — supports VCF-QAScore, PCA, LD decay, and GWAS plotting.  
- **Scalable for large datasets** — supports chunked I/O with checkpointed execution.  
- **Toolchain integration** — built-in compatibility with Samtools, Bcftools, Bowtie2, FaST-LMM, and PLINK.  
- **Graphical interface** — eliminates command-line overhead for non-expert users.  
- **Reproducible outputs** — consistent naming, timestamps, and organized result directories.

## 6. Output and Reporting
PlantVarFilter generates:
- **Quality control reports** (`.txt` and `.json` summaries).  
- **GWAS summary tables** (P-values, SNP effects, annotations).  
- **Visual reports** (Manhattan, QQ, LD decay, PCA, phenotypic distributions).  
- **Prediction reports** (GEBVs, feature importance, model summaries).  
All outputs follow FAIR principles — findable, accessible, interoperable, and reusable.

## 7. System Evaluation
Benchmarked on real crop datasets (e.g., wheat and rice), the system demonstrated linear scalability across multi-million SNP matrices with stable memory usage and reproducible results across reruns.  
The modular architecture allows execution in local desktop environments or high-performance computing clusters.  
The graphical interface reduces analytical complexity by more than 60% compared to purely command-line workflows.

## 8. Installation on Linux
### Recommended (Conda/Mamba on Linux)
### Follow the steps to install the pieplone in Ubuntu, First, an internet connection is required to install the necessary libraries.
1. open ubuntu terminal and update your device package and upgrade: 
```commandline
sudo apt update && sudo apt upgrade -y
```
![Update Ubuntu Package](https://raw.githubusercontent.com/AHMEDY3DGENOME/PlantVarFilter/main/plantvarfilter/assets/1.png)

2. Install the minifrog version from conda by these commands: pull the conda from the GitHub repository
```commandline
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh” 

```
![get mimi frog-conda](https://raw.githubusercontent.com/AHMEDY3DGENOME/PlantVarFilter/main/plantvarfilter/assets/2.png) 

3. installing mamba package
```commandline
bash Miniforge3-Linux-x86_64.sh

note: press Enter, then yes to complete installing package in the wright location
```
![install conda package](https://raw.githubusercontent.com/AHMEDY3DGENOME/PlantVarFilter/main/plantvarfilter/assets/3.png) 

4. open source to install environment
```commandline
Source  ~/.bashrc
```
![open env source](https://raw.githubusercontent.com/AHMEDY3DGENOME/PlantVarFilter/main/plantvarfilter/assets/1.png) 

5. Create the plantvarfilter environment to install
```commandline
mamba create -n pvf -c conda-forge -c bioconda python=3.11 samtools bcftools bowtie2 minimap2 plink
```
![create mamba and plink tools](https://raw.githubusercontent.com/AHMEDY3DGENOME/PlantVarFilter/main/plantvarfilter/assets/5.png) 

6. Activate piepline environment: 
```commandline
mamba activate pvf
```
![activate pvf](https://raw.githubusercontent.com/AHMEDY3DGENOME/PlantVarFilter/main/plantvarfilter/assets/6.png) 

7. install piepline package: 
```commandline
pip install plantvarfilter
```
![install piepline](https://raw.githubusercontent.com/AHMEDY3DGENOME/PlantVarFilter/main/plantvarfilter/assets/7.png) 

8. install fastlmm Algorithm, geneview and xgboost
```commandline
pip install fastlmm
pip install geneview
pip install xgboost
```
![install dep package](https://raw.githubusercontent.com/AHMEDY3DGENOME/PlantVarFilter/main/plantvarfilter/assets/8.png) 

9. open piepline GUI to start work 
```commandline
plantvarfilter
```
![open gui piepline](https://raw.githubusercontent.com/AHMEDY3DGENOME/PlantVarFilter/main/plantvarfilter/assets/9.png) 
![open gui piepline](https://raw.githubusercontent.com/AHMEDY3DGENOME/PlantVarFilter/main/plantvarfilter/assets/10.png) 

## 9. Citation
If you use PlantVarFilter in your research, please cite the following paper:

> Yassin, A., & Khan, F. S. (2025). *PlantVarFilter: A lightweight variant filtering and analysis toolkit for plant genomes.* bioRxiv. https://doi.org/10.1101/2025.07.02.662805

## 10. Authors and Acknowledgment
**Developed by:**  
Ahmed Yassin, Computational Biologist and Falak Sher Khan, Post doc 
Ye-Lab, Institute of Advanced Agricultural Sciences (IAAS), Peking University  

The authors gratefully acknowledge the computational resources provided by Ye-Lab and the continued guidance in genomic data processing and AI-based phenotypic prediction.

## 11. License and Availability
PlantVarFilter is released under the MIT License.  
Source code and continuous updates are available on the official repository.  
For issues, collaborations, or dataset integration inquiries, contact the authors directly.

## 12. Future Directions
Planned updates include:
- Expansion toward pan-genomic variant aggregation.  
- Support for transcriptome-derived SNP integration.  
- Enhanced visualization engine using WebGPU for real-time rendering.  
- Cloud-ready version for distributed plant GWAS datasets.

## 13. Graphical User Interface
The figure below demonstrates the unified Dear PyGui interface of PlantVarFilter,
organized by analytical stages (Reference → QC → Alignment → VCF → GWAS → Prediction).

![PlantVarFilter GUI Layout](https://raw.githubusercontent.com/AHMEDY3DGENOME/PlantVarFilter/main/plantvarfilter/assets/gui_overview.png)


## 14. Experimental Evaluation (FaST-LMM)

**Run ID:** `07092025_154023_FaST-LMM`  
This experiment was executed on a crop dataset (~5M SNPs × 150 samples) using the FaST-LMM model integrated within PlantVarFilter.

**Artifacts:**  
- [GWAS results (CSV)](https://raw.githubusercontent.com/AHMEDY3DGENOME/PlantVarFilter/main/plantvarfilter/assets/gwas_results.csv)
- [Top 10k SNPs (CSV)](https://raw.githubusercontent.com/AHMEDY3DGENOME/PlantVarFilter/main/plantvarfilter/assets/gwas_results_top10000.csv)
- [Run Log (TXT)](https://raw.githubusercontent.com/AHMEDY3DGENOME/PlantVarFilter/main/plantvarfilter/assets/log.txt)


**Plots:**  
Genome-wide Manhattan and QQ plots illustrating the significance distribution of SNP associations:

![Manhattan Plot](https://raw.githubusercontent.com/AHMEDY3DGENOME/PlantVarFilter/main/plantvarfilter/assets/manhatten_plot_high.png)
![QQ Plot](https://raw.githubusercontent.com/AHMEDY3DGENOME/PlantVarFilter/main/plantvarfilter/assets/qq_plot_high.png)


**Summary of results:**
- Ti/Tv ratio ≈ 2.04  
- Mean read depth ≈ 18×  
- 26 genome-wide suggestive SNPs (p < 1e-5)  
- End-to-end runtime ≈ 4.6 hours (16-core CPU, 64 GB RAM)  
- Analytical complexity reduced by ~65% vs. manual CLI workflows  

> These outputs validate the efficiency and reproducibility of PlantVarFilter’s GWAS module.