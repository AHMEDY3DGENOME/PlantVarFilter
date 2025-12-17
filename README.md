# PlantVarFilter: An Integrated GWAS and Genomic Prediction Pipeline for Plant Genomes

## Developers & Contributors

| Developer                         | Expertise                     |
|----------------------------------|-------------------------------|
| **Ahmed Yassin**                 | Computational Biologist       |
| **Falak Sher Khan**              | Computational Biologist       |
| **Plantvarfilter Software (Affiliation)** | Ye-Lab (PKU-IAAS)        |

------------------------------------------------------------------------------------------------------------------------
## Developed By: 

- Ahmed Yassin, Computational Biologist
- Falak Sher Khan, Computational Biologist ( Peking University Institute of Advanced Agricultural Science, PKU-IAAS)


## Acknowledgment

<div style="border: 1px solid #74c0fc; padding: 14px; border-radius: 8px; background: #e7f5ff; color: #0b7285;">
<b>The authors gratefully acknowledge the computational resources provided by Prof. Wenxiu Ye (Ye-Lab) ( Peking University Institute of Advanced Agricultural Science, PKU-IAAS) and the continued guidance in genomic data processing, phenotypic prediction and through support to complete the pipeline.</b>
</div>

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
![open env source](https://raw.githubusercontent.com/AHMEDY3DGENOME/PlantVarFilter/main/plantvarfilter/assets/4.png) 

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

## If you use PlantVarFilter in your research, please cite the following paper:  Manuscript are under write up.


## 10. Authors and Acknowledgment
**Developed by:**  
Ahmed Yassin, Computational Biologist and Falak Sher Khan, Post doc 
Ye-Lab, Institute of Advanced Agricultural Sciences (IAAS), Peking University  

The authors gratefully acknowledge the computational resources provided by Ye-Lab and the continued guidance in genomic data processing and AI-based phenotypic prediction.

## 11. License and Availability
PlantVarFilter is licensed under the **Ye-Lab (Peking University Institute of Advanced Agricultural Science, PKU-IAAS)** 
Copyright © 2025  
 

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

## 14. Full Test from Piepline 
This is a description of the entire experience, starting from building indexing to GWAS analysis and Genomics prediction. 

- From the beginning, we can build indexing from the reference and readings file. These files are raw files in the format (FASTQ)


![Building indexing](https://raw.githubusercontent.com/AHMEDY3DGENOME/PlantVarFilter/main/plantvarfilter/assets/11.png)

- the result for Building indexing 

![Result indexing](https://raw.githubusercontent.com/AHMEDY3DGENOME/PlantVarFilter/main/plantvarfilter/assets/12.png)
------------------------------------------------------------------------------------------------------------------------
## Pangenome Support (Graph-Based Reference Framework)

PlantVarFilter incorporates native support for **pangenome-based genomic analysis** alongside conventional linear reference workflows. This capability enables the construction, management, and utilization of **graph-based pangenome references**, allowing a more comprehensive representation of genomic diversity across multiple individuals or assemblies.

### Conceptual Framework

The pangenome functionality in PlantVarFilter follows a **reference-agnostic design philosophy**, ensuring flexibility without compromising compatibility with established analysis pipelines:

- Researchers possessing a **linear reference genome (FASTA)** may use the software directly without modification.
- Researchers with an existing **pangenome graph (GFA)** may register and manage it as a primary reference.
- Researchers who only have a linear reference may **construct a pangenome internally** using the integrated Pangenome Builder module.

This approach enables seamless transition between linear and graph-based representations depending on analytical requirements.

---

### Supported Reference Representations

PlantVarFilter supports multiple reference genome formats:

- **Linear Reference Genome (FASTA)**  
  A traditional single-reference genome used in alignment, variant calling, and genome-wide association studies (GWAS).

- **Graph-Based Pangenome (GFA)**  
  A non-linear reference structure integrating multiple genomes or consensus assemblies, generated using graph-based genome construction algorithms (e.g., *minigraph*).

- **Linearized Pangenome Reference (optional)**  
  A FASTA-compatible representation derived from a pangenome graph to maintain compatibility with tools requiring linear coordinate systems.

---

## Pangenome Builder Module

The **Pangenome Builder** module enables the construction of a **graph-based pangenome reference (GFA)** by integrating:

- A base linear reference genome (FASTA)
- One or more genome assemblies or consensus FASTA sequences

## Input Requirements
- Base reference genome (FASTA)
- Assemblies or consensus genomes provided as:
  - A single multi-FASTA file, or
  - A directory containing multiple FASTA files
- Output directory
- Build configuration parameters (subset or full dataset mode, minimum contig length)

## Output Products
- `pangenome.gfa` — a graph-based pangenome reference
- Build reports and log files documenting parameters, inputs, and execution details

The resulting GFA file represents a **true pangenome graph**, composed of sequence segments and graph links that capture alternative genomic paths and structural variation.

---

## Reference Utilization and Compatibility

Once generated, the pangenome graph may be registered as an active reference within PlantVarFilter.  
Graph-based references are treated as **first-class reference entities**, with automatic validation of compatibility prior to downstream analysis.

Modules that require a linear coordinate system (e.g., conventional GWAS pipelines) are explicitly validated, and users are guided to select or export a compatible linear reference representation when necessary.

---

## GWAS Considerations

Standard GWAS methodologies operate on linear genomic coordinates (VCF/PLINK formats). Therefore:

- Graph-based pangenome references (GFA) are not directly consumed by GWAS tools.
- PlantVarFilter maintains GWAS compatibility by enabling the use of **linear reference representations derived from the pangenome**.

This ensures methodological correctness while preserving the advantages of pangenome-based preprocessing.

---

## Design Rationale

The pangenome framework in PlantVarFilter is intentionally modular and non-prescriptive:

- Linear reference workflows remain fully supported.
- Pangenome-based references are optional and user-driven.
- Advanced graph-based analyses can be introduced incrementally without disrupting existing pipelines.

This design allows PlantVarFilter to serve both traditional GWAS studies and emerging pangenome-centered research paradigms.

## Pangenome module interface 

![PanGenome_module](https://raw.githubusercontent.com/AHMEDY3DGENOME/PlantVarFilter/main/plantvarfilter/assets/3352623.png)

### Test Result

### Core Files Used for Pangenome Generation

The pangenome reference in this example was generated directly from the following core input files:

- `VHP.hap1.fna`  
  Linear reference genome used as the base reference.

- `SRR23801158_1.fastq.gz`  
  Sequencing reads used to derive a consensus genome incorporated into the pangenome.

#### Output File

- `pangenome.gfa`  
  Graph-based pangenome reference generated directly from the base reference and the sequencing reads.

### OutPut Files for generate Pangenome from PlantVarFilter

![PanGenome_input](https://raw.githubusercontent.com/AHMEDY3DGENOME/PlantVarFilter/main/plantvarfilter/assets/123123.png)

### After the operation is complete, the reference is output in a separate file.

![PanGenome_output](https://raw.githubusercontent.com/AHMEDY3DGENOME/PlantVarFilter/main/plantvarfilter/assets/678678.png)


------------------------------------------------------------------------------------------------------------------------
### Reference Validation (Minimal Example)

The generated pangenome reference was validated using basic structural integrity checks to ensure compliance with accepted bioinformatics standards.

```aiignore
grep -c '^S' pangenome.gfa
grep -c '^L' pangenome.gfa

```
![PanGenome_img](https://raw.githubusercontent.com/AHMEDY3DGENOME/PlantVarFilter/main/plantvarfilter/assets/000111.png)


### Interpretation of Pangenome Graph Statistics
 - Basic inspection of the generated pangenome graph produced the following results:
The pangenome graph contains 5,001 sequence segments (S) and 7,153 graph links (L).

Sequence segments represent distinct genomic regions derived from the base reference and the incorporated sequencing data. Graph links describe the structural relationships between these segments, including alternative genomic paths and branching points.

The presence of a greater number of links than segments indicates a non-linear genome structure, confirming that the output is a true graph-based pangenome rather than a linear reference or a simple sequence concatenation.

These statistics validate the structural integrity of the generated pangenome and demonstrate successful integration of genomic variation beyond a single linear reference

-----------------------------------------------------------------------------------------------------------------------
## Reference Usage Scenarios in PlantVarFilter

PlantVarFilter is designed with a reference-agnostic architecture that supports multiple genomic analysis scenarios. Depending on the available data, users may follow one of the three workflows outlined below.

---

### Scenario 1: Using an Existing Linear Reference Genome

If the user already possesses a standard linear reference genome (FASTA format), PlantVarFilter can be used directly without any additional preprocessing.

In this scenario:
- The linear reference genome is registered in the Reference Manager.
- All downstream modules (alignment, variant calling, GWAS, PCA, and genomic prediction) operate using conventional linear-coordinate workflows.

This mode ensures full compatibility with established GWAS methodologies and external tools such as PLINK.

---

### Scenario 2: Generating a Pangenome Reference from Raw Data

If the user does not have a pangenome reference but possesses:
- A base linear reference genome (FASTA), and
- Sequencing reads or consensus assemblies,

PlantVarFilter enables internal construction of a **graph-based pangenome reference** using the integrated Pangenome Builder module.

In this scenario:
- The pangenome is generated directly from the provided data.
- The resulting graph-based reference (GFA) captures genomic variation beyond the base reference.
- The generated pangenome may be registered and managed within the software as a reference entity.

This workflow supports exploratory and population-aware genomic analysis while preserving downstream compatibility.

---

### Scenario 3: Using an Existing Pangenome Reference

If the user already has a previously generated pangenome reference (GFA format), it can be directly imported into PlantVarFilter.

In this scenario:
- The pangenome reference is validated and registered without reconstruction.
- The software treats the pangenome as a first-class reference object.
- Users may proceed with compatible downstream analyses or derive linear representations when required.

---

### Design Philosophy

These three scenarios ensure that PlantVarFilter accommodates both traditional linear-reference studies and emerging pangenome-driven research, allowing users to transition seamlessly between workflows without altering their experimental design.


-----------------------------------------------------------------------------------------------------------------------
## Alignment Stage
- At this stage, alignment is made between the reference and the raw reading.

![Alignment result](https://raw.githubusercontent.com/AHMEDY3DGENOME/PlantVarFilter/main/plantvarfilter/assets/222.png)


- The alignment result is displayed via the pipeline terminal.

![Alignment](https://raw.githubusercontent.com/AHMEDY3DGENOME/PlantVarFilter/main/plantvarfilter/assets/13.png)

After this stage, a SAM File. 

-----------------------------------------------------------------------------------------------------------------------
## VCF Stage
- At this stage, after the file is produced VCF File from the pipeline, we check its quality via the pipeline.

![Vcf qc](https://raw.githubusercontent.com/AHMEDY3DGENOME/PlantVarFilter/main/plantvarfilter/assets/14.png)

Convert VCF File to plink. 

![plink](https://raw.githubusercontent.com/AHMEDY3DGENOME/PlantVarFilter/main/plantvarfilter/assets/15.png)

Plink Result ths files. 

![plink result](https://raw.githubusercontent.com/AHMEDY3DGENOME/PlantVarFilter/main/plantvarfilter/assets/16.png)

------------------------------------------------------------------------------------------------------------------------

## GWAS Stage. 

- At this stage, we upload the resulting files from the VCF after conversion, 
- along with the phenotype file, then start the analysis from the piepline interface. 
- The results after processing will then appear in the results display terminal.

![GWAS result](https://raw.githubusercontent.com/AHMEDY3DGENOME/PlantVarFilter/main/plantvarfilter/assets/17.png)

![GWAS result one](https://raw.githubusercontent.com/AHMEDY3DGENOME/PlantVarFilter/main/plantvarfilter/assets/18.png)

![GWAS result two](https://raw.githubusercontent.com/AHMEDY3DGENOME/PlantVarFilter/main/plantvarfilter/assets/19.png)

![GWAS result three](https://raw.githubusercontent.com/AHMEDY3DGENOME/PlantVarFilter/main/plantvarfilter/assets/20.png)

------------------------------------------------------------------------------------------------------------------------

## LD Analysis Stage
- Through the new pipeline we offer, we can conduct LD Analysis from UI Piepline. 

![Ld analysis](https://raw.githubusercontent.com/AHMEDY3DGENOME/PlantVarFilter/main/plantvarfilter/assets/21.png)

- Then, a simulation of the data is displayed in the interface, which the user can download for use.

![Ld analysis result](https://raw.githubusercontent.com/AHMEDY3DGENOME/PlantVarFilter/main/plantvarfilter/assets/22.png)

![Ld analysis result o](https://raw.githubusercontent.com/AHMEDY3DGENOME/PlantVarFilter/main/plantvarfilter/assets/23.png)

![Ld analysis result t](https://raw.githubusercontent.com/AHMEDY3DGENOME/PlantVarFilter/main/plantvarfilter/assets/24.png)

------------------------------------------------------------------------------------------------------------------------

## PCA Kinship Stage
- We can also conduct PCA/KINSHIP analysis across the pipeline interface.

![pca result t](https://raw.githubusercontent.com/AHMEDY3DGENOME/PlantVarFilter/main/plantvarfilter/assets/25.png)

- The results 

![pca result tu](https://raw.githubusercontent.com/AHMEDY3DGENOME/PlantVarFilter/main/plantvarfilter/assets/26.png)

![pca result tuu](https://raw.githubusercontent.com/AHMEDY3DGENOME/PlantVarFilter/main/plantvarfilter/assets/27.png)

![pca result tuu](https://raw.githubusercontent.com/AHMEDY3DGENOME/PlantVarFilter/main/plantvarfilter/assets/28.png)


------------------------------------------------------------------------------------------------------------------------

## Genomics Prediction Stage
- In this section we can perform Genomics prediction analysis

The results

![genomics result tudu](https://raw.githubusercontent.com/AHMEDY3DGENOME/PlantVarFilter/main/plantvarfilter/assets/29.png)

![genomics result toudu](https://raw.githubusercontent.com/AHMEDY3DGENOME/PlantVarFilter/main/plantvarfilter/assets/30.png)

![genomics results toudu](https://raw.githubusercontent.com/AHMEDY3DGENOME/PlantVarFilter/main/plantvarfilter/assets/31.png)

![genomics resultsd toudu](https://raw.githubusercontent.com/AHMEDY3DGENOME/PlantVarFilter/main/plantvarfilter/assets/32.png)

## 15. Experimental Evaluation (FaST-LMM)

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


> These outputs validate the efficiency and reproducibility of PlantVarFilter’s GWAS module.