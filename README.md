
# PlantVarFilter: An Intelligent and Modular Platform for Genome-Wide Association Studies (GWAS) and Genomic Prediction

## Abstract
PlantVarFilter is a cross-platform research toolkit designed to enhance Genome-Wide Association Studies (GWAS) through a unified integration of classical bioinformatics workflows and AI-driven genomic prediction. It combines standard tools such as **bcftools**, **samtools**, **plink**, and **tabix** with advanced statistical and machine learning models, including **FaST-LMM**, **Ridge Regression**, **Random Forest**, and **XGBoost**, within an automated and user-friendly environment. The platform supports both **CLI** and **GUI** interfaces, enabling reproducible, large-scale association analyses with publication-quality visual outputs.

---

## 1. System Overview
The PlantVarFilter architecture is designed around modularity and reproducibility. It consists of multiple subsystems that correspond to core GWAS workflow stages:

| Subsystem | Core Modules | Description |
|------------|---------------|--------------|
| **Preprocessing & QC** | `bcftools_utils.py`, `samtools_utils.py`, `vcf_quality.py` | Handles normalization, sorting, filtering, and integrity checks for VCF/BED files. |
| **Association Analysis** | `gwas_pipeline.py` | Implements single-trait GWAS using FaST-LMM, Linear Regression, and AI-based models. |
| **Batch Analysis** | `batch_gwas.py` | Automates GWAS across multiple traits, chromosomes, or environments. |
| **Genomic Prediction** | `genomic_prediction_pipeline.py` | Predicts quantitative traits using ML regressors and evaluates model accuracy. |
| **Visualization** | `pipeline_plots.py` | Produces Manhattan, QQ, LD, phenotype, and genotype statistics plots. |
| **User Interface** | `ui/main_ui.py` | Provides a graphical control center built with DearPyGUI for non-programmatic users. |

---

## 2. Methodology

### 2.1 Data Preprocessing
Using **bcftools** and **samtools**, the system performs standardized operations including normalization, left-alignment, sorting, multi-allelic splitting, and filtering based on user-defined thresholds. Each step ensures reproducibility by logging intermediate results.

### 2.2 Quality Control (QC)
`vcf_quality.py` implements site-level and sample-level metrics such as MAF, missingness, heterozygosity, and Hardy-Weinberg Equilibrium checks. Integration with **pysnptools** allows scalable handling of millions of SNPs.

### 2.3 GWAS Core Algorithms
- **FaST-LMM:** Linear mixed model for population-structured datasets.
- **Linear Regression:** Basic SNP-to-trait associations.
- **AI-driven Models:** Random Forest, Ridge Regression, and XGBoost enable nonlinear and ensemble-based feature discovery.
- Each model outputs a standardized results CSV with SNP, chromosome, position, and effect or p-value columns.

### 2.4 Batch GWAS
`batch_gwas.py` allows multi-trait and multi-environmental analysis through parallel execution, handling multiple phenotype columns in a single manifest. This is ideal for multi-trait plant breeding and environmental interaction studies.

### 2.5 Genomic Prediction
`genomic_prediction_pipeline.py` provides cross-validation-based prediction using RRBLUP (via Ridge), Random Forest, and XGBoost regressors. It produces predicted vs. observed trait plots and exports full CSV reports for further meta-analysis.

---

## 3. Implementation Details

### Supported Platforms
- Linux (Ubuntu ≥ 20.04)
- Windows 10/11 (via embedded binaries)
- macOS (Intel/ARM)

### Bundled Binaries
The toolkit ships with embedded versions of `bcftools`, `samtools`, `plink`, `bowtie2`, and `tabix` under `PlantVarFilter/linux` and `PlantVarFilter/windows` directories.

### Dependencies
```bash
python >= 3.10
fastlmm >= 0.6
pysnptools >= 0.5
xgboost >= 1.7
dearpygui >= 1.11
scikit-learn >= 1.3
pandas, numpy, matplotlib, seaborn, geneview
```

### Launching the GUI
```bash
python ui/main_ui.py
```

### Running via CLI
```bash
python gwas_pipeline.py --vcf data/sample.vcf.gz --pheno data/traits.csv --pheno-col yield --covars data/cov.csv --out results/demo
```

---

## 4. Experimental Workflow

### Input Files
| Type | Format | Example |
|------|---------|----------|
| Genotype | VCF/PLINK | `sample.vcf.gz` / `dataset.bed` |
| Phenotype | CSV/TSV | `FID, IID, Trait1, Trait2, ...` |
| Covariates | CSV/TSV | `IID, PC1, PC2, sex, age` |

### Output Artifacts
- Association tables (`*_gwas_results.csv`)
- Top SNPs summary (`gwas_top_snps.csv`)
- Manhattan and QQ plots
- Genomic prediction reports
- Phenotype & genotype statistics PDFs

---

## 5. Visualization and Analytics

The `pipeline_plots.py` module automates visualization of GWAS results. It generates publication-ready **Manhattan**, **QQ**, and **density plots** using Matplotlib and Seaborn. Additionally, `plot_geno_statistics()` and `plot_pheno_statistics()` compute and visualize allele frequency distributions, heterozygosity, PCA, and heritability estimates, exporting comprehensive statistical PDFs.

---

## 6. Reproducibility and Configuration

All parameters are stored in `default_settings.ini`:
```ini
[qc]
maf = 0.05
geno_missing = 0.1
sample_missing = 0.1

[gwas]
model = Linear regression
threads = 8
```
PlantVarFilter enforces full version pinning via `requirements.txt` and logs all runtime parameters, making experiments fully reproducible.

---

## 7. Discussion and Comparative Analysis

Compared with existing platforms such as **GWAStic**, **TASSEL**, and **PLINK 2.0**, PlantVarFilter provides an integrated workflow that spans from raw sequence preprocessing to machine learning-based association and genomic prediction — all in one environment. Its hybrid CLI–GUI design bridges the gap between bioinformaticians and experimental biologists, while its AI module provides competitive accuracy in nonlinear trait prediction.

---

## 8. Conclusion and Future Work

PlantVarFilter demonstrates a unified architecture for scalable GWAS and genomic prediction workflows. Future development includes:
- Integration of deep learning models (CNN/Transformer architectures)
- Cloud-based execution and collaborative data sharing
- Expansion of visualization dashboard with real-time analytics

---

## 9. Citation
If you use this software, please cite as:
```
Ahmed Yassin, Computational Biologist and Dr. Falak sher Khan et al. (2025). PlantVarFilter: An Intelligent and Modular GWAS and Genomic Prediction Toolkit.
Bioinformatics Advances, 2025.
https://github.com/AHMEDY3DGENOME/PlantVarFilter.git
```

---

## 10. License
Released under the **MIT License**.

---

## 11. Acknowledgements
PlantVarFilter builds upon open-source contributions from:
- The **bcftools**, **samtools**, and **plink** teams
- The **pysnptools** and **fastlmm** developers
- The **dearpygui** community for GUI components