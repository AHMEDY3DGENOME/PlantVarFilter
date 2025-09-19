
# PlantVarFilter

PlantVarFilter is an end-to-end genomic data processing and GWAS analysis tool with a GUI built on DearPyGUI. It integrates VCF preprocessing, PLINK conversion, GWAS pipelines (FaST-LMM, Linear Regression, Ridge Regression, Random Forest, XGBoost, GLM, SAIGE), kinship matrix computation, QC plots, and result visualization.

## Features

- **VCF Preprocessing**:
  - Normalize, split multiallelic, left-align indels (bcftools)
  - Sort, filter, set variant IDs, compress and index
- **BAM/SAM Preprocessing**:
  - Sort, mark/remove duplicates, index, compute QC reports (samtools)
- **Variant Calling**:
  - Call variants from BAM/BAM-list + FASTA using bcftools
- **PLINK Conversion**:
  - Convert VCF to PLINK BED/BIM/FAM files
- **GWAS Analysis**:
  - Supported Algorithms:
    - FaST-LMM (mixed model)
    - Linear Regression
    - Ridge Regression
    - Random Forest (AI)
    - XGBoost (AI)
    - GLM (PLINK2)
    - SAIGE (mixed model)
  - Support for covariates and optional kinship matrix
  - Pheno/Geno statistics report generation
  - Manhattan & QQ plots
- **Kinship Analysis**:
  - Load precomputed kinship matrices or compute from BED
- **Batch GWAS Mode**:
  - Automate multiple GWAS runs in one session
- **Genomic Prediction Pipeline**:
  - End-to-end prediction with cross-validation
- **Interactive GUI**:
  - Navigation sidebar
  - Dynamic file dialogs
  - Integrated results window
  - Log window for live process tracking

## Installation

1. Clone the repository:
```bash
git clone https://github.com/yourusername/PlantVarFilter.git
cd PlantVarFilter
```

2. Create a virtual environment and install dependencies:
```bash
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

3. Make sure required binaries are available:
   - **bcftools**, **samtools**, **plink/plink2**, **bgzip**, **tabix**, **SAIGE** must be in PATH or in `PlantVarFilter/linux/`

## Usage

Run the GUI:
```bash
python main_ui.py
```

### Example Workflow

1. Preprocess VCF using **Preprocess (bcftools)** tab
2. Convert VCF to PLINK using **PLINK Conversion**
3. Open **GWAS** page:
   - Select BED/Phenotype/Covariate/Kinship
   - Choose algorithm (FaST-LMM, GLM, SAIGE, etc.)
   - Run GWAS
4. View results:
   - Manhattan plot
   - QQ plot
   - Top SNPs table
   - Downloadable CSV results

## Input / Output

- **Input:** `.vcf/.vcf.gz`, `.bam`, `.bed/.bim/.fam`, phenotype `.txt/.csv/.xlsx`
- **Output:** Processed VCF, PLINK BED/BIM/FAM, QC reports, GWAS CSV results, plots (PNG/PDF)

## Project Structure

```
PlantVarFilter/
├── main_ui.py              # GUI entry point
├── ui_components.py        # UI builders
├── ui_header.py
├── ui_pages.py
├── ui_theme.py
├── watermark.py
├── gwas_pipeline.py        # GWAS backend
├── gwas_AI_model.py        # AI models for GWAS
├── batch_gwas.py           # Batch GWAS runner
├── genomic_prediction_pipeline.py
├── pipeline_plots.py
├── samtools_utils.py
├── bcftools_utils.py
├── variant_caller_utils.py
├── vcf_quality.py
├── fastq_qc.py
├── aligner.py
├── reference_manager.py
├── helpers.py
└── assets/
```

## Troubleshooting

- **Segmentation fault (SIGSEGV)**: Check DearPyGUI version (`pip install dearpygui==1.11.1`)
- **Missing binaries**: Ensure `bcftools`, `samtools`, `plink2`, `SAIGE` are installed
- **No plots generated**: Make sure matplotlib backend is set to `Agg` for headless mode

## License

MIT License (c) 2025 Develop By Ahmed Yassin and Falak Sher Khan