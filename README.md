# PlantVarFilter

**PlantVarFilter** is a Python toolkit designed for efficient filtering and annotation of plant genomic variants, enabling researchers to link genetic variants with phenotypic traits and perform preliminary genome-wide association studies (GWAS). It addresses challenges in handling large variant datasets and supports integrative analysis combining genomic and trait data.

> ⚠️ Requires **Python 3.12+**  
> **Current Version: 0.1.0** — This is the first stable release.  
> Future releases aim to introduce advanced statistical models, automated reports, and interactive visualizations for plant genomics research.


---
# Citations: 
> This tool is described in the following preprint:

**Ahmed Yassin** (2025). *PlantVarFilter: A flexible tool for variant filtering and multi-trait GWAS analysis in plants*. bioRxiv.  
[https://doi.org/10.1101/2025.07.02.662805] 

 Please cite this work if you use PlantVarFilter in your research.



## Features:

- Filter variants by consequence type (e.g., missense_variant, stop_gained, synonymous_variant, frameshift_variant).
- Include or exclude intergenic regions.
- Annotate variants with gene information from GFF3 files.
- Link genes with trait scores from CSV/TSV files.
- Perform basic GWAS analyses using t-tests and multiple linear regression.
- Generate summary plots including variant consequence distribution, variant type proportions, and Manhattan plots.
- Support for compressed input files (`.gz`).
- Configurable output formats: CSV, TSV, JSON, XLSX, Feather.


---

## Project Structure

```
PlantVarFilter/
├── src/
│   └── plantvarfilter/
│       ├── __init__.py
│       ├── annotator.py
│       ├── cli.py
│       ├── filter.py
│       ├── parser.py
│       ├── regression_gwas.py
│       └── visualize.py
├── setup.py
├── README.md
└── LICENSE
```

---

## Installation

```bash
pip install .
```

Make sure you have the following dependencies installed:

`pandas`, `pyarrow`, `scipy`, `seaborn`, `matplotlib`, `numpy`, `scikit-learn`

We recommend using a Python virtual environment:

```bash
python3 -m venv env
source env/bin/activate  # Linux/macOS
env\Scripts\activate   # Windows
pip install .
```
---

## Usage

### Initialize a new analysis project

```bash
plantvarfilter init /path/to/project
```

This creates the following structure:

- `input/` — for your input data files (VCF, GFF3, trait CSV)
- `output/` — for result files and plots
- `config.json` — template configuration file

### Run the full analysis pipeline

```bash
plantvarfilter run --config /path/to/project/config.json
```

Pipeline steps:

- Filter variants based on consequence types
- Annotate variants with genes
- Annotate variants with trait data
- Perform GWAS analysis (if enabled)
- Generate output files and plots

### Generate plots from existing GWAS results

```bash
plantvarfilter plot-only --config /path/to/project/config.json
```

Requires the config file to include:

```json
{
  "plot_only": true,
  "output_dir": "output/",
  "gwas_results": "output/gwas_basic_results.csv"
}
```

> Requires config to include:
>
```json
{
  "plot_only": true,
  "output_dir": "output/",
  "gwas_results": "output/gwas_basic_results.csv"
}
```

---

## Configuration File Example (`config.json`)

```json
{
  "vcf": "input/data.vcf.gz",
  "gff": "input/annotation.gff3.gz",
  "traits": "input/traits.csv",
  "include_intergenic": true,
  "consequence_types": [
    "missense_variant",
    "stop_gained",
    "synonymous_variant"
  ],
  "output_format": "csv",
  "output_dir": "output/",
  "plot": true,
  "gwas": true
}
```

---

## Output Files

- `filtered_variants.csv` — Filtered and annotated variant dataset.
- `gwas_basic_results.csv` — GWAS association results with p-values.
- `plots/` directory contains:
  - `consequence_distribution.png`
  - `variant_type_pie.png`
  - `manhattan_plot.png`
  - `manhattan_plot_from_file.png` (for `plot-only` mode)
- `run.log` — Execution log.


---

##  Example Experiment Walkthrough

```bash
# Step 1: Create the project folder structure
plantvarfilter init ~/Desktop/PlantTestRun

# Step 2: Place your prepared input files in the input folder:
#   - expanded_variants.vcf.gz
#   - expanded_annotations.gff3.gz
#   - expanded_traits.csv

# Step 3: Update the config.json as:

{
  "vcf": "input/expanded_variants.vcf.gz",
  "gff": "input/expanded_annotations.gff3.gz",
  "traits": "input/expanded_traits.csv",
  "include_intergenic": true,
  "consequence_types": ["MODERATE", "HIGH", "LOW", "MODIFIER"],
  "output_format": "csv",
  "output": "output/filtered_variants.csv",
  "plot": true,
  "gwas": true,
  "output_dir": "output/"
}

# Step 4: Run the full pipeline
plantvarfilter run --config ~/Desktop/PlantTestRun/config.json

# Step 5 (optional): If you only want to regenerate the Manhattan Plot from a modified GWAS CSV
{
  "plot_only": true,
  "output_dir": "output/",
  "gwas_results": "output/gwas_basic_results.csv"
}

plantvarfilter plot-only --config ~/Desktop/PlantTestRun/config.json
```

###  Output Example

#### Consequence Distribution
![Consequence Plot](docs/images/consequence_distribution.png)

#### Variant Type Pie
![Variant Pie](docs/images/variant_type_pie.png)

#### Manhattan Plot
![Manhattan Plot](docs/images/manhattan_plot_from_file.png)

---

##  Future Enhancements

- Support for advanced GWAS models
- Auto-generated PDF/HTML reports
- Interactive Streamlit-based UI
- REST API
- Unit testing and test datasets

---

## License

MIT License. See `LICENSE` for details.

---

##  Author

- Ahmed Yassin || Computational Biologist
-  ahmedyassin300@outlook.com
