# PlantVarFilter

**PlantVarFilter** is a Python-based toolkit for filtering, annotating, and analyzing plant genomic variants (VCF), linking them with gene annotations (GFF3) and trait scores for basic GWAS analysis. The package provides an end-to-end command-line pipeline for plant genomics researchers.

---

## Features

- Filter variants by consequence type (missense, stop_gained, etc.)
- Include/exclude intergenic regions
- Annotate variants with genes using GFF3
- Link genes with trait scores
- Perform basic GWAS (t-test based)
- Generate summary plots (count plots, pie charts, Manhattan plot)
- Configurable output format: CSV, TSV, JSON, XLSX, Feather

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

Make sure the following dependencies are installed:
`pandas`, `pyarrow`, `scipy`, `seaborn`, `matplotlib`, `numpy`.

---

## Usage

### Initialize a new analysis project

To create a new project folder anywhere on your machine:

```bash
plantvarfilter init /desired/path/to/project
```

This creates the following structure:
- `input/`: for input data files (VCF, GFF, traits)
- `output/`: for result files and plots
- `config.json`: configuration template

### Run the full pipeline

```bash
plantvarfilter run --config /desired/path/to/project/config.json
```

If the `--config` option is not provided, it defaults to:

```
~/.plantvarfilter_data/config.json
```

This command performs:
- Filtering of variants
- Gene annotation
- Trait annotation
- Basic GWAS analysis (if enabled)
- Output generation and plots

---

## CLI Commands

### `plantvarfilter init <path>`

Initializes a new project folder at the given `<path>`.

### `plantvarfilter run [--config <path>]`

Runs the complete pipeline using the config file at `<path>`. If omitted, defaults to the user's home config.

### View all commands:

```bash
plantvarfilter --help
```

---

## Configuration File (`config.json`)

A typical configuration looks like this:

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

- Make sure the `input/` folder contains all necessary files.
- The `output/` folder will be populated with results and figures.

---

## Output Files

- `filtered_variants.csv`: Main filtered and annotated variant data
- `gwas_basic_results.csv`: GWAS results (p-values)
- `plots/`: Contains:
  - `consequence_distribution.png`
  - `variant_type_pie.png`
  - `manhattan_plot.png`
- `run.log`: Full log of the run

---

## Future Enhancements

- Support for advanced GWAS models
- Auto-generated PDF/HTML reports
- Interactive Streamlit-based UI
- REST API support
- Unit testing and validation datasets

---

## License

MIT License. See `LICENSE` for details.

---

## Authors

- Ahmed Yassin || Computational Biologist


For inquiries: ahmedyassin300@outlook.com
