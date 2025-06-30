```markdown
# 🌿 PlantVarFilter

**PlantVarFilter** is a Python-based toolkit for filtering, annotating, and analyzing plant genomic variants (VCF), linking them with gene annotations (GFF3) and trait scores for basic GWAS analysis.

> ⚠️ Requires Python 3.12+
> 
> 🧪 Current Version: 0.1.0  
> 🚀 Future versions will include interactive UIs, advanced GWAS models, and full automation for agricultural bioinformatics.

---

## 🔧 Features

- Filter variants by consequence types (e.g. missense, stop_gained, synonymous)
- Include or exclude intergenic variants
- Annotate genes using GFF3
- Link gene variants with phenotypic trait scores
- Run basic GWAS analysis using two-sample t-tests
- Visualize results (pie chart for variant types, bar chart for consequence types, Manhattan plot for GWAS)
- Output formats supported: CSV, TSV, JSON, XLSX, Feather
- Supports modular execution: full pipeline or plotting-only mode
- Works well with compressed `.vcf.gz` and `.gff3.gz` files
- Automatically builds local gene databases and logs all actions

---

## 📊 Real Use Case Example

A sample analysis using mock variant and trait data:

### Initialize a project folder:
```bash
plantvarfilter init ~/Desktop/PlantTestRun
```

### Add your input files into the `input/` directory:
- `expanded_variants.vcf.gz`
- `expanded_annotations.gff3.gz`
- `expanded_traits.csv`

### Configure the pipeline by editing `config.json`:
```json
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
  "output_dir": "output"
}
```

### Run the full analysis:
```bash
plantvarfilter run --config ~/Desktop/PlantTestRun/config.json
```

### Sample Output:
- `filtered_variants.csv`: filtered annotated variants
- `gwas_basic_results.csv`: GWAS results
- Plots in `/output/` folder

> You can also run **plotting only** later using:
```bash
plantvarfilter plot-only --config ~/Desktop/PlantTestRun/config.json
```

---

## 🌐 Documentation Pages

- [`usage.md`](./usage.md): Full CLI usage and commands
- [`config.md`](./config.md): All configuration fields explained
- [`gallery.md`](./gallery.md): Sample visualization output
```
