```markdown
# 🚀 Usage Guide

## == Initialize Project ==
```bash
plantvarfilter init ~/Desktop/PlantTestRun
```
Creates a structured folder:
- `input/`: for input data
- `output/`: stores results and plots
- `config.json`: main configuration file

## == Run Full Pipeline ==
```bash
plantvarfilter run --config ~/Desktop/PlantTestRun/config.json
```
Performs:
- VCF parsing and filtering
- GFF3-based annotation
- Trait association
- t-test based GWAS
- Plot generation (optional)

## == Plotting Only Mode ==
If you only want to generate plots from a GWAS file:
```bash
plantvarfilter plot-only --config ~/Desktop/PlantTestRun/config.json
```
Config must include:
```json
{
  "plot_only": true,
  "output_dir": "output/",
  "gwas_results": "output/gwas_basic_results.csv"
}
```
```

---