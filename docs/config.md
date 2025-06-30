```markdown
# ⚙️ Configuration File (config.json)

PlantVarFilter uses a JSON file to control behavior of the pipeline.

## 🌿 Example
```json
{
  "vcf": "input/data.vcf.gz",
  "gff": "input/annotation.gff3.gz",
  "traits": "input/traits.csv",
  "include_intergenic": true,
  "consequence_types": ["MODERATE", "HIGH", "LOW", "MODIFIER"],
  "output_format": "csv",
  "output_dir": "output/",
  "plot": true,
  "gwas": true
}
```

### 🔍 Parameter Breakdown

| Field               | Description                                        |
|--------------------|----------------------------------------------------|
| `vcf`              | Compressed VCF file with variant calls             |
| `gff`              | GFF3 file with gene annotations                   |
| `traits`           | CSV linking genes to trait scores                 |
| `include_intergenic` | Whether to include non-coding variants           |
| `consequence_types` | List of variant consequences to include           |
| `output_format`    | Output type: csv, tsv, json, feather, xlsx        |
| `output_dir`       | Path to save output files                         |
| `plot`             | Whether to auto-generate plots                    |
| `gwas`             | Run basic GWAS t-test                             |
| `output`           | (Optional) Path for filtered result               |
| `plot_only`        | If true, skip filtering and use GWAS file directly|
| `gwas_results`     | Path to a GWAS result CSV file (for plotting only)|
```
