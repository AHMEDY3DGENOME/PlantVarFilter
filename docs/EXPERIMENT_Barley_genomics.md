
#  Multi-Trait GWAS Analysis using PlantVarFilter

##  Objective
To conduct a **multi-trait Genome-Wide Association Study (GWAS)** on barley (`Hordeum vulgare`) variant data using the `PlantVarFilter` package with support for multiple phenotypic traits across different environments.

---

##  Dataset
# Source 
- https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/current/variation/vcf/hordeum_vulgare/
- https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/current/gff3/hordeum_vulgare/

### 1. VCF File
- **Path:** `/home/ahmed/Desktop/barley_project/input/hordeum_vulgare_incl_consequences.vcf.gz`
- **Content:** SNP and INDEL variants with consequence annotations.

### 2. GFF3 File
- **Path:** `/home/ahmed/Desktop/barley_project/input/Hordeum_vulgare.MorexV3_pseudomolecules_assembly.61.gff3.gz`
- **Content:** Gene annotations including coordinates and functional classification.

### 3. Trait Data File
- **Path:** `/home/ahmed/Desktop/barley_project/input/Barley_MultiTrait_Traits.csv`
- **Format:**

| Gene              | Trait_Env1 | Trait_Env2 |
|-------------------|------------|------------|
| HORVU.MOREX.r2... | 4.13       | 4.93       |
| HORVU.MOREX.r2... | 2.33       | 3.23       |
| ...               | ...        | ...        |

---

##  Config File (JSON)

```json
{
  "vcf": "/home/ahmed/Desktop/barley_project/input/hordeum_vulgare_incl_consequences.vcf.gz",
  "gff": "/home/ahmed/Desktop/barley_project/input/Hordeum_vulgare.MorexV3_pseudomolecules_assembly.61.gff3.gz",
  "traits": "/home/ahmed/Desktop/barley_project/input/Barley_MultiTrait_Traits.csv",
  "include_intergenic": true,
  "consequence_types": [
    "missense_variant",
    "stop_gained",
    "synonymous_variant",
    "frameshift_variant"
  ],
  "output_format": "csv",
  "output": "/home/ahmed/Desktop/barley_project/output/filtered_variants.csv",
  "plot": true,
  "gwas": true,
  "multi_trait_gwas": true,
  "multi_trait_output": "/home/ahmed/Desktop/barley_project/output/gwas_multi_trait_results.csv",
  "output_dir": "/home/ahmed/Desktop/barley_project/output"
}
```

---

## Execution Command

```bash
plantvarfilter run --config /home/ahmed/Desktop/barley_project/config.json
```

---

## Output

###  Successful Steps
- Gene DB built with 35826 entries
- 7006 variants annotated with traits
- **Basic GWAS skipped** due to missing column `Trait_Score`
- **Multi-trait GWAS completed using sklearn**

###  Output File:
- `/home/ahmed/Desktop/barley_project/output/gwas_multi_trait_results.csv`

---

##  Multi-Trait GWAS Results (Sample)

| Trait     | Effect_of_Variant | R_squared |
|-----------|-------------------|-----------|
| Trait_Env1 | 0.127             | 0.845     |
| Trait_Env2 | 0.219             | 0.812     |

---

##  Notes
- `LinearRegression` was used to measure the effect of variant presence across multiple traits.
- The analysis was robust to varying trait columns, skipping basic GWAS automatically when `Trait_Score` is not present.





---

© 2025 - PlantVarFilter Project