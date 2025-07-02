# Advanced GWAS Analysis with Covariates on Rice Dataset

**Date:** 2025-07-02

##  Objective

Perform genome-wide association analysis (GWAS) on the rice (Oryza sativa) dataset with the inclusion of covariates such as **Age** and **Environment**, to simulate a more realistic and complex biological scenario.

---

## Input Files

| File | Description |
|------|-------------|
| `oryza_sativa_incl_consequences.vcf.gz` | VCF file with variant annotations |
| `Oryza_sativa.IRGSP-1.0.61.gff3.gz` | Gene annotation GFF3 file |
| `rice_traits_with_covariates.csv` | Trait scores for each gene along with Age and Environment covariates |
| `config.json` | Experiment configuration JSON file |

## Data Source
- https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/current/variation/vcf/oryza_sativa/
- https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/current/gff3/oryza_sativa/

---

## 🧰 Tools & Packages Used

- **Python 3.12**
- `pandas`, `numpy`
- `scikit-learn` for regression with covariates
- `LinearRegression` model from `sklearn`
- `ColumnTransformer` and `OneHotEncoder` for handling categorical covariates
- `matplotlib` for Manhattan plot visualization

---

## ⚙️ Configuration (`config.json`)

```json
{
  "vcf": "/home/ahmed/Desktop/rice_project/input/oryza_sativa_incl_consequences.vcf.gz",
  "gff": "/home/ahmed/Desktop/rice_project/input/Oryza_sativa.IRGSP-1.0.61.gff3.gz",
  "traits": "/home/ahmed/Desktop/rice_project/input/rice_traits_with_covariates.csv",
  "include_intergenic": false,
  "consequence_types": [
    "missense_variant",
    "stop_gained",
    "synonymous_variant",
    "frameshift_variant"
  ],
  "output_format": "csv",
  "output": "/home/ahmed/Desktop/rice_project/output/filtered_rice.csv",
  "plot": true,
  "gwas": true,
  "gwas_results": "/home/ahmed/Desktop/rice_project/output/gwas_basic_results.csv",
  "gwas_with_covariates": true,
  "gwas_covariate_results": "/home/ahmed/Desktop/rice_project/output/gwas_regression_results.csv",
  "output_dir": "/home/ahmed/Desktop/rice_project/output"
}
```

---

## Analysis Steps

1. **VCF Parsing** – Filtered variants based on selected consequences.
2. **Gene Annotation** – Mapped variants to genes using GFF3.
3. **Trait Annotation** – Annotated genes with traits + covariates.
4. **Basic GWAS (t-test)** – Trait comparison based on variant presence.
5. **Regression GWAS (Linear Model)** – `Trait_Score ~ Has_Variant + Age + Environment`.

---

## Output Files

| File | Description |
|------|-------------|
| `gwas_basic_results.csv` | Results of t-test comparison |
| `gwas_regression_results.csv` | Regression model coefficients |
| `manhattan_plot.png` | Visual summary of significance (-log10 p-value) |

---

## ✅ Outcome

 All steps executed successfully, including both basic and covariate-aware GWAS.

 Regression results saved to:
`/home/ahmed/Desktop/rice_project/output/gwas_regression_results.csv`

---

## Notes

- One-hot encoding used for environment types.
- Linear regression model includes numerical and categorical variables.
- Data validated with no missing records.

---

_Auto-generated as part of the PlantVarFilter experiment log._