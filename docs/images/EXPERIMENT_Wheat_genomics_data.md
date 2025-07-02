# Wheat Variant Filtering and GWAS Analysis using PlantVarFilter

This project demonstrates the use of the `PlantVarFilter` Python package to filter genomic variants and perform a basic GWAS (Genome-Wide Association Study) on **Triticum aestivum** (wheat) using real genomic data.

---

## Data Sources

The data used in this experiment was downloaded from Ensembl Plants FTP:

- **VCF Variants**:  
  [`triticum_aestivum` VCF](https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/current/variation/vcf/triticum_aestivum/)

- **Gene Annotations**:  
  [`triticum_aestivum` GFF3](https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/current/gff3/triticum_aestivum/)

- **Traits File**:  
  A synthetic CSV file was generated for testing, containing repeated trait scores per gene.

---

## Scenario Description

We executed the **Wheat Scenario** designed to:
- **Exclude intergenic variants**
- Filter based on functional consequences
- Annotate with genes and traits
- Perform a basic t-test GWAS
- Generate summary plots

---

## Config Used

Below is the `config.json` used in the pipeline:

```json
{
  "vcf": "input/triticum_aestivum_incl_consequences.vcf.gz",
  "gff": "input/Triticum_aestivum.IWGSC.61.gff3.gz",
  "traits": "input/wheat_traits_with_repeats.csv",
  "include_intergenic": false,
  "consequence_types": ["missense_variant", "stop_gained", "synonymous_variant", "frameshift_variant"],
  "output_format": "csv",
  "output": "output/filtered_variants.csv",
  "plot": true,
  "gwas": true,
  "gwas_results": "output/gwas_basic_results.csv",
  "output_dir": "output"
}
```

---

##  Output Summary

After running the pipeline:

- **Variants Annotated**: 3776
- **Intergenic Variants Skipped**: 90M+
- **GWAS P-values**: Successfully calculated
- **Plots Generated**:
  - `variant_type_pie.png`
  - `consequence_distribution.png`
  - `manhattan_plot.png`

---

## Notes

- The traits file was manually synthesized with repeated entries per gene to enable P-value calculation.
- This setup is ideal for testing the GWAS and plotting pipeline functionality.

---

