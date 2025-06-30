# Experiment 01: Variant Filtering and GWAS Analysis on Zea mays Dataset

## 1. Introduction
This experiment demonstrates the application of the PlantVarFilter pipeline on real Zea mays genomic data obtained from Ensembl Plants.
The main objective is to perform variant filtering, gene annotation, trait association, and genome-wide association study (GWAS) analysis to identify genetic variants correlated with phenotypic traits.

## 2. Data Sources
- **Variant Call Format (VCF) file:** `zea_mays_incl_consequences.vcf.gz` (~500 MB)  
  Source: [Ensembl Plants Variation VCF](https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/current/variation/vcf/zea_mays/)  
  Contains high-confidence SNP and indel calls including consequence annotations.

- **Genome Annotation (GFF3) file:** `Zea_mays.Zm-B73-REFERENCE-NAM-5.0.61.gff3.gz` (~12 MB)  
  Source: [Ensembl Plants GFF3](https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/current/gff3/zea_mays/)  
  Provides structural annotation of genes and features for Zea mays reference genome.

- **Trait Data (CSV):** `traits_score.csv` (custom curated)  
  Contains phenotypic trait scores associated with genes, used for GWAS.

## 3. Methodology

### 3.1 Pipeline Initialization
- Created project directory structure with subfolders for `input` and `output`.
- Placed the above data files in the `input` folder.
- Configured the pipeline parameters via a JSON config file specifying input paths, variant consequence types to filter (`missense_variant`, `stop_gained`, `synonymous_variant`, `frameshift_variant`), and output preferences.

### 3.2 Commands and Execution Scenario

The analysis was performed from the terminal using the following sequence of commands:

```bash
python -m plantvarfilter init /home/ahmed/Desktop/Test_01_Zea_mays

# (Edit the config file at /home/ahmed/Desktop/Test_01_Zea_mays/config.json to set correct input/output paths)

python -m plantvarfilter run --config /home/ahmed/Desktop/Test_01_Zea_mays/config.json

python -m plantvarfilter plot-only --config /home/ahmed/Desktop/Test_01_Zea_mays/config_plot.json
```

### The pipeline processed a ~500 MB compressed VCF file and a 12 MB compressed GFF3 annotation file.
### The entire analysis, from variant filtering to GWAS and plot generation, completed in approximately 6 minutes on the test system.

### 3.3 Output
### Filtered variant data saved as CSV at ```/home/ahmed/Desktop/Test_01_Zea_mays/output/filtered_variants.csv.```
### GWAS summary statistics saved as CSV at ```/home/ahmed/Desktop/Test_01_Zea_mays/output/gwas_basic_results.csv.```
### Visualizations including Manhattan plots and variant consequence distributions saved under ```/home/ahmed/Desktop/Test_01_Zea_mays/output/plots/.```

### 4 Results
### The PlantVarFilter pipeline efficiently handles large-scale plant genomic data and produces meaningful variant-trait associations and visualizations.
### This first real-data test validates the workflow's robustness and performance, setting the stage for further biological interpretation and pipeline enhancements.

### Data retrieved from Ensembl Plants FTP server: https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/current/gff3/

### Analysis performed using PlantVarFilter package version 0.1.0

### Prepared by: Ahmed Yassin || Computational Biologist
#### Date: 2025-07-01
