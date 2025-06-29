import argparse
import gzip
import sys
import pandas as pd
import time
import logging
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime
from pathlib import Path
from plantvarfilter.filter import improved_filter_variants
from plantvarfilter.annotator import (
    build_gene_db,
    annotate_variants_with_genes,
    annotate_with_traits,
)
from plantvarfilter.parser import smart_open, read_gene_traits
import os
import pyarrow.feather as feather
from scipy import stats
import numpy as np

log_dir = Path.home() / "Desktop" / "PlantVarFilter_Outputs"
os.makedirs(log_dir, exist_ok=True)
log_file = log_dir / "run.log"

logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(log_file)
    ]
)

def generate_plots(df: pd.DataFrame, output_dir: Path):
    plot_dir = output_dir / "plots"
    plot_dir.mkdir(exist_ok=True)

    plt.figure(figsize=(10, 6))
    sns.countplot(y=df["Consequence"], order=df["Consequence"].value_counts().index)
    plt.title("Distribution of Variant Consequences")
    plt.tight_layout()
    plt.savefig(plot_dir / "consequence_distribution.png")
    plt.close()

    plt.figure(figsize=(6, 6))
    df["Variant_Type"].value_counts().plot.pie(autopct='%1.1f%%')
    plt.title("Variant Type Proportions")
    plt.ylabel("")
    plt.tight_layout()
    plt.savefig(plot_dir / "variant_type_pie.png")
    plt.close()

def run_basic_gwas(df: pd.DataFrame, traits_df: pd.DataFrame, output_dir: Path):
    result_path = output_dir / "gwas_basic_results.csv"
    manhattan_path = output_dir / "plots" / "manhattan_plot.png"
    if "Gene" not in df.columns:
        logging.warning("❗ GWAS skipped: 'Gene' column missing.")
        return

    logging.info("🧬 Running basic GWAS analysis (t-test)...")
    gwas_results = []

    traits_df.columns = traits_df.columns.str.strip()
    all_genes = traits_df["Gene"].unique()
    for gene in all_genes:
        group1 = traits_df[traits_df["Gene"] == gene]["Trait_Score"].astype(float)
        group2 = traits_df[traits_df["Gene"] != gene]["Trait_Score"].astype(float)

        if group1.empty or group2.empty:
            continue

        try:
            stat, p_value = stats.ttest_ind(group1, group2, equal_var=False)
            gwas_results.append({
                "Gene": gene,
                "Mean_Trait_With_Variant": group1.mean(),
                "Mean_Trait_Without": group2.mean(),
                "P_Value": p_value
            })
        except Exception as e:
            logging.warning(f"⚠ Error in GWAS for gene {gene}: {e}")

    results_df = pd.DataFrame(gwas_results)
    results_df.to_csv(result_path, index=False)
    logging.info(f"✅ GWAS results saved to: {result_path}")

    if not results_df.empty:
        results_df = results_df.dropna(subset=["P_Value"])
        if not results_df.empty:
            plt.figure(figsize=(12, 6))
            results_df["-log10(P_Value)"] = -np.log10(results_df["P_Value"])
            plt.scatter(range(len(results_df)), results_df["-log10(P_Value)"])
            plt.title("Manhattan Plot")
            plt.xlabel("Gene Index")
            plt.ylabel("-log10(P-value)")
            plt.tight_layout()
            plt.savefig(manhattan_path)
            plt.close()
            logging.info(f"📈 Manhattan Plot saved to: {manhattan_path}")
        else:
            logging.warning("❗ No valid P-values for Manhattan Plot.")

def main():
    parser = argparse.ArgumentParser(
        description="PlantVarFilter: Command-line variant filtering for plant genomics"
    )
    parser.add_argument("--vcf", required=True, help="Path to VCF or VCF.GZ file")
    parser.add_argument("--gff", required=True, help="Path to GFF3 or GFF3.GZ file")
    parser.add_argument("--traits", required=True, help="Path to gene trait annotation file (CSV/TSV)")
    parser.add_argument("--include-intergenic", action="store_true", help="Include intergenic variants")
    parser.add_argument("--output", default=None, help="Output file name")
    parser.add_argument("--output-format", choices=["csv", "tsv", "json", "feather", "xlsx"], default="csv", help="Output format")
    parser.add_argument("--consequence-types", nargs="*", help="Filter by consequence types (e.g. missense_variant stop_gained)")
    parser.add_argument("--plot", action="store_true", help="Generate plots of results")
    parser.add_argument("--gwas", action="store_true", help="Run basic GWAS analysis")

    args = parser.parse_args()

    logging.info("🚀 STARTING VARIANT FILTERING...")

    for label, path in [("VCF", args.vcf), ("GFF", args.gff), ("Traits", args.traits)]:
        if not os.path.exists(path):
            logging.error(f"❌ File not found: {path}")
            sys.exit(1)

    start_time = time.time()

    logging.info("🔍 Reading VCF...")
    with (gzip.open(args.vcf) if args.vcf.endswith(".gz") else open(args.vcf, "rb")) as vcf_stream:
        feather_path = improved_filter_variants(
            vcf_stream,
            include_intergenic=args.include_intergenic,
            store_as_feather=True,
            consequence_types=args.consequence_types
        )
        variants_df = pd.read_feather(feather_path)

    logging.info(f"✅ Total variants after filtering: {len(variants_df)}")

    if variants_df.empty:
        logging.warning("❌ No variants found after filtering.")
        sys.exit(1)

    logging.info("📖 Building gene database from GFF...")
    with (gzip.open(args.gff) if args.gff.endswith(".gz") else open(args.gff, "rb")) as gff_stream:
        gene_db = build_gene_db(gff_stream)

    logging.info(f"✅ Gene DB loaded: {len(gene_db)} genes")

    logging.info("🧬 Annotating variants with genes...")
    annotated_df = annotate_variants_with_genes(variants_df, gene_db, include_intergenic=args.include_intergenic)

    logging.info("📄 Reading gene traits...")
    with smart_open(args.traits) as traits_file:
        traits_df = read_gene_traits(traits_file)

    logging.info(f"📊 Traits loaded: {len(traits_df)} entries")

    logging.info("🔗 Annotating variants with traits...")
    final_df = annotate_with_traits(annotated_df, traits_df)

    timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    default_name = f"filtered_variants_{timestamp}.{args.output_format}"
    output_filename = args.output or default_name
    output_path = log_dir / output_filename

    logging.info(f"📂 Saving results to: {output_path}")

    if args.output_format == "csv":
        final_df.to_csv(output_path, index=False)
    elif args.output_format == "tsv":
        final_df.to_csv(output_path, sep='\t', index=False)
    elif args.output_format == "json":
        final_df.to_json(output_path, orient="records", lines=True)
    elif args.output_format == "feather":
        feather.write_feather(final_df, output_path)
    elif args.output_format == "xlsx":
        final_df.to_excel(output_path, index=False)

    if args.plot:
        logging.info("📈 Generating plots...")
        generate_plots(final_df, log_dir)

    if args.gwas:
        run_basic_gwas(final_df, traits_df, log_dir)

    elapsed = round(time.time() - start_time, 2)
    logging.info(f"✅ Done in {elapsed} seconds.")
    logging.shutdown()

if __name__ == "__main__":
    main()
