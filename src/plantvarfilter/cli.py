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
from plantvarfilter import run_regression_gwas
import os
import pyarrow.feather as feather
from scipy import stats
import numpy as np
import json
from sklearn.linear_model import LinearRegression

log_file = None

logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout)
    ]
)

def initialize_user_data(path: str):
    base_path = Path(path).expanduser().resolve()
    input_dir = base_path / "input"
    output_dir = base_path / "output"
    config_path = base_path / "config.json"

    input_dir.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)

    sample_config = {
        "vcf": str(input_dir / "your_data.vcf.gz"),
        "gff": str(input_dir / "your_annotation.gff3.gz"),
        "traits": str(input_dir / "your_traits.csv"),
        "include_intergenic": True,
        "consequence_types": ["missense_variant", "stop_gained", "synonymous_variant", "frameshift_variant"],
        "output_format": "csv",
        "output": str(output_dir / "filtered_variants.csv"),
        "plot": True,
        "gwas": True,
        "output_dir": str(output_dir)
    }

    with open(config_path, "w") as f:
        json.dump(sample_config, f, indent=4)

    print(f"✅ Project initialized at {base_path}")
    print(f"📂 - Input folder: {input_dir}")
    print(f"📂 - Output folder: {output_dir}")
    print(f"🗘️ - Config file: {config_path}")

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
    plot_dir = output_dir / "plots"
    plot_dir.mkdir(exist_ok=True)
    manhattan_path = plot_dir / "manhattan_plot.png"

    if "Gene" not in df.columns:
        logging.warning("❗ GWAS skipped: 'Gene' column missing.")
        return

    logging.info("Running basic GWAS analysis (t-test)...")
    gwas_results = []

    traits_df.columns = traits_df.columns.str.strip()
    if "Trait_Score" not in traits_df.columns:
        logging.warning("❗ GWAS skipped: 'Trait_Score' column missing in traits file.")
        return

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
            logging.warning(f"Error in GWAS for gene {gene}: {e}")

    results_df = pd.DataFrame(gwas_results)
    results_df.to_csv(result_path, index=False)
    logging.info(f"GWAS results saved to: {result_path}")

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
            logging.info(f"Manhattan Plot saved to: {manhattan_path}")
        else:
            logging.warning("No valid P-values for Manhattan Plot.")
def run_multi_trait_gwas_sklearn(variant_df: pd.DataFrame, traits_df: pd.DataFrame, output_path: str):
    if "Gene" not in variant_df.columns:
        logging.error("❌ Column 'Gene' missing in variant data.")
        return

    trait_cols = [col for col in traits_df.columns if col != "Gene"]
    if len(trait_cols) < 2:
        logging.error("❌ At least 2 traits required for multi-trait GWAS.")
        return

    variant_genes = set(variant_df["Gene"])
    traits_df["Has_Variant"] = traits_df["Gene"].apply(lambda g: 1 if g in variant_genes else 0)

    results = []
    for trait in trait_cols:
        try:
            X = traits_df[["Has_Variant"]].values
            y = traits_df[trait].values
            model = LinearRegression()
            model.fit(X, y)
            coef = model.coef_[0]
            r2 = model.score(X, y)
            results.append({
                "Trait": trait,
                "Effect_of_Variant": coef,
                "R_squared": r2
            })
        except Exception as e:
            logging.warning(f"⚠️ Failed to process trait {trait}: {e}")

    results_df = pd.DataFrame(results)
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    results_df.to_csv(output_path, index=False)
    logging.info(f"✅ Multi-trait GWAS (sklearn) saved to: {output_path}")


def load_config_file(path):
    if not path:
        return None
    if not os.path.exists(path):
        logging.error(f"Config file not found: {path}")
        sys.exit(1)
    try:
        with open(path, "r") as f:
            return json.load(f)
    except Exception as e:
        logging.error(f"Failed to load config file: {e}")
        sys.exit(1)

def run_pipeline(config):
    def get_value(key):
        return config.get(key)

    try:
        vcf_path = Path(get_value("vcf"))
        gff_path = Path(get_value("gff"))
        traits_path = Path(get_value("traits"))
    except Exception as e:
        logging.error(f"Invalid path in config: {e}")
        sys.exit(1)

    for label, path in [("VCF", vcf_path), ("GFF", gff_path), ("Traits", traits_path)]:
        if not path.exists():
            logging.error(f"{label} file not found at: {path}")
            sys.exit(1)

    output_path = Path(get_value("output"))
    output_dir = Path(get_value("output_dir")) if get_value("output_dir") else output_path.parent
    output_dir.mkdir(parents=True, exist_ok=True)

    global log_file
    log_file = output_dir / "run.log"
    logging.getLogger().handlers.clear()
    logging.basicConfig(
        level=logging.INFO,
        format='[%(asctime)s] %(levelname)s - %(message)s',
        handlers=[
            logging.StreamHandler(sys.stdout),
            logging.FileHandler(log_file)
        ]
    )

    logging.info("Reading VCF...")
    with (gzip.open(vcf_path) if str(vcf_path).endswith(".gz") else open(vcf_path, "rb")) as vcf_stream:
        feather_path = improved_filter_variants(
            vcf_stream,
            include_intergenic=get_value("include_intergenic"),
            store_as_feather=True,
            consequence_types=get_value("consequence_types")
        )
        variants_df = pd.read_feather(feather_path)

    if variants_df.empty:
        logging.warning("No variants found after filtering.")
        sys.exit(1)

    logging.info("Building gene database...")
    with (gzip.open(gff_path) if str(gff_path).endswith(".gz") else open(gff_path, "rb")) as gff_stream:
        gene_db = build_gene_db(gff_stream)

    logging.info("Annotating variants with genes...")
    annotated_df = annotate_variants_with_genes(variants_df, gene_db, include_intergenic=get_value("include_intergenic"))

    logging.info("Reading trait data...")
    traits_df = read_gene_traits(traits_path)

    logging.info("Annotating variants with traits...")
    final_df = annotate_with_traits(annotated_df, traits_df)

    fmt = get_value("output_format")
    if fmt == "csv":
        final_df.to_csv(output_path, index=False)
    elif fmt == "tsv":
        final_df.to_csv(output_path, sep="\t", index=False)
    elif fmt == "json":
        final_df.to_json(output_path, orient="records", lines=True)
    elif fmt == "feather":
        feather.write_feather(final_df, output_path)
    elif fmt == "xlsx":
        final_df.to_excel(output_path, index=False)

    if get_value("plot"):
        logging.info("Generating plots...")
        generate_plots(final_df, output_dir)

    if get_value("gwas"):
        if "Trait_Score" in traits_df.columns:
            run_basic_gwas(final_df, traits_df, output_dir)
            run_regression_gwas(final_df, traits_df, str(output_dir / "gwas_regression_results.csv"))
        else:
            logging.warning("⚠️ Skipping basic GWAS: 'Trait_Score' column not found. Only multi-trait GWAS will run.")

    if get_value("multi_trait_gwas"):
        run_multi_trait_gwas_sklearn(final_df, traits_df, get_value("multi_trait_output"))


def main():
    parser = argparse.ArgumentParser(description="PlantVarFilter - Variant Filtering for Plant Genomics")
    subparsers = parser.add_subparsers(dest="command", required=True)

    init_parser = subparsers.add_parser("init", help="Initialize a project in a custom directory")
    init_parser.add_argument("path", type=str, help="Target directory to create the project in")

    run_parser = subparsers.add_parser("run", help="Run the full analysis pipeline")
    run_parser.add_argument("--config", help="Optional path to config.json")

    plot_parser = subparsers.add_parser("plot-only", help="Generate plots from existing GWAS results CSV")
    plot_parser.add_argument("--config", type=str, help="Path to config file (should contain gwas_results and output_dir)")

    args = parser.parse_args()

    if args.command == "init":
        initialize_user_data(args.path)
    elif args.command == "run":
        config_path = args.config or (Path.home() / ".plantvarfilter_data" / "config.json")
        config = load_config_file(config_path)
        config["start_time"] = time.time()
        run_pipeline(config)
    elif args.command == "plot-only":
        config_path = args.config or (Path.home() / ".plantvarfilter_data" / "config.json")
        config = load_config_file(config_path)
        gwas_file = config.get("gwas_results")
        output_dir = Path(config.get("output_dir", "."))

        if not gwas_file or not os.path.exists(gwas_file):
            logging.error(f"GWAS results file not found: {gwas_file}")
            sys.exit(1)

        df = pd.read_csv(gwas_file)
        df = df.dropna(subset=["P_Value"])
        if df.empty:
            logging.warning("No valid P-values to plot.")
            sys.exit(0)

        plot_dir = output_dir / "plots"
        plot_dir.mkdir(exist_ok=True)

        df["-log10(P_Value)"] = -np.log10(df["P_Value"])
        plt.figure(figsize=(12, 6))
        plt.scatter(range(len(df)), df["-log10(P_Value)"])
        plt.title("Manhattan Plot (from existing file)")
        plt.xlabel("Gene Index")
        plt.ylabel("-log10(P-value)")
        plt.tight_layout()
        plt.savefig(plot_dir / "manhattan_plot_from_file.png")
        plt.close()
        logging.info(f"✅ Manhattan Plot generated from: {gwas_file}")

if __name__ == "__main__":
    main()
