import argparse
import gzip
import sys
import pandas as pd
import time
import logging
from plantvarfilter.filter import improved_filter_variants
from plantvarfilter.annotator import (
    build_gene_db,
    annotate_variants_with_genes,
    annotate_with_traits,
)
from plantvarfilter.parser import smart_open, read_gene_traits
import os
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s] %(levelname)s - %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)]
)

def main():
    parser = argparse.ArgumentParser(
        description="PlantVarFilter: Command-line variant filtering for plant genomics"
    )
    parser.add_argument("--vcf", required=True, help="Path to VCF or VCF.GZ file")
    parser.add_argument("--gff", required=True, help="Path to GFF3 or GFF3.GZ file")
    parser.add_argument("--traits", required=True, help="Path to gene trait annotation file (CSV/TSV)")
    parser.add_argument("--include-intergenic", action="store_true", help="Include intergenic variants")
    parser.add_argument("--output", default="filtered_variants.csv", help="Output CSV file name")

    args = parser.parse_args()

    logging.info("🚀 STARTING VARIANT FILTERING...")

    for path_label, path_value in [("VCF", args.vcf), ("GFF", args.gff), ("Traits", args.traits)]:
        if not os.path.exists(path_value):
            logging.error(f"❌ File not found: {path_value}")
            sys.exit(1)

    start_time = time.time()

    logging.info("🔍 Reading VCF...")
    with (gzip.open(args.vcf) if args.vcf.endswith(".gz") else open(args.vcf, "rb")) as vcf_stream:
        variants_df = improved_filter_variants(vcf_stream, include_intergenic=args.include_intergenic)

    if isinstance(variants_df, str):
        variants_df = pd.read_feather(variants_df)

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

    desktop_output_dir = Path.home() / "Desktop" / "PlantVarFilter_Outputs"
    desktop_output_dir.mkdir(parents=True, exist_ok=True)
    output_path = desktop_output_dir / args.output

    logging.info(f"📂 Saving results to: {output_path}")
    final_df.to_csv(output_path, index=False)

    elapsed = round(time.time() - start_time, 2)
    logging.info(f"✅ Done in {elapsed} seconds.")

if __name__ == "__main__":
    main()
