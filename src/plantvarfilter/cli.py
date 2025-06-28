import argparse
import gzip
import sys
import pandas as pd
from plantvarfilter import (
    improved_filter_variants,
    build_gene_db,
    annotate_variants_with_genes,
    annotate_with_traits,
)
from plantvarfilter.parser import smart_open, read_gene_traits


def main():
    parser = argparse.ArgumentParser(
        description="PlantVarFilter: Command-line variant filtering for plant genomics"
    )
    parser.add_argument("--vcf", required=True, help="Path to VCF or VCF.GZ file")
    parser.add_argument("--gff", required=True, help="Path to GFF3 or GFF3.GZ file")
    parser.add_argument("--traits", required=True, help="Path to gene trait annotation file (CSV/TSV)")
    parser.add_argument("--include-intergenic", action="store_true", help="Include intergenic variants")
    parser.add_argument("--output", default="filtered_variants.csv", help="Output CSV file")

    args = parser.parse_args()

    print("🔍 Reading VCF...")
    with (gzip.open(args.vcf) if args.vcf.endswith(".gz") else open(args.vcf, "rb")) as vcf_stream:
        variants_df = improved_filter_variants(vcf_stream, include_intergenic=args.include_intergenic)

    print(f"✅ Total variants after filtering: {len(variants_df)}")

    if variants_df.empty:
        print("❌ No variants found after filtering.")
        sys.exit(1)

    print("📖 Building gene database from GFF...")
    with (gzip.open(args.gff) if args.gff.endswith(".gz") else open(args.gff, "rb")) as gff_stream:
        gene_db = build_gene_db(gff_stream)

    print(f"✅ Gene DB loaded: {len(gene_db)} genes")

    print("🧬 Annotating variants with genes...")
    annotated_df = annotate_variants_with_genes(variants_df, gene_db, include_intergenic=args.include_intergenic)

    print("📄 Reading gene traits...")
    with smart_open(args.traits) as traits_file:
        traits_df = read_gene_traits(traits_file)

    print(f"📊 Traits loaded: {len(traits_df)} entries")

    print("🔗 Annotating variants with traits...")
    final_df = annotate_with_traits(annotated_df, traits_df)

    print(f"💾 Saving results to: {args.output}")
    final_df.to_csv(args.output, index=False)
    print("✅ Done.")


if __name__ == "__main__":
    main()
