import argparse
import gzip
import sys
import pandas as pd
import re
from typing import Union, TextIO, Optional
from plantvarfilter.annotator import (
    build_gene_db,
    annotate_variants_with_genes,
    annotate_with_traits,
)
from plantvarfilter.parser import smart_open, read_gene_traits
import os

def parse_info_field(info_str: str) -> dict:
    info = {}
    for item in info_str.split(';'):
        if '=' in item:
            key, value = item.split('=', 1)
            info[key] = value
        elif item:  # skip empty strings
            info[item] = True
    return info

def parse_csq_field(info_dict: dict, csq_header: Optional[list[str]] = None) -> list:
    csq_entries = info_dict.get("CSQ", "")
    if not csq_entries:
        return []
    return [entry.split('|') for entry in csq_entries.split(',')]

def improved_filter_variants(vcf_stream: Union[str, TextIO], include_intergenic: bool = False) -> pd.DataFrame:
    variants = []
    csq_header = None

    def read_lines():
        if hasattr(vcf_stream, 'read'):
            for line in vcf_stream:
                line = line.decode('utf-8') if isinstance(line, bytes) else line
                print("📄 LINE:", line.strip())
                yield line
        else:
            open_fn = gzip.open if str(vcf_stream).endswith('.gz') else open
            with open_fn(vcf_stream, 'rt') as f:
                for line in f:
                    print("📄 LINE:", line.strip())
                    yield line

    for line in read_lines():
        if line.startswith('##INFO=<ID=CSQ'):
            print("📌 Found CSQ header line:", line.strip())
            match = re.search(r'Format=([^">]+)', line)
            if match:
                csq_header = match.group(1).split('|')
                print("✅ Extracted CSQ header:", csq_header)
            else:
                print("❌ Failed to extract CSQ format!")
        elif line.startswith('#CHROM'):
            continue
        elif not line.startswith('#'):
            fields = line.strip().split('\t')
            if len(fields) < 8:
                continue
            chrom, pos, vid, ref, alt, _, _, info_str = fields[:8]
            info_dict = parse_info_field(info_str)
            csq_entries = parse_csq_field(info_dict, csq_header)

            print("🧪 INFO string:", info_str)
            print("🧪 Parsed INFO dict:", info_dict)
            print("🧪 Parsed CSQ entries:", csq_entries)

            print(f"🔎 {chrom}:{pos} → Found {len(csq_entries)} CSQ entries")

            for csq in csq_entries:
                csq_data = dict(zip(csq_header[:len(csq)], csq)) if csq_header else {}
                consequence = csq_data.get("Consequence", "") or (csq[1] if len(csq) > 1 else "")
                gene = csq_data.get("Feature", "") or (csq[3] if len(csq) > 3 else "")

                print(f"   → Consequence: {consequence}, Gene: {gene}")

                if not consequence.strip():
                    continue
                if "intergenic_variant" in consequence and not include_intergenic:
                    print("   ⚠ Skipped intergenic variant (not included)")
                    continue

                variants.append({
                    "CHROM": chrom,
                    "POS": int(pos),
                    "ID": vid,
                    "REF": ref,
                    "ALT": alt,
                    "Consequence": consequence,
                    "Gene": gene
                })

    print(f"✅ Total variants after filtering: {len(variants)}")
    return pd.DataFrame(variants)

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

    print("🚀 STARTING VARIANT FILTERING...")

    if not os.path.exists(args.vcf):
        print(f"❌ File not found: {args.vcf}")
        sys.exit(1)

    if not os.path.exists(args.gff):
        print(f"❌ File not found: {args.gff}")
        sys.exit(1)

    if not os.path.exists(args.traits):
        print(f"❌ File not found: {args.traits}")
        sys.exit(1)

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

    print("🧜‍♂️ Annotating variants with genes...")
    annotated_df = annotate_variants_with_genes(variants_df, gene_db, include_intergenic=args.include_intergenic)

    print("📄 Reading gene traits...")
    with smart_open(args.traits) as traits_file:
        traits_df = read_gene_traits(traits_file)

    print(f"📊 Traits loaded: {len(traits_df)} entries")

    print("🔗 Annotating variants with traits...")
    final_df = annotate_with_traits(annotated_df, traits_df)

    print(f"📂 Saving results to: {args.output}")
    final_df.to_csv(args.output, index=False)
    print("✅ Done.")

if __name__ == "__main__":
    main()
