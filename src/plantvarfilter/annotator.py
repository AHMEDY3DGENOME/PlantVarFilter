import pandas as pd
import re

def build_gene_db(gff_stream):
    gene_db = {}
    for line in gff_stream:
        line = line.decode('utf-8') if isinstance(line, bytes) else line
        if line.startswith('#'):
            continue
        parts = line.strip().split('\t')
        if len(parts) != 9 or parts[2] != 'gene':
            continue
        match = re.search(r'ID=([^;]+)', parts[8])
        if match:
            gene_db[match.group(1)] = {
                "seqid": parts[0],
                "start": int(parts[3]),
                "end": int(parts[4]),
                "strand": parts[6]
            }
    print(f"✅ Gene DB built with {len(gene_db)} entries")
    return gene_db

def annotate_variants_with_genes(variants_df: pd.DataFrame, gene_db: dict, include_intergenic: bool = True) -> pd.DataFrame:
    def find_nearest_gene(chrom, pos):
        closest_gene = None
        min_distance = float('inf')

        for gene_id, info in gene_db.items():
            if info["seqid"] != chrom:
                continue
            if pos < info["start"]:
                dist = info["start"] - pos
            elif pos > info["end"]:
                dist = pos - info["end"]
            else:
                return gene_id

            if dist < min_distance:
                min_distance = dist
                closest_gene = gene_id

        return closest_gene if include_intergenic else None

    variants_df["Gene"] = variants_df.apply(
        lambda row: row["Gene"] if pd.notnull(row["Gene"]) else find_nearest_gene(row["CHROM"], row["POS"]),
        axis=1
    )
    print("✅ Annotation with nearest genes completed")
    return variants_df

def annotate_with_traits(variants_df: pd.DataFrame, traits_df: pd.DataFrame, keep_unmatched: bool = True) -> pd.DataFrame:
    result = variants_df.merge(
        traits_df,
        how='left' if keep_unmatched else 'inner',
        on='Gene'
    )
    print(f"✅ Trait annotation complete: {len(result)} variants annotated")
    return result
