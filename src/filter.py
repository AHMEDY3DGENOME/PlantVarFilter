import pandas as pd
import gzip
import io
import re

def parse_info_field(info_str):
    info = {}
    for item in info_str.split(';'):
        if '=' in item:
            key, value = item.split('=', 1)
            info[key] = value
        else:
            info[item] = True
    return info

def parse_csq_field(info_dict, csq_header=None):
    csq_entries = info_dict.get("CSQ", "")
    if not csq_entries:
        return []

    entries = csq_entries.split(',')
    parsed = []
    for entry in entries:
        fields = entry.split('|')
        parsed.append(fields)
    return parsed

def improved_filter_variants(vcf_stream, include_intergenic=False):
    variants = []
    csq_header = None

    if hasattr(vcf_stream, 'read'):
        lines = (line.decode('utf-8') if isinstance(line, bytes) else line for line in vcf_stream)
    else:
        with gzip.open(vcf_stream, 'rt') if str(vcf_stream).endswith('.gz') else open(vcf_stream, 'r') as f:
            lines = f.readlines()

    for line in lines:
        if line.startswith('##INFO=<ID=CSQ'):
            match = re.search(r'Format=([^">]+)', line)
            if match:
                csq_header = match.group(1).split('|')
        elif line.startswith('#CHROM'):
            header = line.strip().split('\t')
        elif not line.startswith('#'):
            fields = line.strip().split('\t')
            if len(fields) < 8:
                continue
            chrom, pos, vid, ref, alt, qual, flt, info_str = fields[:8]
            info_dict = parse_info_field(info_str)
            csq_entries = parse_csq_field(info_dict, csq_header)

            if csq_entries:
                for csq in csq_entries:
                    if csq_header and len(csq) <= len(csq_header):
                        csq_data = dict(zip(csq_header[:len(csq)], csq))
                        consequence = csq_data.get("Consequence", "")
                        gene = csq_data.get("Gene", "")
                        if not consequence.strip():
                            continue
                        if "intergenic_variant" in consequence and not include_intergenic:
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
            else:
                # fallback for basic VCF
                variants.append({
                    "CHROM": chrom,
                    "POS": int(pos),
                    "ID": vid,
                    "REF": ref,
                    "ALT": alt,
                    "Consequence": "unknown",
                    "Gene": None
                })

    return pd.DataFrame(variants)

def build_gene_db(gff_stream):
    gene_db = {}
    for line in gff_stream:
        if isinstance(line, bytes):
            line = line.decode('utf-8')
        if line.startswith('#'):
            continue
        parts = line.strip().split('\t')
        if len(parts) != 9:
            continue
        seqid, source, feature_type, start, end, score, strand, phase, attributes = parts
        if feature_type != 'gene':
            continue
        match = re.search(r'ID=([^;]+)', attributes)
        if match:
            gene_id = match.group(1)
            gene_db[gene_id] = {
                "seqid": seqid,
                "start": int(start),
                "end": int(end),
                "strand": strand
            }
    return gene_db

def read_gene_traits(file_obj):
    df = pd.read_csv(file_obj, sep=None, engine='python')
    if df.columns.duplicated().any():
        df = df.loc[:, ~df.columns.duplicated()]
    if "Gene_ID" in df.columns:
        df.rename(columns={"Gene_ID": "Gene"}, inplace=True)
    return df

def annotate_variants_with_genes(variants_df, gene_db, include_intergenic=True):
    def find_gene(chrom, pos):
        for gene_id, info in gene_db.items():
            if info["seqid"] == chrom and info["start"] <= pos <= info["end"]:
                return gene_id
        return "intergenic" if include_intergenic else None

    if "Gene" not in variants_df.columns:
        variants_df["Gene"] = None

    variants_df["Gene"] = variants_df.apply(
        lambda row: row["Gene"] if row["Gene"] else find_gene(row["CHROM"], row["POS"]),
        axis=1
    )
    return variants_df

def annotate_with_traits(variants_df, traits_df, keep_unmatched=True):
    merged = variants_df.merge(
        traits_df,
        how='left' if keep_unmatched else 'inner',
        left_on='Gene',
        right_on='Gene'
    )
    return merged
