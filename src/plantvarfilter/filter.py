import gzip
import pandas as pd
import re
import tempfile
from typing import Union, TextIO, Optional
import pyarrow.feather as feather
import logging
import os
CHUNK_SIZE = 10000

logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler()
    ]
)

def parse_info_field(info_str: str) -> dict:
    info = {}
    for item in info_str.split(';'):
        if '=' in item:
            key, value = item.split('=', 1)
            info[key] = value
        elif item:
            info[item] = True
    return info

def parse_csq_field(info_dict: dict, csq_header: Optional[list[str]] = None) -> list:
    csq_entries = info_dict.get("CSQ", "")
    if not csq_entries:
        return []
    return [entry.split('|') for entry in csq_entries.split(',')]

def classify_variant_type(ref: str, alt: str) -> str:
    if len(ref) == 1 and len(alt) == 1:
        return "SNV"
    elif len(ref) < len(alt):
        return "Insertion"
    elif len(ref) > len(alt):
        return "Deletion"
    else:
        return "Complex"

def improved_filter_variants(
    vcf_stream: Union[str, TextIO],
    include_intergenic: bool = False,
    store_as_feather: bool = True,
    consequence_types: Optional[list[str]] = None
) -> Union[pd.DataFrame, str]:
    tmp_dir = tempfile.gettempdir()
    feather_path = os.path.join(tmp_dir, "filtered_variants.feather")
    variants = []
    csq_header = None
    count = 0
    written = False
    skipped_intergenic = 0
    skipped_no_consequence = 0

    def read_lines():
        if hasattr(vcf_stream, 'read'):
            for line in vcf_stream:
                line = line.decode('utf-8') if isinstance(line, bytes) else line
                yield line
        else:
            open_fn = gzip.open if str(vcf_stream).endswith('.gz') else open
            with open_fn(vcf_stream, 'rt') as f:
                for line in f:
                    yield line

    for line in read_lines():
        if line.startswith('##INFO=<ID=CSQ'):
            match = re.search(r'Format=([^">]+)', line)
            if match:
                csq_header = match.group(1).split('|')
        elif line.startswith('#CHROM'):
            continue
        elif not line.startswith('#'):
            fields = line.strip().split('\t')
            if len(fields) < 8:
                continue
            chrom, pos, vid, ref, alt, _, _, info_str = fields[:8]
            info_dict = parse_info_field(info_str)
            csq_entries = parse_csq_field(info_dict, csq_header)
            variant_type = classify_variant_type(ref, alt)

            for csq in csq_entries:
                csq_data = dict(zip(csq_header[:len(csq)], csq)) if csq_header else {}
                consequence = csq_data.get("Consequence", "") or (csq[1] if len(csq) > 1 else "")
                gene = csq_data.get("Feature", "") or (csq[3] if len(csq) > 3 else "")

                if not consequence.strip():
                    skipped_no_consequence += 1
                    continue
                if "intergenic_variant" in consequence and not include_intergenic:
                    skipped_intergenic += 1
                    continue
                if consequence_types and all(ct not in consequence for ct in consequence_types):
                    continue

                variants.append({
                    "CHROM": chrom,
                    "POS": int(pos),
                    "ID": vid,
                    "REF": ref,
                    "ALT": alt,
                    "Variant_Type": variant_type,
                    "Consequence": consequence,
                    "Gene": gene
                })

                count += 1
                if count % CHUNK_SIZE == 0:
                    chunk_df = pd.DataFrame(variants)
                    feather.write_feather(chunk_df, feather_path) if not written else feather.write_feather(chunk_df, feather_path, compression='uncompressed')
                    written = True
                    variants = []

    if variants:
        chunk_df = pd.DataFrame(variants)
        feather.write_feather(chunk_df, feather_path) if not written else feather.write_feather(chunk_df, feather_path, compression='uncompressed')

    logging.info(f"\U0001F4CC Skipped variants with no consequence: {skipped_no_consequence}")
    logging.info(f"\U0001F4CC Skipped intergenic variants (excluded): {skipped_intergenic}")

    return feather_path
