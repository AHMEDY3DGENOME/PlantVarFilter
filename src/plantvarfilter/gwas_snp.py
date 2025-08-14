import logging
import numpy as np
import pandas as pd
import gzip
from scipy import stats

def _open(vcf_path):
    return gzip.open(vcf_path, "rt", encoding="utf-8") if str(vcf_path).endswith(".gz") else open(vcf_path, "rt", encoding="utf-8")

def vcf_to_dosage_matrix(vcf_path, max_variants=None):
    """
    Build dosage matrix (rows=SNP 'chr:pos:ref:alt', cols=samples, values in {0,1,2,NaN})
    Requires GT in FORMAT.
    """
    samples = []
    rows = []
    meta = []

    with _open(vcf_path) as f:
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header = line.strip().lstrip("#").split("\t")
                samples = header[9:]
                continue

            parts = line.strip().split("\t")
            if len(parts) < 10:
                continue

            chrom, pos, vid, ref, alt, qual, flt, info, fmt = parts[:9]
            sample_fields = parts[9:]

            fmt_keys = fmt.split(":")
            try:
                gt_idx = fmt_keys.index("GT")
            except ValueError:
                continue

            def gt_to_dosage(gt):
                if gt is None or gt == "." or gt.startswith("."):
                    return np.nan
                g = gt.replace("|", "/").split("/")
                if len(g) != 2:
                    return np.nan
                try:
                    a = int(g[0] if g[0] != "." else -1)
                    b = int(g[1] if g[1] != "." else -1)
                except ValueError:
                    return np.nan
                if a < 0 or b < 0:
                    return np.nan
                return a + b

            dos = []
            for sf in sample_fields:
                fields = sf.split(":")
                gt = fields[gt_idx] if gt_idx < len(fields) else "."
                dos.append(gt_to_dosage(gt))

            var_key = f"{chrom}:{pos}:{ref}:{alt}"
            rows.append(dos)
            meta.append((chrom, int(pos), ref, alt, vid if vid != "." else var_key))

            if max_variants is not None and len(rows) >= max_variants:
                break

    if not rows:
        return pd.DataFrame(), pd.DataFrame()

    dosage_df = pd.DataFrame(rows, columns=samples)
    index = [f"{m[0]}:{m[1]}:{m[2]}:{m[3]}" for m in meta]
    dosage_df.index = index
    meta_df = pd.DataFrame(meta, columns=["CHROM", "POS", "REF", "ALT", "ID"])
    meta_df.index = index
    return dosage_df, meta_df

def run_snp_level_gwas(dosage_df: pd.DataFrame, traits_df: pd.DataFrame, output_csv: str):
    """
    Linear regression per SNP: Trait_Score ~ Dosage.
    traits_df must contain 'Sample' and 'Trait_Score'.
    """
    if "Sample" not in traits_df.columns:
        raise ValueError("Traits file must contain 'Sample'.")
    if "Trait_Score" not in traits_df.columns:
        raise ValueError("Traits file must contain 'Trait_Score'.")

    traits_df = traits_df.dropna(subset=["Trait_Score"])
    common = [s for s in dosage_df.columns if s in traits_df["Sample"].values]
    if not common:
        raise ValueError("No overlapping samples between VCF and traits.")

    y = traits_df.set_index("Sample").loc[common, "Trait_Score"].astype(float).values

    results = []
    for snp, row in dosage_df[common].iterrows():
        x = row.values.astype(float)
        mask = ~np.isnan(x) & ~np.isnan(y)
        if mask.sum() < 3:
            continue
        slope, intercept, r, p, stderr = stats.linregress(x[mask], y[mask])
        results.append({"SNP": snp, "Beta": slope, "R": r, "P": p, "N": int(mask.sum())})

    if not results:
        logging.warning("No valid SNP regressions.")
        pd.DataFrame(columns=["SNP", "Beta", "R", "P", "N"]).to_csv(output_csv, index=False)
        return

    out = pd.DataFrame(results).sort_values("P", ascending=True)
    out.to_csv(output_csv, index=False)
    logging.info(f"SNP-level GWAS saved: {output_csv}")
