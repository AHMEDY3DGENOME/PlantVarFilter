# plantvarfilter/batch_gwas.py
import os
import tempfile
import pandas as pd
from pysnptools.snpreader import Bed, Pheno
import pysnptools.util as pstutil

def _read_pheno_multi(pheno_path):
    try:
        df = pd.read_csv(pheno_path)
        if df.shape[1] < 3:
            raise ValueError
        return df
    except Exception:
        df = pd.read_csv(pheno_path, sep=r"\s+|,", engine="python", header=0)
        if df.shape[1] < 3:
            raise ValueError(f"Phenotype file must have at least FID IID and one trait column: {pheno_path}")
        return df

def run_batch_gwas_for_all_traits(
    gwas, helper, bed_path, pheno_path, cov_path, algorithm, out_dir,
    log_fn, nr_jobs=-1, gb_goal=0, train_size=0.3, estimators=200, model_nr=1,
    max_depth=3, aggregation_method="sum", snp_limit=None
):
    os.makedirs(out_dir, exist_ok=True)
    dfp = _read_pheno_multi(pheno_path)

    cols = list(dfp.columns)
    fid_col, iid_col = cols[0], cols[1]
    trait_cols = cols[2:]
    if len(trait_cols) == 0:
        raise ValueError("No trait columns found. Expecting: FID, IID, <Trait1>, <Trait2>, ...")

    chrom_mapping = helper.replace_with_integers(bed_path.replace(".bed", ".bim"))

    summary_rows = []
    per_trait_outputs = []

    bed = Bed(str(bed_path), count_A1=False, chrom_map=chrom_mapping)

    cov = Pheno(str(cov_path)) if cov_path else None

    with tempfile.TemporaryDirectory() as tmpdir:
        for trait in trait_cols:
            trait_pheno_path = os.path.join(tmpdir, f"pheno_{trait}.txt")
            sub = dfp[[fid_col, iid_col, trait]].copy()
            sub.columns = ["FID", "IID", "Value"]
            sub = sub.dropna(subset=["Value"])
            sub.to_csv(trait_pheno_path, sep="\t", index=False, header=False)

            pheno = Pheno(str(trait_pheno_path))

            bed_i, pheno_i = pstutil.intersect_apply([bed, pheno])
            bed_fixed = gwas.filter_out_missing(bed_i)

            log_fn(f"[Batch-GWAS] Trait '{trait}': N_SNP={bed_fixed.sid_count}, N_samples={pheno_i.iid_count}")

            res_csv = os.path.join(out_dir, f"{trait}_gwas_results.csv")
            man_png = os.path.join(out_dir, f"{trait}_manhattan.png")
            qq_png  = os.path.join(out_dir, f"{trait}_qq.png")

            gwas_df, df_plot = None, None
            if algorithm in ("FaST-LMM", "Linear regression"):
                gwas_df, df_plot = gwas.run_gwas_lmm(
                    bed_fixed, pheno_i, chrom_mapping, log_fn,
                    res_csv, algorithm, bed_path, cov, gb_goal
                )
            elif algorithm == "Random Forest (AI)":
                gwas_df, df_plot = gwas.run_gwas_rf(
                    bed_fixed, pheno_i, bed_path, 1.0 - train_size, estimators,
                    res_csv, chrom_mapping, log_fn, model_nr, nr_jobs, aggregation_method
                )
            elif algorithm == "XGBoost (AI)":
                gwas_df, df_plot = gwas.run_gwas_xg(
                    bed_fixed, pheno_i, bed_path, 1.0 - train_size, estimators,
                    res_csv, chrom_mapping, log_fn, model_nr, max_depth, nr_jobs, aggregation_method
                )
            elif algorithm == "Ridge Regression":
                gwas_df, df_plot = gwas.run_gwas_ridge(
                    bed_fixed, pheno_i, bed_path, 1.0 - train_size, 1.0,
                    res_csv, chrom_mapping, log_fn, model_nr, aggregation_method
                )
            else:
                raise ValueError(f"Unsupported algorithm for batch mode: {algorithm}")

            if gwas_df is None:
                log_fn(f"[Batch-GWAS] Trait '{trait}': no results produced", warn=True)
                continue

            try:
                gwas.plot_gwas(df_plot, snp_limit, algorithm, man_png, qq_png, chrom_mapping)
            except Exception:
                pass

            top_snp, top_score = None, None
            if algorithm in ("FaST-LMM", "Linear regression") and "PValue" in gwas_df.columns:
                gwas_df = gwas_df.sort_values("PValue", ascending=True)
                row0 = gwas_df.iloc[0]
                top_snp, top_score = str(row0.get("SNP")), float(row0.get("PValue"))
            else:
                colname = "SNP effect"
                if colname in gwas_df.columns:
                    gwas_df = gwas_df.sort_values(colname, ascending=False)
                    row0 = gwas_df.iloc[0]
                    top_snp, top_score = str(row0.get("SNP")), float(row0.get(colname))

            summary_rows.append({
                "Trait": trait,
                "Top_hit": top_snp,
                "Score": top_score,
                "ResultsCSV": os.path.basename(res_csv),
                "ManhattanPNG": os.path.basename(man_png) if os.path.exists(man_png) else "",
                "QQPNG": os.path.basename(qq_png) if os.path.exists(qq_png) else "",
            })
            per_trait_outputs.append({
                "trait": trait,
                "csv": res_csv,
                "manhattan": man_png if os.path.exists(man_png) else None,
                "qq": qq_png if os.path.exists(qq_png) else None,
            })
            log_fn(f"[Batch-GWAS] Trait '{trait}': done → {os.path.basename(res_csv)}")

    summary_path = os.path.join(out_dir, "batch_gwas_summary.csv")
    pd.DataFrame(summary_rows).to_csv(summary_path, index=False)
    return {
        "summary_csv": summary_path,
        "per_trait": per_trait_outputs
    }
