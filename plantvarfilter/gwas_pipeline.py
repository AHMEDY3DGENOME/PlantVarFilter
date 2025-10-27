# gwas_pipeline.py
import os
import sys
import time
import glob
import json
import logging
import subprocess
from typing import Optional, Tuple, List

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import Ridge

from fastlmm.association import single_snp, single_snp_linreg
from pysnptools.util import log_in_place
import geneview as gv
import gzip
from bisect import bisect_left
from plantvarfilter.helpers import HELPERS

plt.switch_backend('Agg')


def _read_lines_with_fallback(path: str) -> List[str]:
    encodings = ('utf-8', 'cp1256', 'cp1252', 'latin-1', 'utf-8-sig')
    for enc in encodings:
        try:
            with open(path, 'r', encoding=enc) as f:
                return f.readlines()
        except UnicodeDecodeError:
            continue
    with open(path, 'r', encoding='utf-8', errors='replace') as f:
        return f.readlines()


def _which(exe: str) -> Optional[str]:
    p = shutil.which(exe) if 'shutil' in sys.modules else None
    if p:
        return p
    import shutil as _sh
    return _sh.which(exe)


def _ensure_executable(path: str):
    try:
        if path and os.path.exists(path):
            os.chmod(path, 0o755)
    except Exception:
        pass


def _safe_float(x):
    try:
        return float(x)
    except Exception:
        return np.nan


class GWAS:
    def __init__(self):
        self.helper = HELPERS()

    # ---------------------------
    # VCF -> PLINK BED Conversion
    # ---------------------------
    def vcf_to_bed(self, vcf_file, id_file, file_out, maf, geno):
        script_dir = os.path.dirname(__file__)
        if sys.platform.startswith('win'):
            abs_file_path = os.path.join(script_dir, "windows", "plink")
        elif sys.platform.startswith('linux'):
            abs_file_path = os.path.join(script_dir, "linux", "plink")
            _ensure_executable(abs_file_path)
        elif sys.platform.startswith('darwin'):
            abs_file_path = os.path.join(script_dir, "mac", "plink")
            _ensure_executable(abs_file_path)
        else:
            raise RuntimeError("Unsupported platform for PLINK binaries.")

        if not vcf_file or not os.path.exists(vcf_file):
            raise RuntimeError(f"VCF not found: {vcf_file}")

        threads = "4"
        memory_mb = "16000"

        cmd = [
            abs_file_path, "--vcf", vcf_file,
            "--make-bed", "--out", file_out,
            "--allow-extra-chr", "--set-missing-var-ids", "@:#",
            "--maf", str(maf), "--geno", str(geno),
            "--double-id",
            "--threads", threads,
            "--memory", memory_mb,
        ]
        if id_file:
            cmd += ["--keep", id_file]

        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        stdout, stderr = proc.stdout, proc.stderr

        bed, bim, fam = f"{file_out}.bed", f"{file_out}.bim", f"{file_out}.fam"
        bed_ok = all(os.path.exists(p) for p in (bed, bim, fam))

        if (proc.returncode != 0) or (not bed_ok):
            msg = [
                "[PLINK] conversion failed.",
                f"Command  : {' '.join(cmd)}",
                f"ExitCode : {proc.returncode}",
                "STDOUT   :",
                stdout or "(empty)",
                "STDERR   :",
                stderr or "(empty)",
            ]
            raise RuntimeError("\n".join(msg))

        logging.info("[PLINK] conversion completed.")
        return stdout or "PLINK conversion completed."

    # ---------------------------
    # Utilities
    # ---------------------------
    def filter_out_missing(self, bed):
        sid_batch_size = 1000
        all_nan = []
        with log_in_place("read snp #", logging.INFO) as updater:
            for sid_start in range(0, bed.sid_count, sid_batch_size):
                updater(f"sid {sid_start:,} of {bed.sid_count:,}")
                snp_data_batch = bed[:, sid_start:sid_start + sid_batch_size].read()
                nan_in_batch = np.isnan(snp_data_batch.val).all(axis=0)
                all_nan.append(nan_in_batch)

        all_nan = np.concatenate(all_nan, axis=0)
        logging.info(f"number of all missing columns is {np.sum(all_nan):,}. They are: ")
        logging.info(bed.sid[all_nan])
        bed_fixed = bed[:, ~all_nan]
        return bed_fixed

    def validate_gwas_input_files(self, bed_file: str, pheno_file: str):
        import re, itertools

        def _report(**kw):
            lines = ["[VALIDATE-REPORT]"]
            for k, v in kw.items():
                lines.append(f"{k}: {v}")
            msg = "\n".join(lines)
            print(msg)  # يظهر في الـLog
            return msg

        if not bed_file or not os.path.exists(bed_file):
            msg = _report(status="ERROR", reason="BED missing", bed_file=bed_file)
            return False, msg

        prefix = bed_file.replace(".bed", "")
        fam_path = prefix + ".fam"
        bim_path = prefix + ".bim"
        issues = {}

        for need in ((".bim", bim_path), (".fam", fam_path), ("pheno", pheno_file)):
            if not need[1] or not os.path.exists(need[1]):
                msg = _report(status="ERROR", reason=f"Missing {need[0]}", path=need[1])
                return False, msg

        fam_lines = _read_lines_with_fallback(fam_path)
        fam_pairs, fam_iids = [], set()
        for ln in fam_lines:
            parts = ln.strip().split()
            if len(parts) < 2:
                msg = _report(status="ERROR", reason="FAM not whitespace-delimited", line=ln[:80])
                return False, msg
            fid, iid = parts[0].strip(), parts[1].strip()
            fam_pairs.append((fid, iid))
            fam_iids.add(iid)

        ext = os.path.splitext(pheno_file)[1].lower()
        rows = []
        if ext in (".xlsx", ".xls"):
            try:
                dfp = pd.read_excel(pheno_file, header=0).dropna(how="all")
            except Exception as e:
                msg = _report(status="ERROR", reason=f"Excel read failed: {e}")
                return False, msg
            rows = dfp.values.tolist()
        else:
            raw = _read_lines_with_fallback(pheno_file)
            splitter = re.compile(r"[,\t;|:\s]+")
            for ln in raw:
                s = ln.replace("\u00A0", " ").strip()
                if not s or s.startswith("#"):
                    continue
                parts = [p for p in splitter.split(s) if p]
                if parts:
                    rows.append(parts)

        if not rows:
            msg = _report(status="ERROR", reason="Phenotype empty after parsing")
            return False, msg

        def looks_like_header(rec):
            sample = " ".join(rec[:3]).upper()
            tokens = {"FID", "IID", "PHENO", "TRAIT", "VALUE", "PHENOTYPE"}
            return any(tok in sample for tok in tokens)

        header_used = False
        if looks_like_header(rows[0]):
            header_used = True
            rows = rows[1:]
            if not rows:
                msg = _report(status="ERROR", reason="Phenotype only header")
                return False, msg

        col_counts = [len(r) for r in rows]
        common_cols = max(set(col_counts), key=col_counts.count)
        if common_cols < 2:
            msg = _report(status="ERROR", reason="Phenotype <2 columns", col_counts=col_counts[:10])
            return False, msg

        pheno_pairs = []
        for r in rows:
            if len(r) < common_cols:
                continue
            if common_cols >= 3:
                fid, iid = str(r[0]).strip(), str(r[1]).strip()
            else:
                fid, iid = "", str(r[0]).strip()
            pheno_pairs.append((fid, iid))

        if not pheno_pairs:
            msg = _report(status="ERROR", reason="No usable phenotype records")
            return False, msg

        pheno_iids = {iid for _, iid in pheno_pairs}
        overlap = fam_iids & pheno_iids

        msg = _report(
            status="OK" if overlap else "MISMATCH",
            bed_file=bed_file,
            fam_path=fam_path,
            pheno_file=pheno_file,
            fam_rows=len(fam_pairs),
            pheno_rows=len(pheno_pairs),
            header_detected=header_used,
            common_cols=common_cols,
            fam_IIDs_sample=list(itertools.islice(iter(fam_iids), 5)),
            pheno_IIDs_sample=list(itertools.islice(iter(pheno_iids), 5)),
            overlap_count=len(overlap),
            missing_in_pheno_sample=list(itertools.islice((x for x in fam_iids if x not in pheno_iids), 5)),
        )

        if not overlap:
            return False, msg
        return True, "Input files validated."

    # ---------------------------
    # GWAS (FaST-LMM / Linear)
    # ---------------------------
    def run_gwas_lmm(self, bed_fixed, pheno, chrom_mapping, add_log,
                     gwas_result_name, algorithm, bed_file, cov_file, gb_goal,
                     kinship_path: Optional[str] = None):
        t1 = time.time()
        if gb_goal == 0:
            gb_goal = None

        K_kwargs = {}
        if kinship_path and os.path.exists(kinship_path):
            try:
                if kinship_path.lower().endswith(".npy"):
                    K = np.load(kinship_path)
                    K_kwargs = {"K0": K}
                    add_log("[LMM] Using provided kinship matrix (K0).")
                else:
                    add_log("[LMM] Kinship provided but not .npy; ignoring for fastlmm.", warn=True)
            except Exception as e:
                add_log(f"[LMM] Failed to load kinship: {e}", warn=True)

        if algorithm == 'FaST-LMM':
            if cov_file:
                df_lmm_gwas = single_snp(
                    bed_fixed, pheno, output_file_name=gwas_result_name,
                    covar=cov_file, GB_goal=gb_goal, **K_kwargs
                )
            else:
                df_lmm_gwas = single_snp(
                    bed_fixed, pheno, output_file_name=gwas_result_name, GB_goal=gb_goal, **K_kwargs
                )

        elif algorithm == 'Linear regression':
            if cov_file:
                df_lmm_gwas = single_snp_linreg(
                    test_snps=bed_fixed, pheno=pheno,
                    output_file_name=gwas_result_name, covar=cov_file, GB_goal=gb_goal
                )
            else:
                df_lmm_gwas = single_snp_linreg(
                    test_snps=bed_fixed, pheno=pheno,
                    output_file_name=gwas_result_name, GB_goal=gb_goal
                )
        else:
            raise ValueError(f"Unknown algorithm: {algorithm}")

        df_lmm_gwas.dropna(subset=['PValue'], inplace=True)
        df_plot = df_lmm_gwas.copy(deep=True)

        if len(chrom_mapping) > 0:
            reversed_chrom_map = {value: key for key, value in chrom_mapping.items()}
            df_lmm_gwas["Chr"] = df_lmm_gwas["Chr"].apply(lambda x: reversed_chrom_map[x])

        df_lmm_gwas.to_csv(gwas_result_name, index=False)
        TOP_N = 50
        if "PValue" in df_lmm_gwas.columns:
            top_snps = df_lmm_gwas.nsmallest(TOP_N, "PValue")
            add_log(f"Top {TOP_N} SNPs Saved to gwas_top_snps.csv")
        t3 = round((time.time() - t1) / 60, 2)
        add_log('Final run time (minutes): ' + str(t3))
        return df_lmm_gwas, df_plot

    # ---------------------------
    # GWAS (XGBoost / RF / Ridge)
    # ---------------------------
    def run_gwas_xg(self, bed_fixed, pheno, bed_file, test_size, estimators,
                    gwas_result_name, chrom_mapping, add_log, model_nr, max_dep_set, nr_jobs, method):
        t1 = time.time()
        dataframes = []

        df_bim = pd.read_csv(bed_file.replace('bed', 'bim'), sep=r'\s+', header=None, engine='python')
        df_bim.columns = ['Chr', 'SNP', 'NA1', 'ChrPos', 'NA2', 'NA3']

        snp_data = bed_fixed.read().val
        snp_data[np.isnan(snp_data)] = -1

        import xgboost as xgb  # lazy import to avoid hard dependency on environments without xgboost

        for i in range(int(model_nr)):
            add_log('Model Iteration: ' + str(i + 1))

            X_train, X_test, y_train, y_test = train_test_split(
                snp_data, pheno.read().val, test_size=test_size
            )

            scaler = StandardScaler()
            X_train = scaler.fit_transform(X_train)

            xg_model = xgb.XGBRegressor(
                n_estimators=estimators, learning_rate=0.1,
                max_depth=max_dep_set, nthread=nr_jobs
            )
            xg_model.fit(X_train, y_train)

            data = []
            snp_ids = df_bim.iloc[:, 1].tolist()
            for col, score in zip(snp_ids, xg_model.feature_importances_):
                data.append((col, score))
            df = pd.DataFrame(data, columns=['SNP', 'PValue'])
            df = pd.merge(df, df_bim[['Chr', 'SNP', 'ChrPos']], on='SNP', how='left')
            df['Chr'] = df['Chr'].replace(chrom_mapping)
            dataframes.append(df)

        df_all_xg = self.helper.merge_models(dataframes, method)
        df_all_xg = df_all_xg.sort_values(by='PValue', ascending=False)

        df_plot = df_all_xg.copy(deep=True)
        if len(chrom_mapping) > 0:
            reversed_chrom_map = {value: key for key, value in chrom_mapping.items()}
            df_all_xg["Chr"] = df_all_xg["Chr"].apply(lambda x: reversed_chrom_map[x])

        df_all_xg.columns = df_all_xg.columns.str.replace('PValue', 'SNP effect')
        df_all_xg.to_csv(gwas_result_name, index=False)

        t3 = round((time.time() - t1) / 60, 2)
        add_log('Final run time (minutes): ' + str(t3))
        return df_all_xg, df_plot

    def run_gwas_rf(self, bed_fixed, pheno, bed_file, test_size, estimators,
                    gwas_result_name, chrom_mapping, add_log, model_nr, nr_jobs, method):
        t1 = time.time()
        dataframes = []

        df_bim = pd.read_csv(bed_file.replace('bed', 'bim'), sep=r'\s+', header=None, engine='python')
        df_bim.columns = ['Chr', 'SNP', 'NA1', 'ChrPos', 'NA2', 'NA3']

        snp_data = bed_fixed.read().val
        snp_data[np.isnan(snp_data)] = -1

        for i in range(int(model_nr)):
            add_log('Model Iteration: ' + str(i + 1))

            X_train, X_test, y_train, y_test = train_test_split(
                snp_data, pheno.read().val, test_size=test_size
            )
            scaler = StandardScaler()
            X_train = scaler.fit_transform(X_train)

            rf_model = RandomForestRegressor(n_estimators=estimators, n_jobs=nr_jobs)
            rf_model.fit(X_train, y_train.ravel())

            data = []
            snp_ids = df_bim.iloc[:, 1].tolist()
            for col, score in zip(snp_ids, rf_model.feature_importances_):
                data.append((col, score))

            df = pd.DataFrame(data, columns=['SNP', 'PValue'])
            df = pd.merge(df, df_bim[['Chr', 'SNP', 'ChrPos']], on='SNP', how='left')
            df['Chr'] = df['Chr'].replace(chrom_mapping)
            dataframes.append(df)

        df_all_rf = self.helper.merge_models(dataframes, method)
        df_all_rf = df_all_rf.sort_values(by='PValue', ascending=False)

        df_plot = df_all_rf.copy(deep=True)
        if len(chrom_mapping) > 0:
            reversed_chrom_map = {value: key for key, value in chrom_mapping.items()}
            df_all_rf["Chr"] = df_all_rf["Chr"].apply(lambda x: reversed_chrom_map[x])

        df_all_rf.columns = df_all_rf.columns.str.replace('PValue', 'SNP effect')
        df_all_rf.to_csv(gwas_result_name, index=False)

        t3 = round((time.time() - t1) / 60, 2)
        add_log('Final run time (minutes): ' + str(t3))
        return df_all_rf, df_plot

    def run_gwas_ridge(self, bed_fixed, pheno, bed_file, test_size, alpha,
                       gwas_result_name, chrom_mapping, add_log, model_nr, method):
        t1 = time.time()
        dataframes = []

        df_bim = pd.read_csv(bed_file.replace('bed', 'bim'), sep=r'\s+', header=None, engine='python')
        df_bim.columns = ['Chr', 'SNP', 'NA1', 'ChrPos', 'NA2', 'NA3']

        snp_data = bed_fixed.read().val
        snp_data[np.isnan(snp_data)] = -1

        for i in range(int(model_nr)):
            add_log('Model Iteration: ' + str(i + 1))

            X_train, X_test, y_train, y_test = train_test_split(
                snp_data, pheno.read().val, test_size=test_size
            )
            scaler = StandardScaler()
            X_train = scaler.fit_transform(X_train)

            ridge_model = Ridge(alpha=alpha)
            ridge_model.fit(X_train, y_train.ravel())

            data = []
            snp_ids = df_bim.iloc[:, 1].tolist()
            coefs = ridge_model.coef_.ravel()
            for col, coef in zip(snp_ids, coefs):
                data.append((col, float(coef)))

            df = pd.DataFrame(data, columns=['SNP', 'PValue'])
            df = pd.merge(df, df_bim[['Chr', 'SNP', 'ChrPos']], on='SNP', how='left')
            df['Chr'] = df['Chr'].replace(chrom_mapping)
            dataframes.append(df)

        df_all_ridge = self.helper.merge_models(dataframes, method)
        df_all_ridge = df_all_ridge.sort_values(by='PValue', ascending=False)

        df_plot = df_all_ridge.copy(deep=True)
        if len(chrom_mapping) > 0:
            reversed_chrom_map = {value: key for key, value in chrom_mapping.items()}
            df_all_ridge["Chr"] = df_all_ridge["Chr"].apply(lambda x: reversed_chrom_map[x])

        df_all_ridge.columns = df_all_ridge.columns.str.replace('PValue', 'SNP effect')
        df_all_ridge.to_csv(gwas_result_name, index=False)

        t3 = round((time.time() - t1) / 60, 2)
        add_log('Final run time (minutes): ' + str(t3))
        return df_all_ridge, df_plot

    # ---------------------------
    # GWAS (GLM via PLINK2)
    # ---------------------------
    def run_gwas_glm_plink2(self,
                             bed_file: str,
                             pheno_file: str,
                             cov_file: Optional[str],
                             out_csv: str,
                             chrom_mapping: dict,
                             add_log,
                             plink2_bin: str = "plink2",
                             glm_model: Optional[str] = None) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Run PLINK2 --glm. Auto-detects linear/logistic by 'glm_model' or lets PLINK2 decide.
        Expects bed_file prefix (*.bed, *.bim, *.fam) and phenotype file with FID IID Value.
        """
        prefix = bed_file.replace(".bed", "")
        out_prefix = os.path.splitext(out_csv)[0]
        out_prefix = out_prefix + ".plink2"

        bin_path = _which(plink2_bin) or plink2_bin
        cmd = [
            bin_path,
            "--bfile", prefix,
            "--pheno", pheno_file,
            "--allow-extra-chr",
            "--glm", "hide-covar", "omit-ref", "no-x-sex",
            "--out", out_prefix
        ]
        if cov_file:
            cmd += ["--covar", cov_file]
        if glm_model:
            cmd += ["--glm", glm_model]

        add_log("$ " + " ".join(cmd))
        proc = subprocess.run(cmd, text=True, capture_output=True)
        if proc.stdout:
            for ln in proc.stdout.splitlines():
                add_log(ln)
        if proc.stderr:
            for ln in proc.stderr.splitlines():
                add_log(ln, warn=True)
        if proc.returncode != 0:
            raise RuntimeError("plink2 --glm failed")

        candidates = glob.glob(out_prefix + "*.glm*")
        if not candidates:
            raise RuntimeError("No GLM output files produced by plink2.")

        use_file = None
        for pat in (".glm.linear", ".glm.logistic.hybrid", ".glm.logistic"):
            picks = [c for c in candidates if c.endswith(pat)]
            if picks:
                use_file = picks[0]
                break
        if not use_file:
            use_file = candidates[0]

        df_glm = pd.read_csv(use_file, sep=r"\s+|,", engine="python")
        col_chr = next((c for c in df_glm.columns if c.upper() in ("#CHROM", "CHROM", "CHR")), None)
        col_pos = next((c for c in df_glm.columns if c.upper() in ("POS", "BP", "BP_HG19", "BP_B37")), None)
        col_id = next((c for c in df_glm.columns if c.upper() in ("ID", "SNP", "RSID", "VARIANT_ID")), None)
        col_p = next((c for c in df_glm.columns if c.upper() in ("P", "PVAL", "PVALUE")), None)

        if not all([col_chr, col_pos, col_id, col_p]):
            raise RuntimeError(f"Unexpected GLM columns in: {use_file}")

        df = pd.DataFrame({
            "Chr": df_glm[col_chr],
            "ChrPos": df_glm[col_pos],
            "SNP": df_glm[col_id],
            "PValue": df_glm[col_p].apply(_safe_float)
        })
        df = df.dropna(subset=["PValue"])
        df = df.sort_values(["Chr", "ChrPos"]).reset_index(drop=True)

        df_plot = df.copy(deep=True)
        if len(chrom_mapping) > 0:
            reversed_chrom_map = {value: key for key, value in chrom_mapping.items()}
            df["Chr"] = df["Chr"].apply(lambda x: reversed_chrom_map.get(x, x))

        df.to_csv(out_csv, index=False)
        add_log(f"[PLINK2-GLM] Saved results: {out_csv}")
        return df, df_plot

    # ---------------------------
    # GWAS (SAIGE)
    # ---------------------------
    def run_gwas_saige(self,
                       bed_file: str,
                       pheno_file: str,
                       out_csv: str,
                       chrom_mapping: dict,
                       add_log,
                       cov_file: Optional[str] = None,
                       kinship_path: Optional[str] = None,
                       trait_type: str = "quantitative",
                       saige_step1: str = "step1_fitNULLGLMM.R",
                       saige_step2: str = "step2_SPAtests.R",
                       saige_create_grm: str = "createSparseGRM.R",
                       threads: int = 4) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        SAIGE two-step pipeline. Assumes SAIGE R scripts are available in PATH.
        If kinship_path (.rda/.rds) not provided, tries to build sparse GRM from PLINK bed.
        Outputs merged association table converted to standard columns.
        """
        prefix = bed_file.replace(".bed", "")
        out_prefix = os.path.splitext(out_csv)[0] + ".saige"
        work_dir = os.path.dirname(out_prefix) or "."

        fam = prefix + ".fam"
        bim = prefix + ".bim"
        bed = prefix + ".bed"
        for p in (fam, bim, bed, pheno_file):
            if not os.path.exists(p):
                raise RuntimeError(f"[SAIGE] Missing required file: {p}")

        add_log("[SAIGE] Preparing phenotype file (FID IID y)...")
        pheno_df = self._load_pheno_3col(pheno_file)
        pheno_out = os.path.join(work_dir, "saige.pheno.txt")
        pheno_df.to_csv(pheno_out, sep="\t", index=False, header=True)

        cov_args = []
        if cov_file and os.path.exists(cov_file):
            add_log("[SAIGE] Using covariates file (FID IID cov1 cov2 ...)")
            cov_df = self._load_covariates(cov_file)
            cov_out = os.path.join(work_dir, "saige.cov.txt")
            cov_df.to_csv(cov_out, sep="\t", index=False, header=True)
            cov_args = ["--covarColList", ",".join([c for c in cov_df.columns if c not in ("FID", "IID")])]
        else:
            cov_out = None

        null_model_rda = os.path.join(work_dir, "saige.null.model.rda")
        sample_id_col = "IID"

        if kinship_path and os.path.exists(kinship_path) and kinship_path.lower().endswith((".rda", ".rds")):
            add_log("[SAIGE] Using provided GRM/kinship file for null model.")
            grm_prefix = kinship_path
            sparse_grm_mtx = None
        else:
            add_log("[SAIGE] Building sparse GRM from PLINK BED...")
            grm_prefix = os.path.join(work_dir, "saige.sparseGRM")
            cmd_grm = [
                "Rscript", saige_create_grm,
                "--bfile", prefix,
                "--outputPrefix", grm_prefix,
                "--relatednessCutoff", "0.05",
                "--numRandomMarkerforSparseKin", "2000"
            ]
            add_log("$ " + " ".join(cmd_grm))
            p = subprocess.run(cmd_grm, text=True, capture_output=True)
            if p.stdout:
                for ln in p.stdout.splitlines():
                    add_log(ln)
            if p.stderr:
                for ln in p.stderr.splitlines():
                    add_log(ln, warn=True)
            if p.returncode != 0:
                raise RuntimeError("[SAIGE] createSparseGRM failed.")

            sparse_grm_mtx = grm_prefix + ".relatednessCutoff_0.05_2000_randomMarkersUsed.sparseGRM.mtx"
            sparse_grm_col = grm_prefix + ".relatednessCutoff_0.05_2000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt"
            if not os.path.exists(sparse_grm_mtx):
                raise RuntimeError("[SAIGE] Sparse GRM files not found.")

        add_log("[SAIGE] Step1: fit NULL GLMM...")
        cmd_step1 = [
            "Rscript", saige_step1,
            "--plinkFile", prefix,
            "--phenoFile", pheno_out,
            "--phenoCol", "y",
            "--covarColList", ",".join([]),
            "--sampleIDColinphenoFile", sample_id_col,
            "--traitType", "quantitative" if trait_type == "quantitative" else "binary",
            "--outputPrefix", out_prefix + ".null",
            "--nThreads", str(threads),
        ]
        if cov_args:
            cmd_step1[cmd_step1.index("--covarColList") + 1] = cov_args[-1].split()[-1]
        if kinship_path and kinship_path.lower().endswith((".rda", ".rds")):
            cmd_step1 += ["--useSparseGRMtoFitNULL", "TRUE", "--sparseGRMFile", kinship_path]
        elif 'sparse_grm_mtx' in locals() and sparse_grm_mtx:
            cmd_step1 += ["--useSparseGRMtoFitNULL", "TRUE", "--sparseGRMFile", sparse_grm_mtx]

        add_log("$ " + " ".join(cmd_step1))
        p1 = subprocess.run(cmd_step1, text=True, capture_output=True)
        if p1.stdout:
            for ln in p1.stdout.splitlines():
                add_log(ln)
        if p1.stderr:
            for ln in p1.stderr.splitlines():
                add_log(ln, warn=True)
        if p1.returncode != 0:
            raise RuntimeError("[SAIGE] step1_fitNULLGLMM failed.")
        if not os.path.exists(out_prefix + ".null.rda"):
            add_log("[SAIGE] Warning: null model rda not found at expected name; trying generic name.", warn=True)

        add_log("[SAIGE] Step2: association testing...")
        out_assoc = out_prefix + ".assoc.txt"
        cmd_step2 = [
            "Rscript", saige_step2,
            "--bgenFile", "",  # we are using PLINK bed; leave blank for bgen route
            "--plinkFile", prefix,
            "--vcfFile", "",
            "--vcfFileIndex", "",
            "--minMAF", "0.0",
            "--minMAC", "1",
            "--GMMATmodelFile", out_prefix + ".null.rda",
            "--varianceRatioFile", out_prefix + ".null.varianceRatio.txt",
            "--SAIGEOutputFile", out_assoc,
            "--numLinesOutput", "2"
        ]
        add_log("$ " + " ".join([c for c in cmd_step2 if c]))
        p2 = subprocess.run([c for c in cmd_step2 if c], text=True, capture_output=True)
        if p2.stdout:
            for ln in p2.stdout.splitlines():
                add_log(ln)
        if p2.stderr:
            for ln in p2.stderr.splitlines():
                add_log(ln, warn=True)
        if p2.returncode != 0:
            raise RuntimeError("[SAIGE] step2_SPAtests failed.")

        if not os.path.exists(out_assoc):
            raise RuntimeError("[SAIGE] Association output not found.")

        df_s = pd.read_csv(out_assoc, sep=r"\s+|,", engine="python")
        col_chr = next((c for c in df_s.columns if c.upper() in ("CHR", "CHROM", "#CHROM")), None)
        col_pos = next((c for c in df_s.columns if c.upper() in ("POS", "BP")), None)
        col_id = next((c for c in df_s.columns if "SNPID" in c.upper() or c.upper() in ("MARKERID", "SNP", "ID")), None)
        col_p = next((c for c in df_s.columns if c.upper() in ("PVAL", "PVALUE", "P", "SPA.PVAL", "P.VALUE")), None)

        if not all([col_chr, col_pos, col_id, col_p]):
            raise RuntimeError("[SAIGE] Unexpected columns in association output.")

        df = pd.DataFrame({
            "Chr": df_s[col_chr],
            "ChrPos": df_s[col_pos],
            "SNP": df_s[col_id],
            "PValue": df_s[col_p].apply(_safe_float)
        })
        df = df.dropna(subset=["PValue"])
        df = df.sort_values(["Chr", "ChrPos"]).reset_index(drop=True)

        df_plot = df.copy(deep=True)
        if len(chrom_mapping) > 0:
            reversed_chrom_map = {value: key for key, value in chrom_mapping.items()}
            df["Chr"] = df["Chr"].apply(lambda x: reversed_chrom_map.get(x, x))

        df.to_csv(out_csv, index=False)
        add_log(f"[SAIGE] Saved results: {out_csv}")
        return df, df_plot

    # ---------------------------
    # Plotting
    # ---------------------------
    def plot_gwas(
            self,
            df,
            limit,
            algorithm,
            manhatten_plot_name,
            qq_plot_name,
            chrom_mapping,
            region: str = None,
            region_only_csv: str = None,
            title_suffix: str = None,
    ):
        import math
        import re

        def _parse_region(s: str):
            if not s or not str(s).strip():
                return None, None, None
            s = str(s).strip()
            s = s.replace(" ", "")
            m = re.match(r"^(chr)?([A-Za-z0-9]+)(?::([0-9eE\+\-]+)-([0-9eE\+\-]+))?$", s)
            if not m:
                return None, None, None
            chrom_token = m.group(2)
            try:
                start = int(float(m.group(3))) if m.group(3) else None
                end = int(float(m.group(4))) if m.group(4) else None
            except Exception:
                start, end = None, None
            if start is not None and end is not None and end < start:
                start, end = end, start
            return chrom_token, start, end

        def _apply_region_filter(df_in, region_str):
            if not region_str:
                return df_in.copy(), None
            chrom_token, start, end = _parse_region(region_str)
            if chrom_token is None:
                return df_in.copy(), None

            df_loc = df_in.copy()
            if "Chr" not in df_loc.columns or "ChrPos" not in df_loc.columns:
                return df_in.copy(), None

            valid_ints = set()
            try:
                as_int = int(chrom_token)
                valid_ints.add(as_int)
            except Exception:
                pass
            for key, mapped in (chrom_mapping or {}).items():
                kt = key.lower().replace("chr", "")
                if kt == chrom_token.lower().replace("chr", ""):
                    try:
                        valid_ints.add(int(mapped))
                    except Exception:
                        pass

            if valid_ints:
                mask = df_loc["Chr"].astype(int).isin(valid_ints)
            else:
                mask = pd.Series([True] * len(df_loc), index=df_loc.index)

            if start is not None:
                mask &= (df_loc["ChrPos"].astype(int) >= int(start))
            if end is not None:
                mask &= (df_loc["ChrPos"].astype(int) <= int(end))

            df_f = df_loc.loc[mask].copy()
            desc = f"{chrom_token}"
            if start is not None and end is not None:
                desc += f":{start}-{end}"
            return df_f, desc

        df_work = df.copy()
        if limit not in ("", None):
            try:
                lim = int(limit)
                if lim > 0:
                    df_work = df_work.head(lim)
            except Exception:
                pass
        if "ChrPos" in df_work.columns:
            df_work = df_work.dropna(subset=["ChrPos"])

        region_desc = None
        if region:
            df_work, region_desc = _apply_region_filter(df_work, region)

        if region and region_only_csv is None:
            base = os.path.splitext(manhatten_plot_name)[0]
            suffix = region_desc.replace(":", "_").replace("-", "_") if region_desc else "region"
            region_only_csv = f"{base}.{suffix}.csv"

        try:
            if region_only_csv and len(df_work) > 0:
                out_df = df_work.copy()
                flipped = {v: k for k, v in (chrom_mapping or {}).items()}
                try:
                    out_df["Chr"] = out_df["Chr"].astype(float).replace(flipped)
                except Exception:
                    pass
                out_df.to_csv(region_only_csv, index=False)
        except Exception:
            pass

        if algorithm in ("FaST-LMM", "Linear regression"):
            if "PValue" in df_work.columns:
                df_work = df_work[df_work["PValue"] != 0.0]
            df_work = df_work.sort_values(by=["Chr", "ChrPos"])
            df_work["Chr"] = df_work["Chr"].astype(int)
            df_work["ChrPos"] = df_work["ChrPos"].astype(int)

            plt_params = {
                "font.sans-serif": "Arial",
                "legend.fontsize": 14,
                "axes.titlesize": 18,
                "axes.labelsize": 16,
                "xtick.labelsize": 14,
                "ytick.labelsize": 14,
            }
            plt.rcParams.update(plt_params)

            sugg_line = 1.0 / max(1, len(df_work.get("SNP", [])))
            gen_line = 0.05 / max(1, len(df_work.get("SNP", [])))

            f, ax = plt.subplots(figsize=(12, 5), facecolor="w", edgecolor="k")
            flipped_dict = {value: key for key, value in (chrom_mapping or {}).items()}
            try:
                df_plot = df_work.copy()
                df_plot["Chr"] = df_plot["Chr"].astype(float).replace(flipped_dict)
            except Exception:
                df_plot = df_work

            plot_title = f"Manhattan Plot ({algorithm})"
            if title_suffix:
                plot_title += f" — {title_suffix}"
            if region_desc:
                plot_title += f" [{region_desc}]"

            _ = gv.manhattanplot(
                data=df_plot,
                chrom="Chr",
                pos="ChrPos",
                pv="PValue",
                snp="SNP",
                marker=".",
                color=["#4297d8", "#eec03c", "#423496", "#495227", "#d50b6f", "#e76519", "#d580b7", "#84d3ac"],
                sign_marker_color="r",
                title=plot_title,
                xlabel="Chromosome",
                ylabel=r"$-log_{10}{(P)}$",
                sign_line_cols=["#D62728", "#2CA02C"],
                hline_kws={"linestyle": "--", "lw": 1.3},
                text_kws={"fontsize": 12, "arrowprops": dict(arrowstyle="-", color="k", alpha=0.6)},
                logp=True,
                ax=ax,
                xticklabel_kws={"rotation": "vertical"},
                suggestiveline=sugg_line,
                genomewideline=gen_line,
            )
            plt.tight_layout(pad=1)
            manh_out = manhatten_plot_name if not region_desc else manhatten_plot_name.replace(
                ".png", f".{region_desc.replace(':', '_').replace('-', '_')}.png"
            )
            plt.savefig(manh_out)
            plt.savefig(manh_out.replace("manhatten_plot", "manhatten_plot_high"), dpi=300)

            if "PValue" in df_work.columns and len(df_work) > 0:
                f, ax = plt.subplots(figsize=(6, 6), facecolor="w", edgecolor="k")
                qq_title = f"QQ Plot ({algorithm})"
                if title_suffix:
                    qq_title += f" — {title_suffix}"
                if region_desc:
                    qq_title += f" [{region_desc}]"
                _ = gv.qqplot(
                    data=df_work["PValue"],
                    marker="o",
                    title=qq_title,
                    xlabel=r"Expected $-log_{10}{(P)}$",
                    ylabel=r"Observed $-log_{10}{(P)}$",
                    ax=ax,
                )
                plt.tight_layout(pad=1)
                qq_out = qq_plot_name if not region_desc else qq_plot_name.replace(
                    ".png", f".{region_desc.replace(':', '_').replace('-', '_')}.png"
                )
                plt.savefig(qq_out)
                plt.savefig(qq_out.replace("qq_plot", "qq_plot_high"), dpi=300)
            else:
                qq_out = None

        else:
            plt_params = {
                "font.sans-serif": "Arial",
                "legend.fontsize": 10,
                "axes.titlesize": 14,
                "axes.labelsize": 12,
                "xtick.labelsize": 10,
                "ytick.labelsize": 10,
            }
            plt.rcParams.update(plt_params)

            df_work = df_work.sort_values(by=["Chr", "ChrPos"])
            df_work["Chr"] = df_work["Chr"].astype(int)
            df_work["ChrPos"] = df_work["ChrPos"].astype(int)

            flipped_dict = {value: key for key, value in (chrom_mapping or {}).items()}
            try:
                df_plot = df_work.copy()
                df_plot["Chr"] = df_plot["Chr"].astype(float).replace(flipped_dict)
            except Exception:
                df_plot = df_work

            f, ax = plt.subplots(figsize=(12, 6), facecolor="w", edgecolor="k")
            algorithm2 = algorithm.replace(" (AI)", "")
            plot_title = f"Manhattan Plot ({algorithm2})"
            if title_suffix:
                plot_title += f" — {title_suffix}"
            if region_desc:
                plot_title += f" [{region_desc}]"

            _ = gv.manhattanplot(
                data=df_plot,
                chrom="Chr",
                pos="ChrPos",
                pv="PValue",
                snp="SNP",
                logp=False,
                title=plot_title,
                color=["#4297d8", "#eec03c", "#423496", "#495227", "#d50b6f", "#e76519", "#d580b7", "#84d3ac"],
                xlabel="Chromosome",
                ylabel=r"SNP effect",
                xticklabel_kws={"rotation": "vertical"},
            )
            plt.tight_layout(pad=1)
            manh_out = manhatten_plot_name if not region_desc else manhatten_plot_name.replace(
                ".png", f".{region_desc.replace(':', '_').replace('-', '_')}.png"
            )
            plt.savefig(manh_out)
            plt.savefig(manh_out.replace("manhatten_plot", "manhatten_plot_high"), dpi=300)
            qq_out = None

        return {
            "manhattan_png": manh_out,
            "qq_png": qq_out,
            "region_csv": region_only_csv if (region and len(df_work) > 0) else None,
            "filtered_count": int(len(df_work)),
            "region_desc": region_desc,
        }

    #------------------------------
     #Annotations
    #------------------------------
    def _smart_open(self, path):
        return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "r")

    def _parse_gtf_attributes(self, attr_str: str):
        out = {}
        # Handles GTF/GFF attribute styles
        for part in attr_str.strip().split(";"):
            part = part.strip()
            if not part:
                continue
            if "=" in part:  # GFF style
                k, v = part.split("=", 1)
            else:  # GTF style key "value"
                seg = part.split(" ", 1)
                if len(seg) == 1:
                    k, v = seg[0], ""
                else:
                    k, v = seg[0], seg[1].strip()
            v = v.strip().strip('"').strip("'")
            out[k] = v
        return out

    def build_gtf_index(self, gtf_path: str, log_fn=None):
        """
        Build a minimal index of TSS positions from a GTF/GFF:
          returns: dict {chrom: [(tss_pos, gene_id, gene_name, feature_type), ...]} sorted by tss_pos
        """
        if log_fn: log_fn(f"[ANNOT] Building GTF/GFF index: {gtf_path}")
        idx = {}
        n_lines = 0
        with self._smart_open(gtf_path) as fh:
            for ln in fh:
                if not ln or ln.startswith("#"):
                    continue
                n_lines += 1
                parts = ln.rstrip("\n").split("\t")
                if len(parts) < 9:
                    continue
                chrom, source, feature, start, end, score, strand, frame, attrs = parts
                if feature not in ("gene", "transcript", "mRNA"):
                    # We only need a representative TSS anchor; genes/transcripts are fine.
                    continue
                try:
                    start_i = int(start)
                    end_i = int(end)
                except Exception:
                    continue
                tss = start_i if strand == "+" else end_i
                a = self._parse_gtf_attributes(attrs)
                gene_id = a.get("gene_id") or a.get("ID") or ""
                gene_name = a.get("gene_name") or a.get("Name") or gene_id
                if not gene_id and not gene_name:
                    continue
                lst = idx.setdefault(chrom, [])
                lst.append((tss, gene_id, gene_name, feature))
        # sort
        for chrom in list(idx.keys()):
            idx[chrom].sort(key=lambda x: x[0])
        if log_fn:
            tot = sum(len(v) for v in idx.values())
            log_fn(f"[ANNOT] Indexed {tot:,} TSS anchors across {len(idx)} chromosomes.")
        return idx

    def _nearest_within_kb(self, sorted_list, pos, window_bp):
        """
        Given sorted_list = [(tss, gene_id, gene_name, feat), ...] (sorted by tss),
        return nearest entry and distance if within window_bp. Else (None, None).
        """
        if not sorted_list:
            return None, None
        tss_only = [t[0] for t in sorted_list]
        i = bisect_left(tss_only, pos)
        cand = []
        if i < len(sorted_list):
            cand.append(sorted_list[i])
        if i > 0:
            cand.append(sorted_list[i - 1])
        best = None
        best_dist = None
        for t in cand:
            d = abs(t[0] - pos)
            if best is None or d < best_dist:
                best = t
                best_dist = d
        if best is not None and best_dist <= window_bp:
            return best, best_dist
        return None, None

    def annotate_gwas_results(self, gwas_csv: str, gtf_path: str, out_csv: str,
                              window_kb: int = 50, log_fn=None):
        """
        Annotate GWAS table (expects columns: Chr, ChrPos, SNP or similar)
        with nearest gene within +/- window_kb around TSS.
        """
        if log_fn: log_fn(f"[ANNOT] Loading GWAS table: {gwas_csv}")
        df = pd.read_csv(gwas_csv)
        # Normalize expected columns
        if "Chr" not in df.columns or "ChrPos" not in df.columns:
            # Try alternative spellings
            if "Chr" not in df.columns:
                raise ValueError("GWAS table must contain 'Chr' column")
            if "ChrPos" not in df.columns:
                if "BP" in df.columns:
                    df["ChrPos"] = df["BP"]
                elif "Position" in df.columns:
                    df["ChrPos"] = df["Position"]
                else:
                    raise ValueError("GWAS table must contain 'ChrPos' (or BP/Position)")

        # Read index
        idx = self.build_gtf_index(gtf_path, log_fn=log_fn)
        window_bp = int(window_kb) * 1000

        ann_gene = []
        ann_feat = []
        ann_dist = []
        ann_hit = []

        if log_fn: log_fn(f"[ANNOT] Annotating with ±{window_kb} kb around TSS")
        for _, r in df.iterrows():
            chrom = str(r["Chr"])
            try:
                pos = int(r["ChrPos"])
            except Exception:
                pos = None
            if pos is None or chrom not in idx:
                ann_gene.append("")
                ann_feat.append("")
                ann_dist.append(np.nan)
                ann_hit.append(False)
                continue
            best, dist = self._nearest_within_kb(idx[chrom], pos, window_bp)
            if best is None:
                ann_gene.append("")
                ann_feat.append("")
                ann_dist.append(np.nan)
                ann_hit.append(False)
            else:
                tss, gid, gname, feat = best
                ann_gene.append(gname if gname else gid)
                ann_feat.append(feat)
                ann_dist.append(int(dist))
                ann_hit.append(True)

        df["NearestGene"] = ann_gene
        df["NearestFeature"] = ann_feat
        df["DistanceToTSS"] = ann_dist
        df["WithinWindow"] = ann_hit

        df.to_csv(out_csv, index=False)
        if log_fn: log_fn(f"[ANNOT] Saved: {out_csv}")
        return out_csv

    # ---------------------------
    # Helpers (SAIGE)
    # ---------------------------
    def _load_pheno_3col(self, pheno_file: str) -> pd.DataFrame:
        ext = os.path.splitext(pheno_file)[1].lower()
        if ext in (".xlsx", ".xls"):
            df = pd.read_excel(pheno_file, header=None)
            df = df.iloc[:, :3]
            df.columns = ["FID", "IID", "y"]
            return df
        raw = _read_lines_with_fallback(pheno_file)
        rows = []
        for ln in raw:
            ln = ln.strip()
            if not ln or ln.startswith("#"):
                continue
            parts = ln.replace("\u00A0", " ").split()
            if len(parts) < 3:
                parts = [p for p in ln.replace(",", " ").replace(";", " ").replace("|", " ").split(" ") if p]
            if len(parts) >= 3:
                rows.append(parts[:3])
        df = pd.DataFrame(rows, columns=["FID", "IID", "y"])
        return df

    def _load_covariates(self, cov_file: str) -> pd.DataFrame:
        raw = _read_lines_with_fallback(cov_file)
        rows = []
        for ln in raw:
            ln = ln.strip()
            if not ln or ln.startswith("#"):
                continue
            parts = ln.replace("\u00A0", " ").split()
            if len(parts) < 2:
                parts = [p for p in ln.replace(",", " ").replace(";", " ").replace("|", " ").split(" ") if p]
            rows.append(parts)
        max_cols = max(len(r) for r in rows)
        cols = ["FID", "IID"] + [f"cov{i}" for i in range(1, max_cols - 1)]
        df = pd.DataFrame([r + [""] * (max_cols - len(r)) for r in rows], columns=cols)
        return df
