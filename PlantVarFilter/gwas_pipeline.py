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

from PlantVarFilter.helpers import HELPERS

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

    def validate_gwas_input_files(self, bed_file, pheno_file):
        fam_path = bed_file.replace('.bed', '.fam')
        if not os.path.exists(fam_path):
            return False, f"FAM file not found: {fam_path}"
        if not os.path.exists(pheno_file):
            return False, f"Phenotypic file not found: {pheno_file}"

        fam_data = _read_lines_with_fallback(fam_path)
        fam_ids = []
        for line in fam_data:
            parts = line.strip().split()
            if len(parts) > 1:
                fam_ids.append(parts[0])
            else:
                return False, "FAM file is not whitespace-delimited."

        ext = os.path.splitext(pheno_file)[1].lower()
        pheno_ids = []
        columns_seen = None

        if ext in (".xlsx", ".xls"):
            try:
                dfp = pd.read_excel(pheno_file, header=None)
            except Exception as e:
                return False, f"Failed to read Excel phenotype: {e}"

            dfp = dfp.dropna(how="all")
            if dfp.shape[1] < 3:
                return False, "Invalid phenotypic Excel file. Needs 3 columns (FID IID Value)."

            dfp = dfp.iloc[:, :3]
            columns_seen = 3
            pheno_ids = [str(x) for x in dfp.iloc[:, 0].astype(str).tolist()]
        else:
            raw_lines = _read_lines_with_fallback(pheno_file)

            def split_smart(s: str):
                s = s.replace("\u00A0", " ").strip()
                if not s or s.startswith("#"):
                    return None
                candidates = []
                candidates.append(s.split())
                candidates.append(s.split(","))
                candidates.append(s.split(";"))
                candidates.append(s.split("\t"))
                candidates.append(s.split("|"))
                candidates.append(s.split(":"))
                cleaned = []
                for parts in candidates:
                    parts = [p.strip().strip('"').strip("'") for p in parts if p.strip() != ""]
                    cleaned.append(parts)
                best = max(cleaned, key=lambda x: len(x))
                return best if len(best) > 0 else None

            rows = []
            for ln in raw_lines:
                parts = split_smart(ln)
                if parts is None:
                    continue
                rows.append(parts)
                if columns_seen is None:
                    columns_seen = len(parts)

            if not rows:
                return False, "Phenotypic file appears empty."

            from collections import Counter
            cnt = Counter(len(r) for r in rows)
            common_cols = cnt.most_common(1)[0][0]

            if common_cols < 3:
                return False, "Invalid phenotypic file. File should contain 3 columns (FID IID Value)."

            columns_seen = common_cols
            pheno_ids = [r[0] for r in rows if len(r) >= common_cols]

        check_ids = any(x in pheno_ids for x in fam_ids)

        if check_ids and columns_seen == 3:
            return True, "Input files validated."
        elif columns_seen != 3:
            return False, "Invalid phenotypic file. File should contain 3 columns (FID IID Value)."
        else:
            return False, "Phenotypic IDs do not match .fam file IDs."

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
    def plot_gwas(self, df, limit, algorithm, manhatten_plot_name, qq_plot_name, chrom_mapping):
        if algorithm in ('FaST-LMM', 'Linear regression'):
            if limit != '':
                df = df.head(int(limit))
            df = df.dropna(subset=['ChrPos'])

            sugg_line = 1 / len(df['SNP'])
            gen_line = 0.05 / len(df['SNP'])

            df = df.sort_values(by=['Chr', 'ChrPos'])
            df = df[df['PValue'] != 0.0]
            df['Chr'] = df['Chr'].astype(int)
            df['ChrPos'] = df['ChrPos'].astype(int)

            plt_params = {
                "font.sans-serif": "Arial",
                "legend.fontsize": 14,
                "axes.titlesize": 18,
                "axes.labelsize": 16,
                "xtick.labelsize": 14,
                "ytick.labelsize": 14
            }
            plt.rcParams.update(plt_params)

            f, ax = plt.subplots(figsize=(12, 5), facecolor="w", edgecolor="k")
            flipped_dict = {value: key for key, value in chrom_mapping.items()}
            df['Chr'] = df['Chr'].astype(float).replace(flipped_dict)

            _ = gv.manhattanplot(
                data=df, chrom='Chr', pos="ChrPos", pv="PValue", snp="SNP", marker=".",
                color=['#4297d8', '#eec03c', '#423496', '#495227', '#d50b6f', '#e76519', '#d580b7', '#84d3ac'],
                sign_marker_color="r", title=f"Manhattan Plot ({algorithm})",
                xlabel="Chromosome", ylabel=r"$-log_{10}{(P)}$",
                sign_line_cols=["#D62728", "#2CA02C"],
                hline_kws={"linestyle": "--", "lw": 1.3},
                text_kws={"fontsize": 12, "arrowprops": dict(arrowstyle="-", color="k", alpha=0.6)},
                logp=True, ax=ax, xticklabel_kws={"rotation": "vertical"},
                suggestiveline=sugg_line, genomewideline=gen_line
            )
            plt.tight_layout(pad=1)
            plt.savefig(manhatten_plot_name)
            plt.savefig(manhatten_plot_name.replace('manhatten_plot', 'manhatten_plot_high'), dpi=300)

            f, ax = plt.subplots(figsize=(6, 6), facecolor="w", edgecolor="k")
            _ = gv.qqplot(
                data=df["PValue"], marker="o", title=f"QQ Plot ({algorithm})",
                xlabel=r"Expected $-log_{10}{(P)}$", ylabel=r"Observed $-log_{10}{(P)}$", ax=ax
            )
            plt.tight_layout(pad=1)
            plt.savefig(qq_plot_name)
            plt.savefig(qq_plot_name.replace('qq_plot', 'qq_plot_high'), dpi=300)

        else:
            plt_params = {
                "font.sans-serif": "Arial",
                "legend.fontsize": 10,
                "axes.titlesize": 14,
                "axes.labelsize": 12,
                "xtick.labelsize": 10,
                "ytick.labelsize": 10
            }
            plt.rcParams.update(plt_params)

            df = df.sort_values(by=['Chr', 'ChrPos'])
            df['Chr'] = df['Chr'].astype(int)
            df['ChrPos'] = df['ChrPos'].astype(int)

            flipped_dict = {value: key for key, value in chrom_mapping.items()}
            df['Chr'] = df['Chr'].astype(float).replace(flipped_dict)

            f, ax = plt.subplots(figsize=(12, 6), facecolor="w", edgecolor="k")
            algorithm2 = algorithm.replace(' (AI)', '')
            _ = gv.manhattanplot(
                data=df, chrom='Chr', pos="ChrPos", pv="PValue", snp="SNP", logp=False,
                title=f"Manhatten Plot ({algorithm2})",
                color=['#4297d8', '#eec03c', '#423496', '#495227', '#d50b6f', '#e76519', '#d580b7', '#84d3ac'],
                xlabel="Chromosome", ylabel=r"SNP effect", xticklabel_kws={"rotation": "vertical"}
            )
            plt.tight_layout(pad=1)
            plt.savefig(manhatten_plot_name)
            plt.savefig(manhatten_plot_name.replace('manhatten_plot', 'manhatten_plot_high'), dpi=300)

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
