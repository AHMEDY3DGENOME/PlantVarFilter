# gwas_pipeline.py
import os
import sys
import time
import logging
import subprocess

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import Ridge, ridge_regression  # noqa: F401
import xgboost as xgb

from fastlmm.association import single_snp, single_snp_linreg, single_snp_scale  # noqa: F401
from pysnptools.util import log_in_place
import geneview as gv

from PlantVarFilter.helpers import HELPERS

# Use non-interactive backend for plotting
plt.switch_backend('Agg')


def _read_lines_with_fallback(path):
    """
    Read text lines with robust encoding fallback.
    Tries: utf-8 -> cp1256 -> cp1252 -> latin-1 -> utf-8-sig -> (utf-8 with errors='replace').
    """
    encodings = ('utf-8', 'cp1256', 'cp1252', 'latin-1', 'utf-8-sig')
    for enc in encodings:
        try:
            with open(path, 'r', encoding=enc) as f:
                return f.readlines()
        except UnicodeDecodeError:
            continue
    with open(path, 'r', encoding='utf-8', errors='replace') as f:
        return f.readlines()


class GWAS:
    """GWAS class."""

    def __init__(self):
        # self.gwas_ai = GWASAI()
        self.helper = HELPERS()

    # ---------------------------
    # VCF -> PLINK BED Conversion
    # ---------------------------
    def vcf_to_bed(self, vcf_file, id_file, file_out, maf, geno):
        """
        Converts VCF (prefer .vcf.gz + index) to PLINK BED/BIM/FAM
        - Same signature for UI compatibility.
        - Strengthened execution: capture stdout/stderr and check returncode + file creation.
        - On failure: raise with full STDERR so UI can display it.
        """
        script_dir = os.path.dirname(__file__)  # absolute dir the script is in

        # Pick bundled PLINK binary per platform
        if sys.platform.startswith('win'):
            abs_file_path = os.path.join(script_dir, "windows", "plink")
        elif sys.platform.startswith('linux'):
            abs_file_path = os.path.join(script_dir, "linux", "plink")
            try:
                os.chmod(abs_file_path, 0o755)
            except Exception:
                pass
        elif sys.platform.startswith('darwin'):
            abs_file_path = os.path.join(script_dir, "mac", "plink")
            try:
                os.chmod(abs_file_path, 0o755)
            except Exception:
                pass
        else:
            raise RuntimeError("Unsupported platform for PLINK binaries.")

        # Ensure input exists
        if not vcf_file or not os.path.exists(vcf_file):
            raise RuntimeError(f"VCF not found: {vcf_file}")

        # Reasonable defaults (adjust if needed)
        threads = "4"
        memory_mb = "16000"  # 16 GB

        # Build PLINK command
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

        # Run and capture output
        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        stdout, stderr = proc.stdout, proc.stderr

        # Verify output files exist
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
            full_msg = "\n".join(msg)
            logging.error(full_msg)
            raise RuntimeError(full_msg)

        logging.info("[PLINK] conversion completed.")
        return stdout or "PLINK conversion completed."

    # ---------------------------
    # Utilities
    # ---------------------------
    def filter_out_missing(self, bed):
        """Drop SNPs that are entirely NaN in batches to reduce memory."""
        sid_batch_size = 1000  # reduce if memory is tight
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
        """
        Validate inputs:
        - Phenotypic must have 3 columns (FID IID Value) separated by whitespace or a common delimiter
          (',', ';', tab, '|', ':') or be an Excel file (.xlsx/.xls).
        - IDs in .fam must intersect with IDs in phenotype.
        """
        fam_path = bed_file.replace('.bed', '.fam')
        if not os.path.exists(fam_path):
            return False, f"FAM file not found: {fam_path}"
        if not os.path.exists(pheno_file):
            return False, f"Phenotypic file not found: {pheno_file}"

        # --- FAM: always whitespace ---
        fam_data = _read_lines_with_fallback(fam_path)
        fam_ids = []
        for line in fam_data:
            parts = line.strip().split()
            if len(parts) > 1:
                fam_ids.append(parts[0])
            else:
                return False, "FAM file is not whitespace-delimited."

        # --- PHENO: Excel or text ---
        ext = os.path.splitext(pheno_file)[1].lower()

        pheno_ids = []
        columns_seen = None

        if ext in (".xlsx", ".xls"):
            # Excel phenotype: take first 3 columns as FID IID Value
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
            # Text path: robust per-line parsing with multiple delimiters
            raw_lines = _read_lines_with_fallback(pheno_file)

            def split_smart(s: str):
                # normalize NBSP to space
                s = s.replace("\u00A0", " ").strip()
                if not s or s.startswith("#"):
                    return None
                candidates = []
                candidates.append(s.split())     # whitespace (space/tab)
                candidates.append(s.split(","))  # comma
                candidates.append(s.split(";"))  # semicolon
                candidates.append(s.split("\t")) # explicit tab
                candidates.append(s.split("|"))  # pipe
                candidates.append(s.split(":"))  # colon
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

            # choose the most common column count across rows (tolerate mixed delimiters)
            from collections import Counter
            cnt = Counter(len(r) for r in rows)
            common_cols = cnt.most_common(1)[0][0]

            if common_cols < 3:
                return False, "Invalid phenotypic file. File should contain 3 columns (FID IID Value)."

            columns_seen = common_cols
            pheno_ids = [r[0] for r in rows if len(r) >= common_cols]

        # Any shared IDs?
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
                     gwas_result_name, algorithm, bed_file, cov_file, gb_goal):
        """GWAS using LMM and linear regression methods from fast-lmm library."""
        t1 = time.time()
        if gb_goal == 0:
            gb_goal = None

        if algorithm == 'FaST-LMM':
            if cov_file:
                df_lmm_gwas = single_snp(
                    bed_fixed, pheno, output_file_name=gwas_result_name,
                    covar=cov_file, GB_goal=gb_goal
                )
            else:
                df_lmm_gwas = single_snp(
                    bed_fixed, pheno, output_file_name=gwas_result_name, GB_goal=gb_goal
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

        # Copy for plotting (may need numeric Chr)
        df_plot = df_lmm_gwas.copy(deep=True)

        if len(chrom_mapping) > 0:
            reversed_chrom_map = {value: key for key, value in chrom_mapping.items()}
            df_lmm_gwas["Chr"] = df_lmm_gwas["Chr"].apply(lambda x: reversed_chrom_map[x])

        df_lmm_gwas.to_csv(gwas_result_name, index=False)
        TOP_N = 50
        if "PValue" in df_lmm_gwas.columns:
            top_snps = df_lmm_gwas.nsmallest(TOP_N, "PValue")
            # top_snps.to_csv("gwas_top_snp.csv", index=False)
            add_log(f"Top {TOP_N} SNPs Saved to gwas_top_snps.csv")
        t3 = round((time.time() - t1) / 60, 2)
        add_log('Final run time (minutes): ' + str(t3))
        return df_lmm_gwas, df_plot

    # ---------------------------
    # GWAS (XGBoost / RF / Ridge)
    # ---------------------------
    def run_gwas_xg(self, bed_fixed, pheno, bed_file, test_size, estimators,
                    gwas_result_name, chrom_mapping, add_log, model_nr, max_dep_set, nr_jobs, method):
        """GWAS using XGBoost with cross validation (feature importance used as SNP effect)."""
        t1 = time.time()
        dataframes = []

        # More robust BIM read (tabs/spaces)
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
        """GWAS using Random Forest with cross validation."""
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
        """GWAS using Ridge Regression with cross validation."""
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
            X_test = scaler.transform(X_test)  # noqa: F841

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
    # Plotting
    # ---------------------------
    def plot_gwas(self, df, limit, algorithm, manhatten_plot_name, qq_plot_name, chrom_mapping):
        """Manhattan & QQ plots."""
        # Extract the topN SNPs if limit set
        if algorithm in ('FaST-LMM', 'Linear regression'):
            if limit != '':
                df = df.head(int(limit))
            df = df.dropna(subset=['ChrPos'])

            # Thresholds
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

            # Manhattan
            f, ax = plt.subplots(figsize=(12, 5), facecolor="w", edgecolor="k")  # noqa: F841
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

            # QQ
            f, ax = plt.subplots(figsize=(6, 6), facecolor="w", edgecolor="k")  # noqa: F841
            _ = gv.qqplot(
                data=df["PValue"], marker="o", title=f"QQ Plot ({algorithm})",
                xlabel=r"Expected $-log_{10}{(P)}$", ylabel=r"Observed $-log_{10}{(P)}$", ax=ax
            )
            plt.tight_layout(pad=1)
            plt.savefig(qq_plot_name)
            plt.savefig(qq_plot_name.replace('qq_plot', 'qq_plot_high'), dpi=300)

        else:
            # Effect maps (no -log10)
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

            f, ax = plt.subplots(figsize=(12, 6), facecolor="w", edgecolor="k")  # noqa: F841
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
