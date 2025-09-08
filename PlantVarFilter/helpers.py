import pandas as pd
import numpy as np
from datetime import datetime
import os
import shutil
import configparser


class HELPERS:
    # def duplicate_column(self, input_file, output_file):
    #     """Depreciated."""
    #     df = pd.read_csv(input_file, header=None)
    #     df[1] = df[0]
    #     df.to_csv(output_file, index=False, header=False, sep=' ')

    def replace_with_integers(self, input_file):
        """Replace string chromosome names with integers."""
        mapping = {}
        current_integer = 1.0
        with open(input_file, 'r', errors="ignore") as infile:
            for line in infile:
                parts = line.strip().split('\t')
                col1_value = parts[0]
                try:
                    col1_value = int(col1_value)
                except ValueError:
                    if col1_value in mapping:
                        parts[0] = str(mapping[col1_value])
                    else:
                        mapping[col1_value] = current_integer
                        parts[0] = str(current_integer)
                        current_integer += 1
        return mapping

    def get_timestamp(self):
        """Get timestamp."""
        now = datetime.now()
        dt_string = now.strftime("%d%m%Y_%H%M%S")
        return dt_string

    def save_raw_data(self, bed, pheno):
        """Save the SNP matrix."""
        np.save('snp', bed.read().val)
        np.savez_compressed('snp.npz', bed.read().val)
        np.save('pheno', pheno.read().val)

    def save_settings(self, default_path):
        """Save the user settings."""
        config = configparser.ConfigParser()
        config['DefaultSettings'] = {'path': default_path, 'algorithm': 'FaST-LMM'}
        with open('settings.ini', 'w') as configfile:
            config.write(configfile)

    def get_settings(self, setting):
        """Get the user settings."""
        config = configparser.ConfigParser()
        config.read('settings.ini')
        return config['DefaultSettings'][setting]

    # ---------- helpers ----------
    def _safe_copy(self, src_path, dst_dir, add_log):
        """Copy a file if it exists. Returns True if copied."""
        if not src_path:
            return False
        if os.path.exists(src_path):
            try:
                os.makedirs(dst_dir, exist_ok=True)
                dst_path = os.path.join(dst_dir, os.path.basename(src_path))
                shutil.copy(src_path, dst_path)
                add_log(f"File saved: {os.path.basename(src_path)}")
                return True
            except Exception as e:
                add_log(f"[SAVE] Could not copy {src_path}: {e}", error=True)
        return False

    def _make_top_snps_csv(self, gwas_result_path, out_path, add_log, top_n=100):
        """Create a concise CSV of top SNPs according to available columns."""
        try:
            if not (gwas_result_path and os.path.exists(gwas_result_path)):
                return
            df = pd.read_csv(gwas_result_path)
            if 'PValue' in df.columns:
                df2 = df.dropna(subset=['PValue']).sort_values('PValue', ascending=True).head(top_n)
            elif 'SNP effect' in df.columns:
                df2 = df.dropna(subset=['SNP effect']).sort_values('SNP effect', ascending=False).head(top_n)
            else:
                df2 = df.head(top_n)
            df2.to_csv(out_path, index=False)
            add_log(f"File saved: {os.path.basename(out_path)}")
        except Exception as e:
            add_log(f"[TopSNPs] Failed to create top SNPs: {e}", warn=True)

    def save_results(
            self,
            current_dir,
            save_dir,
            gwas_result_name,
            gwas_result_name_top,
            manhatten_plot_name,
            qq_plot_name,
            algorithm,
            genomic_predict_name,
            gp_plot_name,
            gp_plot_name_scatter,
            add_log,
            settings_lst,
            pheno_stats_name,
            geno_stats_name
    ):
        """Store analysis files under a writable folder; never drop files in '/' or code dir."""

        def _is_writable(path: str) -> bool:
            try:
                if not path:
                    return False
                path = os.path.abspath(path)
                return os.path.isdir(path) and os.access(path, os.W_OK)
            except Exception:
                return False

        # —— choose a safe base dir ——
        base_dir_candidates = [
            save_dir,
            current_dir,
            os.path.expanduser("~/Desktop"),
            os.path.expanduser("~")
        ]
        base_dir = None
        for cand in base_dir_candidates:
            if _is_writable(cand) and os.path.abspath(cand) != "/":
                base_dir = os.path.abspath(cand)
                break
        if base_dir is None:
            # last resort: process working dir (may still fail, but we tried)
            base_dir = os.getcwd()
            add_log(f"[SAVE] No writable directory found, falling back to {base_dir}", warn=True)
        else:
            if save_dir and os.path.abspath(save_dir) != base_dir:
                add_log(f"[SAVE] Selected folder not writable; using {base_dir} instead.", warn=True)

        # —— create destination folder ——
        safe_algo = str(algorithm or "Run").replace(" ", "_")
        ts = self.get_timestamp() + "_" + safe_algo
        dest_dir = os.path.join(base_dir, ts)
        try:
            os.makedirs(dest_dir, exist_ok=True)
        except Exception as e:
            add_log(f"Can not create folder '{dest_dir}': {e}", error=True)
            return base_dir  # bail out with something sensible

        # —— build helper for resolving relative inputs (from current_dir) ——
        def _abs_or_join(p):
            if not p:
                return None
            return p if os.path.isabs(p) else os.path.join(current_dir, os.path.basename(p))

        # —— create Top10k inside dest_dir (not next to code) ——
        try:
            if gwas_result_name and os.path.exists(gwas_result_name):
                df = pd.read_csv(gwas_result_name)
                df.head(10000).to_csv(os.path.join(dest_dir, "gwas_results_top10000.csv"), index=False)
        except Exception as e:
            add_log(f"[SAVE] Could not create top 10000 GWAS CSV: {e}", warn=True)

        # —— create TopSNPs inside dest_dir ——
        try:
            top_snps_path = os.path.join(dest_dir, "gwas_top_snps.csv")
            self._make_top_snps_csv(gwas_result_name, top_snps_path, add_log, top_n=100)
        except Exception as e:
            add_log(f"[TopSNPs] Could not generate: {e}", warn=True)

        # —— resolve optional “high-res / validation” companion outputs ——
        manh_high = _abs_or_join(str(manhatten_plot_name or "").replace('manhatten_plot', 'manhatten_plot_high'))
        qq_high = _abs_or_join(str(qq_plot_name or "").replace('qq_plot', 'qq_plot_high'))
        gp_scatter_high = _abs_or_join(
            str(gp_plot_name_scatter or "").replace('GP_scatter_plot', 'GP_scatter_plot_high'))
        gp_ba_high = _abs_or_join(str(gp_plot_name or "").replace('Bland_Altman_plot', 'Bland_Altman_plot_high'))
        gp_validation = _abs_or_join(str(genomic_predict_name or "").replace('.csv', '_valdation.csv'))

        # —— files to copy into dest_dir ——
        src_files = [
            _abs_or_join(gwas_result_name),
            _abs_or_join(genomic_predict_name),
            _abs_or_join(manhatten_plot_name),
            _abs_or_join(qq_plot_name),
            _abs_or_join(gp_plot_name),
            _abs_or_join(gp_plot_name_scatter),
            manh_high, qq_high, gp_scatter_high, gp_ba_high, gp_validation,
            _abs_or_join(pheno_stats_name),
            _abs_or_join(geno_stats_name),
        ]

        for src in src_files:
            self._safe_copy(src, dest_dir, add_log)

        # —— write run log inside dest_dir ——
        try:
            log_path = os.path.join(dest_dir, 'log.txt')
            with open(log_path, 'w', encoding='utf-8') as log_file:
                def _get(i, default=""):
                    try:
                        return settings_lst[i]
                    except Exception:
                        return default

                log_file.write('Algorithm: ' + str(_get(0, 'N/A')))
                log_file.write('\nBed file used: ' + str(_get(1, 'N/A')))
                log_file.write('\nPheno file used: ' + str(_get(2, 'N/A')))
                try:
                    tr = 100 - (100 * float(_get(3, 0.3)))
                except Exception:
                    tr = 'N/A'
                log_file.write('\nTraining size: ' + str(tr))
                log_file.write('\nNr of trees: ' + str(_get(4, 'N/A')))
                log_file.write('\nNr of models: ' + str(_get(5, 'N/A')))
                log_file.write('\nMax depth: ' + str(_get(6, 'N/A')))
        except Exception as e:
            add_log(f"[SAVE] Could not write log.txt: {e}", warn=True)

        add_log(f"[INFO] Results saved in: {dest_dir}")
        return dest_dir

    def merge_gp_models(self, dataframes):
        """Merge all GP ML models."""
        df_combined = pd.concat(dataframes)
        df_result = df_combined.groupby(['ID1', 'BED_ID2'])['Predicted_Value'].mean().reset_index()
        df_result = pd.merge(dataframes[0], df_result, on='ID1', how='left')
        df_result = df_result.drop(['Predicted_Value_x', 'BED_ID2_y', 'Pheno_ID2'], axis=1)
        df_result = df_result.rename(columns={'Predicted_Value_y': 'Mean_Predicted_Value'})
        df_result['Difference'] = (df_result['Pheno_Value'] - df_result['Mean_Predicted_Value']).abs()
        df_result['Difference'] = df_result['Difference'].round(decimals=3)
        df_result['Mean_Predicted_Value'] = df_result['Mean_Predicted_Value'].round(decimals=3)
        return df_result

    def merge_models(self, dataframes, method):
        """Merge all GWAS ML models."""
        df_combined = pd.concat(dataframes)
        df_combined = df_combined[df_combined['PValue'] > 0]
        df_grouped = df_combined.groupby(['SNP', 'Chr', 'ChrPos'])['PValue'].agg([method, 'std']).reset_index()
        df_grouped = df_grouped.rename(columns={method: 'PValue_x', 'std': 'PValue_sd'})
        df_result = df_grouped.rename(columns={'PValue_x': 'PValue'})
        df_result['Chr'] = df_result['Chr'].astype(int)
        df_result['ChrPos'] = df_result['ChrPos'].astype(int)
        df_result = df_result.sort_values(by=['Chr', 'ChrPos'])
        return df_result
