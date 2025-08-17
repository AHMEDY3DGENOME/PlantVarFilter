# main_gui.py  — PlantVarFilter (updated)

import os
import dearpygui.dearpygui as dpg
from dearpygui_ext import logger

from PlantVarFilter.gwas_pipeline import GWAS
from PlantVarFilter.genomic_prediction_pipeline import GenomicPrediction
from PlantVarFilter.helpers import HELPERS
from PlantVarFilter.pipeline_plots import Plot

from pysnptools.snpreader import Bed, Pheno
import pysnptools.util as pstutil


def main():
    """Entry point """
    dpg.create_context()
    app = GWASApp()
    app.run()
    dpg.destroy_context()


class GWASApp:
    def __init__(self):
        self.gwas = GWAS()
        self.helper = HELPERS()
        self.genomic_predict_class = GenomicPrediction()
        self.plot_class = Plot()

        # ---------------- Fonts ----------------
        with dpg.font_registry():
            script_dir = os.path.dirname(__file__)
            font_path = os.path.join(script_dir, "test.ttf")
            if os.path.exists(font_path):
                self.font = dpg.add_font(font_path, 40, tag="ttf-font")  # 20*2
            else:
                self.font = None

        # ---------------- Theme (عام) ----------------
        with dpg.theme() as self.our_theme:
            with dpg.theme_component(dpg.mvAll):
                dpg.add_theme_color(dpg.mvThemeCol_WindowBg, (18, 18, 22))
                dpg.add_theme_color(dpg.mvThemeCol_ChildBg, (18, 18, 22))
                dpg.add_theme_color(dpg.mvThemeCol_FrameBg, (34, 34, 40))
                dpg.add_theme_color(dpg.mvThemeCol_Text, (235, 235, 235))
                dpg.add_theme_color(dpg.mvThemeCol_Button, (52, 88, 180))
                dpg.add_theme_color(dpg.mvThemeCol_ButtonHovered, (72, 108, 200))
                dpg.add_theme_color(dpg.mvThemeCol_ButtonActive, (62, 98, 190))
                dpg.add_theme_style(dpg.mvStyleVar_WindowRounding, 8)
                dpg.add_theme_style(dpg.mvStyleVar_FrameRounding, 6)

        # ---------------- State ----------------
        self.vcf_file = ''
        self.pheno_file = ''
        self.vcf_app_data = None
        self.variants_app_data = None  # IDs file (PLINK --keep)
        self.results_directory = None
        self.bed_app_data = None
        self.pheno_app_data = None
        self.cov_app_data = None
        self.default_path = '.'

        # filenames
        self.gwas_result_name = "gwas_results.csv"
        self.gwas_result_name_top = "gwas_results_top10000.csv"
        self.genomic_predict_name = "genomic_prediction_results.csv"
        self.manhatten_plot_name = "manhatten_plot.png"
        self.qq_plot_name = "qq_plot.png"
        self.gp_plot_name = "Bland_Altman_plot.png"
        self.gp_plot_name_scatter = "GP_scatter_plot.png"
        self.pheno_stats_name = 'pheno_statistics.pdf'
        self.geno_stats_name = 'geno_statistics.pdf'

        self.snp_limit = None

        # Log
        self.log_win = dpg.add_window(label="Log", pos=(0, 635), width=1000, height=500, horizontal_scrollbar=True)
        self.logz = logger.mvLogger(self.log_win)

        # Build UI
        self.setup_gui()

    # ---------------- UI ----------------
    def setup_gui(self):
        dpg.create_viewport(title='PlantVarFilter', width=2000, height=1200, resizable=True)

        # File dialogs
        with dpg.file_dialog(directory_selector=False, show=False, callback=self.callback_vcf,
                             file_count=3, tag="file_dialog_vcf", width=700, height=400,
                             default_path=self.default_path):
            dpg.add_file_extension("Source files (*.vcf *.gz){.vcf,.gz}", color=(255, 255, 0, 255))

        with dpg.file_dialog(directory_selector=False, show=False, callback=self.callback_variants,
                             file_count=3, tag="file_dialog_variants", width=700, height=400,
                             default_path=self.default_path):
            dpg.add_file_extension("Text files (*.txt *.csv){.txt,.csv}", color=(255, 255, 0, 255))
            dpg.add_file_extension(".*")

        with dpg.file_dialog(directory_selector=False, show=False, callback=self.callback_pheno,
                             file_count=3, tag="file_dialog_pheno", width=700, height=400,
                             default_path=self.default_path):
            dpg.add_file_extension("Text files (*.txt *.csv){.txt,.csv}", color=(255, 255, 0, 255))
            dpg.add_file_extension(".*")

        with dpg.file_dialog(directory_selector=False, show=False, callback=self.callback_cov,
                             file_count=3, tag="file_dialog_cov", width=700, height=400,
                             default_path=self.default_path):
            dpg.add_file_extension("Text files (*.txt *.csv){.txt,.csv}", color=(255, 255, 0, 255))
            dpg.add_file_extension(".*")

        with dpg.file_dialog(directory_selector=False, show=False, callback=self.callback_bed,
                             file_count=3, tag="file_dialog_bed", width=700, height=400,
                             default_path=self.default_path):
            dpg.add_file_extension(".bed", color=(255, 150, 150, 255))
            dpg.add_file_extension(".*")

        dpg.add_file_dialog(directory_selector=True, show=False, callback=self.callback_save_results,
                            tag="select_directory", cancel_callback=self.cancel_callback_directory,
                            width=700, height=400, default_path=self.default_path)

        # Main window
        with dpg.window(label="PlantVarFilter", width=1000, height=600, no_close=True, horizontal_scrollbar=True):
            # Header
            dpg.add_text("PlantVarFilter", color=(230, 230, 230))
            dpg.add_text("Variant filtering • GWAS • Genomic Prediction", color=(160, 160, 160))
            dpg.add_separator()

            with dpg.tab_bar(label='tabbar'):

                # -------- GWAS --------
                with dpg.tab(label='GWAS Analysis'):
                    dpg.add_text("\nStart GWAS Analysis", indent=50)
                    dpg.add_spacer(height=20)
                    geno = dpg.add_button(label="Choose a BED file", callback=lambda: dpg.show_item("file_dialog_bed"),
                                          indent=50, tag='tooltip_bed')
                    dpg.add_spacer(height=5)
                    pheno = dpg.add_button(label="Choose a phenotype file",
                                           callback=lambda: dpg.show_item("file_dialog_pheno"),
                                           indent=50, tag='tooltip_pheno')
                    dpg.add_spacer(height=5)
                    cov_file = dpg.add_button(label="Choose a covariate file (optional)",
                                              callback=lambda: dpg.show_item("file_dialog_cov"),
                                              indent=50, tag='tooltip_cov')
                    dpg.add_spacer(height=20)
                    self.gwas_combo = dpg.add_combo(label="Select algorithm",
                                                    items=["FaST-LMM", "Linear regression", "Ridge Regression",
                                                           "Random Forest (AI)", "XGBoost (AI)"],
                                                    indent=50, width=200, default_value="FaST-LMM",
                                                    tag='tooltip_algorithm')
                    dpg.add_spacer(height=20)
                    gwas_btn = dpg.add_button(label="Run GWAS", callback=self.run_gwas,
                                              user_data=[geno, pheno], indent=50)
                    dpg.bind_item_theme(gwas_btn, self.our_theme)

                # -------- GP --------
                with dpg.tab(label='Genomic Prediction'):
                    dpg.add_text("\nStart Genomic Prediction", indent=50)
                    dpg.add_spacer(height=20)
                    geno = dpg.add_button(label="Choose a BED file", callback=lambda: dpg.show_item("file_dialog_bed"),
                                          indent=50)
                    pheno = dpg.add_button(label="Choose a phenotype file",
                                           callback=lambda: dpg.show_item("file_dialog_pheno"), indent=50)
                    dpg.add_spacer(height=20)
                    self.gwas_gp = dpg.add_combo(label="Select Algorithm",
                                                 items=["XGBoost (AI)", "Random Forest (AI)",
                                                        "Ridge Regression", 'GP_LMM', 'val'],
                                                 indent=50, width=200, default_value="XGBoost (AI)")
                    dpg.add_spacer(height=20)
                    gwas_btn = dpg.add_button(label="Run Genomic Prediction",
                                              callback=self.run_genomic_prediction,
                                              user_data=[geno, pheno], indent=50)
                    dpg.bind_item_theme(gwas_btn, self.our_theme)

                # -------- Convert VCF --------
                with dpg.tab(label='Convert to PLINK'):
                    dpg.add_text("\nConvert a VCF file into PLINK BED and apply MAF/missing genotype filters.",
                                 indent=50)
                    dpg.add_spacer(height=20)
                    dpg.add_text("Select files:", indent=50)
                    vcf = dpg.add_button(label="Choose a VCF file",
                                         callback=lambda: dpg.show_item("file_dialog_vcf"),
                                         indent=50, tag='tooltip_vcf')

                    # IDs file (PLINK --keep)
                    variant_ids = dpg.add_button(label="Choose an IDs file (optional)",
                                                 tag='tooltip_variant',
                                                 callback=lambda: dpg.show_item("file_dialog_variants"),
                                                 indent=50)
                    dpg.add_spacer(height=20)
                    dpg.add_text("Apply filters:", indent=50)
                    maf_input = dpg.add_input_float(label="Minor allele frequency (MAF)",
                                                    width=150, default_value=0.05, step=0.005,
                                                    indent=50, tag='tooltip_maf')
                    geno_input = dpg.add_input_float(label="Missing genotype rate",
                                                     width=150, default_value=0.1, step=0.005,
                                                     indent=50, tag='tooltip_missing')
                    dpg.add_spacer(height=20)
                    convert_btn = dpg.add_button(label="Convert VCF", callback=self.convert_vcf,
                                                 user_data=[maf_input, geno_input, vcf, variant_ids], indent=50)
                    dpg.bind_item_theme(convert_btn, self.our_theme)

                # -------- Settings --------
                with dpg.tab(label='Settings'):
                    dpg.add_spacer(height=10)
                    dpg.add_text("General Settings", indent=50, color=(72, 138, 199))
                    dpg.add_spacer(height=7)
                    self.nr_jobs = dpg.add_input_int(label="Number of jobs to run", width=150, default_value=-1, step=1,
                                                     indent=50, min_value=-1, max_value=50, min_clamped=True,
                                                     max_clamped=True, tag='tooltip_nr_jobs')
                    dpg.add_spacer(height=7)
                    self.gb_goal = dpg.add_input_int(label="Gigabytes of memory per run", width=150, default_value=0,
                                                     step=4, indent=50, min_value=0, max_value=512, min_clamped=True,
                                                     max_clamped=True, tag='tooltip_gb_goal')
                    dpg.add_spacer(height=7)
                    self.snp_limit = dpg.add_input_text(label="SNP limit", indent=50, width=150, default_value='',
                                                        tag="tooltip_limit")
                    dpg.add_spacer(height=7)
                    self.plot_stats = dpg.add_checkbox(label="Advanced Plotting", indent=50, default_value=False,
                                                       tag="tooltip_stats")

                    dpg.add_spacer(height=20)
                    dpg.add_text("Machine Learning Settings", indent=50, color=(72, 138, 199))
                    dpg.add_spacer(height=10)
                    self.train_size_set = dpg.add_input_int(label="Training size", width=150, default_value=70, step=10,
                                                            indent=50, min_value=0, max_value=100, min_clamped=True,
                                                            max_clamped=True, tag='tooltip_training')
                    dpg.add_spacer(height=7)
                    self.model_nr = dpg.add_input_int(label="Nr. of models", width=150, default_value=1, step=1,
                                                      indent=50, min_value=1, max_value=50, min_clamped=True,
                                                      max_clamped=True, tag='tooltip_model')
                    dpg.add_spacer(height=7)
                    self.aggregation_method = dpg.add_combo(("sum", "median", "mean"),
                                                            label="Aggregation Method", indent=50, width=150,
                                                            default_value='sum', tag='tooltip_aggr')
                    dpg.add_spacer(height=7)
                    self.estim_set = dpg.add_input_int(label="Number of trees", width=150, default_value=200, step=10,
                                                       indent=50, min_value=1, min_clamped=True, tag='tooltip_trees')
                    dpg.add_spacer(height=7)
                    self.max_dep_set = dpg.add_input_int(label="Max depth", width=150, default_value=3, step=10,
                                                         indent=50, min_value=0, max_value=100, min_clamped=True,
                                                         max_clamped=True, tag='tooltip_depth')

            # ---------------- Tooltips ----------------
            with dpg.tooltip("tooltip_vcf"):
                dpg.add_text("Select a Variant Call Format file (.vcf or .vcf.gz).", color=[79, 128, 226])

            with dpg.tooltip("tooltip_variant"):
                dpg.add_text(
                    "Optional sample IDs list (PLINK --keep):\n"
                    "two columns, space-separated → FID IID\n"
                    "Use to subset individuals present in the VCF.",
                    color=[79, 128, 226]
                )

            with dpg.tooltip("tooltip_maf"):
                dpg.add_text(
                    "Minor Allele Frequency threshold.\nVariants with MAF below this value are removed.",
                    color=[79, 128, 226]
                )

            with dpg.tooltip("tooltip_missing"):
                dpg.add_text(
                    "Maximum allowed missing genotype rate per variant (PLINK --geno).",
                    color=[79, 128, 226]
                )

            with dpg.tooltip("tooltip_bed"):
                dpg.add_text(
                    "Select a PLINK .bed file (requires matching .bim and .fam in the same folder).",
                    color=[79, 128, 226]
                )

            with dpg.tooltip("tooltip_pheno"):
                dpg.add_text(
                    "Phenotype file (space-separated) with 3 columns:\n"
                    "FID IID Value\n"
                    "IDs must match those in the .fam file.\n"
                    "No header line.",
                    color=[79, 128, 226]
                )

            with dpg.tooltip("tooltip_cov"):
                dpg.add_text(
                    "Optional covariates file (space-separated):\n"
                    "FID IID <cov1> <cov2> ...\n"
                    "IDs must match those in the .fam file.",
                    color=[79, 128, 226]
                )

            with dpg.tooltip("tooltip_algorithm"):
                dpg.add_text("Select the algorithm to use for analysis.", color=[79, 128, 226])

            with dpg.tooltip("tooltip_training"):
                dpg.add_text("Percentage of data used for training (rest is testing).", color=[79, 128, 226])

            with dpg.tooltip("tooltip_trees"):
                dpg.add_text("Number of trees (RF/XGB). Higher can improve accuracy but costs time.", color=[79, 128, 226])

            with dpg.tooltip("tooltip_model"):
                dpg.add_text("Number of models to train (used for model aggregation in GWAS-ML).", color=[79, 128, 226])

            with dpg.tooltip("tooltip_depth"):
                dpg.add_text("Maximum tree depth (complexity).", color=[79, 128, 226])

            with dpg.tooltip("tooltip_nr_jobs"):
                dpg.add_text("Number of CPU cores ( -1 = use all cores ).", color=[79, 128, 226])

            with dpg.tooltip("tooltip_gb_goal"):
                dpg.add_text(
                    "Target GB of RAM per run. If 0, uses block-wise reading (memory efficient).",
                    color=[79, 128, 226]
                )

            with dpg.tooltip("tooltip_limit"):
                dpg.add_text(
                    "Limit number of SNPs plotted in Manhattan/QQ (speed up plotting on huge datasets).\n"
                    "Empty = use all.",
                    color=[79, 128, 226]
                )

            with dpg.tooltip("tooltip_stats"):
                dpg.add_text(
                    "Enable advanced PDF plots for phenotype/genotype statistics.\n"
                    "Useful for publication; increases runtime on large data.",
                    color=[79, 128, 226]
                )

            with dpg.tooltip("tooltip_aggr"):
                dpg.add_text(
                    "Aggregate SNP effects from multiple models (GWAS-ML):\n"
                    "sum / median / mean.",
                    color=[79, 128, 226]
                )

            # Bind font & theme
            if self.font:
                dpg.bind_font(self.font)
            dpg.set_global_font_scale(0.6)
            dpg.bind_theme(self.our_theme)

    # ---------------- Callbacks ----------------
    def callback_vcf(self, sender, app_data):
        self.vcf_app_data = app_data
        vcf_path, current_path = self.get_selection_path(self.vcf_app_data)
        dpg.configure_item("file_dialog_variants", default_path=current_path)
        self.add_log('VCF Selected: ' + vcf_path)

    def callback_bed(self, sender, app_data):
        self.bed_app_data = app_data
        try:
            bed_path, current_path = self.get_selection_path(self.bed_app_data)
            dpg.configure_item("file_dialog_cov", default_path=current_path)
            dpg.configure_item("file_dialog_pheno", default_path=current_path)
            self.add_log('BED file Selected: ' + bed_path)
        except TypeError:
            self.add_log('Invalid BED file Selected', error=True)

    def callback_variants(self, sender, app_data):
        """IDs list file for PLINK --keep (اختياري)."""
        self.variants_app_data = app_data
        variants_path, current_path = self.get_selection_path(self.variants_app_data)
        dpg.configure_item("file_dialog_vcf", default_path=current_path)
        self.add_log('IDs file Selected: ' + variants_path)

    def callback_pheno(self, sender, app_data):
        self.pheno_app_data = app_data
        try:
            pheno_path, current_path = self.get_selection_path(self.pheno_app_data)
            dpg.configure_item("file_dialog_cov", default_path=current_path)
            dpg.configure_item("file_dialog_bed", default_path=current_path)
            self.add_log('Pheno File Selected: ' + pheno_path)
        except TypeError:
            self.add_log('Wrong Pheno File Selected', error=True)

    def callback_cov(self, sender, app_data):
        self.cov_app_data = app_data
        try:
            cov_path, current_path = self.get_selection_path(self.cov_app_data)
            dpg.configure_item("file_dialog_bed", default_path=current_path)
            dpg.configure_item("file_dialog_pheno", default_path=current_path)
            self.add_log('Covariates File Selected: ' + cov_path)
        except TypeError:
            self.add_log('Wrong Covariates File Selected', error=True)

    def get_selection_path(self, app_data):
        current_path = app_data['current_path'] + '/'
        k = app_data['selections']
        try:
            for _, value in k.items():
                file_path = value
            return file_path, current_path
        except UnboundLocalError:
            pass

    def callback_save_results(self, sender, app_data):
        self.results_directory = app_data
        results_path, current_path = self.get_selection_path(self.results_directory)
        save_dir = self.helper.save_results(
            os.getcwd(), current_path,
            self.gwas_result_name, self.gwas_result_name_top,
            self.manhatten_plot_name, self.qq_plot_name, self.algorithm,
            self.genomic_predict_name, self.gp_plot_name, self.gp_plot_name_scatter,
            self.add_log, self.settings_lst, self.pheno_stats_name, self.geno_stats_name
        )
        self.add_log('Results saved in: ' + save_dir)

    def cancel_callback_directory(self, sender, app_data):
        self.add_log('Process Canceled')

    # ---------------- Helpers ----------------
    def delete_files(self):
        for tag in ["manhatten_image", "manhatten_tag", "qq_image", "qq_tag",
                    "table_gwas", "table_gp", "ba_tag", "ba_tag2", "ba_image", "ba_image2"]:
            if dpg.does_item_exist(tag):
                dpg.delete_item(tag)

        file_names = [
            self.gwas_result_name, self.gwas_result_name_top, self.genomic_predict_name, self.gp_plot_name,
            self.manhatten_plot_name, self.qq_plot_name, self.gp_plot_name_scatter,
            self.manhatten_plot_name.replace('manhatten_plot', 'manhatten_plot_high'),
            self.qq_plot_name.replace('qq_plot', 'qq_plot_high'),
            self.gp_plot_name_scatter.replace('GP_scatter_plot', 'GP_scatter_plot_high'),
            self.gp_plot_name.replace('Bland_Altman_plot', 'Bland_Altman_plot_high'),
            self.genomic_predict_name.replace('.csv', '_valdation.csv'),
            self.pheno_stats_name, self.geno_stats_name
        ]
        for f in file_names:
            if os.path.exists(f):
                os.remove(f)

    def add_log(self, message, warn=False, error=False):
        if warn:
            self.logz.log_warning(message)
        elif error:
            self.logz.log_error(message)
        else:
            self.logz.log_info(message)

    # ---------------- Actions ----------------
    def convert_vcf(self, sender, data, user_data):
        maf = str(dpg.get_value(user_data[0]))
        geno = str(dpg.get_value(user_data[1]))
        vcf_path, _ = self.get_selection_path(self.vcf_app_data)

        if self.variants_app_data is None:
            variants_path = None  # IDs file not provided
        else:
            variants_path, _ = self.get_selection_path(self.variants_app_data)

        self.add_log('Start converting VCF to BED...')
        out_prefix = f"{vcf_path.split('.')[0]}_maf{round(float(maf), 2)}_geno{round(float(geno), 2)}"
        plink_log = self.gwas.vcf_to_bed(vcf_path, variants_path, out_prefix, maf, geno)
        self.add_log(plink_log)

    def run_gwas(self, sender, data, user_data):
        self.delete_files()

        # Read settings
        train_size_set = (100 - dpg.get_value(self.train_size_set)) / 100
        estimators = dpg.get_value(self.estim_set)
        model_nr = dpg.get_value(self.model_nr)
        snp_limit = dpg.get_value(self.snp_limit)
        nr_jobs = int(dpg.get_value(self.nr_jobs)) or -1
        gb_goal = int(dpg.get_value(self.gb_goal))
        max_dep_set = dpg.get_value(self.max_dep_set)
        self.algorithm = dpg.get_value(self.gwas_combo)
        aggregation_method = str(dpg.get_value(self.aggregation_method))

        try:
            self.add_log('Reading files...')
            bed_path, _ = self.get_selection_path(self.bed_app_data)
            pheno_path, _ = self.get_selection_path(self.pheno_app_data)
            try:
                cov_path, _ = self.get_selection_path(self.cov_app_data)
            except Exception:
                cov_path = None

            self.add_log('Validating files...')
            check_input_data = self.gwas.validate_gwas_input_files(bed_path, pheno_path)

            chrom_mapping = self.helper.replace_with_integers(bed_path.replace('.bed', '.bim'))
            self.settings_lst = [self.algorithm, bed_path, pheno_path, train_size_set, estimators, model_nr, max_dep_set]

            if check_input_data[0]:
                bed = Bed(str(bed_path), count_A1=False, chrom_map=chrom_mapping)
                pheno = Pheno(str(pheno_path))
                cov = Pheno(str(cov_path)) if cov_path else None

                bed, pheno = pstutil.intersect_apply([bed, pheno])
                bed_fixed = self.gwas.filter_out_missing(bed)

                s3 = f"Dataset after intersection: SNPs: {bed.sid_count}  Pheno IDs: {pheno.iid_count}"
                self.add_log(s3, warn=True)
                self.add_log('Starting Analysis, this might take a while...')

                if self.algorithm in ('FaST-LMM', 'Linear regression'):
                    gwas_df, df_plot = self.gwas.run_gwas_lmm(
                        bed_fixed, pheno, chrom_mapping, self.add_log,
                        self.gwas_result_name, self.algorithm, bed_path, cov, gb_goal
                    )
                elif self.algorithm == 'Random Forest (AI)':
                    gwas_df, df_plot = self.gwas.run_gwas_rf(
                        bed_fixed, pheno, bed_path, train_size_set, estimators,
                        self.gwas_result_name, chrom_mapping, self.add_log,
                        model_nr, nr_jobs, aggregation_method
                    )
                elif self.algorithm == 'XGBoost (AI)':
                    gwas_df, df_plot = self.gwas.run_gwas_xg(
                        bed_fixed, pheno, bed_path, train_size_set, estimators,
                        self.gwas_result_name, chrom_mapping, self.add_log,
                        model_nr, max_dep_set, nr_jobs, aggregation_method
                    )
                elif self.algorithm == 'Ridge Regression':
                    gwas_df, df_plot = self.gwas.run_gwas_ridge(
                        bed_fixed, pheno, bed_path, train_size_set, 1.0,
                        self.gwas_result_name, chrom_mapping, self.add_log,
                        model_nr, aggregation_method
                    )
            else:
                self.add_log(check_input_data[1], error=True)
                return

            if gwas_df is not None:
                self.add_log('GWAS Analysis done.')
                self.add_log('GWAS Results Plotting...')
                if dpg.get_value(self.plot_stats):
                    self.plot_class.plot_pheno_statistics(pheno_path, self.pheno_stats_name)
                    self.plot_class.plot_geno_statistics(bed_fixed, pheno, self.geno_stats_name)

                self.gwas.plot_gwas(df_plot, snp_limit, self.algorithm,
                                    self.manhatten_plot_name, self.qq_plot_name, chrom_mapping)
                self.add_log('Done...')
                self.show_results_window(gwas_df, self.algorithm, genomic_predict=False)

                # clear selections
                self.bed_app_data = None
                self.pheno_app_data = None
                self.cov_app_data = None
            else:
                self.add_log('Error, GWAS Analysis could not be started.', error=True)

        except TypeError:
            self.add_log('Please select a phenotype and genotype file. ', error=True)

    def run_genomic_prediction(self, sender, data, user_data):
        self.delete_files()
        self.add_log('Reading Bed file...')

        self.algorithm = dpg.get_value(self.gwas_gp)
        test_size = (100 - dpg.get_value(self.train_size_set)) / 100
        estimators = dpg.get_value(self.estim_set)
        max_dep_set = dpg.get_value(self.max_dep_set)
        model_nr = dpg.get_value(self.model_nr)
        nr_jobs = int(dpg.get_value(self.nr_jobs)) or -1

        try:
            self.add_log('Reading files...')
            bed_path, _ = self.get_selection_path(self.bed_app_data)
            pheno_path, _ = self.get_selection_path(self.pheno_app_data)

            self.add_log('Validating files...')
            check_input_data = self.gwas.validate_gwas_input_files(bed_path, pheno_path)
            chrom_mapping = self.helper.replace_with_integers(bed_path.replace('.bed', '.bim'))
            self.settings_lst = [self.algorithm, bed_path, pheno_path, test_size, estimators, model_nr, max_dep_set]

            if check_input_data[0]:
                bed = Bed(str(bed_path), count_A1=False, chrom_map=chrom_mapping)
                pheno = Pheno(str(pheno_path))
                bed, pheno = pstutil.intersect_apply([bed, pheno])
                bed_fixed = self.gwas.filter_out_missing(bed)

                s3 = f"Dataset after intersection: SNPs: {bed.sid_count}  Pheno IDs: {pheno.iid_count}"
                self.add_log(s3, warn=True)
                self.add_log('Starting Analysis, this might take a while...')

                if self.algorithm == 'GP_LMM':
                    gp_df = self.genomic_predict_class.run_lmm_gp(
                        bed_fixed, pheno, self.genomic_predict_name, model_nr, self.add_log,
                        bed_path, chrom_mapping
                    )
                elif self.algorithm == 'Random Forest (AI)':
                    gp_df = self.genomic_predict_class.run_gp_rf(
                        bed_fixed, pheno, bed_path, test_size, estimators,
                        self.genomic_predict_name, chrom_mapping, self.add_log,
                        model_nr, nr_jobs
                    )
                elif self.algorithm == 'XGBoost (AI)':
                    gp_df = self.genomic_predict_class.run_gp_xg(
                        bed_fixed, pheno, bed_path, test_size, estimators,
                        self.genomic_predict_name, chrom_mapping, self.add_log,
                        model_nr, max_dep_set, nr_jobs
                    )
                elif self.algorithm == 'Ridge Regression':
                    gp_df = self.genomic_predict_class.run_gp_ridge(
                        bed_fixed, pheno, bed_path, test_size, 1.0,
                        self.genomic_predict_name, chrom_mapping, self.add_log,
                        model_nr
                    )
                else:
                    self.genomic_predict_class.model_validation(
                        bed_fixed, pheno, bed_path, test_size, estimators,
                        self.genomic_predict_name, chrom_mapping, self.add_log,
                        model_nr, max_dep_set, validation_size=0.1
                    )
            else:
                self.add_log(check_input_data[1], error=True)
                return

            if gp_df is not None:
                self.add_log('Genomic Prediction done.')
                self.add_log('Genomic Prediction Plotting...')
                self.genomic_predict_class.plot_gp(gp_df, self.gp_plot_name, self.algorithm)
                self.genomic_predict_class.plot_gp_scatter(gp_df, self.gp_plot_name_scatter, self.algorithm)
                self.add_log('Done...')
                self.show_results_window(gp_df, self.algorithm, genomic_predict=True)

                self.bed_app_data = None
                self.pheno_app_data = None
            else:
                self.add_log('Error, GWAS Analysis could not be started.', error=True)

        except TypeError:
            self.add_log('Please select a phenotype and genotype file. ', error=True)

    # ---------------- Results Window ----------------
    def show_results_window(self, df, algorithm, genomic_predict):
        with dpg.window(label="Results", width=975, height=600, horizontal_scrollbar=True, pos=(1000, 35)):
            dpg.add_button(label="Export Results", pos=(400, 40), callback=lambda: dpg.show_item("select_directory"))
            dpg.add_spacer(height=60)

            if genomic_predict:
                width, height, channels, data = dpg.load_image(self.gp_plot_name)
                with dpg.texture_registry(show=False):
                    dpg.add_static_texture(width=width, height=height, default_value=data, tag="ba_tag")

                width, height, channels, data = dpg.load_image(self.gp_plot_name_scatter)
                with dpg.texture_registry(show=False):
                    dpg.add_static_texture(width=width, height=height, default_value=data, tag="ba_tag2")

                with dpg.tab_bar(label='tabbar'):
                    with dpg.tab(label="Genomic Prediction Results"):
                        df = df[['ID1', 'BED_ID2_x', 'Mean_Predicted_Value', 'Pheno_Value', 'Difference']]
                        df.columns = ['FID', 'IID', 'Predicted_Value', 'Pheno_Value', 'Difference']
                        with dpg.table(label='DatasetTable2', row_background=True, borders_innerH=True,
                                       borders_outerH=True, borders_innerV=True, borders_outerV=True, tag='table_gp',
                                       sortable=True, resizable=True):
                            for i in range(df.shape[1]):
                                dpg.add_table_column(label=df.columns[i], parent='table_gp')
                            for i in range(len(df)):
                                with dpg.table_row():
                                    for j in range(df.shape[1]):
                                        value = df.iloc[i, j]
                                        formatted_value = f"{value:.2f}" if isinstance(value, float) else str(value)
                                        dpg.add_text(formatted_value)
                    with dpg.tab(label="Correlation Plot (Predicted vs. Phenotype)"):
                        dpg.add_image(texture_tag="ba_tag2", tag="ba_image2", width=750, height=450)
                    with dpg.tab(label="Bland-Altman Plot (Model Accuracy)"):
                        dpg.add_image(texture_tag="ba_tag", tag="ba_image", width=750, height=450)

            else:
                width, height, channels, data = dpg.load_image(self.manhatten_plot_name)
                with dpg.texture_registry(show=False):
                    dpg.add_static_texture(width=width, height=height, default_value=data, tag="manhatten_tag")
                with dpg.tab_bar(label='tabbar'):
                    with dpg.tab(label="Manhattan Plot"):
                        if algorithm in ("FaST-LMM", "Linear regression"):
                            dpg.add_image(texture_tag="manhatten_tag", tag="manhatten_image", width=950, height=400)
                        else:
                            dpg.add_image(texture_tag="manhatten_tag", tag="manhatten_image", width=900, height=300)
                    if algorithm in ("FaST-LMM", "Linear regression"):
                        width2, height2, channels2, data2 = dpg.load_image(self.qq_plot_name)
                        with dpg.texture_registry(show=False):
                            dpg.add_static_texture(width=width2, height=height2, default_value=data2, tag="qq_tag")
                        df = df.sort_values(by=['PValue'], ascending=True)
                        with dpg.tab(label="QQ-Plot"):
                            dpg.add_image(texture_tag="qq_tag", tag="qq_image", height=450, width=450)
                    else:
                        df = df.sort_values(by=['SNP effect'], ascending=False)

                    with dpg.tab(label="GWAS Results (Top 500)"):
                        if algorithm in ("FaST-LMM", "Linear regression"):
                            df = df[['SNP', 'Chr', 'ChrPos', 'PValue']]
                        else:
                            df.columns = df.columns.str.replace('SNP effect_sd', 'SNP effect SD')

                        max_rows = min(len(df), 500)
                        with dpg.table(label='DatasetTable', row_background=True, borders_innerH=True,
                                       borders_outerH=True, borders_innerV=True, borders_outerV=True, tag='table_gwas',
                                       sortable=True):
                            for i in range(df.shape[1]):
                                dpg.add_table_column(label=df.columns[i], parent='table_gwas')
                            for i in range(max_rows):
                                with dpg.table_row():
                                    for j in range(df.shape[1]):
                                        value = df.iloc[i, j]
                                        dpg.add_text(value)

    # ---------------- Runner ----------------
    def run(self):
        dpg.setup_dearpygui()
        dpg.show_viewport()
        dpg.start_dearpygui()


if __name__ == "__main__":
    main()
