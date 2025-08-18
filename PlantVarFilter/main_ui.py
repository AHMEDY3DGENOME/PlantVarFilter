# main_gui.py — PlantVarFilter (plant-themed UI + docking + improved dialogs + header icon)

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
    dpg.create_context()
    dpg.configure_app(docking=True, docking_space=True)
    app = GWASApp()
    app.run()
    dpg.destroy_context()


class GWASApp:
    def __init__(self):
        self.gwas = GWAS()
        self.helper = HELPERS()
        self.genomic_predict_class = GenomicPrediction()
        self.plot_class = Plot()

        # ---------- Fonts ----------
        with dpg.font_registry():
            script_dir = os.path.dirname(__file__)
            font_path = os.path.join(script_dir, "test.ttf")
            self.font = dpg.add_font(font_path, 40, tag="ttf-font") if os.path.exists(font_path) else None

        # ---------- Base window themes (Dark/Light plant palette) ----------
        with dpg.theme() as self.theme_plant_dark:
            with dpg.theme_component(dpg.mvAll):
                dpg.add_theme_color(dpg.mvThemeCol_WindowBg, (18, 24, 20))
                dpg.add_theme_color(dpg.mvThemeCol_ChildBg, (24, 32, 26))
                dpg.add_theme_color(dpg.mvThemeCol_PopupBg, (26, 34, 28))
                dpg.add_theme_color(dpg.mvThemeCol_Text, (220, 235, 220))
                dpg.add_theme_color(dpg.mvThemeCol_TextDisabled, (150, 165, 150))
                dpg.add_theme_color(dpg.mvThemeCol_Border, (70, 80, 75))
                dpg.add_theme_color(dpg.mvThemeCol_Separator, (70, 80, 75))
                dpg.add_theme_color(dpg.mvThemeCol_Tab, (28, 36, 30))
                dpg.add_theme_color(dpg.mvThemeCol_TabHovered, (38, 50, 42))
                dpg.add_theme_color(dpg.mvThemeCol_TabActive, (36, 48, 40))
                dpg.add_theme_color(dpg.mvThemeCol_TitleBg, (28, 40, 34))
                dpg.add_theme_color(dpg.mvThemeCol_TitleBgActive, (46, 84, 66))
                dpg.add_theme_style(dpg.mvStyleVar_WindowPadding, 14, 14)
                dpg.add_theme_style(dpg.mvStyleVar_FramePadding, 10, 8)
                dpg.add_theme_style(dpg.mvStyleVar_ItemSpacing, 10, 8)
                dpg.add_theme_style(dpg.mvStyleVar_WindowRounding, 10)
                dpg.add_theme_style(dpg.mvStyleVar_FrameRounding, 8)
                dpg.add_theme_style(dpg.mvStyleVar_ChildRounding, 10)
                dpg.add_theme_style(dpg.mvStyleVar_GrabRounding, 8)
                dpg.add_theme_style(dpg.mvStyleVar_TabRounding, 8)

        with dpg.theme() as self.theme_plant_light:
            with dpg.theme_component(dpg.mvAll):
                dpg.add_theme_color(dpg.mvThemeCol_WindowBg, (240, 245, 238))
                dpg.add_theme_color(dpg.mvThemeCol_ChildBg, (235, 242, 232))
                dpg.add_theme_color(dpg.mvThemeCol_PopupBg, (230, 238, 226))
                dpg.add_theme_color(dpg.mvThemeCol_Text, (30, 45, 35))
                dpg.add_theme_color(dpg.mvThemeCol_TextDisabled, (110, 125, 115))
                dpg.add_theme_color(dpg.mvThemeCol_Border, (180, 190, 185))
                dpg.add_theme_color(dpg.mvThemeCol_Separator, (180, 190, 185))
                dpg.add_theme_color(dpg.mvThemeCol_Tab, (222, 235, 224))
                dpg.add_theme_color(dpg.mvThemeCol_TabHovered, (208, 228, 214))
                dpg.add_theme_color(dpg.mvThemeCol_TabActive, (206, 226, 210))
                dpg.add_theme_color(dpg.mvThemeCol_TitleBg, (206, 226, 210))
                dpg.add_theme_color(dpg.mvThemeCol_TitleBgActive, (186, 214, 194))
                dpg.add_theme_style(dpg.mvStyleVar_WindowPadding, 14, 14)
                dpg.add_theme_style(dpg.mvStyleVar_FramePadding, 10, 8)
                dpg.add_theme_style(dpg.mvStyleVar_ItemSpacing, 10, 8)
                dpg.add_theme_style(dpg.mvStyleVar_WindowRounding, 10)
                dpg.add_theme_style(dpg.mvStyleVar_FrameRounding, 8)
                dpg.add_theme_style(dpg.mvStyleVar_ChildRounding, 10)
                dpg.add_theme_style(dpg.mvStyleVar_GrabRounding, 8)
                dpg.add_theme_style(dpg.mvStyleVar_TabRounding, 8)

        # ---------- Component themes ----------
        with dpg.theme() as self.theme_primary_btn_dark:
            with dpg.theme_component(dpg.mvButton):
                dpg.add_theme_color(dpg.mvThemeCol_Button, (60, 179, 113))
                dpg.add_theme_color(dpg.mvThemeCol_ButtonHovered, (80, 200, 140))
                dpg.add_theme_color(dpg.mvThemeCol_ButtonActive, (50, 160, 100))
                dpg.add_theme_style(dpg.mvStyleVar_FrameRounding, 10)
                dpg.add_theme_style(dpg.mvStyleVar_FramePadding, 12, 8)

        with dpg.theme() as self.theme_primary_btn_light:
            with dpg.theme_component(dpg.mvButton):
                dpg.add_theme_color(dpg.mvThemeCol_Button, (90, 190, 120))
                dpg.add_theme_color(dpg.mvThemeCol_ButtonHovered, (110, 205, 140))
                dpg.add_theme_color(dpg.mvThemeCol_ButtonActive, (75, 175, 110))
                dpg.add_theme_style(dpg.mvStyleVar_FrameRounding, 10)
                dpg.add_theme_style(dpg.mvStyleVar_FramePadding, 12, 8)

        with dpg.theme() as self.theme_secondary_btn_dark:
            with dpg.theme_component(dpg.mvButton):
                dpg.add_theme_color(dpg.mvThemeCol_Button, (45, 55, 50))
                dpg.add_theme_color(dpg.mvThemeCol_ButtonHovered, (60, 75, 65))
                dpg.add_theme_color(dpg.mvThemeCol_ButtonActive, (40, 50, 45))
                dpg.add_theme_style(dpg.mvStyleVar_FrameRounding, 8)
                dpg.add_theme_style(dpg.mvStyleVar_FramePadding, 10, 6)

        with dpg.theme() as self.theme_secondary_btn_light:
            with dpg.theme_component(dpg.mvButton):
                dpg.add_theme_color(dpg.mvThemeCol_Button, (210, 232, 216))
                dpg.add_theme_color(dpg.mvThemeCol_ButtonHovered, (198, 224, 206))
                dpg.add_theme_color(dpg.mvThemeCol_ButtonActive, (188, 214, 198))
                dpg.add_theme_style(dpg.mvStyleVar_FrameRounding, 8)
                dpg.add_theme_style(dpg.mvStyleVar_FramePadding, 10, 6)

        with dpg.theme() as self.theme_inputs_dark:
            with dpg.theme_component(dpg.mvInputFloat):
                dpg.add_theme_color(dpg.mvThemeCol_FrameBg, (30, 40, 32))
                dpg.add_theme_color(dpg.mvThemeCol_Text, (220, 235, 220))
                dpg.add_theme_style(dpg.mvStyleVar_FrameRounding, 8)
            with dpg.theme_component(dpg.mvInputInt):
                dpg.add_theme_color(dpg.mvThemeCol_FrameBg, (30, 40, 32))
                dpg.add_theme_style(dpg.mvStyleVar_FrameRounding, 8)
            with dpg.theme_component(dpg.mvInputText):
                dpg.add_theme_color(dpg.mvThemeCol_FrameBg, (30, 40, 32))
                dpg.add_theme_style(dpg.mvStyleVar_FrameRounding, 8)
            with dpg.theme_component(dpg.mvCombo):
                dpg.add_theme_color(dpg.mvThemeCol_FrameBg, (30, 40, 32))
                dpg.add_theme_style(dpg.mvStyleVar_FrameRounding, 8)

        with dpg.theme() as self.theme_inputs_light:
            light_bg = (215, 228, 215)
            with dpg.theme_component(dpg.mvInputFloat):
                dpg.add_theme_color(dpg.mvThemeCol_FrameBg, light_bg)
                dpg.add_theme_color(dpg.mvThemeCol_Text, (30, 45, 35))
                dpg.add_theme_style(dpg.mvStyleVar_FrameRounding, 8)
            with dpg.theme_component(dpg.mvInputInt):
                dpg.add_theme_color(dpg.mvThemeCol_FrameBg, light_bg)
                dpg.add_theme_style(dpg.mvStyleVar_FrameRounding, 8)
            with dpg.theme_component(dpg.mvInputText):
                dpg.add_theme_color(dpg.mvThemeCol_FrameBg, light_bg)
                dpg.add_theme_style(dpg.mvStyleVar_FrameRounding, 8)
            with dpg.theme_component(dpg.mvCombo):
                dpg.add_theme_color(dpg.mvThemeCol_FrameBg, light_bg)
                dpg.add_theme_style(dpg.mvStyleVar_FrameRounding, 8)

        # ---------- File Dialog themes ----------
        with dpg.theme() as self.theme_fd_dark:
            with dpg.theme_component(dpg.mvFileDialog):
                dpg.add_theme_color(dpg.mvThemeCol_WindowBg, (22, 28, 24))
                dpg.add_theme_color(dpg.mvThemeCol_TitleBg, (28, 40, 34))
                dpg.add_theme_color(dpg.mvThemeCol_TitleBgActive, (46, 84, 66))
                dpg.add_theme_color(dpg.mvThemeCol_Header, (40, 56, 48))
                dpg.add_theme_color(dpg.mvThemeCol_HeaderHovered, (54, 76, 62))
                dpg.add_theme_color(dpg.mvThemeCol_HeaderActive, (60, 96, 76))
                dpg.add_theme_color(dpg.mvThemeCol_TableHeaderBg, (34, 48, 40))
                dpg.add_theme_color(dpg.mvThemeCol_FrameBg, (30, 40, 32))
                dpg.add_theme_color(dpg.mvThemeCol_Button, (60, 179, 113))
                dpg.add_theme_color(dpg.mvThemeCol_ButtonHovered, (80, 200, 140))
                dpg.add_theme_color(dpg.mvThemeCol_ButtonActive, (50, 160, 100))
                dpg.add_theme_color(dpg.mvThemeCol_Text, (220, 235, 220))
                dpg.add_theme_style(dpg.mvStyleVar_WindowRounding, 10)
                dpg.add_theme_style(dpg.mvStyleVar_FrameRounding, 8)

        with dpg.theme() as self.theme_fd_light:
            with dpg.theme_component(dpg.mvFileDialog):
                dpg.add_theme_color(dpg.mvThemeCol_WindowBg, (244, 248, 242))
                dpg.add_theme_color(dpg.mvThemeCol_TitleBg, (206, 226, 210))
                dpg.add_theme_color(dpg.mvThemeCol_TitleBgActive, (186, 214, 194))
                dpg.add_theme_color(dpg.mvThemeCol_Header, (208, 228, 214))
                dpg.add_theme_color(dpg.mvThemeCol_HeaderHovered, (198, 222, 206))
                dpg.add_theme_color(dpg.mvThemeCol_HeaderActive, (188, 214, 198))
                dpg.add_theme_color(dpg.mvThemeCol_TableHeaderBg, (210, 230, 214))
                dpg.add_theme_color(dpg.mvThemeCol_FrameBg, (215, 228, 215))
                dpg.add_theme_color(dpg.mvThemeCol_Button, (90, 190, 120))
                dpg.add_theme_color(dpg.mvThemeCol_ButtonHovered, (110, 205, 140))
                dpg.add_theme_color(dpg.mvThemeCol_ButtonActive, (75, 175, 110))
                dpg.add_theme_color(dpg.mvThemeCol_Text, (30, 45, 35))
                dpg.add_theme_style(dpg.mvStyleVar_WindowRounding, 10)
                dpg.add_theme_style(dpg.mvStyleVar_FrameRounding, 8)

        # ---------- Track items to re-bind themes on toggle ----------
        self._primary_buttons = []
        self._secondary_buttons = []
        self._inputs = []
        self._file_dialogs = []
        self.night_mode = True

        # Header tags to update on theme toggle
        self.header_group = None
        self.brand_title = None
        self.brand_tagline = None
        self.logo_draw = None

        # ---------- State ----------
        self.vcf_app_data = None
        self.variants_app_data = None
        self.results_directory = None
        self.bed_app_data = None
        self.pheno_app_data = None
        self.cov_app_data = None
        self.default_path = os.path.expanduser("~")

        # Filenames
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

        self.setup_gui()

    # ---------- UI ----------
    def setup_gui(self):
        dpg.create_viewport(title='PlantVarFilter', width=2000, height=1200, resizable=True)

        # Root dockspace
        with dpg.window(tag="PrimaryWindow", no_close=True, no_move=True, no_resize=True,
                        no_title_bar=True, no_background=True):
            pass
        dpg.set_primary_window("PrimaryWindow", True)

        # File dialogs (modal + larger)
        with dpg.file_dialog(directory_selector=False, show=False, callback=self.callback_vcf,
                             file_count=3, tag="file_dialog_vcf", width=920, height=560,
                             default_path=self.default_path, modal=True):
            dpg.add_file_extension("Source files (*.vcf *.gz){.vcf,.gz}", color=(255, 255, 0, 255))
        self._file_dialogs.append("file_dialog_vcf")

        with dpg.file_dialog(directory_selector=False, show=False, callback=self.callback_variants,
                             file_count=3, tag="file_dialog_variants", width=920, height=560,
                             default_path=self.default_path, modal=True):
            dpg.add_file_extension("Text files (*.txt *.csv){.txt,.csv}", color=(255, 255, 0, 255))
            dpg.add_file_extension(".*")
        self._file_dialogs.append("file_dialog_variants")

        with dpg.file_dialog(directory_selector=False, show=False, callback=self.callback_pheno,
                             file_count=3, tag="file_dialog_pheno", width=920, height=560,
                             default_path=self.default_path, modal=True):
            dpg.add_file_extension("Text files (*.txt *.csv){.txt,.csv}", color=(255, 255, 0, 255))
            dpg.add_file_extension(".*")
        self._file_dialogs.append("file_dialog_pheno")

        with dpg.file_dialog(directory_selector=False, show=False, callback=self.callback_cov,
                             file_count=3, tag="file_dialog_cov", width=920, height=560,
                             default_path=self.default_path, modal=True):
            dpg.add_file_extension("Text files (*.txt *.csv){.txt,.csv}", color=(255, 255, 0, 255))
            dpg.add_file_extension(".*")
        self._file_dialogs.append("file_dialog_cov")

        with dpg.file_dialog(directory_selector=False, show=False, callback=self.callback_bed,
                             file_count=3, tag="file_dialog_bed", width=920, height=560,
                             default_path=self.default_path, modal=True):
            dpg.add_file_extension(".bed", color=(255, 150, 150, 255))
            dpg.add_file_extension(".*")
        self._file_dialogs.append("file_dialog_bed")

        dpg.add_file_dialog(directory_selector=True, show=False, callback=self.callback_save_results,
                            tag="select_directory", cancel_callback=self.cancel_callback_directory,
                            width=900, height=540, default_path=self.default_path, modal=True)
        self._file_dialogs.append("select_directory")

        # Main window (dockable)
        with dpg.window(label="PlantVarFilter", tag="MainWindow", no_close=True, width=1100, height=700):
            # Header (icon + title + tagline)
            self.build_header(parent="MainWindow")

            dpg.add_spacer(height=4)
            dpg.add_separator()
            dpg.add_spacer(height=6)

            dpg.add_tab_bar(tag="main_tabbar")

            # -------- GWAS --------
            with dpg.tab(label='GWAS Analysis', parent="main_tabbar"):
                dpg.add_text("\nStart GWAS Analysis", indent=10)
                dpg.add_spacer(height=10)
                with dpg.group(horizontal=True, horizontal_spacing=60):
                    with dpg.group():
                        geno = dpg.add_button(label="Choose a BED file",
                                              callback=lambda: dpg.show_item("file_dialog_bed"),
                                              width=220, tag='tooltip_bed')
                        self._secondary_buttons.append(geno)

                        dpg.add_spacer(height=6)
                        pheno = dpg.add_button(label="Choose a phenotype file",
                                               callback=lambda: dpg.show_item("file_dialog_pheno"),
                                               width=220, tag='tooltip_pheno')
                        self._secondary_buttons.append(pheno)

                        dpg.add_spacer(height=6)
                        cov_file = dpg.add_button(label="Choose covariate file (option)",
                                                  callback=lambda: dpg.show_item("file_dialog_cov"),
                                                  width=220, tag='tooltip_cov')
                        self._secondary_buttons.append(cov_file)

                    with dpg.group():
                        self.gwas_combo = dpg.add_combo(
                            label="Analysis Algorithms",
                            items=["FaST-LMM", "Linear regression", "Ridge Regression",
                                   "Random Forest (AI)", "XGBoost (AI)"],
                            width=240, default_value="FaST-LMM", tag='tooltip_algorithm'
                        )
                        self._inputs.append(self.gwas_combo)

                        dpg.add_spacer(height=14)
                        gwas_btn = dpg.add_button(
                            label="Run GWAS", callback=self.run_gwas,
                            user_data=[geno, pheno], width=180, height=36
                        )
                        self._primary_buttons.append(gwas_btn)

            # -------- Genomic Prediction --------
            with dpg.tab(label='Genomic Prediction', parent="main_tabbar"):
                dpg.add_text("\nStart Genomic Prediction", indent=10)
                dpg.add_spacer(height=10)
                with dpg.group(horizontal=True, horizontal_spacing=60):
                    with dpg.group():
                        geno = dpg.add_button(label="Choose a BED file",
                                              callback=lambda: dpg.show_item("file_dialog_bed"),
                                              width=220, tag='tooltip_bed_gp')
                        self._secondary_buttons.append(geno)

                        dpg.add_spacer(height=6)
                        pheno = dpg.add_button(label="Choose a phenotype file",
                                               callback=lambda: dpg.show_item("file_dialog_pheno"),
                                               width=220, tag='tooltip_pheno_gp')
                        self._secondary_buttons.append(pheno)

                    with dpg.group():
                        self.gwas_gp = dpg.add_combo(
                            label="Analysis Algorithms",
                            items=["XGBoost (AI)", "Random Forest (AI)",
                                   "Ridge Regression", "GP_LMM", "val"],
                            width=240, default_value="XGBoost (AI)", tag='tooltip_algorithm_gp'
                        )
                        self._inputs.append(self.gwas_gp)

                        dpg.add_spacer(height=14)
                        gp_btn = dpg.add_button(
                            label="Run Genomic Prediction", callback=self.run_genomic_prediction,
                            user_data=[geno, pheno], width=220, height=36
                        )
                        self._primary_buttons.append(gp_btn)

            # -------- Convert to PLINK --------
            with dpg.tab(label='Convert to PLINK', parent="main_tabbar"):
                dpg.add_text("\nConvert a VCF file into PLINK BED and apply MAF/missing genotype filters.",
                             indent=10)
                dpg.add_spacer(height=10)
                with dpg.group(horizontal=True, horizontal_spacing=60):
                    with dpg.group():
                        dpg.add_text("Select files:", indent=0)
                        vcf = dpg.add_button(
                            label="Choose a VCF file",
                            callback=lambda: dpg.show_item("file_dialog_vcf"),
                            width=220, tag='tooltip_vcf'
                        )
                        self._secondary_buttons.append(vcf)

                        dpg.add_spacer(height=6)
                        variant_ids = dpg.add_button(
                            label="Choose IDs file (option)",
                            callback=lambda: dpg.show_item("file_dialog_variants"),
                            width=220, tag='tooltip_variant'
                        )
                        self._secondary_buttons.append(variant_ids)

                    with dpg.group():
                        dpg.add_text("Apply filters:", indent=0)
                        maf_input = dpg.add_input_float(
                            label="Minor allele frequency (MAF)",
                            width=220, default_value=0.05, step=0.005, tag='tooltip_maf'
                        )
                        self._inputs.append(maf_input)

                        dpg.add_spacer(height=6)
                        geno_input = dpg.add_input_float(
                            label="Missing genotype rate",
                            width=220, default_value=0.10, step=0.005, tag='tooltip_missing'
                        )
                        self._inputs.append(geno_input)

                        dpg.add_spacer(height=14)
                        convert_btn = dpg.add_button(
                            label="Convert VCF",
                            callback=self.convert_vcf,
                            user_data=[maf_input, geno_input, vcf, variant_ids],
                            width=160, height=36
                        )
                        self._primary_buttons.append(convert_btn)

            # -------- Settings --------
            with dpg.tab(label='Settings', parent="main_tabbar"):
                dpg.add_spacer(height=10)
                with dpg.table(header_row=False, borders_innerH=False, borders_outerH=False,
                               borders_innerV=False, borders_outerV=False, resizable=False):
                    dpg.add_table_column()
                    dpg.add_table_column()
                    with dpg.table_row():
                        with dpg.group():
                            dpg.add_text("General Settings", color=(200, 180, 90))
                            dpg.add_spacer(height=8)
                            self.night_toggle = dpg.add_checkbox(
                                label="Night Mode (Dark)", default_value=True, callback=self.toggle_theme
                            )
                            dpg.add_spacer(height=6)
                            self.nr_jobs = dpg.add_input_int(
                                label="Number of jobs to run", width=220,
                                default_value=-1, step=1, min_value=-1, max_value=50,
                                min_clamped=True, max_clamped=True, tag="tooltip_nr_jobs"
                            )
                            self._inputs.append(self.nr_jobs)
                            self.gb_goal = dpg.add_input_int(
                                label="Gigabytes of memory per run", width=220,
                                default_value=0, step=4, min_value=0, max_value=512,
                                min_clamped=True, max_clamped=True, tag="tooltip_gb_goal"
                            )
                            self._inputs.append(self.gb_goal)
                            self.plot_stats = dpg.add_checkbox(
                                label="Advanced Plotting", default_value=False, tag="tooltip_stats"
                            )
                            self.snp_limit = dpg.add_input_text(
                                label="SNP limit", width=220, default_value='', tag="tooltip_limit"
                            )
                            self._inputs.append(self.snp_limit)

                        with dpg.group():
                            dpg.add_text("Machine Learning Settings", color=(200, 180, 90))
                            dpg.add_spacer(height=8)
                            self.train_size_set = dpg.add_input_int(
                                label="Training size %", width=220, default_value=70, step=10,
                                min_value=0, max_value=100, min_clamped=True, max_clamped=True,
                                tag="tooltip_training"
                            )
                            self._inputs.append(self.train_size_set)
                            self.estim_set = dpg.add_input_int(
                                label="Number of trees", width=220, default_value=200, step=10,
                                min_value=1, min_clamped=True, tag="tooltip_trees"
                            )
                            self._inputs.append(self.estim_set)
                            self.max_dep_set = dpg.add_input_int(
                                label="Max depth", width=220, default_value=3, step=10,
                                min_value=0, max_value=100, min_clamped=True, max_clamped=True, tag="tooltip_depth"
                            )
                            self._inputs.append(self.max_dep_set)
                            self.model_nr = dpg.add_input_int(
                                label="Nr. of models", width=220, default_value=1, step=1,
                                min_value=1, max_value=50, min_clamped=True, max_clamped=True, tag="tooltip_model"
                            )
                            self._inputs.append(self.model_nr)
                            self.aggregation_method = dpg.add_combo(
                                ("sum", "median", "mean"),
                                label="Aggregation Method", width=220, default_value='sum', tag="tooltip_aggr"
                            )
                            self._inputs.append(self.aggregation_method)

        # Log (dockable)
        self.log_win = dpg.add_window(label="Log", tag="LogWindow", width=800, height=400)
        self.logz = logger.mvLogger(self.log_win)

        # Results (dockable, initially hidden)
        dpg.add_window(label="Results", tag="ResultsWindow", width=1000, height=600, show=False)

        # Tooltips
        with dpg.tooltip("tooltip_vcf"):
            dpg.add_text("Select a Variant Call Format file (.vcf or .vcf.gz).", color=[79, 128, 90])
        with dpg.tooltip("tooltip_variant"):
            dpg.add_text("Optional sample IDs list (PLINK --keep): FID IID (space-separated).", color=[79, 128, 90])
        with dpg.tooltip("tooltip_maf"):
            dpg.add_text("Minor Allele Frequency threshold.", color=[79, 128, 90])
        with dpg.tooltip("tooltip_missing"):
            dpg.add_text("Maximum allowed missing genotype rate per variant.", color=[79, 128, 90])
        with dpg.tooltip("tooltip_bed"):
            dpg.add_text("Select a PLINK .bed file (needs .bim and .fam).", color=[79, 128, 90])
        with dpg.tooltip("tooltip_pheno"):
            dpg.add_text("Phenotype file: FID IID Value (no header).", color=[79, 128, 90])
        with dpg.tooltip("tooltip_cov"):
            dpg.add_text("Covariates file: FID IID <cov1> <cov2> ...", color=[79, 128, 90])
        with dpg.tooltip("tooltip_algorithm"):
            dpg.add_text("Select the algorithm to use for analysis.", color=[79, 128, 90])
        with dpg.tooltip("tooltip_algorithm_gp"):
            dpg.add_text("Select the algorithm for genomic prediction.", color=[79, 128, 90])
        with dpg.tooltip("tooltip_training"):
            dpg.add_text("Percent of data used for training.", color=[79, 128, 90])
        with dpg.tooltip("tooltip_trees"):
            dpg.add_text("Number of trees (RF/XGB).", color=[79, 128, 90])
        with dpg.tooltip("tooltip_model"):
            dpg.add_text("Number of models for aggregation.", color=[79, 128, 90])
        with dpg.tooltip("tooltip_depth"):
            dpg.add_text("Maximum tree depth.", color=[79, 128, 90])
        with dpg.tooltip("tooltip_nr_jobs"):
            dpg.add_text("CPU cores (-1 = use all cores).", color=[79, 128, 90])
        with dpg.tooltip("tooltip_gb_goal"):
            dpg.add_text("Target GB of RAM per run. 0 = block-wise reading.", color=[79, 128, 90])
        with dpg.tooltip("tooltip_limit"):
            dpg.add_text("Limit SNPs in plots on huge datasets. Empty = all.", color=[79, 128, 90])
        with dpg.tooltip("tooltip_stats"):
            dpg.add_text("Enable advanced PDF plots for pheno/geno stats.", color=[79, 128, 90])
        with dpg.tooltip("tooltip_aggr"):
            dpg.add_text("Aggregate SNP effects: sum / median / mean.", color=[79, 128, 90])

        # Initial theme + per-item bindings
        self.apply_component_themes()

    # ---------- Header (icon + texts) ----------
    def build_header(self, parent):
        if self.header_group and dpg.does_item_exist(self.header_group):
            dpg.delete_item(self.header_group)

        with dpg.group(horizontal=True, parent=parent) as self.header_group:
            self.logo_draw = dpg.add_drawlist(width=48, height=48)
            self.brand_title = dpg.add_text("PlantVarFilter")
            self.brand_tagline = dpg.add_text("Software for GWAS Analysis")

        self.redraw_header_icon_and_colors()

    def header_palette(self):
        if self.night_mode:
            return {
                "ring": (34, 70, 52, 255),
                "leaf": (76, 175, 110, 255),
                "vein": (230, 245, 235, 220),
                "title": (210, 230, 210),
                "tagline": (190, 175, 95)
            }
        else:
            return {
                "ring": (190, 220, 200, 255),
                "leaf": (72, 160, 100, 255),
                "vein": (25, 50, 35, 200),
                "title": (30, 45, 35),
                "tagline": (90, 120, 70)
            }

    def redraw_header_icon_and_colors(self):
        pal = self.header_palette()

        # clear and draw simple leaf icon
        if self.logo_draw and dpg.does_item_exist(self.logo_draw):
            dpg.delete_item(self.logo_draw, children_only=True)
            # outer soft ring
            dpg.draw_circle(center=(24, 24), radius=22, color=pal["ring"], fill=(0, 0, 0, 0), parent=self.logo_draw, thickness=2)
            # leaf polygon
            leaf_points = [(24, 10), (38, 22), (22, 38), (14, 26)]
            dpg.draw_polygon(points=leaf_points, color=(0, 0, 0, 0), fill=pal["leaf"], parent=self.logo_draw)
            # vein
            dpg.draw_line(p1=(24, 12), p2=(24, 34), color=pal["vein"], thickness=2, parent=self.logo_draw)
            dpg.draw_line(p1=(24, 22), p2=(34, 22), color=pal["vein"], thickness=2, parent=self.logo_draw)

        # text colors
        if self.brand_title and dpg.does_item_exist(self.brand_title):
            dpg.configure_item(self.brand_title, color=pal["title"])
        if self.brand_tagline and dpg.does_item_exist(self.brand_tagline):
            dpg.configure_item(self.brand_tagline, color=pal["tagline"])

    # ---------- Theme binding helpers ----------
    def apply_component_themes(self):
        dpg.bind_theme(self.theme_plant_dark if self.night_mode else self.theme_plant_light)

        primary = self.theme_primary_btn_dark if self.night_mode else self.theme_primary_btn_light
        secondary = self.theme_secondary_btn_dark if self.night_mode else self.theme_secondary_btn_light
        inputs = self.theme_inputs_dark if self.night_mode else self.theme_inputs_light
        fd_theme = self.theme_fd_dark if self.night_mode else self.theme_fd_light

        for b in self._primary_buttons:
            if dpg.does_item_exist(b):
                dpg.bind_item_theme(b, primary)
        for b in self._secondary_buttons:
            if dpg.does_item_exist(b):
                dpg.bind_item_theme(b, secondary)
        for it in self._inputs:
            if dpg.does_item_exist(it):
                dpg.bind_item_theme(it, inputs)
        for t in self._file_dialogs:
            if dpg.does_item_exist(t):
                dpg.bind_item_theme(t, fd_theme)

        # update header colors
        self.redraw_header_icon_and_colors()

        if self.font:
            dpg.bind_font(self.font)
        dpg.set_global_font_scale(0.62)

    # ---------- Callbacks (files) ----------
    def callback_vcf(self, s, app_data):
        self.vcf_app_data = app_data
        vcf_path, current_path = self.get_selection_path(self.vcf_app_data)
        dpg.configure_item("file_dialog_variants", default_path=current_path)
        self.add_log('VCF Selected: ' + vcf_path)

    def callback_bed(self, s, app_data):
        self.bed_app_data = app_data
        try:
            bed_path, current_path = self.get_selection_path(self.bed_app_data)
            dpg.configure_item("file_dialog_cov", default_path=current_path)
            dpg.configure_item("file_dialog_pheno", default_path=current_path)
            self.add_log('BED file Selected: ' + bed_path)
        except TypeError:
            self.add_log('Invalid BED file Selected', error=True)

    def callback_variants(self, s, app_data):
        self.variants_app_data = app_data
        variants_path, current_path = self.get_selection_path(self.variants_app_data)
        dpg.configure_item("file_dialog_vcf", default_path=current_path)
        self.add_log('IDs file Selected: ' + variants_path)

    def callback_pheno(self, s, app_data):
        self.pheno_app_data = app_data
        try:
            pheno_path, current_path = self.get_selection_path(self.pheno_app_data)
            dpg.configure_item("file_dialog_cov", default_path=current_path)
            dpg.configure_item("file_dialog_bed", default_path=current_path)
            self.add_log('Pheno File Selected: ' + pheno_path)
        except TypeError:
            self.add_log('Wrong Pheno File Selected', error=True)

    def callback_cov(self, s, app_data):
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

    def callback_save_results(self, s, app_data):
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

    def cancel_callback_directory(self, s, app_data):
        self.add_log('Process Canceled')

    # ---------- Actions ----------
    def convert_vcf(self, s, data, user_data):
        maf = str(dpg.get_value(user_data[0]))
        geno = str(dpg.get_value(user_data[1]))
        vcf_path, _ = self.get_selection_path(self.vcf_app_data)
        variants_path = None if self.variants_app_data is None else self.get_selection_path(self.variants_app_data)[0]
        self.add_log('Start converting VCF to BED...')
        out_prefix = f"{vcf_path.split('.')[0]}_maf{round(float(maf), 2)}_geno{round(float(geno), 2)}"
        plink_log = self.gwas.vcf_to_bed(vcf_path, variants_path, out_prefix, maf, geno)
        self.add_log(plink_log)

    def run_gwas(self, s, data, user_data):
        self.delete_files()
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
            bed_path = self.get_selection_path(self.bed_app_data)[0]
            pheno_path = self.get_selection_path(self.pheno_app_data)[0]
            cov_path = None
            try:
                cov_path = self.get_selection_path(self.cov_app_data)[0]
            except Exception:
                pass

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

                self.add_log(f"Dataset after intersection: SNPs: {bed.sid_count}  Pheno IDs: {pheno.iid_count}", warn=True)
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

                self.bed_app_data = None
                self.pheno_app_data = None
                self.cov_app_data = None
            else:
                self.add_log('Error, GWAS Analysis could not be started.', error=True)

        except TypeError:
            self.add_log('Please select a phenotype and genotype file.', error=True)

    def run_genomic_prediction(self, s, data, user_data):
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
            bed_path = self.get_selection_path(self.bed_app_data)[0]
            pheno_path = self.get_selection_path(self.pheno_app_data)[0]

            self.add_log('Validating files...')
            check_input_data = self.gwas.validate_gwas_input_files(bed_path, pheno_path)
            chrom_mapping = self.helper.replace_with_integers(bed_path.replace('.bed', '.bim'))
            self.settings_lst = [self.algorithm, bed_path, pheno_path, test_size, estimators, model_nr, max_dep_set]

            if check_input_data[0]:
                bed = Bed(str(bed_path), count_A1=False, chrom_map=chrom_mapping)
                pheno = Pheno(str(pheno_path))
                bed, pheno = pstutil.intersect_apply([bed, pheno])
                bed_fixed = self.gwas.filter_out_missing(bed)

                self.add_log(f"Dataset after intersection: SNPs: {bed.sid_count}  Pheno IDs: {pheno.iid_count}", warn=True)
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
            self.add_log('Please select a phenotype and genotype file.', error=True)

    # ---------- Results window ----------
    def show_results_window(self, df, algorithm, genomic_predict):
        dpg.configure_item("ResultsWindow", show=True)
        dpg.delete_item("ResultsWindow", children_only=True)

        dpg.add_button(label="Export Results", parent="ResultsWindow",
                       callback=lambda: dpg.show_item("select_directory"))
        dpg.add_spacer(height=10, parent="ResultsWindow")

        if genomic_predict:
            w1, h1, c1, data1 = dpg.load_image(self.gp_plot_name)
            with dpg.texture_registry(show=False):
                dpg.add_static_texture(width=w1, height=h1, default_value=data1, tag="ba_tag")
            w2, h2, c2, data2 = dpg.load_image(self.gp_plot_name_scatter)
            with dpg.texture_registry(show=False):
                dpg.add_static_texture(width=w2, height=h2, default_value=data2, tag="ba_tag2")

            with dpg.tab_bar(label='tabbar_results_gp', parent="ResultsWindow"):
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
                                    v = df.iloc[i, j]
                                    dpg.add_text(f"{v:.2f}" if isinstance(v, float) else str(v))
                with dpg.tab(label="Correlation Plot (Predicted vs. Phenotype)"):
                    dpg.add_image(texture_tag="ba_tag2", tag="ba_image2", width=750, height=450)
                with dpg.tab(label="Bland-Altman Plot (Model Accuracy)"):
                    dpg.add_image(texture_tag="ba_tag", tag="ba_image", width=750, height=450)

        else:
            w, h, c, data = dpg.load_image(self.manhatten_plot_name)
            with dpg.texture_registry(show=False):
                dpg.add_static_texture(width=w, height=h, default_value=data, tag="manhatten_tag")

            with dpg.tab_bar(label='tabbar_results_gwas', parent="ResultsWindow"):
                with dpg.tab(label="Manhattan Plot"):
                    if algorithm in ("FaST-LMM", "Linear regression"):
                        dpg.add_image(texture_tag="manhatten_tag", tag="manhatten_image", width=950, height=400)
                    else:
                        dpg.add_image(texture_tag="manhatten_tag", tag="manhatten_image", width=900, height=300)

                if algorithm in ("FaST-LMM", "Linear regression"):
                    w2, h2, c2, data2 = dpg.load_image(self.qq_plot_name)
                    with dpg.texture_registry(show=False):
                        dpg.add_static_texture(width=w2, height=h2, default_value=data2, tag="qq_tag")
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
                                    dpg.add_text(df.iloc[i, j])

    # ---------- Theme toggle ----------
    def toggle_theme(self, sender, app_data):
        self.night_mode = bool(app_data)
        self.apply_component_themes()
        self.add_log("Dark mode enabled" if self.night_mode else "Light mode enabled")

    # ---------- Utils ----------
    def delete_files(self):
        for tag in ["manhatten_image", "manhatten_tag", "qq_image", "qq_tag",
                    "table_gwas", "table_gp", "ba_tag", "ba_tag2", "ba_image", "ba_image2"]:
            if dpg.does_item_exist(tag):
                dpg.delete_item(tag)
        for f in [
            self.gwas_result_name, self.gwas_result_name_top, self.genomic_predict_name, self.gp_plot_name,
            self.manhatten_plot_name, self.qq_plot_name, self.gp_plot_name_scatter,
            self.manhatten_plot_name.replace('manhatten_plot', 'manhatten_plot_high'),
            self.qq_plot_name.replace('qq_plot', 'qq_plot_high'),
            self.gp_plot_name_scatter.replace('GP_scatter_plot', 'GP_scatter_plot_high'),
            self.gp_plot_name.replace('Bland_Altman_plot', 'Bland_Altman_plot_high'),
            self.genomic_predict_name.replace('.csv', '_valdation.csv'),
            self.pheno_stats_name, self.geno_stats_name
        ]:
            if os.path.exists(f):
                try:
                    os.remove(f)
                except Exception:
                    pass

    def add_log(self, message, warn=False, error=False):
        if warn:
            self.logz.log_warning(message)
        elif error:
            self.logz.log_error(message)
        else:
            self.logz.log_info(message)

    # ---------- Runner ----------
    def run(self):
        self.apply_component_themes()
        dpg.setup_dearpygui()
        dpg.show_viewport()
        dpg.start_dearpygui()


if __name__ == "__main__":
    main()
