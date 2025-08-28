# main_gui.py — PlantVarFilter
# (plant-themed UI + docking + improved dialogs + header icon + VCF QC tab + QC plots
#  + file labels + samtools preprocess + bcftools preprocess + variant calling BAM→VCF)

import os
import shutil
import subprocess
import dearpygui.dearpygui as dpg
from dearpygui_ext import logger

from vcf_quality import VCFQualityChecker

from PlantVarFilter.gwas_pipeline import GWAS
from PlantVarFilter.genomic_prediction_pipeline import GenomicPrediction
from PlantVarFilter.helpers import HELPERS
from PlantVarFilter.pipeline_plots import Plot

from pysnptools.snpreader import Bed, Pheno
import pysnptools.util as pstutil

from bcftools_utils import BCFtools, BCFtoolsError
from samtools_utils import Samtools, SamtoolsError

# —— variant calling helper (new) ——
try:
    from variant_caller_utils import VariantCaller, VariantCallerError
except Exception:  # graceful fallback if file missing
    VariantCaller = None

# Prefer bundled binaries under PlantVarFilter/linux, else PATH
try:
    from PlantVarFilter.linux import resolve_tool
except Exception:
    def resolve_tool(name: str) -> str:
        return shutil.which(name) or name


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

        self.vcf_qc_checker = VCFQualityChecker(max_sites_scan=200_000, min_sites_required=200)
        self.blacklist_app_data = None

        # bcftools helpers
        self.bcft = BCFtools(
            bcftools_bin=resolve_tool("bcftools"),
            bgzip_bin=resolve_tool("bgzip"),
            tabix_bin=resolve_tool("tabix")
        )
        self.fasta_app_data = None
        self.bcf_out_last = None

        # samtools helper
        self.sam = Samtools(exe=resolve_tool("samtools"))
        self.bam_app_data = None
        self.sam_out_last = None

        # variant caller (bcftools mpileup+call)
        self.vcaller = VariantCaller(
            bcftools_bin=resolve_tool("bcftools"),
            bgzip_bin=resolve_tool("bgzip"),
            tabix_bin=resolve_tool("tabix")
        ) if VariantCaller else None
        self.bam_vc_app_data = None        # single BAM for calling
        self.bamlist_app_data = None       # .list of BAMs for calling
        self.vc_out_last = None            # produced VCF

        with dpg.font_registry():
            script_dir = os.path.dirname(__file__)
            font_path = os.path.join(script_dir, "test.ttf")
            self.font = dpg.add_font(font_path, 40, tag="ttf-font") if os.path.exists(font_path) else None

        # ——— Themes ———
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

        self._primary_buttons = []
        self._secondary_buttons = []
        self._inputs = []
        self._file_dialogs = []
        self.night_mode = True

        self.header_group = None
        self.brand_title = None
        self.brand_tagline = None
        self.logo_draw = None

        self.vcf_app_data = None
        self.variants_app_data = None
        self.results_directory = None
        self.bed_app_data = None
        self.pheno_app_data = None
        self.cov_app_data = None
        self.default_path = os.path.expanduser("~")

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
        self._check_cli_versions()

    # ——— UI ———
    def setup_gui(self):
        dpg.create_viewport(title='PlantVarFilter', width=2000, height=1200, resizable=True)

        with dpg.window(tag="PrimaryWindow", no_close=True, no_move=True, no_resize=True,
                        no_title_bar=True, no_background=True):
            pass
        dpg.set_primary_window("PrimaryWindow", True)

        # File dialogs
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

        with dpg.file_dialog(directory_selector=False, show=False, callback=self.callback_blacklist,
                             file_count=1, tag="file_dialog_blacklist", width=920, height=560,
                             default_path=self.default_path, modal=True):
            dpg.add_file_extension(".bed", color=(255, 200, 150, 255))
            dpg.add_file_extension(".*")
        self._file_dialogs.append("file_dialog_blacklist")

        # FASTA for bcftools/variant calling
        with dpg.file_dialog(directory_selector=False, show=False, callback=self.callback_fasta,
                             file_count=1, tag="file_dialog_fasta", width=920, height=560,
                             default_path=self.default_path, modal=True):
            dpg.add_file_extension("Reference (*.fa *.fasta *.fa.gz){.fa,.fasta,.fa.gz}", color=(150, 200, 255, 255))
        self._file_dialogs.append("file_dialog_fasta")

        # BAM for preprocess
        with dpg.file_dialog(directory_selector=False, show=False, callback=self.callback_bam,
                             file_count=1, tag="file_dialog_bam", width=920, height=560,
                             default_path=self.default_path, modal=True):
            dpg.add_file_extension(".bam", color=(255, 180, 120, 255))
            dpg.add_file_extension(".*")
        self._file_dialogs.append("file_dialog_bam")

        # BAM (single) for variant calling
        with dpg.file_dialog(directory_selector=False, show=False, callback=self.callback_bam_vc,
                             file_count=1, tag="file_dialog_bam_vc", width=920, height=560,
                             default_path=self.default_path, modal=True):
            dpg.add_file_extension(".bam", color=(255, 180, 120, 255))
            dpg.add_file_extension(".*")
        self._file_dialogs.append("file_dialog_bam_vc")

        # BAM list for variant calling
        with dpg.file_dialog(directory_selector=False, show=False, callback=self.callback_bamlist_vc,
                             file_count=1, tag="file_dialog_bamlist", width=920, height=560,
                             default_path=self.default_path, modal=True):
            dpg.add_file_extension("BAM list (*.list){.list}", color=(255, 230, 140, 255))
            dpg.add_file_extension(".*")
        self._file_dialogs.append("file_dialog_bamlist")

        dpg.add_file_dialog(directory_selector=True, show=False, callback=self.callback_save_results,
                            tag="select_directory", cancel_callback=self.cancel_callback_directory,
                            width=900, height=540, default_path=self.default_path, modal=True)
        self._file_dialogs.append("select_directory")

        # Main window
        with dpg.window(label="PlantVarFilter", tag="MainWindow", no_close=True, width=1100, height=700):
            self.build_header(parent="MainWindow")

            dpg.add_spacer(height=4)
            dpg.add_separator()
            dpg.add_spacer(height=6)

            dpg.add_tab_bar(tag="main_tabbar")

            # —— Preprocess (samtools) ——
            with dpg.tab(label='Preprocess (samtools)', parent="main_tabbar"):
                dpg.add_text("\nClean BAM: sort / fixmate / markdup / index + QC reports", indent=10)
                dpg.add_spacer(height=10)
                with dpg.group(horizontal=True, horizontal_spacing=60):
                    with dpg.group():
                        bam_btn = dpg.add_button(
                            label="Choose a BAM file",
                            callback=lambda: dpg.show_item("file_dialog_bam"),
                            width=220, tag='tooltip_bam_sam'
                        )
                        self._secondary_buttons.append(bam_btn)
                        dpg.add_text("", tag="sam_bam_path_lbl", wrap=500)

                    with dpg.group():
                        self.sam_threads = dpg.add_input_int(
                            label="Threads", width=220, default_value=4, min_value=1, min_clamped=True
                        )
                        self._inputs.append(self.sam_threads)
                        self.sam_remove_dups = dpg.add_checkbox(
                            label="Remove duplicates (instead of marking)", default_value=False
                        )
                        self._inputs.append(self.sam_remove_dups)
                        self.sam_compute_stats = dpg.add_checkbox(
                            label="Compute QC reports (flagstat/stats/idxstats/depth)", default_value=True
                        )
                        self._inputs.append(self.sam_compute_stats)

                        dpg.add_spacer(height=6)
                        self.sam_out_prefix = dpg.add_input_text(
                            label="Output prefix (optional)",
                            hint="Leave empty to auto-generate next to the input file",
                            width=320
                        )
                        self._inputs.append(self.sam_out_prefix)

                        dpg.add_spacer(height=12)
                        run_sam = dpg.add_button(
                            label="Run samtools preprocess",
                            callback=self.run_samtools_preprocess,
                            width=240, height=38
                        )
                        self._primary_buttons.append(run_sam)

            # —— Variant Calling (BAM → VCF) ——
            with dpg.tab(label='Variant Calling (BAM→VCF)', parent="main_tabbar"):
                dpg.add_text("\nCall variants with bcftools mpileup + call", indent=10)
                dpg.add_spacer(height=10)
                with dpg.group(horizontal=True, horizontal_spacing=60):

                    with dpg.group():
                        b1 = dpg.add_button(
                            label="Choose BAM (single)",
                            callback=lambda: dpg.show_item("file_dialog_bam_vc"),
                            width=220, tag='tooltip_bam_vc'
                        )
                        self._secondary_buttons.append(b1)
                        dpg.add_text("", tag="vc_bam_path_lbl", wrap=500)

                        dpg.add_spacer(height=6)
                        b2 = dpg.add_button(
                            label="Choose BAM-list (.list)",
                            callback=lambda: dpg.show_item("file_dialog_bamlist"),
                            width=220, tag='tooltip_bamlist_vc'
                        )
                        self._secondary_buttons.append(b2)
                        dpg.add_text("", tag="vc_bamlist_path_lbl", wrap=500)

                        dpg.add_spacer(height=6)
                        fa_btn = dpg.add_button(
                            label="Choose reference FASTA",
                            callback=lambda: dpg.show_item("file_dialog_fasta"),
                            width=220, tag='tooltip_fa_vc'
                        )
                        self._secondary_buttons.append(fa_btn)
                        dpg.add_text("", tag="vc_ref_path_lbl", wrap=500)

                        dpg.add_spacer(height=6)
                        reg_btn2 = dpg.add_button(
                            label="Choose regions BED (optional)",
                            callback=lambda: dpg.show_item("file_dialog_blacklist"),
                            width=220, tag='tooltip_reg_vc'
                        )
                        self._secondary_buttons.append(reg_btn2)
                        dpg.add_text("", tag="vc_regions_path_lbl", wrap=500)

                    with dpg.group():
                        self.vc_threads = dpg.add_input_int(label="Threads", width=220, default_value=4, min_value=1, min_clamped=True)
                        self._inputs.append(self.vc_threads)
                        self.vc_ploidy = dpg.add_input_int(label="Ploidy", width=220, default_value=2, min_value=1, min_clamped=True, tag='tooltip_ploidy')
                        self._inputs.append(self.vc_ploidy)
                        self.vc_min_bq = dpg.add_input_int(label="Min BaseQ", width=220, default_value=20, min_value=0, min_clamped=True, tag='tooltip_bq')
                        self._inputs.append(self.vc_min_bq)
                        self.vc_min_mq = dpg.add_input_int(label="Min MapQ", width=220, default_value=20, min_value=0, min_clamped=True, tag='tooltip_mq')
                        self._inputs.append(self.vc_min_mq)

                        dpg.add_spacer(height=6)
                        self.vc_out_prefix = dpg.add_input_text(
                            label="Output prefix (optional)",
                            hint="Leave empty to auto-generate next to the BAM",
                            width=320
                        )
                        self._inputs.append(self.vc_out_prefix)

                        dpg.add_spacer(height=12)
                        run_vc = dpg.add_button(
                            label="Call variants (bcftools)",
                            callback=self.run_variant_calling,
                            width=240, height=38
                        )
                        self._primary_buttons.append(run_vc)

            # —— Preprocess (bcftools) ——
            with dpg.tab(label='Preprocess (bcftools)', parent="main_tabbar"):
                dpg.add_text("\nNormalize / split multiallelic / sort / filter / set IDs (bcftools)", indent=10)
                dpg.add_spacer(height=10)
                with dpg.group(horizontal=True, horizontal_spacing=60):

                    with dpg.group():
                        vcf_btn_bcf = dpg.add_button(
                            label="Choose a VCF file",
                            callback=lambda: dpg.show_item("file_dialog_vcf"),
                            width=220, tag='tooltip_vcf_bcf'
                        )
                        self._secondary_buttons.append(vcf_btn_bcf)
                        dpg.add_text("", tag="bcf_vcf_path_lbl", wrap=500)

                        dpg.add_spacer(height=6)
                        fasta_btn = dpg.add_button(
                            label="Choose reference FASTA (for left-align)",
                            callback=lambda: dpg.show_item("file_dialog_fasta"),
                            width=220, tag='tooltip_fa_bcf'
                        )
                        self._secondary_buttons.append(fasta_btn)
                        dpg.add_text("", tag="bcf_ref_path_lbl", wrap=500)

                        dpg.add_spacer(height=6)
                        reg_btn = dpg.add_button(
                            label="Choose regions BED (optional)",
                            callback=lambda: dpg.show_item("file_dialog_blacklist"),
                            width=220, tag='tooltip_reg_bcf'
                        )
                        self._secondary_buttons.append(reg_btn)
                        dpg.add_text("", tag="bcf_regions_path_lbl", wrap=500)

                    with dpg.group():
                        self.bcf_split = dpg.add_checkbox(label="Split multiallelic", default_value=True)
                        self._inputs.append(self.bcf_split)
                        self.bcf_left  = dpg.add_checkbox(label="Left-align indels (needs FASTA)", default_value=True)
                        self._inputs.append(self.bcf_left)
                        self.bcf_sort  = dpg.add_checkbox(label="Sort", default_value=True)
                        self._inputs.append(self.bcf_sort)
                        self.bcf_setid = dpg.add_checkbox(label="Set ID to CHR:POS:REF:ALT", default_value=True)
                        self._inputs.append(self.bcf_setid)
                        self.bcf_compr = dpg.add_checkbox(label="Compress output (.vcf.gz)", default_value=True)
                        self._inputs.append(self.bcf_compr)
                        self.bcf_index = dpg.add_checkbox(label="Index output (tabix)", default_value=True)
                        self._inputs.append(self.bcf_index)
                        self.bcf_rmflt = dpg.add_checkbox(label="Keep only PASS (remove filtered)", default_value=False)
                        self._inputs.append(self.bcf_rmflt)

                        dpg.add_spacer(height=6)
                        self.bcf_filter_expr = dpg.add_input_text(
                            label="bcftools filter expression (optional)",
                            hint="Example: QUAL>=30 && INFO/DP>=10",
                            width=320
                        )
                        self._inputs.append(self.bcf_filter_expr)

                        dpg.add_spacer(height=6)
                        self.bcf_out_prefix = dpg.add_input_text(
                            label="Output prefix (optional)",
                            hint="Leave empty to auto-generate next to the input file",
                            width=320
                        )
                        self._inputs.append(self.bcf_out_prefix)

                        dpg.add_spacer(height=12)
                        run_bcf = dpg.add_button(
                            label="Run bcftools preprocess",
                            callback=self.run_bcftools_preprocess,
                            width=240, height=38
                        )
                        self._primary_buttons.append(run_bcf)

            # —— Check VCF File (QC) ——
            with dpg.tab(label='Check VCF File', parent="main_tabbar"):
                dpg.add_text("\nCheck VCF quality before conversion/analysis", indent=10)
                dpg.add_spacer(height=10)
                with dpg.group(horizontal=True, horizontal_spacing=60):
                    with dpg.group():
                        vcf_btn_qc = dpg.add_button(
                            label="Choose a VCF file",
                            callback=lambda: dpg.show_item("file_dialog_vcf"),
                            width=220, tag='tooltip_vcf_qc'
                        )
                        self._secondary_buttons.append(vcf_btn_qc)
                        dpg.add_text("", tag="qc_vcf_path_lbl", wrap=500)

                        dpg.add_spacer(height=6)
                        bl_btn = dpg.add_button(
                            label="Choose blacklist BED (optional)",
                            callback=lambda: dpg.show_item("file_dialog_blacklist"),
                            width=220, tag='tooltip_bl_qc'
                        )
                        self._secondary_buttons.append(bl_btn)
                        dpg.add_text("", tag="qc_bl_path_lbl", wrap=500)

                    with dpg.group():
                        self.deep_scan = dpg.add_checkbox(label="Deep scan", default_value=False)
                        self._inputs.append(self.deep_scan)
                        dpg.add_spacer(height=12)
                        run_qc_btn = dpg.add_button(
                            label="Run Quality Check",
                            callback=self.run_vcf_qc,
                            width=200, height=36
                        )
                        self._primary_buttons.append(run_qc_btn)

            # —— Convert to PLINK ——
            with dpg.tab(label='Convert to PLINK', parent="main_tabbar"):
                dpg.add_text("\nConvert a VCF file into PLINK BED and apply MAF/missing genotype filters.", indent=10)
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
                        dpg.add_text("", tag="conv_vcf_path_lbl", wrap=500)

                        dpg.add_spacer(height=6)
                        variant_ids = dpg.add_button(
                            label="Choose IDs file (option)",
                            callback=lambda: dpg.show_item("file_dialog_variants"),
                            width=220, tag='tooltip_variant'
                        )
                        self._secondary_buttons.append(variant_ids)
                        dpg.add_text("", tag="ids_path_lbl", wrap=500)

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

            # —— GWAS Analysis ——
            with dpg.tab(label='GWAS Analysis', parent="main_tabbar"):
                dpg.add_text("\nStart GWAS Analysis", indent=10)
                dpg.add_spacer(height=10)
                with dpg.group(horizontal=True, horizontal_spacing=60):
                    with dpg.group():
                        geno = dpg.add_button(label="Choose a BED file",
                                              callback=lambda: dpg.show_item("file_dialog_bed"),
                                              width=220, tag='tooltip_bed')
                        self._secondary_buttons.append(geno)
                        dpg.add_text("", tag="gwas_bed_path_lbl", wrap=500)

                        dpg.add_spacer(height=6)
                        pheno = dpg.add_button(label="Choose a phenotype file",
                                               callback=lambda: dpg.show_item("file_dialog_pheno"),
                                               width=220, tag='tooltip_pheno')
                        self._secondary_buttons.append(pheno)
                        dpg.add_text("", tag="gwas_pheno_path_lbl", wrap=500)

                        dpg.add_spacer(height=6)
                        cov_file = dpg.add_button(label="Choose covariate file (option)",
                                                  callback=lambda: dpg.show_item("file_dialog_cov"),
                                                  width=220, tag='tooltip_cov')
                        self._secondary_buttons.append(cov_file)
                        dpg.add_text("", tag="gwas_cov_path_lbl", wrap=500)

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

            # —— Genomic Prediction ——
            with dpg.tab(label='Genomic Prediction', parent="main_tabbar"):
                dpg.add_text("\nStart Genomic Prediction", indent=10)
                dpg.add_spacer(height=10)
                with dpg.group(horizontal=True, horizontal_spacing=60):
                    with dpg.group():
                        geno = dpg.add_button(label="Choose a BED file",
                                              callback=lambda: dpg.show_item("file_dialog_bed"),
                                              width=220, tag='tooltip_bed_gp')
                        self._secondary_buttons.append(geno)
                        dpg.add_text("", tag="gp_bed_path_lbl", wrap=500)

                        dpg.add_spacer(height=6)
                        pheno = dpg.add_button(label="Choose a phenotype file",
                                               callback=lambda: dpg.show_item("file_dialog_pheno"),
                                               width=220, tag='tooltip_pheno_gp')
                        self._secondary_buttons.append(pheno)
                        dpg.add_text("", tag="gp_pheno_path_lbl", wrap=500)

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

            # —— Settings ——
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

        # Log + Results windows
        self.log_win = dpg.add_window(label="Log", tag="LogWindow", width=800, height=400)
        self.logz = logger.mvLogger(self.log_win)
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

        # samtools tooltip
        with dpg.tooltip("tooltip_bam_sam"):
            dpg.add_text("Select input BAM to clean and index.", color=[79, 128, 90])

        # variant calling tooltips
        with dpg.tooltip("tooltip_bam_vc"):
            dpg.add_text("Single sample BAM for calling.", color=[79, 128, 90])
        with dpg.tooltip("tooltip_bamlist_vc"):
            dpg.add_text("Text file with one BAM per line (-b list) for joint calling.", color=[79, 128, 90])
        with dpg.tooltip("tooltip_fa_vc"):
            dpg.add_text("Reference FASTA (must be indexed .fai).", color=[79, 128, 90])
        with dpg.tooltip("tooltip_reg_vc"):
            dpg.add_text("Optional BED of regions to restrict calling.", color=[79, 128, 90])
        with dpg.tooltip("tooltip_ploidy"):
            dpg.add_text("Assumed ploidy used by bcftools call.", color=[79, 128, 90])
        with dpg.tooltip("tooltip_bq"):
            dpg.add_text("Minimum base quality for mpileup.", color=[79, 128, 90])
        with dpg.tooltip("tooltip_mq"):
            dpg.add_text("Minimum mapping quality for mpileup.", color=[79, 128, 90])

        # bcftools tooltips
        with dpg.tooltip("tooltip_vcf_bcf"):
            dpg.add_text("Select input VCF/VCF.GZ to preprocess.", color=[79, 128, 90])
        with dpg.tooltip("tooltip_fa_bcf"):
            dpg.add_text("Reference FASTA is required for accurate left alignment in bcftools norm.", color=[79, 128, 90])
        with dpg.tooltip("tooltip_reg_bcf"):
            dpg.add_text("Regions BED to include (bcftools view -R). Optional.", color=[79, 128, 90])

        with dpg.tooltip("tooltip_vcf_qc"):
            dpg.add_text("Select a VCF/VCF.GZ file to evaluate.", color=[79, 128, 90])
        with dpg.tooltip("tooltip_bl_qc"):
            dpg.add_text("Optional BED of low-mappability/blacklist regions.", color=[79, 128, 90])

        self.apply_component_themes()

    # ——— Header ———
    def build_header(self, parent):
        if self.header_group and dpg.does_item_exist(self.header_group):
            dpg.delete_item(self.header_group)

        with dpg.group(horizontal=True, parent=parent) as self.header_group:
            logo_path = os.path.join(os.path.dirname(__file__), "assets", "logo.png")
            if os.path.exists(logo_path):
                w, h, c, data = dpg.load_image(logo_path)
                if not dpg.does_item_exist("plant_logo_tex"):
                    with dpg.texture_registry(show=False):
                        dpg.add_static_texture(width=w, height=h, default_value=data, tag="plant_logo_tex")
                dpg.add_image("plant_logo_tex", width=48, height=48)
            else:
                self.logo_draw = dpg.add_drawlist(width=48, height=48)

            self.brand_title = dpg.add_text("PlantVarFilter")
            self.brand_tagline = dpg.add_text("Software for GWAS Analysis")

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

        if self.logo_draw and dpg.does_item_exist(self.logo_draw):
            dpg.delete_item(self.logo_draw, children_only=True)
            dpg.draw_circle(center=(24, 24), radius=22, color=pal["ring"], fill=(0, 0, 0, 0), parent=self.logo_draw, thickness=2)

        if self.brand_title and dpg.does_item_exist(self.brand_title):
            dpg.configure_item(self.brand_title, color=pal["title"])
        if self.brand_tagline and dpg.does_item_exist(self.brand_tagline):
            dpg.configure_item(self.brand_tagline, color=pal["tagline"])

    # ——— Helpers ———
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

        self.redraw_header_icon_and_colors()

        if self.font:
            dpg.bind_font(self.font)
        dpg.set_global_font_scale(0.62)

    def _set_text(self, tag: str, value: str):
        if dpg.does_item_exist(tag):
            dpg.set_value(tag, value)

    def _fmt_name(self, path: str) -> str:
        try:
            return os.path.basename(path) or path
        except Exception:
            return str(path)

    def _get_appdata_path_safe(self, app_data):
        try:
            return self.get_selection_path(app_data)[0]
        except Exception:
            return None

    def _set_virtual_selection(self, attr_name: str, file_path: str):
        setattr(self, attr_name, {
            'current_path': os.path.dirname(file_path) + '/',
            'selections': {'file': file_path}
        })

    # ——— Callbacks (files) ———
    def callback_vcf(self, s, app_data):
        self.vcf_app_data = app_data
        vcf_path, current_path = self.get_selection_path(self.vcf_app_data)
        if not vcf_path:
            return
        dpg.configure_item("file_dialog_variants", default_path=current_path or self.default_path)
        self.add_log('VCF Selected: ' + vcf_path)
        name = self._fmt_name(vcf_path)
        self._set_text("qc_vcf_path_lbl", name)
        self._set_text("conv_vcf_path_lbl", name)
        self._set_text("bcf_vcf_path_lbl", name)

    def callback_bed(self, s, app_data):
        self.bed_app_data = app_data
        try:
            bed_path, current_path = self.get_selection_path(self.bed_app_data)
            if not bed_path:
                return
            dpg.configure_item("file_dialog_cov", default_path=current_path or self.default_path)
            dpg.configure_item("file_dialog_pheno", default_path=current_path or self.default_path)
            self.add_log('BED file Selected: ' + bed_path)
            name = self._fmt_name(bed_path)
            self._set_text("gwas_bed_path_lbl", name)
            self._set_text("gp_bed_path_lbl", name)
        except TypeError:
            self.add_log('Invalid BED file Selected', error=True)

    def callback_variants(self, s, app_data):
        self.variants_app_data = app_data
        variants_path, current_path = self.get_selection_path(self.variants_app_data)
        if not variants_path:
            return
        dpg.configure_item("file_dialog_vcf", default_path=current_path or self.default_path)
        self.add_log('IDs file Selected: ' + variants_path)
        self._set_text("ids_path_lbl", self._fmt_name(variants_path))

    def callback_pheno(self, s, app_data):
        self.pheno_app_data = app_data
        try:
            pheno_path, current_path = self.get_selection_path(self.pheno_app_data)
            if not pheno_path:
                return
            dpg.configure_item("file_dialog_cov", default_path=current_path or self.default_path)
            dpg.configure_item("file_dialog_bed", default_path=current_path or self.default_path)
            self.add_log('Pheno File Selected: ' + pheno_path)
            name = self._fmt_name(pheno_path)
            self._set_text("gwas_pheno_path_lbl", name)
            self._set_text("gp_pheno_path_lbl", name)
        except TypeError:
            self.add_log('Wrong Pheno File Selected', error=True)

    def callback_cov(self, s, app_data):
        self.cov_app_data = app_data
        try:
            cov_path, current_path = self.get_selection_path(self.cov_app_data)
            if not cov_path:
                return
            dpg.configure_item("file_dialog_bed", default_path=current_path or self.default_path)
            dpg.configure_item("file_dialog_pheno", default_path=current_path or self.default_path)
            self.add_log('Covariates File Selected: ' + cov_path)
            self._set_text("gwas_cov_path_lbl", self._fmt_name(cov_path))
        except TypeError:
            self.add_log('Wrong Covariates File Selected', error=True)

    def callback_blacklist(self, s, app_data):
        self.blacklist_app_data = app_data
        try:
            bl_path, _ = self.get_selection_path(self.blacklist_app_data)
            if not bl_path:
                return
            self.add_log('Blacklist/Regions BED Selected: ' + bl_path)
            # update all places that show a regions BED
            for tag in ("qc_bl_path_lbl", "bcf_regions_path_lbl", "vc_regions_path_lbl"):
                self._set_text(tag, self._fmt_name(bl_path))
        except TypeError:
            self.add_log('Invalid blacklist file', error=True)

    def callback_fasta(self, s, app_data):
        self.fasta_app_data = app_data
        try:
            fasta_path, _ = self.get_selection_path(self.fasta_app_data)
            if not fasta_path:
                return
            self.add_log('FASTA Selected: ' + fasta_path)
            # update both bcftools + variant calling labels if present
            for tag in ("bcf_ref_path_lbl", "vc_ref_path_lbl"):
                self._set_text(tag, self._fmt_name(fasta_path))
        except TypeError:
            self.add_log('Invalid FASTA file', error=True)

    def callback_bam(self, s, app_data):
        self.bam_app_data = app_data
        try:
            bam_path, _ = self.get_selection_path(self.bam_app_data)
            if not bam_path:
                return
            self.add_log('BAM Selected: ' + bam_path)
            self._set_text("sam_bam_path_lbl", self._fmt_name(bam_path))
        except Exception:
            self.add_log('Invalid BAM file', error=True)

    # variant calling callbacks
    def callback_bam_vc(self, s, app_data):
        self.bam_vc_app_data = app_data
        try:
            bam_path, _ = self.get_selection_path(self.bam_vc_app_data)
            if not bam_path:
                return
            self.add_log('VC-BAM Selected: ' + bam_path)
            self._set_text("vc_bam_path_lbl", self._fmt_name(bam_path))
        except Exception:
            self.add_log('Invalid BAM file', error=True)

    def callback_bamlist_vc(self, s, app_data):
        self.bamlist_app_data = app_data
        try:
            lst_path, _ = self.get_selection_path(self.bamlist_app_data)
            if not lst_path:
                return
            self.add_log('BAM-list Selected: ' + lst_path)
            self._set_text("vc_bamlist_path_lbl", self._fmt_name(lst_path))
        except Exception:
            self.add_log('Invalid BAM-list file', error=True)

    def get_selection_path(self, app_data):
        if not app_data:
            return None, None
        current_path = (app_data.get('current_path') or '') + '/'
        sels = app_data.get('selections') or {}
        file_path = None
        for _, value in sels.items():
            file_path = value
            break
        return file_path, current_path

    def callback_save_results(self, s, app_data):
        self.results_directory = app_data
        results_path, current_path = self.get_selection_path(self.results_directory)
        if not current_path:
            self.add_log('No directory selected.', error=True)
            return
        save_dir = self.helper.save_results(
            os.getcwd(), current_path,
            self.gwas_result_name, self.gwas_result_name_top,
            self.manhatten_plot_name, self.qq_plot_name, getattr(self, "algorithm", "Unknown"),
            self.genomic_predict_name, self.gp_plot_name, self.gp_plot_name_scatter,
            self.add_log, getattr(self, "settings_lst", []), self.pheno_stats_name, self.geno_stats_name
        )
        self.add_log('Results saved in: ' + save_dir)

    def cancel_callback_directory(self, s, app_data):
        self.add_log('Process Canceled')

    # ——— QC plotting helpers ———
    def _add_hist_plot(self, parent, title, values, bins=30, x_range=None,
                       p10=None, p50=None, p90=None):
        if not values:
            with dpg.tab(label=title, parent=parent):
                dpg.add_text("No data")
            return

        mn = min(values); mx = max(values)
        if x_range:
            mn, mx = x_range
        if mx <= mn:
            mx = mn + 1e-6

        width = (mx - mn) / bins
        counts = [0] * bins
        for v in values:
            if v < mn or v > mx:
                continue
            idx = min(bins - 1, int((v - mn) / (mx - mn) * bins))
            counts[idx] += 1
        centers = [mn + (i + 0.5) * width for i in range(bins)]

        with dpg.tab(label=title, parent=parent):
            with dpg.plot(label=title, height=380, width=720, anti_aliased=True):
                x_axis = dpg.add_plot_axis(dpg.mvXAxis, label="Value")
                y_axis = dpg.add_plot_axis(dpg.mvYAxis, label="Count")
                dpg.add_bar_series(centers, counts, parent=y_axis, label="Histogram")

                for val, lab in [(p10, "P10"), (p50, "P50"), (p90, "P90")]:
                    if val is not None:
                        try:
                            dpg.add_vline_series([float(val)], parent=y_axis, label=lab)
                        except Exception:
                            pass

    # ——— Actions ———
    def run_samtools_preprocess(self, s, data):
        try:
            bam_in = self.get_selection_path(self.bam_app_data)[0]
        except Exception:
            self.add_log('Please select a BAM file first (Preprocess samtools).', error=True)
            return

        if not bam_in:
            self.add_log('Please select a BAM file first (Preprocess samtools).', error=True)
            return

        threads = max(1, int(dpg.get_value(self.sam_threads)))
        remove_dups = bool(dpg.get_value(self.sam_remove_dups))
        compute_stats = bool(dpg.get_value(self.sam_compute_stats))
        out_prefix = (dpg.get_value(self.sam_out_prefix) or None)

        self.add_log("Running samtools preprocess...")
        try:
            outs = self.sam.preprocess(
                input_bam=bam_in,
                out_prefix=out_prefix,
                threads=threads,
                remove_dups=remove_dups,
                compute_stats=compute_stats,
                log=self.add_log,
                keep_temps=False
            )
            self.sam_out_last = outs.final_bam
            self.add_log(f"samtools preprocess: Done. Final BAM: {outs.final_bam}")
            if outs.bai:
                self.add_log(f"Index: {outs.bai}")
            for k, p in outs.stats_files.items():
                self.add_log(f"{k} report: {p}")
        except SamtoolsError as e:
            self.add_log(f"samtools error: {e}", error=True)
        except Exception as e:
            self.add_log(f"Unexpected error: {e}", error=True)

    def run_variant_calling(self, s, data):
        if not self.vcaller:
            self.add_log("variant_caller_utils not found. Please add it to the project.", error=True)
            return

        # Prefer BAM-list if provided, else single BAM
        bamlist = self._get_appdata_path_safe(self.bamlist_app_data)
        bam_single = self._get_appdata_path_safe(self.bam_vc_app_data)
        bam_spec = bamlist or bam_single
        if not bam_spec:
            self.add_log('Please select a BAM (or BAM-list) and FASTA first.', error=True)
            return

        fasta_path = self._get_appdata_path_safe(self.fasta_app_data)
        if not fasta_path:
            self.add_log('Please select a reference FASTA first.', error=True)
            return

        regions_path = self._get_appdata_path_safe(self.blacklist_app_data)
        threads = max(1, int(dpg.get_value(self.vc_threads)))
        ploidy = max(1, int(dpg.get_value(self.vc_ploidy)))
        min_bq = max(0, int(dpg.get_value(self.vc_min_bq)))
        min_mq = max(0, int(dpg.get_value(self.vc_min_mq)))
        outpfx = dpg.get_value(self.vc_out_prefix) or None

        self.add_log("Running variant calling...")
        try:
            # API expected from variant_caller_utils
            vcf_gz, tbi = self.vcaller.call_variants(
                input_path=bam_spec,
                reference_fasta=fasta_path,
                out_prefix=outpfx,
                threads=threads,
                regions_bed=regions_path,
                min_baseq=min_bq,
                min_mapq=min_mq,
                ploidy=ploidy,
                log=self.add_log
            )
            self.vc_out_last = vcf_gz

            # Make downstream tabs pick it up automatically
            self._set_virtual_selection('vcf_app_data', vcf_gz)
            name = self._fmt_name(vcf_gz)
            for tag in ("bcf_vcf_path_lbl", "qc_vcf_path_lbl", "conv_vcf_path_lbl"):
                self._set_text(tag, name)

            self.add_log(f"Variant calling: Done → {vcf_gz}")
            if tbi and os.path.exists(tbi):
                self.add_log(f"Index: {tbi}")
        except VariantCallerError as e:
            self.add_log(f"variant calling error: {e}", error=True)
        except TypeError:
            # Fallback if your function uses another signature (older/newer)
            try:
                vcf_gz, tbi = self.vcaller.call_variants(
                    bam_spec, fasta_path, out_prefix=outpfx, threads=threads,
                    regions_bed=regions_path, min_baseq=min_bq, min_mapq=min_mq,
                    ploidy=ploidy, log=self.add_log
                )
                self.vc_out_last = vcf_gz
                self._set_virtual_selection('vcf_app_data', vcf_gz)
                name = self._fmt_name(vcf_gz)
                for tag in ("bcf_vcf_path_lbl", "qc_vcf_path_lbl", "conv_vcf_path_lbl"):
                    self._set_text(tag, name)
                self.add_log(f"Variant calling: Done → {vcf_gz}")
            except Exception as e2:
                self.add_log(f"variant calling error: {e2}", error=True)
        except Exception as e:
            self.add_log(f"Unexpected error: {e}", error=True)

    def run_bcftools_preprocess(self, s, data):
        vcf_path = self._get_appdata_path_safe(self.vcf_app_data)
        if not vcf_path:
            self.add_log('Please select a VCF file first (Preprocess tab).', error=True)
            return

        fasta_path = self._get_appdata_path_safe(self.fasta_app_data)
        regions_path = self._get_appdata_path_safe(self.blacklist_app_data)

        split_m = bool(dpg.get_value(self.bcf_split))
        left_al = bool(dpg.get_value(self.bcf_left))
        do_sort = bool(dpg.get_value(self.bcf_sort))
        set_id  = bool(dpg.get_value(self.bcf_setid))
        compr   = bool(dpg.get_value(self.bcf_compr))
        index   = bool(dpg.get_value(self.bcf_index))
        rmflt   = bool(dpg.get_value(self.bcf_rmflt))
        filt    = dpg.get_value(self.bcf_filter_expr) or None
        outpfx  = dpg.get_value(self.bcf_out_prefix) or None

        self.add_log("Running bcftools preprocess...")
        try:
            final_vcf, stats = self.bcft.preprocess(
                input_vcf=vcf_path,
                out_prefix=outpfx,
                log=self.add_log,
                ref_fasta=fasta_path,
                regions_bed=regions_path,
                split_multiallelic=split_m,
                left_align=left_al,
                do_sort=do_sort,
                set_id_from_fields=set_id,
                filter_expr=filt,
                remove_filtered=rmflt,
                compress_output=compr,
                index_output=index,
                keep_temps=False
            )
            self.bcf_out_last = final_vcf

            self._set_virtual_selection('vcf_app_data', final_vcf)
            name = self._fmt_name(final_vcf)
            self._set_text("bcf_vcf_path_lbl", name)
            self._set_text("qc_vcf_path_lbl", name)
            self._set_text("conv_vcf_path_lbl", name)

            self.add_log("bcftools preprocess: Done.")
            if stats and os.path.exists(stats):
                self.add_log(f"bcftools stats saved: {stats}")
        except BCFtoolsError as e:
            self.add_log(f"bcftools error: {e}", error=True)
        except Exception as e:
            self.add_log(f"Unexpected error: {e}", error=True)

    def run_vcf_qc(self, s, data):
        try:
            vcf_path = self.get_selection_path(self.vcf_app_data)[0]
        except Exception:
            self.add_log('Please select a VCF file first.', error=True)
            return

        if not vcf_path:
            self.add_log('Please select a VCF file first.', error=True)
            return

        try:
            _ = self.get_selection_path(self.blacklist_app_data)[0]
        except Exception:
            pass

        _ = bool(dpg.get_value(self.deep_scan))
        self.add_log('Running VCF Quality Check...')

        report = self.vcf_qc_checker.evaluate(vcf_path, log_fn=self.add_log)

        dpg.configure_item("ResultsWindow", show=True)
        dpg.delete_item("ResultsWindow", children_only=True)

        dpg.add_text(f"VCF-QAScore: {report.score:.1f}   |   Verdict: {report.verdict}", parent="ResultsWindow")
        if report.hard_fail_reasons:
            dpg.add_spacer(height=4, parent="ResultsWindow")
            dpg.add_text("Hard fail reasons:", parent="ResultsWindow")
            for r in report.hard_fail_reasons:
                dpg.add_text(f"- {r}", parent="ResultsWindow")
            return

        dpg.add_spacer(height=10, parent="ResultsWindow")
        dpg.add_text("Recommendations:", parent="ResultsWindow")
        if report.recommendations:
            for r in report.recommendations:
                dpg.add_text(f"- {r}", parent="ResultsWindow")
        else:
            dpg.add_text("- No specific recommendations.", parent="ResultsWindow")

        dpg.add_spacer(height=10, parent="ResultsWindow")
        dpg.add_text("Metrics:", parent="ResultsWindow")
        with dpg.table(row_background=True, borders_innerH=True, borders_outerH=True,
                       borders_innerV=True, borders_outerV=True, parent="ResultsWindow"):
            dpg.add_table_column(label="Metric")
            dpg.add_table_column(label="Value")
            for k, v in sorted(report.metrics.items()):
                with dpg.table_row():
                    dpg.add_text(str(k))
                    dpg.add_text(f"{v:.6g}" if isinstance(v, (int, float)) else str(v))

        dpg.add_spacer(height=10, parent="ResultsWindow")
        dpg.add_text("QC Plots:", parent="ResultsWindow")

        with dpg.tab_bar(parent="ResultsWindow") as qc_tabbar:
            d = report.dists or {}
            self._add_hist_plot(
                qc_tabbar, "Depth (DP)",
                d.get("dp", []), bins=30,
                p10=report.metrics.get("dp_p10"),
                p50=report.metrics.get("dp_median"),
                p90=report.metrics.get("dp_p90"),
            )
            self._add_hist_plot(
                qc_tabbar, "Genotype Quality (GQ)",
                d.get("gq", []), bins=30,
                p10=report.metrics.get("gq_p10"),
                p50=report.metrics.get("gq_median"),
                p90=report.metrics.get("gq_p90"),
            )
            self._add_hist_plot(
                qc_tabbar, "Allele Balance |AB - 0.5| (hets)",
                d.get("ab_dev", []), bins=30, x_range=(0.0, 1.0),
                p90=report.metrics.get("ab_dev_p90"),
            )
            self._add_hist_plot(
                qc_tabbar, "Site Missingness",
                d.get("site_missing", []), bins=30, x_range=(0.0, 1.0),
                p90=report.metrics.get("site_missing_p90"),
            )

    def convert_vcf(self, s, data, user_data):
        maf = str(dpg.get_value(user_data[0]))
        geno = str(dpg.get_value(user_data[1]))
        vcf_path, _ = self.get_selection_path(self.vcf_app_data)
        variants_path = None if self.variants_app_data is None else self.get_selection_path(self.variants_app_data)[0]
        if not vcf_path:
            self.add_log('Please select a VCF file first.', error=True)
            return
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

    # ——— Results window ———
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

    # ——— Theme toggle ———
    def toggle_theme(self, sender, app_data):
        self.night_mode = bool(app_data)
        self.apply_component_themes()
        self.add_log("Dark mode enabled" if self.night_mode else "Light mode enabled")

    # ——— Utils ———
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

    def _check_cli_versions(self):
        try:
            self.add_log(self.sam.version())
        except Exception as e:
            self.add_log(f"samtools check failed: {e}", error=True)
        try:
            out = subprocess.run([resolve_tool("bcftools"), "--version"], capture_output=True, text=True)
            line0 = out.stdout.splitlines()[0] if out.stdout else "bcftools (unknown)"
            self.add_log(line0)
        except Exception as e:
            self.add_log(f"bcftools check failed: {e}", error=True)
        try:
            out = subprocess.run([resolve_tool("bgzip"), "--version"], capture_output=True, text=True)
            self.add_log(out.stdout.splitlines()[0] if out.stdout else "bgzip (unknown)")
        except Exception as e:
            self.add_log(f"bgzip check failed: {e}", error=True)
        try:
            out = subprocess.run([resolve_tool("tabix"), "--version"], capture_output=True, text=True)
            self.add_log(out.stdout.splitlines()[0] if out.stdout else "tabix (unknown)")
        except Exception as e:
            self.add_log(f"tabix check failed: {e}", error=True)

    # ——— Runner ———
    def run(self):
        self.apply_component_themes()
        dpg.setup_dearpygui()
        dpg.show_viewport()
        dpg.start_dearpygui()


if __name__ == "__main__":
    main()
