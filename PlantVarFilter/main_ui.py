# main_ui.py
import os
import shutil
import subprocess
import dearpygui.dearpygui as dpg

# optional logger
try:
    from dearpygui_ext import logger
except Exception:
    logger = None

# project backends
from vcf_quality import VCFQualityChecker
from PlantVarFilter.gwas_pipeline import GWAS
from PlantVarFilter.genomic_prediction_pipeline import GenomicPrediction
from PlantVarFilter.helpers import HELPERS
from PlantVarFilter.pipeline_plots import Plot
from pysnptools.snpreader import Bed, Pheno
import pysnptools.util as pstutil
from bcftools_utils import BCFtools, BCFtoolsError
from samtools_utils import Samtools, SamtoolsError

try:
    from variant_caller_utils import VariantCaller, VariantCallerError
except Exception:
    VariantCaller = None

# resolve tool path
try:
    from PlantVarFilter.linux import resolve_tool
except Exception:
    def resolve_tool(name: str) -> str:
        return shutil.which(name) or name

# UI modules
from ui.ui_theme import (
    setup_app_chrome,
    build_dark_theme,
    build_light_theme,
    build_component_themes,
    apply_theme,
    set_font_scale,
    set_accent_color,  # accent updater (runtime tags)
    get_primary_button_theme_tag  # runtime-aware theme getter
)

from ui.ui_pages import build_pages


class GWASApp:
    def __init__(self):
        # common state
        self._inputs = []
        self._primary_buttons = []
        self._secondary_buttons = []
        self._file_dialogs = []
        self.default_path = os.path.expanduser("~")
        self.night_mode = True

        # watermark params
        self._wm_alpha = 48
        self._wm_scale = 0.35

        # backends
        self.gwas = GWAS()
        self.helper = HELPERS()
        self.genomic_predict_class = GenomicPrediction()
        self.plot_class = Plot()
        self.vcf_qc_checker = VCFQualityChecker(max_sites_scan=200_000, min_sites_required=200)
        self.bcft = BCFtools(
            bcftools_bin=resolve_tool("bcftools"),
            bgzip_bin=resolve_tool("bgzip"),
            tabix_bin=resolve_tool("tabix")
        )
        self.sam = Samtools(exe=resolve_tool("samtools"))
        self.vcaller = VariantCaller(
            bcftools_bin=resolve_tool("bcftools"),
            bgzip_bin=resolve_tool("bgzip"),
            tabix_bin=resolve_tool("tabix")
        ) if VariantCaller else None

        # selections
        self.vcf_app_data = None
        self.variants_app_data = None
        self.results_directory = None
        self.bed_app_data = None
        self.pheno_app_data = None
        self.cov_app_data = None
        self.blacklist_app_data = None
        self.fasta_app_data = None
        self.bam_app_data = None
        self.bam_vc_app_data = None
        self.bamlist_app_data = None

        # last outputs
        self.bcf_out_last = None
        self.sam_out_last = None
        self.vc_out_last = None

        # plot/result filenames
        self.gwas_result_name = "gwas_results.csv"
        self.gwas_result_name_top = "gwas_results_top10000.csv"
        self.genomic_predict_name = "genomic_prediction_results.csv"
        self.manhatten_plot_name = "manhatten_plot.png"
        self.qq_plot_name = "qq_plot.png"
        self.gp_plot_name = "Bland_Altman_plot.png"
        self.gp_plot_name_scatter = "GP_scatter_plot.png"
        self.pheno_stats_name = 'pheno_statistics.pdf'
        self.geno_stats_name = 'geno_statistics.pdf'

        # UI handles
        self.log_win = None
        self.logz = None

        # handles created by pages
        self.sam_threads = None
        self.sam_remove_dups = None
        self.sam_compute_stats = None
        self.sam_out_prefix = None

        self.vc_threads = None
        self.vc_ploidy = None
        self.vc_min_bq = None
        self.vc_min_mq = None
        self.vc_out_prefix = None

        self.bcf_split = None
        self.bcf_left = None
        self.bcf_sort = None
        self.bcf_setid = None
        self.bcf_compr = None
        self.bcf_index = None
        self.bcf_rmflt = None
        self.bcf_filter_expr = None
        self.bcf_out_prefix = None

        self.deep_scan = None

        self.gwas_combo = None
        self.gwas_gp = None

        self.nr_jobs = None
        self.gb_goal = None
        self.plot_stats = None
        self.snp_limit = None
        self.train_size_set = None
        self.estim_set = None
        self.max_dep_set = None
        self.model_nr = None
        self.aggregation_method = None

        # navigation state
        self._nav_buttons = {}  # key -> button id
        self._pages = {}  # key -> page tag (alias)
        self._active_key = None

        # nav order and labels
        self._nav_items = [
            ("pre_sam", "Preprocess (samtools)"),
            ("vc", "Variant Calling (BAM/VCF)"),
            ("pre_bcf", "Preprocess (bcftools)"),
            ("check_vcf", "Check VCF File"),
            ("plink", "Convert to PLINK"),
            ("gwas", "GWAS Analysis"),
            ("gp", "Genomic Prediction"),
            ("settings", "Settings"),
        ]

    # ---------------- UI build ----------------
    def setup_gui(self):
        setup_app_chrome(base_size=18)
        build_dark_theme()
        build_light_theme()
        build_component_themes()
        apply_theme(dark=True)

        with dpg.window(
                tag="PrimaryWindow",
                no_title_bar=True, no_move=True, no_resize=True, no_close=True, no_collapse=True,
                pos=(0, 0)
        ):
            pass
        dpg.set_primary_window("PrimaryWindow", True)

        self._build_file_dialogs()

        with dpg.window(label="Workspace", tag="WorkspaceWindow", width=1600, height=1000, pos=(10, 10)):
            with dpg.group(horizontal=True, horizontal_spacing=12):

                # Sidebar
                with dpg.child_window(tag="Sidebar", width=240, border=True):
                    dpg.add_spacer(height=8)
                    self._build_header(parent="Sidebar")
                    dpg.add_spacer(height=6)
                    dpg.add_separator()
                    dpg.add_spacer(height=8)

                    for key, label in self._nav_items:
                        btn = dpg.add_button(
                            label=label,
                            width=-1,
                            height=36,
                            callback=self._nav_click,
                            user_data=key
                        )
                        self._nav_buttons[key] = btn
                        self._bind_nav_button_theme(btn, active=False)
                        dpg.add_spacer(height=6)

                # Content area
                with dpg.child_window(tag="content_area", width=-1, border=True):
                    self._build_header(parent="content_area", big=True)
                    dpg.add_spacer(height=6)
                    dpg.add_separator()
                    dpg.add_spacer(height=8)

                    built = build_pages(self, parent="content_area")
                    if isinstance(built, dict):
                        self._pages.update(built)

        # index pages from UI
        self._index_pages_from_ui()
        self.add_log(f"[UI] Pages discovered: {list(self._pages.keys()) or 'NONE'}", warn=not bool(self._pages))

        # settings callbacks
        settings_cbs = {
            "settings_dark_toggle": self.toggle_theme,
            "settings_font_scale": self._on_font_scale_change,
            "settings_accent_combo": self._on_accent_change,
        }
        for tag, cb in settings_cbs.items():
            if dpg.does_item_exist(tag):
                dpg.set_item_callback(tag, cb)

        self._build_tooltips()

        self.log_win = dpg.add_window(label="Log", tag="LogWindow", width=800, height=400, show=True, pos=(20, 720))
        if logger:
            self.logz = logger.mvLogger(self.log_win)
        dpg.add_window(label="Results", tag="ResultsWindow", width=1000, height=600, show=False)

        cb_map = {
            "file_dialog_vcf": self.callback_vcf,
            "file_dialog_variants": self.callback_variants,
            "file_dialog_pheno": self.callback_pheno,
            "file_dialog_cov": self.callback_cov,
            "file_dialog_bed": self.callback_bed,
            "file_dialog_blacklist": self.callback_blacklist,
            "file_dialog_fasta": self.callback_fasta,
            "file_dialog_bam": self.callback_bam,
            "file_dialog_bam_vc": self.callback_bam_vc,
            "file_dialog_bamlist": self.callback_bamlist_vc,
            "select_directory": self.callback_save_results,
        }
        for tag, fn in cb_map.items():
            if dpg.does_item_exist(tag):
                dpg.set_item_callback(tag, fn)

        # hide all pages first
        for tag in self._pages.values():
            if dpg.does_item_exist(tag):
                dpg.configure_item(tag, show=False)

        # default page
        self.show_page("pre_sam")

        self.apply_component_themes()
        self._check_cli_versions()

    # ----------- nav helpers -----------
    def _nav_click(self, sender, app_data, user_data):
        key = user_data
        self.add_log(f"[NAV] Button clicked -> '{key}'")
        self.show_page(key)

    def _get_alias(self, item_id):
        try:
            return dpg.get_item_alias(item_id)
        except Exception:
            try:
                info = dpg.get_item_info(item_id)
                return info.get("alias")
            except Exception:
                return None

    def _index_pages_from_ui(self):
        try:
            children = dpg.get_item_children("content_area", 1) or []
        except Exception:
            children = []
        found = {}
        for ch in children:
            alias = self._get_alias(ch)
            if not alias:
                continue
            if alias.startswith("page_"):
                key = alias[5:]
                found[key] = alias
        for k, v in found.items():
            self._pages.setdefault(k, v)

    # ---------------- dialogs ----------------
    def _build_file_dialogs(self):
        with dpg.file_dialog(directory_selector=False, show=False, callback=None,
                             file_count=3, tag="file_dialog_vcf", width=920, height=560,
                             default_path=self.default_path, modal=True):
            dpg.add_file_extension("Source files (*.vcf *.gz){.vcf,.gz}", color=(255, 255, 0, 255))
        self._file_dialogs.append("file_dialog_vcf")

        with dpg.file_dialog(directory_selector=False, show=False, callback=None,
                             file_count=3, tag="file_dialog_variants", width=920, height=560,
                             default_path=self.default_path, modal=True):
            dpg.add_file_extension("Text files (*.txt *.csv){.txt,.csv}", color=(255, 255, 0, 255))
            dpg.add_file_extension(".*")
        self._file_dialogs.append("file_dialog_variants")

        with dpg.file_dialog(directory_selector=False, show=False, callback=None,
                             file_count=3, tag="file_dialog_pheno", width=920, height=560,
                             default_path=self.default_path, modal=True):
            dpg.add_file_extension("Text files (*.txt *.csv){.txt,.csv}", color=(255, 255, 0, 255))
            dpg.add_file_extension(".*")
        self._file_dialogs.append("file_dialog_pheno")

        with dpg.file_dialog(directory_selector=False, show=False, callback=None,
                             file_count=3, tag="file_dialog_cov", width=920, height=560,
                             default_path=self.default_path, modal=True):
            dpg.add_file_extension("Text files (*.txt *.csv){.txt,.csv}", color=(255, 255, 0, 255))
            dpg.add_file_extension(".*")
        self._file_dialogs.append("file_dialog_cov")

        with dpg.file_dialog(directory_selector=False, show=False, callback=None,
                             file_count=3, tag="file_dialog_bed", width=920, height=560,
                             default_path=self.default_path, modal=True):
            dpg.add_file_extension(".bed", color=(255, 150, 150, 255))
            dpg.add_file_extension(".*")
        self._file_dialogs.append("file_dialog_bed")

        with dpg.file_dialog(directory_selector=False, show=False, callback=None,
                             file_count=1, tag="file_dialog_blacklist", width=920, height=560,
                             default_path=self.default_path, modal=True):
            dpg.add_file_extension(".bed", color=(255, 200, 150, 255))
            dpg.add_file_extension(".*")
        self._file_dialogs.append("file_dialog_blacklist")

        with dpg.file_dialog(directory_selector=False, show=False, callback=None,
                             file_count=1, tag="file_dialog_fasta", width=920, height=560,
                             default_path=self.default_path, modal=True):
            dpg.add_file_extension("Reference (*.fa *.fasta *.fa.gz){.fa,.fasta,.fa.gz}", color=(150, 200, 255, 255))
        self._file_dialogs.append("file_dialog_fasta")

        with dpg.file_dialog(directory_selector=False, show=False, callback=None,
                             file_count=1, tag="file_dialog_bam", width=920, height=560,
                             default_path=self.default_path, modal=True):
            dpg.add_file_extension(".bam", color=(255, 180, 120, 255))
            dpg.add_file_extension(".*")
        self._file_dialogs.append("file_dialog_bam")

        with dpg.file_dialog(directory_selector=False, show=False, callback=None,
                             file_count=1, tag="file_dialog_bam_vc", width=920, height=560,
                             default_path=self.default_path, modal=True):
            dpg.add_file_extension(".bam", color=(255, 180, 120, 255))
            dpg.add_file_extension(".*")
        self._file_dialogs.append("file_dialog_bam_vc")

        with dpg.file_dialog(directory_selector=False, show=False, callback=None,
                             file_count=1, tag="file_dialog_bamlist", width=920, height=560,
                             default_path=self.default_path, modal=True):
            dpg.add_file_extension("BAM list (*.list){.list}", color=(255, 230, 140, 255))
            dpg.add_file_extension(".*")
        self._file_dialogs.append("file_dialog_bamlist")

        dpg.add_file_dialog(directory_selector=True, show=False, callback=None,
                            tag="select_directory", cancel_callback=self.cancel_callback_directory,
                            width=900, height=540, default_path=self.default_path, modal=True)
        self._file_dialogs.append("select_directory")

    # ---------------- header ----------------
    def _build_header(self, parent, big: bool = False):
        with dpg.group(parent=parent, horizontal=True, horizontal_spacing=8):
            logo_path = os.path.join(os.path.dirname(__file__), "assets", "logo.png")
            if os.path.exists(logo_path):
                w, h, c, data = dpg.load_image(logo_path)
                if not dpg.does_item_exist("plant_logo_tex"):
                    with dpg.texture_registry(show=False):
                        dpg.add_static_texture(width=w, height=h, default_value=data, tag="plant_logo_tex")
                dpg.add_image("plant_logo_tex", width=40 if not big else 52, height=40 if not big else 52)
            else:
                dl = dpg.add_drawlist(width=40 if not big else 52, height=40 if not big else 52)
                dpg.draw_circle(center=(20 if not big else 26, 20 if not big else 26),
                                radius=18 if not big else 24,
                                color=(76, 175, 110, 255), thickness=2, parent=dl)

            dpg.add_text("PlantVarFilter", color=(210, 230, 210) if self.night_mode else (30, 45, 35))
            if big:
                dpg.add_spacer(width=10)
                dpg.add_text("GWAS, Genomic Prediction, and VCF Utilities",
                             color=(190, 175, 95) if self.night_mode else (90, 120, 70))

    # ---------------- navigation ----------------
    def _bind_nav_button_theme(self, btn, active: bool):
        try:
            dpg.bind_item_theme(btn, get_primary_button_theme_tag() if active else "theme_button_secondary")
        except Exception:
            pass

    def show_page(self, key: str):
        if not self._pages or key not in self._pages:
            self._index_pages_from_ui()
            if key not in self._pages:
                fallback = f"page_{key}"
                if dpg.does_item_exist(fallback):
                    self._pages[key] = fallback
                    self.add_log(f"[NAV] Using fallback for '{key}' -> '{fallback}'")
                else:
                    self.add_log(f"[NAV] Page '{key}' not found. Available: {list(self._pages.keys())}", error=True)
                    return

        self.add_log(f"[NAV] Switching to page '{key}'")
        for k, page_tag in self._pages.items():
            if dpg.does_item_exist(page_tag):
                dpg.configure_item(page_tag, show=(k == key))
        self._active_key = key
        for k, btn in self._nav_buttons.items():
            self._bind_nav_button_theme(btn, active=(k == key))

        # refresh watermark after page switch
        self._refresh_watermark()

    # ---------------- tooltips ----------------
    def _build_tooltips(self):
        tooltip_pairs = [
            ("tooltip_vcf", "Select a Variant Call Format file (.vcf or .vcf.gz)."),
            ("tooltip_variant", "Optional sample IDs list (PLINK --keep): FID IID (space-separated)."),
            ("tooltip_maf", "Minor Allele Frequency threshold."),
            ("tooltip_missing", "Maximum allowed missing genotype rate per variant."),
            ("tooltip_bed", "Select a PLINK .bed file (needs .bim and .fam)."),
            ("tooltip_pheno", "Phenotype file: FID IID Value (no header)."),
            ("tooltip_cov", "Covariates file: FID IID <cov1> <cov2> ..."),
            ("tooltip_algorithm", "Select the algorithm to use for analysis."),
            ("tooltip_algorithm_gp", "Select the algorithm for genomic prediction."),
            ("tooltip_training", "Percent of data used for training."),
            ("tooltip_trees", "Number of trees (RF/XGB)."),
            ("tooltip_model", "Number of models for aggregation."),
            ("tooltip_depth", "Maximum tree depth."),
            ("tooltip_nr_jobs", "CPU cores (-1 = use all cores)."),
            ("tooltip_gb_goal", "Target GB of RAM per run. 0 = block-wise reading."),
            ("tooltip_limit", "Limit SNPs in plots on huge datasets. Empty = all."),
            ("tooltip_stats", "Enable advanced PDF plots for pheno/geno stats."),
            ("tooltip_bam_sam", "Select input BAM to clean and index."),
            ("tooltip_bam_vc", "Single sample BAM for calling."),
            ("tooltip_bamlist_vc", "Text file with one BAM per line (-b list) for joint calling."),
            ("tooltip_fa_vc", "Reference FASTA (must be indexed .fai)."),
            ("tooltip_reg_vc", "Optional BED of regions to restrict calling."),
            ("tooltip_vcf_bcf", "Select input VCF/VCF.GZ to preprocess."),
            ("tooltip_fa_bcf", "Reference FASTA is required for accurate left alignment in bcftools norm."),
            ("tooltip_reg_bcf", "Regions BED to include (bcftools view -R). Optional."),
            ("tooltip_vcf_qc", "Select a VCF/VCF.GZ file to evaluate."),
            ("tooltip_bl_qc", "Optional BED of low-mappability/blacklist regions."),
        ]
        for tag, txt in tooltip_pairs:
            if dpg.does_item_exist(tag):
                with dpg.tooltip(tag):
                    dpg.add_text(txt, color=[79, 128, 90])

    # ---------------- lifecycle ----------------
    def run(self):
        dpg.create_context()
        try:
            self.setup_gui()
            dpg.create_viewport(title="PlantVarFilter", width=1600, height=1000, resizable=True)
            dpg.setup_dearpygui()
            dpg.show_viewport()

            # initial watermark + ensure refresh after first frame
            self._refresh_watermark()
            try:
                dpg.set_frame_callback(1, lambda: self._refresh_watermark())
            except Exception:
                pass

            self._hook_viewport_resize()
            self._on_viewport_resize(None, None)
            dpg.start_dearpygui()
        finally:
            dpg.destroy_context()

    def _hook_viewport_resize(self):
        if hasattr(dpg, "add_viewport_resize_handler"):
            with dpg.handler_registry():
                dpg.add_viewport_resize_handler(callback=self._on_viewport_resize)
        elif hasattr(dpg, "set_viewport_resize_callback"):
            dpg.set_viewport_resize_callback(self._on_viewport_resize)

    def _on_viewport_resize(self, sender, app_data):
        if dpg.does_item_exist("PrimaryWindow"):
            dpg.set_item_width("PrimaryWindow", dpg.get_viewport_client_width())
            dpg.set_item_height("PrimaryWindow", dpg.get_viewport_client_height())
        self._refresh_watermark()

    # ---- watermark helper ----
    def _refresh_watermark(self):
        try:
            if not dpg.does_item_exist("content_area"):
                return

            from ui.watermark import setup as setup_watermark

            # draw after one frame so the layout/rects are valid
            def _do():
                try:
                    setup_watermark(
                        alpha=self._wm_alpha,
                        scale=self._wm_scale,
                        target_window_tag="content_area",
                        front=True,  # set False if you want it behind everything
                    )
                except Exception as ex:
                    self.add_log(f"[wm] draw failed: {ex}", warn=True)

            dpg.set_frame_callback(dpg.get_frame_count() + 1, _do)

        except Exception as e:
            self.add_log(f"[wm] refresh failed: {e}", warn=True)

    # -------- settings callbacks --------
    def _on_font_scale_change(self, sender, value):
        try:
            val = float(value)
        except Exception:
            val = 1.0
        set_font_scale(val)
        self.add_log(f"[Theme] Font scale set to {val:.2f}")

    def _on_accent_change(self, sender, value):
        presets = {
            "Evergreen (Green)": ((46, 125, 50), (67, 160, 71), (27, 94, 32)),
            "Teal": ((0, 121, 107), (0, 150, 136), (0, 105, 97)),
            "Blue": ((33, 105, 170), (54, 134, 204), (22, 78, 140)),
            "Amber": ((221, 140, 20), (240, 170, 40), (190, 120, 10)),
            "Purple": ((121, 82, 179), (150, 110, 210), (98, 60, 160)),
        }
        base, hov, act = presets.get(value, presets["Evergreen (Green)"])
        set_accent_color(base, hov, act)
        self.apply_component_themes()
        self.add_log(f"[Theme] Accent set to: {value}")

    # ---------------- misc helpers ----------------
    def add_log(self, message, warn=False, error=False):
        if self.logz:
            if warn:
                self.logz.log_warning(message)
            elif error:
                self.logz.log_error(message)
            else:
                self.logz.log_info(message)
        else:
            prefix = "[WARN]" if warn else "[ERROR]" if error else "[INFO]"
            print(prefix, message)

    def _fmt_name(self, path: str) -> str:
        try:
            return os.path.basename(path) or path
        except Exception:
            return str(path)

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

    def _get_appdata_path_safe(self, app_data):
        try:
            return self.get_selection_path(app_data)[0]
        except Exception:
            return None

    def _set_text(self, tag: str, value: str):
        if dpg.does_item_exist(tag):
            dpg.set_value(tag, value)

    def _set_virtual_selection(self, attr_name: str, file_path: str):
        setattr(self, attr_name, {
            'current_path': os.path.dirname(file_path) + '/',
            'selections': {'file': file_path}
        })

    # ---------------- callbacks: file dialogs ----------------
    def callback_vcf(self, s, app_data):
        self.vcf_app_data = app_data
        vcf_path, current_path = self.get_selection_path(self.vcf_app_data)
        if not vcf_path:
            return
        # enable Convert VCF button once VCF is chosen
        if dpg.does_item_exist("convert_vcf_btn"):
            dpg.configure_item("convert_vcf_btn", enabled=True)

        dpg.configure_item("file_dialog_variants", default_path=current_path or self.default_path)
        self.add_log('VCF Selected: ' + vcf_path)
        name = self._fmt_name(vcf_path)
        for tag in ("qc_vcf_path_lbl", "conv_vcf_path_lbl", "bcf_vcf_path_lbl"):
            self._set_text(tag, name)

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
            for tag in ("gwas_bed_path_lbl", "gp_bed_path_lbl"):
                self._set_text(tag, name)
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
            for tag in ("gwas_pheno_path_lbl", "gp_pheno_path_lbl"):
                self._set_text(tag, name)
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

    # ---------------- actions ----------------
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
            self._set_virtual_selection('vcf_app_data', vcf_gz)
            name = self._fmt_name(vcf_gz)
            for tag in ("bcf_vcf_path_lbl", "qc_vcf_path_lbl", "conv_vcf_path_lbl"):
                self._set_text(tag, name)
            self.add_log(f"Variant calling: Done → {vcf_gz}")
            if tbi and os.path.exists(tbi):
                self.add_log(f"Index: {tbi}")
        except VariantCallerError as e:
            self.add_log(f"variant calling error: {e}", error=True)
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
        set_id = bool(dpg.get_value(self.bcf_setid))
        compr = bool(dpg.get_value(self.bcf_compr))
        index = bool(dpg.get_value(self.bcf_index))
        rmflt = bool(dpg.get_value(self.bcf_rmflt))
        filt = dpg.get_value(self.bcf_filter_expr) or None
        outpfx = dpg.get_value(self.bcf_out_prefix) or None

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
            for tag in ("bcf_vcf_path_lbl", "qc_vcf_path_lbl", "conv_vcf_path_lbl"):
                self._set_text(tag, name)
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
            self._add_hist_plot(qc_tabbar, "Depth (DP)", d.get("dp", []), bins=30,
                                p10=report.metrics.get("dp_p10"),
                                p50=report.metrics.get("dp_median"),
                                p90=report.metrics.get("dp_p90"))
            self._add_hist_plot(qc_tabbar, "Genotype Quality (GQ)", d.get("gq", []), bins=30,
                                p10=report.metrics.get("gq_p10"),
                                p50=report.metrics.get("gq_median"),
                                p90=report.metrics.get("gq_p90"))
            self._add_hist_plot(qc_tabbar, "Allele Balance |AB - 0.5| (hets)",
                                d.get("ab_dev", []), bins=30, x_range=(0.0, 1.0),
                                p90=report.metrics.get("ab_dev_p90"))
            self._add_hist_plot(qc_tabbar, "Site Missingness",
                                d.get("site_missing", []), bins=30, x_range=(0.0, 1.0),
                                p90=report.metrics.get("site_missing_p90"))

    def _add_hist_plot(self, parent, title, values, bins=30, x_range=None,
                       p10=None, p50=None, p90=None):
        if not values:
            with dpg.tab(label=title, parent=parent):
                dpg.add_text("No data")
            return
        mn = min(values)
        mx = max(values)
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

    def convert_vcf(self, s, data, user_data):
        """
        Robust VCF -> PLINK conversion:
        - Reads MAF/GENO from user_data dict or fallback tags.
        - Uses selected VCF/IDs from file dialogs (self.vcf_app_data / self.variants_app_data).
        - Logs command output and verifies .bed/.bim/.fam.
        """
        import os

        # 1) Read MAF / GENO safely
        try:
            if isinstance(user_data, dict):
                maf_item = user_data.get("maf", "tooltip_maf")
                geno_item = user_data.get("geno", "tooltip_missing")
            elif isinstance(user_data, (list, tuple)) and len(user_data) >= 2:
                maf_item, geno_item = user_data[0], user_data[1]
            else:
                maf_item, geno_item = "tooltip_maf", "tooltip_missing"

            maf = (str(dpg.get_value(maf_item)) or "0.01").strip()
            geno = (str(dpg.get_value(geno_item)) or "0.1").strip()
        except Exception as e:
            self.add_log(f"[WARN] Using fallback for MAF/GENO due to: {e}")
            maf, geno = "0.01", "0.1"

        # 2) Fetch paths
        vcf_path, _ = self.get_selection_path(self.vcf_app_data)
        variants_path = None if self.variants_app_data is None else self.get_selection_path(self.variants_app_data)[0]

        if not vcf_path or not os.path.exists(vcf_path):
            self.add_log('[ERROR] Please select a valid VCF file first.', error=True)
            return

        # 3) Build output prefix
        base_noext = os.path.splitext(vcf_path)[0]
        try:
            out_prefix = f"{base_noext}_maf{round(float(maf), 3)}_geno{round(float(geno), 3)}"
        except Exception:
            out_prefix = f"{base_noext}_maf{maf}_geno{geno}"

        # 4) Log inputs
        self.add_log("[INFO] Converting VCF → PLINK")
        self.add_log(f"[INFO] VCF: {vcf_path}")
        if variants_path:
            self.add_log(f"[INFO] IDs: {variants_path}")
        self.add_log(f"[INFO] MAF={maf}, GENO={geno}")
        self.add_log(f"[INFO] OUT={out_prefix}")

        # 5) Run conversion and log output
        try:
            plink_log = self.gwas.vcf_to_bed(vcf_path, variants_path, out_prefix, maf, geno)
            self.add_log(plink_log)
        except Exception as e:
            self.add_log(f"[ERROR] PLINK conversion failed: {e}", error=True)
            return

        # 6) Verify outputs and prime next step
        bed, bim, fam = f"{out_prefix}.bed", f"{out_prefix}.bim", f"{out_prefix}.fam"
        if all(os.path.exists(p) for p in (bed, bim, fam)):
            self.add_log("[OK] BED/BIM/FAM generated successfully.")
            # set selection for downstream steps
            self._set_virtual_selection('bed_app_data', bed)
            # also show filename in GWAS/GP labels if present
            for tag in ("gwas_bed_path_lbl", "gp_bed_path_lbl"):
                self._set_text(tag, self._fmt_name(bed))
        else:
            self.add_log("[ERROR] Conversion did not produce BED/BIM/FAM files.", error=True)

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
            self.settings_lst = [self.algorithm, bed_path, pheno_path, train_size_set, estimators, model_nr,
                                 max_dep_set]

            if check_input_data[0]:
                bed = Bed(str(bed_path), count_A1=False, chrom_map=chrom_mapping)
                pheno = Pheno(str(pheno_path))
                cov = Pheno(str(cov_path)) if cov_path else None

                bed, pheno = pstutil.intersect_apply([bed, pheno])
                bed_fixed = self.gwas.filter_out_missing(bed)

                self.add_log(f"Dataset after intersection: SNPs: {bed.sid_count}  Pheno IDs: {pheno.iid_count}",
                             warn=True)
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
            max_dep_set = int(max_dep_set)
        except Exception:
            max_dep_set = 0
        try:
            model_nr = int(model_nr)
        except Exception:
            model_nr = 1

        gp_df = None

        try:
            bed_path = self._get_appdata_path_safe(self.bed_app_data)
            pheno_path = self._get_appdata_path_safe(self.pheno_app_data)
            if not bed_path or not pheno_path:
                self.add_log('Please select a phenotype and genotype file.', error=True)
                return

            self.add_log('Reading files...')
            self.add_log('Validating files...')
            ok, msg = self.gwas.validate_gwas_input_files(bed_path, pheno_path)
            if not ok:
                self.add_log(msg, error=True)
                return

            chrom_mapping = self.helper.replace_with_integers(bed_path.replace('.bed', '.bim'))
            self.settings_lst = [self.algorithm, bed_path, pheno_path, test_size, estimators, model_nr, max_dep_set]

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
                self.add_log("[VAL] Running internal cross-validation (model_validation)...")
                try:
                    _ = self.genomic_predict_class.model_validation(
                        bed_fixed, pheno, bed_path, test_size, estimators,
                        self.genomic_predict_name, chrom_mapping, self.add_log,
                        model_nr, max_dep_set, validation_size=0.1
                    )
                    self.add_log("[VAL] Validation completed.")
                except Exception as e:
                    self.add_log(f"[VAL] Validation failed: {e}", error=True)
                    return

                # Show validation correlation plot(s) if produced by the backend
                possible_plots = [
                    "gp_validation_plot.png",
                    "correlation_plots.png",
                    os.path.join(os.getcwd(), "gp_validation_plot.png"),
                    os.path.join(os.getcwd(), "correlation_plots.png"),
                ]
                plot_path = next((p for p in possible_plots if os.path.exists(p)), None)
                if plot_path:
                    dpg.configure_item("ResultsWindow", show=True)
                    dpg.delete_item("ResultsWindow", children_only=True)
                    try:
                        w, h, c, data = dpg.load_image(plot_path)
                        with dpg.texture_registry(show=False):
                            dpg.add_static_texture(width=w, height=h, default_value=data, tag="gp_val_plot_tag")
                        dpg.add_text("Validation Correlation Plots", parent="ResultsWindow")
                        dpg.add_spacer(height=6, parent="ResultsWindow")
                        dpg.add_image(texture_tag="gp_val_plot_tag", parent="ResultsWindow", width=min(1200, w),
                                      height=int(h * (min(1200, w) / w)))
                    except Exception as ex:
                        self.add_log(f"[VAL] Could not load validation plot: {ex}", warn=True)
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
                self.add_log('Error, Genomic Prediction could not be started.', error=True)

        except Exception as e:
            self.add_log(f'Unexpected error in Genomic Prediction: {e}', error=True)

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

    def toggle_theme(self, sender, app_data):
        self.night_mode = bool(app_data)
        self.apply_component_themes()
        self.add_log("Dark mode enabled" if self.night_mode else "Light mode enabled")

    def apply_component_themes(self):
        apply_theme(dark=self.night_mode)

        primary_tag = get_primary_button_theme_tag()

        for b in self._primary_buttons:
            if dpg.does_item_exist(b):
                dpg.bind_item_theme(b, primary_tag)
        for b in self._secondary_buttons:
            if dpg.does_item_exist(b):
                dpg.bind_item_theme(b, "theme_button_secondary")
        for it in self._inputs:
            if dpg.does_item_exist(it):
                dpg.bind_item_theme(it, "theme_input")
        for t in self._file_dialogs:
            if dpg.does_item_exist(t):
                dpg.bind_item_theme(t, "theme_dialog")
        for k, btn in self._nav_buttons.items():
            self._bind_nav_button_theme(btn, active=(k == self._active_key))

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


def main():
    app = GWASApp()
    app.run()


if __name__ == "__main__":
    main()
