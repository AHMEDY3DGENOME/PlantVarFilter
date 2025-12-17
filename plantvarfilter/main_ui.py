# main_ui.py
import threading
import queue
import traceback
import os
import shutil
import subprocess
import time
import dearpygui.dearpygui as dpg
import sys
from pathlib import Path

# dearpygui-ext is optional; keep the GUI working even if it's missing.
try:
    from dearpygui_ext import logger  # type: ignore
except Exception:
    logger = None


# ----------------------------
# Internal imports (USE RELATIVE IMPORTS INSIDE THE PACKAGE)
# ----------------------------

# VCF quality checker (make sure plantvarfilter/vcf_quality.py exists).
try:
    from .vcf_quality import VCFQualityChecker
except ImportError:
    from plantvarfilter.vcf_quality import VCFQualityChecker

# Core pipelines and utilities (all relative to the package root).
from .gwas_pipeline import GWAS
from .genomic_prediction_pipeline import GenomicPrediction
from .helpers import HELPERS
from .pipeline_plots import Plot
from plantvarfilter.pangenome_builder import build_pangenome_graph

# pysnptools may not be available on all Py3.12 envs; keep the GUI booting anyway.
try:
    from pysnptools.snpreader import Bed, Pheno  # type: ignore
    import pysnptools.util as pstutil  # type: ignore
except Exception:
    Bed = None
    Pheno = None
    pstutil = None

from .bcftools_utils import BCFtools, BCFtoolsError
from .samtools_utils import Samtools, SamtoolsError
from .batch_gwas import run_batch_gwas_for_all_traits
from .annotation_utils import Annotator



# Variant calling is optional; don't block GUI if module is missing.
try:
    from .variant_caller_utils import VariantCaller, VariantCallerError  # type: ignore
except Exception:
    VariantCaller = None
    VariantCallerError = None

# Cross-platform resolve_tool():
# Prefer PATH (conda/bioconda tools). Fall back to platform helpers if present.
try:
    from .linux import resolve_tool as _resolve_tool_linux  # type: ignore
except Exception:
    _resolve_tool_linux = None

try:
    from .windows import resolve_tool as _resolve_tool_windows  # type: ignore
except Exception:
    _resolve_tool_windows = None



def resolve_tool(name: str) -> str:
    """
    Resolve executable path. Prefer system/conda PATH; otherwise try platform-specific helpers.
    Return the name itself if nothing found, so callers can still attempt to run it.
    """
    p = shutil.which(name)
    if p:
        return p
    if _resolve_tool_linux:
        try:
            return _resolve_tool_linux(name)
        except Exception:
            pass
    if _resolve_tool_windows:
        try:
            return _resolve_tool_windows(name)
        except Exception:
            pass
    return name

# GUI theming & pages (relative imports from the ui/ package)
from .ui.ui_theme import (
    setup_app_chrome,
    build_dark_theme,
    build_light_theme,
    build_component_themes,
    apply_theme,
    set_font_scale,
    set_accent_color,
    get_primary_button_theme_tag,
)

from .ui.ui_pages import build_pages  # (remove duplicate import if existed)

class GWASApp:
    def __init__(self):
        self._ui_events = queue.Queue()
        self._active_tasks = {}
        self._task_lock = threading.Lock()
        self._inputs = []
        self._primary_buttons = []
        self._secondary_buttons = []
        self._file_dialogs = []
        self.default_path = self._compute_default_path()
        self.night_mode = True

        self._wm_alpha = 48
        self._wm_scale = 0.35

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

        self.vcf_app_data = None
        self.vcf2_app_data = None
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
        self._chrom_mapping_last = None
        self.bcf_out_last = None
        self.sam_out_last = None
        self.vc_out_last = None
        self.gtf_app_data = None

        self._workspace_dir = os.getcwd()
        self.gwas_result_name = "gwas_results.csv"
        self.gwas_result_name_top = "gwas_results_top10000.csv"
        self.genomic_predict_name = "genomic_prediction_results.csv"
        self.manhatten_plot_name = "manhatten_plot.png"
        self.qq_plot_name = "qq_plot.png"
        self.gp_plot_name = "Bland_Altman_plot.png"
        self.gp_plot_name_scatter = "GP_scatter_plot.png"
        self.pheno_stats_name = "pheno_statistics.pdf"
        self.geno_stats_name = "geno_statistics.pdf"
        self.gwas_top_snps_name = "gwas_top_snps.csv"
        self._rebind_output_paths()

        self.logz = None
        self._log_fallback_tag = None
        self._results_body = None

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
        self.kinship_app_data = None
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

        self._nav_buttons = {}
        self._pages = {}
        self._active_key = None

        self._nav_items = [
            ("ref_manager", "Reference Manager"),
            ("pangenome", "Pangenome Builder"),
            ("fastq_qc", "FASTQ QC"),
            ("alignment", "Alignment"),
            ("pre_sam", "Preprocess (samtools)"),
            ("vc", "Variant Calling (BAM/VCF)"),
            ("pre_bcf", "Preprocess (bcftools)"),
            ("check_vcf", "Check VCF File"),
            ("plink", "Convert to PLINK"),
            ("ld", "LD Analysis"),
            ("gwas", "GWAS Analysis"),
            ("pca", "PCA / Kinship"),
            ("gp", "Genomic Prediction"),
            ("batch", "Batch GWAS"),
            ("settings", "Settings"),
        ]

        self._gui_ready = False

        self._ref_info = {"fasta": None, "out_dir": None, "prefix": None}

    @property
    def workspace_dir(self) -> str:
        return self._workspace_dir

    def _compute_default_path(self) -> str:
        env = os.environ.get("PVF_START_DIR")
        if env and os.path.exists(env):
            return env
        ws = getattr(self, "_workspace_dir", "") or os.getcwd()
        if ws and os.path.exists(ws):
            return ws
        home = os.path.expanduser("~")
        desktop = os.path.join(home, "Desktop")
        return desktop if os.path.exists(desktop) else home

    def _refresh_dialog_default_paths(self):
        self.default_path = self._compute_default_path()
        for tag in self._file_dialogs:
            if dpg.does_item_exist(tag):
                dpg.configure_item(tag, default_path=self.default_path)

    def _set_workspace_dir(self, directory: str):
        try:
            if not directory:
                directory = os.getcwd()
            os.makedirs(directory, exist_ok=True)
            self._workspace_dir = directory
            self._rebind_output_paths()
            self._refresh_dialog_default_paths()
            if getattr(self, "_gui_ready", False):
                self.add_log(f"[WS] Workspace set to: {directory}")
        except Exception as e:
            if getattr(self, "_gui_ready", False):
                self.add_log(f"[WS] Failed to set workspace dir: {e}", warn=True)

    def _rebind_output_paths(self):
        def _mk(name): return os.path.join(self._workspace_dir, name)
        self.gwas_result_name = _mk("gwas_results.csv")
        self.gwas_result_name_top = _mk("gwas_results_top10000.csv")
        self.genomic_predict_name = _mk("genomic_prediction_results.csv")
        self.manhatten_plot_name = _mk("manhatten_plot.png")
        self.qq_plot_name = _mk("qq_plot.png")
        self.gp_plot_name = _mk("Bland_Altman_plot.png")
        self.gp_plot_name_scatter = _mk("GP_scatter_plot.png")
        self.pheno_stats_name = _mk("pheno_statistics.pdf")
        self.geno_stats_name = _mk("geno_statistics.pdf")
        self.gwas_top_snps_name = _mk("gwas_top_snps.csv")

    def _safe_load_image(self, path, texture_tag):
        try:
            if not path or not os.path.exists(path):
                if getattr(self, "_gui_ready", False):
                    self.add_log(f"[IMG] Missing image: {path}", warn=True)
                return None
            w, h, c, data = dpg.load_image(path)
            if not dpg.does_item_exist(texture_tag):
                with dpg.texture_registry(show=False):
                    dpg.add_static_texture(width=w, height=h, default_value=data, tag=texture_tag)
            return (w, h)
        except Exception as e:
            if getattr(self, "_gui_ready", False):
                self.add_log(f"[IMG] Failed to load image '{path}': {e}", error=True)
            return None

    def ensure_log_window(self, show: bool = True) -> None:
        window_tag = "LogWindow"

        if not dpg.does_item_exist(window_tag):
            with dpg.window(
                    label="Log",
                    tag=window_tag,
                    width=700,
                    height=420,
                    pos=(80, 80),
                    no_saved_settings=False,
                    modal=False,
                    no_close=False,
            ):
                if logger:
                    self.logz = logger.mvLogger(parent=window_tag)
                else:
                    self._log_fallback_tag = dpg.add_input_text(
                        multiline=True,
                        readonly=True,
                        width=-1,
                        height=-1,
                        tag="log_fallback_win",
                    )

        if show:
            dpg.configure_item(window_tag, show=True)
            try:
                dpg.focus_item(window_tag)
            except Exception:
                pass

    def ensure_results_window(self, show=True, title="Results"):
        window_tag = "ResultsWindow"

        if not dpg.does_item_exist(window_tag):
            with dpg.window(
                    label=title,
                    tag=window_tag,
                    width=1100,
                    height=700,
                    pos=(300, 120),
                    no_saved_settings=False,
                    modal=False,
                    no_close=False,
            ):
                self._results_body = dpg.add_child_window(
                    tag="results_body",
                    autosize_x=True,
                    autosize_y=True,
                    border=False,
                )
        else:
            if title:
                dpg.configure_item(window_tag, label=title)
            if not dpg.does_item_exist("results_body"):
                self._results_body = dpg.add_child_window(
                    tag="results_body",
                    autosize_x=True,
                    autosize_y=True,
                    border=False,
                    parent=window_tag,
                )
            else:
                self._results_body = "results_body"

        if show:
            dpg.show_item(window_tag)
            if hasattr(dpg, "bring_item_to_front"):
                dpg.bring_item_to_front(window_tag)
            elif hasattr(dpg, "focus_item"):
                dpg.focus_item(window_tag)

    def setup_gui(self):
        setup_app_chrome(base_size=20)
        build_dark_theme()
        build_light_theme()
        build_component_themes()
        apply_theme(dark=True)

        self._font_title = None
        self._font_subtitle = None
        try:
            fonts_dir = os.path.join(os.path.dirname(__file__), "assets", "fonts")
            title_font_path = os.path.join(fonts_dir, "Inter-SemiBold.ttf")
            subtitle_font_path = os.path.join(fonts_dir, "Inter-Medium.ttf")
            with dpg.font_registry():
                if os.path.exists(title_font_path):
                    self._font_title = dpg.add_font(title_font_path, 24)
                if os.path.exists(subtitle_font_path):
                    self._font_subtitle = dpg.add_font(subtitle_font_path, 16)
        except Exception:
            self._font_title = None
            self._font_subtitle = None

        with dpg.window(
                tag="PrimaryWindow",
                no_title_bar=True, no_move=True, no_resize=True, no_close=True, no_collapse=True, pos=(0, 0)
        ):
            pass
        dpg.set_primary_window("PrimaryWindow", True)

        self._build_file_dialogs()
        self._refresh_dialog_default_paths()

        with dpg.window(label="Workspace", tag="WorkspaceWindow", width=1600, height=1000, pos=(10, 10)):
            with dpg.group(horizontal=True, horizontal_spacing=12):

                with dpg.child_window(tag="Sidebar", width=240, border=True):
                    dpg.add_spacer(height=8)
                    self._build_header(parent="Sidebar")
                    dpg.add_spacer(height=6)
                    dpg.add_separator()
                    dpg.add_spacer(height=8)

                    for key, label in self._nav_items:
                        btn = dpg.add_button(
                            label=label, width=-1, height=36,
                            callback=self._nav_click, user_data=key
                        )
                        self._nav_buttons[key] = btn
                        self._bind_nav_button_theme(btn, active=False)
                        dpg.add_spacer(height=6)

                    dpg.add_separator()
                    dpg.add_spacer(height=8)
                    dpg.add_text("Windows", color=(200, 180, 90))
                    dpg.add_spacer(height=4)
                    dpg.add_button(label="Open Log Window", width=-1,
                                   callback=lambda: self.ensure_log_window(show=True))
                    dpg.add_spacer(height=4)
                    dpg.add_button(label="Open Results Window", width=-1,
                                   callback=lambda: self.ensure_results_window(show=True))

                with dpg.child_window(tag="content_area", width=-1, border=True):
                    self._build_header(parent="content_area", big=True)
                    dpg.add_spacer(height=6)
                    dpg.add_separator()
                    dpg.add_spacer(height=8)

                    built = build_pages(self, parent="content_area")
                    if isinstance(built, dict):
                        self._pages.update(built)

        self._index_pages_from_ui()
        self.add_log(f"[UI] Pages discovered: {list(self._pages.keys()) or 'NONE'}", warn=not bool(self._pages))

        settings_cbs = {
            "settings_dark_toggle": self.toggle_theme,
            "settings_font_scale": self._on_font_scale_change,
            "settings_accent_combo": self._on_accent_change,
        }
        for tag, cb in settings_cbs.items():
            if dpg.does_item_exist(tag):
                dpg.set_item_callback(tag, cb)

        self._build_tooltips()

        cb_map = {
            "file_dialog_vcf": self.callback_vcf,
            "file_dialog_vcf2": self.callback_vcf2,
            "file_dialog_variants": self.callback_variants,
            "file_dialog_pheno": self.callback_pheno,
            "file_dialog_cov": self.callback_cov,
            "file_dialog_bed": self.callback_bed,
            "file_dialog_blacklist": self.callback_blacklist,
            "file_dialog_fasta": self.callback_fasta,
            "file_dialog_bam": self.callback_bam,
            "file_dialog_bam_vc": self.callback_bam_vc,
            "file_dialog_kinship": self.callback_kinship,
            "file_dialog_bamlist": self.callback_bamlist_vc,
            "file_dialog_gtf": self.callback_gtf,
            "select_directory": self.callback_save_results,
        }
        for tag, fn in cb_map.items():
            if dpg.does_item_exist(tag):
                dpg.set_item_callback(tag, fn)

        for tag in self._pages.values():
            if dpg.does_item_exist(tag):
                dpg.configure_item(tag, show=False)

        self.show_page("ref_manager")

        self.apply_component_themes()
        self._check_cli_versions()

        self.ensure_log_window(show=False)
        self.ensure_results_window(show=False)

        self._gui_ready = True

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

    def _build_file_dialogs(self):
        def _new_dialog(tag: str, file_count: int, exts: list, directory_selector: bool = False):
            with dpg.file_dialog(
                    directory_selector=directory_selector,
                    show=False,
                    callback=None,
                    file_count=file_count,
                    tag=tag,
                    width=980,
                    height=600,
                    default_path=self.default_path,
                    modal=True
            ):
                for label, color in exts:
                    dpg.add_file_extension(label, color=color)
            self._file_dialogs.append(tag)
            try:
                dpg.bind_item_theme(tag, "theme_dialog")
            except Exception:
                pass

        _new_dialog(
            tag="file_dialog_vcf",
            file_count=3,
            exts=[("Source files (*.vcf *.gz){.vcf,.gz}", (255, 255, 0, 255))]
        )
        _new_dialog(
            tag="file_dialog_vcf2",
            file_count=3,
            exts=[("Source files (*.vcf *.gz){.vcf,.gz}", (255, 255, 0, 255))]
        )

        _new_dialog(
            tag="file_dialog_variants",
            file_count=3,
            exts=[
                ("Text files (*.txt *.csv){.txt,.csv}", (255, 255, 0, 255)),
                (".*", (200, 200, 200, 255)),
            ]
        )
        _new_dialog(
            tag="file_dialog_pheno",
            file_count=3,
            exts=[
                ("Text files (*.txt *.csv){.txt,.csv}", (255, 255, 0, 255)),
                (".*", (200, 200, 200, 255)),
            ]
        )
        _new_dialog(
            tag="file_dialog_cov",
            file_count=3,
            exts=[
                ("Text files (*.txt *.csv){.txt,.csv}", (255, 255, 0, 255)),
                (".*", (200, 200, 200, 255)),
            ]
        )
        _new_dialog(
            tag="file_dialog_bed",
            file_count=3,
            exts=[
                (".bed", (255, 150, 150, 255)),
                (".*", (200, 200, 200, 255)),
            ]
        )
        _new_dialog(
            tag="file_dialog_blacklist",
            file_count=1,
            exts=[
                (".bed", (255, 200, 150, 255)),
                (".*", (200, 200, 200, 255)),
            ]
        )
        _new_dialog(
            tag="file_dialog_fasta",
            file_count=1,
            exts=[("Reference (*.fa *.fasta *.fa.gz){.fa,.fasta,.fa.gz}", (150, 200, 255, 255))]
        )
        _new_dialog(
            tag="file_dialog_bam",
            file_count=1,
            exts=[
                (".sam", (180, 180, 180, 255)),
                (".bam", (255, 180, 120, 255)),
                (".*", (200, 200, 200, 255)),
            ]
        )
        _new_dialog(
            tag="file_dialog_kinship",
            file_count=1,
            exts=[
                ("Kinship (*.npy *.csv *.tsv *.txt){.npy,.csv,.tsv,.txt}", (180, 220, 255, 255)),
                (".*", (200, 200, 200, 255)),
            ],
        )

        _new_dialog(
            tag="file_dialog_bam_vc",
            file_count=1,
            exts=[
                (".bam", (255, 180, 120, 255)),
                (".*", (200, 200, 200, 255)),
            ]
        )
        _new_dialog(
            tag="file_dialog_bamlist",
            file_count=1,
            exts=[
                ("BAM list (*.list){.list}", (255, 230, 140, 255)),
                (".*", (200, 200, 200, 255)),
            ]
        )
        _new_dialog(
            tag="file_dialog_gtf",
            file_count=1,
            exts=[
                ("Annotation (*.gtf *.gff *.gff3 *.gtf.gz *.gff.gz){.gtf,.gff,.gff3,.gtf.gz,.gff.gz}",
                 (180, 255, 180, 255)),
                (".*", (200, 200, 200, 255)),
            ]
        )

        with dpg.file_dialog(
                directory_selector=True,
                show=False,
                callback=None,
                tag="select_directory",
                cancel_callback=self.cancel_callback_directory,
                width=980,
                height=600,
                default_path=self.default_path,
                modal=True
        ):
            pass
        self._file_dialogs.append("select_directory")
        try:
            dpg.bind_item_theme("select_directory", "theme_dialog")
        except Exception:
            pass

    def _build_header(self, parent, big: bool = False):
        with dpg.group(parent=parent, horizontal=True, horizontal_spacing=8):
            logo_path = os.path.join(os.path.dirname(__file__), "assets", "logo.png")
            if self._safe_load_image(logo_path, "plant_logo_tex"):
                dpg.add_image("plant_logo_tex", width=40 if not big else 52, height=40 if not big else 52)
            else:
                dl = dpg.add_drawlist(width=40 if not big else 52, height=40 if not big else 52)
                dpg.draw_circle(center=(20 if not big else 26, 20 if not big else 26),
                                radius=18 if not big else 24,
                                color=(76, 175, 110, 255), thickness=2, parent=dl)

            dpg.add_text("plantvarfilter", color=(210, 230, 210) if self.night_mode else (30, 45, 35))
            if big:
                dpg.add_spacer(width=10)
                dpg.add_text("", color=(220, 200, 120) if self.night_mode else (40, 90, 40))

    def run_ld_analysis(self, s, data):
        bed_path = self._get_appdata_path_safe(self.bed_app_data)
        vcf_path = self._get_appdata_path_safe(self.vcf_app_data)
        if not bed_path and not vcf_path:
            self.add_log("Please select a PLINK BED or a VCF file.", error=True)
            return

        try:
            window_kb = int(dpg.get_value(self.ld_window_kb)) if dpg.does_item_exist(self.ld_window_kb) else 500
            window_snps = int(dpg.get_value(self.ld_window_snp)) if dpg.does_item_exist(self.ld_window_snp) else 5000
            max_kb = int(dpg.get_value(self.ld_max_kb)) if dpg.does_item_exist(self.ld_max_kb) else 1000
            min_r2 = float(dpg.get_value(self.ld_min_r2)) if dpg.does_item_exist(self.ld_min_r2) else 0.1
            region = dpg.get_value(self.ld_region).strip() if dpg.does_item_exist(self.ld_region) else ""
            do_decay = bool(dpg.get_value(self.ld_do_decay)) if dpg.does_item_exist(self.ld_do_decay) else True
            do_heat = bool(dpg.get_value(self.ld_do_heatmap)) if dpg.does_item_exist(self.ld_do_heatmap) else True
            do_div = bool(dpg.get_value(self.ld_do_div)) if dpg.does_item_exist(self.ld_do_div) else True
        except Exception as e:
            self.add_log(f"[LD] Invalid inputs: {e}", error=True)
            return

        input_for_ws = bed_path if bed_path else vcf_path
        self._set_workspace_dir(os.path.dirname(input_for_ws))

        try:
            from plantvarfilter.ld_utils import LDAnalyzer, LDAnalysisError
        except Exception as e:
            self.add_log(f"[LD] Could not import ld_utils: {e}", error=True)

        analyzer = LDAnalyzer()
        base_name = os.path.splitext(os.path.basename(input_for_ws))[0]
        out_base = os.path.join(self._workspace_dir, f"{base_name}.ld")

        results = {}
        try:
            if do_decay:
                self.add_log("[LD] Computing LD decay...")
                res = analyzer.ld_decay(
                    out_prefix=out_base + ".decay",
                    bfile_prefix=bed_path.replace(".bed", "") if bed_path and bed_path.endswith(".bed") else None,
                    vcf_path=vcf_path if vcf_path else None,
                    window_kb=window_kb,
                    window_snps=window_snps,
                    max_dist_kb=max_kb,
                    min_r2=min_r2,
                    region=region or None,
                )
                results.update(res)

            if do_heat:
                self.add_log("[LD] Building LD heatmap...")
                res = analyzer.ld_heatmap(
                    out_prefix=out_base + ".heat",
                    bfile_prefix=bed_path.replace(".bed", "") if bed_path and bed_path.endswith(".bed") else None,
                    vcf_path=vcf_path if vcf_path else None,
                    region=region or None,
                    window_snps=min(window_snps, 2000),
                    min_r2=min_r2,
                )
                results.update(res)

            if do_div:
                self.add_log("[LD] Computing diversity metrics...")
                res = analyzer.diversity(
                    out_prefix=out_base + ".div",
                    bfile_prefix=bed_path.replace(".bed", "") if bed_path and bed_path.endswith(".bed") else None,
                    vcf_path=vcf_path if vcf_path else None,
                    region=region or None,
                )
                results.update(res)
        except Exception as e:
            self.add_log(f"[LD] Analysis failed: {e}", error=True)
            return

        self.ensure_results_window(show=True, title="LD Analysis")
        dpg.delete_item(self._results_body, children_only=True)

        dpg.add_button(label="Export Results (copy files)", parent=self._results_body,
                       callback=lambda: dpg.show_item("select_directory"))
        dpg.add_spacer(height=8, parent=self._results_body)

        if "decay_png" in results and results["decay_png"] and os.path.exists(results["decay_png"]):
            tag = "ld_decay_img"
            wh = self._safe_load_image(results["decay_png"], tag)
            if wh:
                w, h = wh
                dpg.add_text("LD decay", parent=self._results_body)
                dpg.add_image(tag, parent=self._results_body, width=min(900, w), height=int(h * min(900, w) / w))
                dpg.add_text(f"CSV: {results.get('decay_csv', '')}", parent=self._results_body)
                dpg.add_spacer(height=6, parent=self._results_body)

        if "heat_png" in results and results["heat_png"] and os.path.exists(results["heat_png"]):
            tag = "ld_heat_img"
            wh = self._safe_load_image(results["heat_png"], tag)
            if wh:
                w, h = wh
                dpg.add_text("LD heatmap", parent=self._results_body)
                dpg.add_image(tag, parent=self._results_body, width=min(900, w), height=int(h * min(900, w) / w))
                dpg.add_spacer(height=6, parent=self._results_body)

        if "summary_csv" in results and os.path.exists(results["summary_csv"]):
            import pandas as pd
            try:
                df = pd.read_csv(results["summary_csv"])
                dpg.add_text("Diversity metrics", parent=self._results_body)
                with dpg.table(parent=self._results_body, row_background=True,
                               borders_innerH=True, borders_outerH=True,
                               borders_innerV=True, borders_outerV=True):
                    dpg.add_table_column(label="Metric")
                    dpg.add_table_column(label="Value")
                    for _, r in df.iterrows():
                        with dpg.table_row():
                            dpg.add_text(str(r["metric"]))
                            try:
                                dpg.add_text(f"{float(r['value']):.6g}")
                            except Exception:
                                dpg.add_text(str(r["value"]))
            except Exception:
                dpg.add_text(f"Summary: {results['summary_csv']}", parent=self._results_body)

        self.add_log("[LD] Done.")

    def compute_kinship_from_bed(self, s=None, a=None):
        import os
        import numpy as np
        import pandas as pd
        from pysnptools.snpreader import Bed

        bed_path = self._get_appdata_path_safe(self.bed_app_data)
        if not bed_path or not os.path.exists(bed_path):
            self.add_log("Please select a valid BED file first.", error=True)
            return

        try:
            self._set_workspace_dir(os.path.dirname(bed_path))
            bim_path = bed_path.replace(".bed", ".bim")
            chrom_mapping = self.helper.replace_with_integers(bim_path)

            self.add_log("[KIN] Loading genotype matrix...")
            bed = Bed(str(bed_path), count_A1=False, chrom_map=chrom_mapping)
            X = bed.read().val.astype(float)
            iid = bed.iid

            self.add_log("[KIN] Computing allele frequencies...")
            p = np.nanmean(X, axis=0) / 2.0
            valid = np.isfinite(p) & (p > 0.0) & (p < 1.0)
            if not np.any(valid):
                self.add_log("[KIN] No informative SNPs for kinship (all monomorphic/missing).", error=True)
                return

            Xv = X[:, valid]
            pv = p[valid]
            self.add_log("[KIN] Standardizing genotypes...")
            Xv_centered = Xv.copy()
            nan_mask = ~np.isfinite(Xv_centered)
            if np.any(nan_mask):
                Xv_centered[nan_mask] = np.take(2.0 * pv, np.where(nan_mask)[1])
            Xv_centered -= (2.0 * pv)

            denom = np.sqrt(2.0 * pv * (1.0 - pv))
            Z = Xv_centered / denom

            self.add_log("[KIN] Building GRM...")
            M = Z.shape[1]
            G = (Z @ Z.T) / float(M)

            sample_ids = [f"{fid}_{iid_}" for fid, iid_ in iid.tolist()]
            df_grm = pd.DataFrame(G, index=sample_ids, columns=sample_ids)

            base = os.path.splitext(os.path.basename(bed_path))[0]
            out_csv = os.path.join(self._workspace_dir, f"{base}_kinship.csv")
            df_grm.to_csv(out_csv, index=True)

            self.kinship_path_last = out_csv
            try:
                if hasattr(self, "_set_virtual_selection"):
                    self._set_virtual_selection('kinship_app_data', out_csv)
            except Exception:
                pass

            if dpg.does_item_exist("gwas_kinship_path_lbl"):
                dpg.set_value("gwas_kinship_path_lbl", os.path.basename(out_csv))

            self.add_log(f"[KIN] Kinship matrix saved: {out_csv}")
        except Exception as e:
            self.add_log(f"[KIN] Failed to compute kinship: {e}", error=True)

    def run_pca_module(self, s=None, a=None):
        bed_path = self._get_appdata_path_safe(self.bed_app_data)
        if not bed_path or not os.path.exists(bed_path):
            self.add_log("[PCA] Please select a PLINK .bed file first.", error=True)
            return

        self._set_workspace_dir(os.path.dirname(bed_path))

        if bed_path.lower().endswith(".bed"):
            bfile_prefix = bed_path[:-4]
        else:
            bfile_prefix = os.path.splitext(bed_path)[0]

        try:
            npcs = int(dpg.get_value(self.pca_npcs)) if (
                    self.pca_npcs and dpg.does_item_exist(self.pca_npcs)
            ) else 10
        except Exception:
            npcs = 10

        do_kin = bool(dpg.get_value(self.pca_kinship)) if (
                self.pca_kinship and dpg.does_item_exist(self.pca_kinship)
        ) else True

        outpfx = ""
        if self.pca_out_prefix and dpg.does_item_exist(self.pca_out_prefix):
            outpfx = (dpg.get_value(self.pca_out_prefix) or "").strip()
        if not outpfx:
            base = os.path.basename(bfile_prefix)
            outpfx = os.path.join(self._workspace_dir, f"{base}.pca")

        plink2 = resolve_tool("plink2")
        plink19 = resolve_tool("plink")
        use_plink2 = (shutil.which(plink2) is not None) if plink2 else False
        use_plink19 = (shutil.which(plink19) is not None) if plink19 else False
        if not use_plink2 and not use_plink19:
            self.add_log("[PCA] Neither plink2 nor plink (1.9) is available on PATH.", error=True)
            return

        try:
            if use_plink2:
                cmd = [
                    plink2,
                    "--bfile",
                    bfile_prefix,
                    "--pca",
                    str(npcs),
                    "biallelic-var-wts",
                    "approx",
                    "--out",
                    outpfx,
                ]
            else:
                cmd = [plink19, "--bfile", bfile_prefix, "--pca", str(npcs), "--out", outpfx]

            self._run_cmd(cmd)

            eigvec = f"{outpfx}.eigenvec"
            eigval = f"{outpfx}.eigenval"
            if not (os.path.exists(eigvec) and os.path.exists(eigval)):
                self.add_log("[PCA] PCA output files not found (.eigenvec/.eigenval).", error=True)
                return

            kin_out = None
            if do_kin:
                kin_out = outpfx + ".rel"
                if use_plink2:
                    cmd = [plink2, "--bfile", bfile_prefix, "--make-rel", "square", "--out", outpfx]
                    self._run_cmd(cmd)
                    if not os.path.exists(kin_out) and os.path.exists(kin_out + ".gz"):
                        kin_out = kin_out + ".gz"
                elif use_plink19:
                    cmd = [plink19, "--bfile", bfile_prefix, "--make-rel", "square", "--out", outpfx]
                    self._run_cmd(cmd)
                    if not os.path.exists(kin_out) and os.path.exists(kin_out + ".gz"):
                        kin_out = kin_out + ".gz"
                else:
                    kin_out = None

            self.ensure_results_window(show=True, title="PCA / Kinship Results")
            dpg.delete_item(self._results_body, children_only=True)

            dpg.add_text(f"Output prefix: {outpfx}", parent=self._results_body)
            dpg.add_spacer(height=4, parent=self._results_body)
            dpg.add_text(
                f"PCA files: {os.path.basename(eigvec)}, {os.path.basename(eigval)}",
                parent=self._results_body,
            )
            if kin_out:
                dpg.add_text(f"Kinship file: {os.path.basename(kin_out)}", parent=self._results_body)

            self._render_pca_plots(eigvec, eigval)

            dpg.add_spacer(height=10, parent=self._results_body)
            dpg.add_button(
                label="Export Results",
                parent=self._results_body,
                callback=lambda: dpg.show_item("select_directory"),
            )

            self.add_log("[PCA] Done.")
        except Exception as e:
            self.add_log(f"[PCA] Failed: {e}", error=True)

    def _render_pca_plots(self, eigvec_path: str, eigval_path: str) -> None:
        try:
            xs_pc1_pc2_x, xs_pc1_pc2_y = [], []
            xs_pc1_pc3_x, xs_pc1_pc3_y = [], []

            # Read eigenvectors
            with open(eigvec_path, "r", encoding="utf-8") as f:
                for line in f:
                    parts = line.strip().split()
                    if len(parts) < 4:
                        continue
                    try:
                        pc1 = float(parts[2])
                        pc2 = float(parts[3])
                    except ValueError:
                        continue
                    xs_pc1_pc2_x.append(pc1)
                    xs_pc1_pc2_y.append(pc2)

                    if len(parts) >= 5:
                        try:
                            pc3 = float(parts[4])
                        except ValueError:
                            continue
                        xs_pc1_pc3_x.append(pc1)
                        xs_pc1_pc3_y.append(pc3)

            # Read eigenvalues
            eigenvalues = []
            try:
                with open(eigval_path, "r", encoding="utf-8") as f:
                    for line in f:
                        line = line.strip()
                        if not line:
                            continue
                        try:
                            eigenvalues.append(float(line.split()[0]))
                        except ValueError:
                            continue
            except FileNotFoundError:
                eigenvalues = []

            parent = self._results_body
            dpg.add_spacer(height=10, parent=parent)

            with dpg.tab_bar(parent=parent):
                if xs_pc1_pc2_x and xs_pc1_pc2_y:
                    with dpg.tab(label="PC1 vs PC2"):
                        with dpg.plot(height=420, width=720):
                            x_axis = dpg.add_plot_axis(dpg.mvXAxis, label="PC1")
                            y_axis = dpg.add_plot_axis(dpg.mvYAxis, label="PC2")
                            dpg.add_scatter_series(xs_pc1_pc2_x, xs_pc1_pc2_y, label="samples", parent=y_axis)

                if xs_pc1_pc3_x and xs_pc1_pc3_y:
                    with dpg.tab(label="PC1 vs PC3"):
                        with dpg.plot(height=420, width=720):
                            x_axis = dpg.add_plot_axis(dpg.mvXAxis, label="PC1")
                            y_axis = dpg.add_plot_axis(dpg.mvYAxis, label="PC3")
                            dpg.add_scatter_series(xs_pc1_pc3_x, xs_pc1_pc3_y, label="samples", parent=y_axis)

                if eigenvalues:
                    with dpg.tab(label="Scree plot"):
                        pcs_idx = list(range(1, len(eigenvalues) + 1))
                        with dpg.plot(height=420, width=720):
                            x_axis = dpg.add_plot_axis(dpg.mvXAxis, label="PC index")
                            y_axis = dpg.add_plot_axis(dpg.mvYAxis, label="Eigenvalue")
                            dpg.add_line_series(pcs_idx, eigenvalues, label="eigenvalues", parent=y_axis)
        except Exception as ex:
            self.add_log(f"[PCA] Could not render PCA plots: {ex}", warn=True)

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
        self._refresh_watermark()

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

    def run(self):
        from pathlib import Path
        import dearpygui.dearpygui as dpg

        dpg.create_context()
        try:
            self.setup_gui()

            icon_path = Path(__file__).resolve().parent / "assets" / "app_icon.png"

            kwargs = dict(
                title="PlantVarFilter",
                width=1400,
                height=900,
                resizable=True,
            )

            if icon_path.exists():
                kwargs["small_icon"] = str(icon_path)
                kwargs["large_icon"] = str(icon_path)

            dpg.create_viewport(**kwargs)
            dpg.setup_dearpygui()
            dpg.show_viewport()

            self._refresh_watermark()
            try:
                dpg.set_frame_callback(1, lambda: self._refresh_watermark())
            except Exception:
                pass

            if hasattr(self, "_install_ui_poller"):
                try:
                    self._install_ui_poller()
                except Exception:
                    pass

            dpg.start_dearpygui()
        finally:
            dpg.destroy_context()

    def ui_emit(self, fn, *args, **kwargs):
        try:
            self._ui_events.put((fn, args, kwargs))
        except Exception:
            pass

    def _ui_poll(self):
        drained = 0
        while drained < 200:
            try:
                fn, args, kwargs = self._ui_events.get_nowait()
            except Exception:
                break
            try:
                fn(*args, **kwargs)
            except Exception:
                pass
            drained += 1

    def run_task(self, key: str, target, on_done=None, on_error=None):
        with self._task_lock:
            t = self._active_tasks.get(key)
            if t and t.is_alive():
                self.add_log(f"[TASK] '{key}' is already running.", warn=True)
                return

        def _runner():
            try:
                res = target()
                if on_done:
                    self.ui_emit(on_done, res)
            except Exception as e:
                if on_error:
                    self.ui_emit(on_error, e)
                else:
                    tb = traceback.format_exc()
                    self.ui_emit(self.add_log, f"[TASK] '{key}' failed: {e}\n{tb}", True, False, True)

        th = threading.Thread(target=_runner, daemon=True)
        with self._task_lock:
            self._active_tasks[key] = th
        th.start()

    def _install_ui_poller(self):
        if hasattr(dpg, "set_render_callback"):
            dpg.set_render_callback(lambda: self._ui_poll())
            return

        def _tick():
            self._ui_poll()
            try:
                dpg.set_frame_callback(1, lambda: _tick())
            except Exception:
                pass

        try:
            dpg.set_frame_callback(1, lambda: _tick())
        except Exception:
            pass


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

    def _refresh_watermark(self):
        try:
            if not dpg.does_item_exist("content_area"):
                return

            from .ui.watermark import setup as setup_watermark, place_signature as setup_lab_signature

            def _place():
                try:
                    setup_watermark(
                        alpha=self._wm_alpha,
                        scale=self._wm_scale,
                        target_window_tag="content_area",
                        front=True,
                    )
                    setup_lab_signature(
                        target_window_tag="content_area",
                        image_name="logo_lab.png",
                        width=240,
                        margin=(16, 16),
                    )
                except Exception as ex:
                    self.add_log(f"[wm] draw failed: {ex}", warn=True)

                dpg.set_frame_callback(dpg.get_frame_count() + 10, _place)

            dpg.set_frame_callback(dpg.get_frame_count() + 1, _place)

        except Exception as e:
            self.add_log(f"[wm] refresh failed: {e}", warn=True)

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

    def add_log(self, message, warn=False, error=False):
        self.ensure_log_window(show=False)
        if self.logz:
            if warn:
                self.logz.log_warning(message)
            elif error:
                self.logz.log_error(message)
            else:
                self.logz.log_info(message)
        else:
            prefix = "[WARN]" if warn else "[ERROR]" if error else "[INFO]"
            txt = f"{prefix} {message}\n"
            if self._log_fallback_tag and dpg.does_item_exist(self._log_fallback_tag):
                old = dpg.get_value(self._log_fallback_tag) or ""
                dpg.set_value(self._log_fallback_tag, old + txt)
            else:
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

    def _val(self, *tags):
        for t in tags:
            if dpg.does_item_exist(t):
                v = dpg.get_value(t)
                if isinstance(v, str) and v.strip():
                    return v.strip()
        return ""

    def _run_cmd(self, cmd, stdout_path=None):
        self.add_log("$ " + " ".join(cmd))
        if stdout_path:
            with open(stdout_path, "w", encoding="utf-8") as fh:
                p = subprocess.run(cmd, text=True, stdout=fh, stderr=subprocess.PIPE)
            if p.stderr:
                for ln in p.stderr.splitlines():
                    self.add_log(ln, warn=True)
        else:
            p = subprocess.run(cmd, text=True, capture_output=True)
            if p.stdout:
                for ln in p.stdout.splitlines():
                    self.add_log(ln)
            if p.stderr:
                for ln in p.stderr.splitlines():
                    self.add_log(ln, warn=True)
        if p.returncode != 0:
            raise RuntimeError(f"command failed: {' '.join(cmd)}")
        return p.returncode

    def callback_vcf(self, s, app_data):
        self.vcf_app_data = app_data
        vcf_path, current_path = self.get_selection_path(self.vcf_app_data)
        if not vcf_path:
            return
        self._set_workspace_dir(os.path.dirname(vcf_path))
        dpg.configure_item("file_dialog_variants", default_path=current_path or self.default_path)
        self.add_log('VCF Selected: ' + vcf_path)
        name = self._fmt_name(vcf_path)
        for tag in ("qc_vcf_path_lbl", "conv_vcf_path_lbl", "bcf_vcf_path_lbl"):
            self._set_text(tag, name)

    def callback_vcf2(self, s, app_data):
        self.vcf2_app_data = app_data
        vcf_path, current_path = self.get_selection_path(self.vcf2_app_data)
        if not vcf_path:
            return
        self._set_workspace_dir(os.path.dirname(vcf_path))
        dpg.configure_item("file_dialog_variants", default_path=current_path or self.default_path)
        self.add_log('VCF#2 Selected: ' + vcf_path)
        self._set_text("qc_vcf2_path_lbl", self._fmt_name(vcf_path))

    def callback_bed(self, s, app_data):
        self.bed_app_data = app_data
        try:
            bed_path, current_path = self.get_selection_path(self.bed_app_data)
            if not bed_path:
                return

            self._set_workspace_dir(os.path.dirname(bed_path))
            dpg.configure_item("file_dialog_cov", default_path=current_path or self.default_path)
            dpg.configure_item("file_dialog_pheno", default_path=current_path or self.default_path)
            self.add_log("BED file Selected: " + bed_path)

            name = self._fmt_name(bed_path)

            # GWAS / GP labels
            for tag in ("gwas_bed_path_lbl", "gp_bed_path_lbl"):
                self._set_text(tag, name)

            # LD page
            if dpg.does_item_exist("ld_btn_bed"):
                dpg.configure_item("ld_btn_bed", label=name)
            if dpg.does_item_exist("ld_bed_path_lbl"):
                self._set_text("ld_bed_path_lbl", name)

            # PCA page
            if dpg.does_item_exist("pca_btn_bed"):
                dpg.configure_item("pca_btn_bed", label=name)
            if dpg.does_item_exist("pca_bed_path_lbl"):
                self._set_text("pca_bed_path_lbl", name)

        except TypeError:
            self.add_log("Invalid BED file Selected", error=True)

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

    def callback_gtf(self, s, app_data):
        self.gtf_app_data = app_data
        try:
            gtf_path, current_path = self.get_selection_path(self.gtf_app_data)
            if not gtf_path:
                return
            self._set_workspace_dir(os.path.dirname(gtf_path))
            self.add_log('Annotation (GTF/GFF) Selected:' + gtf_path)
            if dpg.does_item_exist("gwas_gtf_path_lbl"):
                self._set_text("gwas_gtf_path_lbl", self._fmt_name(gtf_path))
        except Exception:
            self.add_log('Invalid GTF/GFF annotation file', error=True)

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
            if dpg.does_item_exist("ref_fasta_path_inp"):
                dpg.set_value("ref_fasta_path_inp", fasta_path)
        except TypeError:
            self.add_log('Invalid FASTA file', error=True)

    def callback_bam(self, s, app_data):
        self.bam_app_data = app_data
        try:
            path, _ = self.get_selection_path(self.bam_app_data)
            if not path:
                return
            self._set_workspace_dir(os.path.dirname(path))
            self.add_log('SAM/BAM Selected: ' + path)
            self._set_text("sam_bam_path_lbl", self._fmt_name(path))
        except Exception:
            self.add_log('Invalid SAM/BAM file', error=True)

    def callback_bam_vc(self, s, app_data):
        self.bam_vc_app_data = app_data
        try:
            bam_path, _ = self.get_selection_path(self.bam_vc_app_data)
            if not bam_path:
                return
            self._set_workspace_dir(os.path.dirname(bam_path))
            self.add_log('VC-BAM Selected: ' + bam_path)
            self._set_text("vc_bam_path_lbl", self._fmt_name(bam_path))
            if dpg.does_item_exist(self.vc_out_prefix) and not (dpg.get_value(self.vc_out_prefix) or "").strip():
                base = os.path.splitext(os.path.basename(bam_path))[0]
                dpg.set_value(self.vc_out_prefix, os.path.join(self._workspace_dir, f"{base}.raw"))
        except Exception:
            self.add_log('Invalid BAM file', error=True)

    def callback_bamlist_vc(self, s, app_data):
        self.bamlist_app_data = app_data
        try:
            lst_path, _ = self.get_selection_path(self.bamlist_app_data)
            if not lst_path:
                return
            self._set_workspace_dir(os.path.dirname(lst_path))
            self.add_log('BAM-list Selected: ' + lst_path)
            self._set_text("vc_bamlist_path_lbl", self._fmt_name(lst_path))
            if dpg.does_item_exist(self.vc_out_prefix) and not (dpg.get_value(self.vc_out_prefix) or "").strip():
                base = os.path.splitext(os.path.basename(lst_path))[0]
                dpg.set_value(self.vc_out_prefix, os.path.join(self._workspace_dir, f"{base}.raw"))
        except Exception:
            self.add_log('Invalid BAM-list file', error=True)

    def callback_kinship(self, s, app_data):
        self.kinship_app_data = app_data
        kin_path, _ = self.get_selection_path(self.kinship_app_data)
        if not kin_path:
            return
        self._set_workspace_dir(os.path.dirname(kin_path))
        self.add_log("Kinship file selected: " + kin_path)
        self._set_text("gwas_kinship_path_lbl", self._fmt_name(kin_path))

    def callback_save_results(self, s, app_data):
        sel_path, cur_dir = self.get_selection_path(app_data)
        base = sel_path or cur_dir or self._workspace_dir
        if not base:
            self.add_log('No directory selected.', error=True)
            return
        if not os.path.isdir(base):
            base = os.path.dirname(base)
        base = os.path.abspath(base)
        save_dir = self.helper.save_results(
            current_dir=self._workspace_dir,
            save_dir=base,
            gwas_result_name=self.gwas_result_name,
            gwas_result_name_top=self.gwas_result_name_top,
            manhatten_plot_name=self.manhatten_plot_name,
            qq_plot_name=self.qq_plot_name,
            algorithm=getattr(self, "algorithm", "Unknown"),
            genomic_predict_name=self.genomic_predict_name,
            gp_plot_name=self.gp_plot_name,
            gp_plot_name_scatter=self.gp_plot_name_scatter,
            add_log=self.add_log,
            settings_lst=getattr(self, "settings_lst", []),
            pheno_stats_name=self.pheno_stats_name,
            geno_stats_name=self.geno_stats_name
        )
        self.add_log('Results saved in: ' + save_dir)

    def cancel_callback_directory(self, s, app_data):
        self.add_log('Process Canceled')

    def ref_build_indexes(self, s=None, a=None):
        fasta = self._val("ref_fasta_path_inp", "ref_fasta_path")
        out_dir = self._val("ref_out_dir_inp", "ref_out_dir") or (os.path.dirname(fasta) if fasta else "")
        if not fasta or not os.path.exists(fasta):
            self.add_log("[REF] Please provide a valid FASTA", error=True)
            return
        if not out_dir:
            out_dir = os.path.dirname(fasta)
        os.makedirs(out_dir, exist_ok=True)
        self._set_workspace_dir(out_dir)

        base = os.path.splitext(os.path.basename(fasta))[0]
        bt2_prefix = os.path.join(out_dir, base)
        mmi_path = os.path.join(out_dir, base + ".mmi")
        dict_path = os.path.join(out_dir, base + ".dict")

        try:
            self.add_log("[REF] Building samtools faidx")
            self._run_cmd([resolve_tool("samtools"), "faidx", fasta])

            self.add_log("[REF] Building sequence dictionary")
            self._run_cmd([resolve_tool("samtools"), "dict", fasta, "-o", dict_path])

            self.add_log("[REF] Building minimap2 index (.mmi)")
            self._run_cmd(["minimap2", "-d", mmi_path, fasta])

            self.add_log("[REF] Building bowtie2 index")
            self._run_cmd(["bowtie2-build", fasta, bt2_prefix])

            self._ref_info = {"fasta": fasta, "out_dir": out_dir, "prefix": os.path.join(out_dir, base)}
            self.add_log("[REF] Indexing completed.")
            if dpg.does_item_exist("ref_status"):
                dpg.set_value("ref_status", f"OK\nFASTA: {fasta}\nOUT: {out_dir}")
        except Exception as e:
            self.add_log(f"[REF] Failed: {e}", error=True)

    def run_fastq_qc(self, s=None, a=None):
        r1 = self._val("fq1_inp", "fq1_path_inp", "fq1_path")
        r2 = self._val("fq2_inp", "fq2_path_inp", "fq2_path")
        out_dir = self._val("fq_out_dir_inp", "fq_out_dir") or self.workspace_dir
        if not r1:
            self.add_log("[QC] Please select FASTQ #1", error=True)
            return
        os.makedirs(out_dir, exist_ok=True)
        self._set_workspace_dir(out_dir)

        try:
            self.add_log("[QC] Running FastQC")
            cmd = ["fastqc", "-o", out_dir, "--threads", "4", r1]
            if r2:
                cmd.append(r2)
            self._run_cmd(cmd)

            try:
                self.add_log("[QC] Running MultiQC")
                self._run_cmd(["multiqc", out_dir, "-o", out_dir])
            except FileNotFoundError:
                self.add_log("[QC] multiqc not found; skipping.", warn=True)

            msg = f"Reports in: {out_dir}"
            if dpg.does_item_exist("qc_status"):
                dpg.set_value("qc_status", msg)
            self.add_log("[QC] Done. " + msg)
        except Exception as e:
            self.add_log(f"[QC] Failed: {e}", error=True)

    def run_alignment(self, s=None, a=None):
        platform = (dpg.get_value("align_platform") if dpg.does_item_exist("align_platform") else "illumina").lower()
        r1 = self._val("aln_r1_inp", "aln_r1_path_inp", "aln_r1_path")
        r2 = self._val("aln_r2_inp", "aln_r2_path_inp", "aln_r2_path")
        out_dir = self._val("aln_out_dir_inp", "aln_out_dir") or (self._ref_info.get("out_dir") or self.workspace_dir)

        if not r1:
            self.add_log("[ALN] Please select reads (R1)", error=True)
            return
        if not self._ref_info.get("prefix"):
            self.add_log("[ALN] Please build reference indexes first.", error=True)
            return

        os.makedirs(out_dir, exist_ok=True)
        self._set_workspace_dir(out_dir)

        base = os.path.splitext(os.path.basename(r1))[0]
        sam_out = os.path.join(out_dir, f"{base}.sam")

        try:
            if platform in ("nanopore", "ont", "pacbio", "pb", "long", "longreads"):
                mmi = self._ref_info["prefix"] + ".mmi"
                if not os.path.exists(mmi):
                    self.add_log("[ALN] Missing .mmi index. Build Reference Manager first.", error=True)
                    return
                preset = "map-ont" if platform in ("nanopore", "ont") else "map-pb"
                self.add_log(f"[ALN] minimap2 ({preset}) → {sam_out}")
                cmd = ["minimap2", "-t", "4", "-ax", preset, mmi, r1]
                if r2:
                    cmd.append(r2)
                self._run_cmd(cmd, stdout_path=sam_out)
            else:
                bt2 = self._ref_info["prefix"]
                self.add_log(f"[ALN] bowtie2 → {sam_out}")
                if r2:
                    cmd = ["bowtie2", "-x", bt2, "-1", r1, "-2", r2, "-S", sam_out, "-p", "4"]
                else:
                    cmd = ["bowtie2", "-x", bt2, "-U", r1, "-S", sam_out, "-p", "4"]
                self._run_cmd(cmd)

            self.add_log(f"[ALN] Done. SAM: {sam_out}")
            if dpg.does_item_exist("aln_status"):
                dpg.set_value("aln_status", f"SAM: {sam_out}")
        except Exception as e:
            self.add_log(f"[ALN] Failed: {e}", error=True)

    def run_samtools_preprocess(self, s=None, data=None):
        try:
            in_path = self.get_selection_path(self.bam_app_data)[0]
        except Exception:
            self.add_log('Please select a SAM/BAM file first (Preprocess samtools).', error=True)
            return

        if not in_path or not os.path.exists(in_path):
            self.add_log('Please select a SAM/BAM file first (Preprocess samtools).', error=True)
            return

        self._set_workspace_dir(os.path.dirname(in_path))

        threads = max(1, int(dpg.get_value(self.sam_threads)))
        remove_dups = bool(dpg.get_value(self.sam_remove_dups))
        compute_stats = bool(dpg.get_value(self.sam_compute_stats))
        out_prefix = (dpg.get_value(self.sam_out_prefix) or None)

        ext = os.path.splitext(in_path)[1].lower()
        bam_for_pipeline = in_path
        tmp_created = False

        try:
            if ext == ".sam":
                base = os.path.splitext(os.path.basename(in_path))[0]
                bam_for_pipeline = os.path.join(self._workspace_dir, f"{base}.unsorted.bam")
                self.add_log(f"Converting SAM -> BAM: {bam_for_pipeline}")
                cmd = [
                    resolve_tool("samtools"),
                    "view",
                    "-@", str(threads),
                    "-bS",
                    in_path,
                    "-o",
                    bam_for_pipeline,
                ]
                p = subprocess.run(cmd, text=True, capture_output=True)
                if p.stdout:
                    for ln in p.stdout.splitlines():
                        self.add_log(ln)
                if p.stderr:
                    for ln in p.stderr.splitlines():
                        self.add_log(ln, warn=True)
                if p.returncode != 0 or (not os.path.exists(bam_for_pipeline)):
                    raise RuntimeError("SAM to BAM conversion failed")
                tmp_created = True

            self.add_log("Running samtools preprocess...")
            outs = self.sam.preprocess(
                bam_for_pipeline,  # أول برامتر *positional* فقط
                out_prefix=out_prefix,
                threads=threads,
                remove_dups=remove_dups,
                compute_stats=compute_stats,
                log=self.add_log,
                keep_temps=False,
            )

            self.sam_out_last = outs.final_bam
            self.add_log(f"samtools preprocess: Done. Final BAM: {outs.final_bam}")
            if outs.bai:
                self.add_log(f"Index: {outs.bai}")
            for k, pth in outs.stats_files.items():
                self.add_log(f"{k} report: {pth}")

        except SamtoolsError as e:
            self.add_log(f"samtools error: {e}", error=True)
        except Exception as e:
            self.add_log(f"Unexpected error: {e}", error=True)
        finally:
            if tmp_created:
                try:
                    os.remove(bam_for_pipeline)
                except Exception:
                    pass

    def run_variant_calling(self, s=None, data=None):
        # Lazy-init VariantCaller
        if getattr(self, "vcaller", None) is None:
            try:
                from variant_caller_utils import VariantCaller

                bcftools_bin = resolve_tool("bcftools") if resolve_tool else None
                bgzip_bin = resolve_tool("bgzip") if resolve_tool else None
                tabix_bin = resolve_tool("tabix") if resolve_tool else None

                self.vcaller = VariantCaller(
                    bcftools_bin=bcftools_bin,
                    bgzip_bin=bgzip_bin,
                    tabix_bin=tabix_bin,
                )
                self.add_log("[VC] VariantCaller initialized.")
            except Exception as e:
                self.vcaller = None
                self.add_log(f"[VC] Could not initialize VariantCaller: {e}", error=True)
                return

        # Resolve BAM / BAM-list
        bamlist = self._get_appdata_path_safe(self.bamlist_app_data)
        bam_single = self._get_appdata_path_safe(self.bam_vc_app_data)

        if bamlist:
            bam_spec = bamlist
            bams_arg = bamlist  # treated as BAM-list file by VariantCaller
        elif bam_single:
            bam_spec = bam_single
            bams_arg = [bam_single]  # single BAM
        else:
            self.add_log("Please select a BAM (or BAM-list) and FASTA first.", error=True)
            return

        # Resolve reference FASTA
        fasta_path = self._get_appdata_path_safe(self.fasta_app_data)
        if not fasta_path:
            self.add_log("Please select a reference FASTA first.", error=True)
            return

        # Workspace and indexes
        self._set_workspace_dir(os.path.dirname(bam_spec))

        if bam_single and os.path.isfile(bam_single) and not os.path.exists(bam_single + ".bai"):
            self.add_log("[VC] .bai not found; indexing BAM...")
            subprocess.run(
                [resolve_tool("samtools"), "index", bam_single],
                check=False,
                text=True,
                capture_output=True,
            )

        if not os.path.exists(fasta_path + ".fai"):
            self.add_log("[VC] FASTA .fai not found; building faidx...")
            subprocess.run(
                [resolve_tool("samtools"), "faidx", fasta_path],
                check=False,
                text=True,
                capture_output=True,
            )

        regions_path = self._get_appdata_path_safe(self.blacklist_app_data)

        threads = max(1, int(dpg.get_value(self.vc_threads)))
        ploidy = max(1, int(dpg.get_value(self.vc_ploidy)))
        min_bq = max(0, int(dpg.get_value(self.vc_min_bq)))
        min_mq = max(0, int(dpg.get_value(self.vc_min_mq)))
        outpfx = dpg.get_value(self.vc_out_prefix) or None

        if not outpfx:
            base = os.path.basename(bam_spec)
            base = os.path.splitext(base)[0]
            outpfx = os.path.join(self._workspace_dir, f"{base}.raw")
            dpg.set_value(self.vc_out_prefix, outpfx)

        split_after = bool(dpg.get_value("vc_split_after_calling")) if dpg.does_item_exist(
            "vc_split_after_calling"
        ) else False

        self.add_log("Running variant calling...")
        try:
            vcf_gz, tbi = self.vcaller.call_bcftools(
                bams=bams_arg,
                ref_fasta=fasta_path,
                out_prefix=outpfx,
                regions_bed=regions_path,
                threads=threads,
                min_baseq=min_bq,
                min_mapq=min_mq,
                ploidy=ploidy,
                log=self.add_log,
                split_after_calling=split_after,
            )

            self.vc_out_last = vcf_gz
            self._set_virtual_selection("vcf_app_data", vcf_gz)
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

    def run_bcftools_preprocess(self, s=None, data=None):
        vcf_path = self._get_appdata_path_safe(self.vcf_app_data)
        if not vcf_path:
            self.add_log('Please select a VCF file first (Preprocess tab).', error=True)
            return

        self._set_workspace_dir(os.path.dirname(vcf_path))

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

        filltags = bool(dpg.get_value(self.bcf_filltags)) if dpg.does_item_exist(self.bcf_filltags) else False
        make_snps = bool(dpg.get_value(self.bcf_make_snps)) if dpg.does_item_exist(self.bcf_make_snps) else False
        make_svs = bool(dpg.get_value(self.bcf_make_svs)) if dpg.does_item_exist(self.bcf_make_svs) else False

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
                keep_temps=False,
                fill_tags=filltags,
            )

            self.bcf_out_last = final_vcf
            self._set_virtual_selection('vcf_app_data', final_vcf)
            name = self._fmt_name(final_vcf)
            for tag in ("bcf_vcf_path_lbl", "qc_vcf_path_lbl", "conv_vcf_path_lbl"):
                self._set_text(tag, name)

            self.add_log("bcftools preprocess: Done.")
            if stats and os.path.exists(stats):
                self.add_log(f"bcftools stats saved: {stats}")

            stem = final_vcf
            if stem.endswith(".vcf.gz"):
                stem = stem[:-7]
            elif stem.endswith(".vcf"):
                stem = stem[:-4]

            if make_snps:
                snps_out = f"{stem}.snps.vcf.gz"
                cmd = [resolve_tool("bcftools"), "view", "-v", "snps", final_vcf, "-Oz", "-o", snps_out]
                self.add_log("$ " + " ".join(cmd))
                subprocess.run(cmd, check=False, text=True, capture_output=True)
                subprocess.run([resolve_tool("tabix"), "-f", "-p", "vcf", snps_out],
                               check=False, text=True, capture_output=True)
                self.add_log(f"SNP-only: {snps_out}")

            if make_svs:
                sv_out = f"{stem}.sv.vcf.gz"
                cmd = [resolve_tool("bcftools"), "view", "-V", "snps,indels", final_vcf, "-Oz", "-o", sv_out]
                self.add_log("$ " + " ".join(cmd))
                subprocess.run(cmd, check=False, text=True, capture_output=True)
                subprocess.run([resolve_tool("tabix"), "-f", "-p", "vcf", sv_out],
                               check=False, text=True, capture_output=True)
                self.add_log(f"SV-only: {sv_out}")

        except BCFtoolsError as e:
            self.add_log(f"bcftools error: {e}", error=True)
        except Exception as e:
            self.add_log(f"Unexpected error: {e}", error=True)

    def run_vcf_qc(self, s=None, data=None):
        try:
            vcf1 = self.get_selection_path(self.vcf_app_data)[0]
        except Exception:
            vcf1 = None
        try:
            vcf2 = self.get_selection_path(self.vcf2_app_data)[0]
        except Exception:
            vcf2 = None

        if not vcf1 and not vcf2:
            self.add_log('Please select a VCF file first.', error=True)
            return

        vcf_main = vcf1 or vcf2
        self._set_workspace_dir(os.path.dirname(vcf_main))

        self.add_log('Running VCF Quality Check...')
        reports = []
        names = []

        if vcf1:
            r1 = self.vcf_qc_checker.evaluate(vcf1, log_fn=self.add_log)
            reports.append(r1)
            names.append(self._fmt_name(vcf1))

        if vcf2:
            r2 = self.vcf_qc_checker.evaluate(vcf2, log_fn=self.add_log)
            reports.append(r2)
            names.append(self._fmt_name(vcf2))

        self.ensure_results_window(show=True, title="VCF QC Results")
        dpg.delete_item(self._results_body, children_only=True)

        with dpg.group(parent=self._results_body, horizontal=True, horizontal_spacing=8):
            dpg.add_button(label="Export Results (copy files)", callback=lambda: dpg.show_item("select_directory"))
            if vcf1:
                dpg.add_button(label="Save QC Report VCF#1", callback=lambda: self._save_qc_report(vcf1, reports[0]))
            if vcf2:
                idx = 1 if vcf1 else 0
                dpg.add_button(label="Save QC Report VCF#2", callback=lambda: self._save_qc_report(vcf2, reports[idx]))

        dpg.add_spacer(height=10, parent=self._results_body)

        def _render_one(tab_parent, name, report):
            dpg.add_text(f"VCF: {name}", parent=tab_parent)
            dpg.add_text(f"VCF-QAScore: {report.score:.1f}   |   Verdict: {report.verdict}", parent=tab_parent)
            dpg.add_text(f"Samples: {int(report.metrics.get('samples', 0))}", parent=tab_parent)

            try:
                dt = getattr(report, "data_type", None)
                if not dt:
                    s_val = float(report.metrics.get("snps", 0.0) or 0.0)
                    i_val = float(report.metrics.get("indels", 0.0) or 0.0)
                    dt = "SNPs" if (s_val > 0 and i_val == 0) else "SVs"
                dpg.add_text(f"Data type: {dt}", parent=tab_parent)
            except Exception:
                dpg.add_text("Data type: SVs", parent=tab_parent)

            if report.hard_fail_reasons:
                dpg.add_spacer(height=4, parent=tab_parent)
                dpg.add_text("Hard fail reasons:", parent=tab_parent)
                for r in report.hard_fail_reasons:
                    dpg.add_text(f"- {r}", parent=tab_parent)
                return

            dpg.add_spacer(height=10, parent=tab_parent)
            dpg.add_text("Recommendations:", parent=tab_parent)
            if report.recommendations:
                for r in report.recommendations:
                    dpg.add_text(f"- {r}", parent=tab_parent)
            else:
                dpg.add_text("- No specific recommendations.", parent=tab_parent)

            dpg.add_spacer(height=10, parent=tab_parent)
            dpg.add_text("Metrics:", parent=tab_parent)
            with dpg.table(
                    row_background=True,
                    borders_innerH=True,
                    borders_outerH=True,
                    borders_innerV=True,
                    borders_outerV=True,
                    parent=tab_parent,
            ):
                dpg.add_table_column(label="Metric")
                dpg.add_table_column(label="Value")
                for k, v in sorted(report.metrics.items()):
                    with dpg.table_row():
                        dpg.add_text(str(k))
                        dpg.add_text(f"{v:.6g}" if isinstance(v, (int, float)) else str(v))

            dpg.add_spacer(height=10, parent=tab_parent)
            dpg.add_text("QC Plots:", parent=tab_parent)
            with dpg.tab_bar(parent=tab_parent) as qc_tabbar:
                d = report.dists or {}
                self._add_hist_plot(
                    qc_tabbar,
                    "Depth (DP)",
                    d.get("dp", []),
                    bins=30,
                    p10=report.metrics.get("dp_p10"),
                    p50=report.metrics.get("dp_median"),
                    p90=report.metrics.get("dp_p90"),
                )
                self._add_hist_plot(
                    qc_tabbar,
                    "Genotype Quality (GQ)",
                    d.get("gq", []),
                    bins=30,
                    p10=report.metrics.get("gq_p10"),
                    p50=report.metrics.get("gq_median"),
                    p90=report.metrics.get("gq_p90"),
                )
                self._add_hist_plot(
                    qc_tabbar,
                    "Allele Balance |AB - 0.5| (hets)",
                    d.get("ab_dev", []),
                    bins=30,
                    x_range=(0.0, 1.0),
                    p90=report.metrics.get("ab_dev_p90"),
                )
                self._add_hist_plot(
                    qc_tabbar,
                    "Site Missingness",
                    d.get("site_missing", []),
                    bins=30,
                    x_range=(0.0, 1.0),
                    p90=report.metrics.get("site_missing_p90"),
                )

        with dpg.tab_bar(parent=self._results_body):
            for name, rep in zip(names, reports):
                with dpg.tab(label=name):
                    _render_one(tab_parent=dpg.last_item(), name=name, report=rep)

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
            # skip values outside the plotting range
            if v < mn or v > mx:
                continue

            idx = min(bins - 1, int((v - mn) / (mx - mn) * bins))
            counts[idx] += 1

        centers = [mn + (i + 0.5) * width for i in range(bins)]

        with dpg.tab(label=title, parent=parent):
            with dpg.plot(label=title, height=380, width=720) as plot_tag:
                x_axis = dpg.add_plot_axis(dpg.mvXAxis, label="Value")
                y_axis = dpg.add_plot_axis(dpg.mvYAxis, label="Count")

                dpg.add_bar_series(centers, counts, parent=y_axis, label="Histogram")

                for val, lab in [(p10, "P10"), (p50, "P50"), (p90, "P90")]:
                    if val is not None:
                        try:
                            dpg.add_drag_line(
                                parent=plot_tag,
                                default_value=float(val),
                                vertical=True,
                                label=lab,
                            )
                        except Exception:
                            pass

    def convert_vcf(self, sender=None, app_data=None, user_data=None):
        # Try to show spinner and disable button
        try:
            dpg.configure_item("plink_convert_spinner", show=True)
            dpg.configure_item("convert_vcf_btn", enabled=False)
        except Exception:
            pass

        try:
            # 1) Get VCF path from the existing file dialog (vcf_app_data is already used elsewhere)
            vcf_path = self._get_appdata_path_safe(getattr(self, "vcf_app_data", None))
            if not vcf_path:
                self.add_log("[PLINK] Please select a VCF file first.", error=True)
                return

            # keep workspace next to VCF
            self._set_workspace_dir(os.path.dirname(vcf_path))

            # 2) optional IDs file from the variants dialog
            ids_path = self._get_appdata_path_safe(getattr(self, "variants_app_data", None))

            # base name used for both normalized VCF and PLINK output
            base = os.path.splitext(os.path.basename(vcf_path))[0]

            # 3) normalize VCF to split multi-allelic variants (bcftools norm -m -any)
            vcf_to_use = vcf_path
            bcftools_exe = shutil.which("bcftools")
            if bcftools_exe:
                normalized_vcf = os.path.join(self.workspace_dir, f"{base}.split.vcf.gz")
                norm_cmd = [
                    bcftools_exe,
                    "norm",
                    "-m", "-any",
                    vcf_path,
                    "-Oz",
                    "-o", normalized_vcf,
                ]
                index_cmd = [bcftools_exe, "index", normalized_vcf]

                self.add_log(f"[PLINK] Normalizing VCF with bcftools: {' '.join(norm_cmd)}")
                proc_norm = subprocess.run(
                    norm_cmd,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    text=True,
                )
                self.add_log(proc_norm.stdout)

                if proc_norm.returncode != 0:
                    raise RuntimeError(f"bcftools norm failed with exit code {proc_norm.returncode}")

                proc_index = subprocess.run(
                    index_cmd,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    text=True,
                )
                self.add_log(proc_index.stdout)

                if proc_index.returncode != 0:
                    raise RuntimeError(f"bcftools index failed with exit code {proc_index.returncode}")

                vcf_to_use = normalized_vcf
                self.add_log(f"[PLINK] Using normalized VCF: {vcf_to_use}")
            else:
                self.add_log(
                    "[PLINK] 'bcftools' not found in PATH; using original VCF "
                    "(multiallelic variants may cause PLINK to fail).",
                    error=True,
                )

            # 4) read filters from UI by tag (defined in page_convert_plink)
            try:
                maf = float(dpg.get_value("plink_maf_input") or 0.05)
            except Exception:
                maf = 0.05

            try:
                geno = float(dpg.get_value("plink_missing_input") or 0.10)
            except Exception:
                geno = 0.10

            mind = geno  # reuse geno for mind if there is no separate input

            # 5) choose output prefix next to the VCF
            out_prefix = os.path.join(self.workspace_dir, f"{base}.plink")

            # 6) locate plink/plink2 in PATH
            plink_exe = shutil.which("plink2")
            if not plink_exe:
                self.add_log(
                    "[PLINK] Could not find plink2 in PATH. "
                    "Please install PLINK 2 and make sure the 'plink2' command is available.",
                    error=True,
                )
                return

            # 7) build PLINK command
            cmd = [
                plink_exe,
                "--vcf", vcf_to_use,
                "--make-bed",
                "--maf", str(maf),
                "--geno", str(geno),
                "--mind", str(mind),
                "--out", out_prefix,
            ]
            if ids_path:
                cmd.extend(["--extract", ids_path])

            self.add_log(f"[PLINK] Running: {' '.join(cmd)}")

            # 8) run PLINK and log its output
            proc = subprocess.run(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
            )
            self.add_log(proc.stdout)

            if proc.returncode != 0:
                raise RuntimeError(f"PLINK failed with exit code {proc.returncode}")

            self.add_log(f"[PLINK] Conversion finished. Output prefix: {out_prefix}")
            self.add_log(f"[PLINK] Expected files: {out_prefix}.bed / .bim / .fam")

        except Exception as e:
            self.add_log(f"[PLINK] Error during conversion: {e}", error=True)

        finally:
            try:
                dpg.configure_item("plink_convert_spinner", show=False)
                dpg.configure_item("convert_vcf_btn", enabled=True)
            except Exception:
                pass

    def _build_top_snps_file(self):
        try:
            import pandas as pd
            if not os.path.exists(self.gwas_result_name):
                return
            df = pd.read_csv(self.gwas_result_name)
            if 'PValue' in df.columns:
                out = df.dropna(subset=['PValue']).sort_values('PValue', ascending=True).head(100)
            elif 'SNP effect' in df.columns:
                out = df.dropna(subset=['SNP effect']).sort_values('SNP effect', ascending=False).head(100)
            else:
                out = df.head(100)
            out.to_csv(self.gwas_top_snps_name, index=False)
            self.add_log(f"[TopSNPs] Saved: {self.gwas_top_snps_name}")
        except Exception as e:
            self.add_log(f"[TopSNPs] Failed to build: {e}", warn=True)

    def run_gwas(self, s, data, user_data):
        self.delete_files()
        try:
            train_size_set = (100 - dpg.get_value(self.train_size_set)) / 100
        except Exception:
            train_size_set = 0.7
        try:
            estimators = dpg.get_value(self.estim_set)
        except Exception:
            estimators = 200
        try:
            model_nr = dpg.get_value(self.model_nr)
        except Exception:
            model_nr = 1
        try:
            snp_limit = dpg.get_value(self.snp_limit)
        except Exception:
            snp_limit = None
        try:
            nr_jobs = int(dpg.get_value(self.nr_jobs)) or -1
        except Exception:
            nr_jobs = -1
        try:
            gb_goal = int(dpg.get_value(self.gb_goal))
        except Exception:
            gb_goal = 0
        try:
            max_dep_set = dpg.get_value(self.max_dep_set)
        except Exception:
            max_dep_set = 3
        try:
            self.algorithm = dpg.get_value(self.gwas_combo)
        except Exception:
            self.algorithm = "FaST-LMM"
        try:
            aggregation_method = str(dpg.get_value(self.aggregation_method))
        except Exception:
            aggregation_method = "sum"

        try:
            self.add_log("Reading files...")
            bed_path = self.get_selection_path(self.bed_app_data)[0]
            pheno_path = self.get_selection_path(self.pheno_app_data)[0]

            try:
                kin_path = self._get_appdata_path_safe(getattr(self, "kinship_app_data", None))
            except Exception:
                kin_path = None

            if not bed_path or not pheno_path:
                self.add_log("Please select a phenotype and genotype file.", error=True)
                return

            self._set_workspace_dir(os.path.dirname(bed_path))

            cov_path = None
            try:
                cov_path = self.get_selection_path(self.cov_app_data)[0]
            except Exception:
                pass

            self.add_log("Validating files...")
            check_input_data = self.gwas.validate_gwas_input_files(bed_path, pheno_path)

            chrom_mapping = self.helper.replace_with_integers(bed_path.replace(".bed", ".bim"))
            self._chrom_mapping_last = chrom_mapping
            self.settings_lst = [
                self.algorithm,
                bed_path,
                pheno_path,
                train_size_set,
                estimators,
                model_nr,
                max_dep_set,
            ]

            if not check_input_data[0]:
                self.add_log(check_input_data[1], error=True)
                return

            bed = Bed(str(bed_path), count_A1=False, chrom_map=chrom_mapping)
            pheno = Pheno(str(pheno_path))
            cov = Pheno(str(cov_path)) if cov_path else None

            bed, pheno = pstutil.intersect_apply([bed, pheno])
            bed_fixed = self.gwas.filter_out_missing(bed)

            try:
                import pandas as pd
                import numpy as np
                region_chr = dpg.get_value(self.region_chr).strip() if dpg.does_item_exist(self.region_chr) else ""
                region_start = int(dpg.get_value(self.region_start)) if dpg.does_item_exist(self.region_start) else 0
                region_end = int(dpg.get_value(self.region_end)) if dpg.does_item_exist(self.region_end) else 0
                if region_chr and region_end > region_start >= 0:
                    bim_path = bed_path.replace(".bed", ".bim")
                    df_bim = pd.read_csv(bim_path, sep=r"\s+", header=None, engine="python")
                    df_bim.columns = ["Chr", "SNP", "NA1", "ChrPos", "NA2", "NA3"]
                    chr_vals = df_bim["Chr"].astype(str)
                    mask_chr = (chr_vals == str(region_chr)) | (
                            chr_vals == str(self.helper.replace_with_integers_value(region_chr)
                                            if hasattr(self.helper, "replace_with_integers_value") else region_chr)
                    )
                    mask_pos = (df_bim["ChrPos"].astype(int) >= region_start) & (
                                df_bim["ChrPos"].astype(int) <= region_end)
                    df_sub = df_bim[mask_chr & mask_pos]
                    if not df_sub.empty:
                        snp_keep = set(df_sub["SNP"].astype(str).tolist())
                        snp_mask = np.array([sid in snp_keep for sid in bed.sid])
                        if snp_mask.any():
                            prev = bed_fixed.sid_count
                            bed_fixed = bed_fixed[:, snp_mask]
                            self.add_log(
                                f"[Region] Filter applied: {region_chr}:{region_start}-{region_end} → {bed_fixed.sid_count}/{prev} SNPs kept")
                        else:
                            self.add_log(
                                f"[Region] No SNPs matched {region_chr}:{region_start}-{region_end}. Proceeding without region filter.",
                                warn=True)
                    else:
                        self.add_log(
                            f"[Region] No BIM records matched {region_chr}:{region_start}-{region_end}. Proceeding without region filter.",
                            warn=True)
            except Exception as ex_region:
                self.add_log(f"[Region] Failed to apply region filter: {ex_region}", warn=True)

            self.add_log(f"Dataset after intersection: SNPs: {bed_fixed.sid_count}  Pheno IDs: {pheno.iid_count}",
                         warn=True)
            self.add_log("Starting Analysis, this might take a while...")

            if self.algorithm in ("FaST-LMM", "Linear regression"):
                gwas_df, df_plot = self.gwas.run_gwas_lmm(
                    bed_fixed=bed_fixed,
                    pheno=pheno,
                    chrom_mapping=chrom_mapping,
                    add_log=self.add_log,
                    gwas_result_name=self.gwas_result_name,
                    algorithm=self.algorithm,
                    bed_file=bed_path,
                    cov_file=cov,
                    gb_goal=gb_goal,
                    kinship_path=kin_path,
                )
            elif self.algorithm == "Random Forest (AI)":
                gwas_df, df_plot = self.gwas.run_gwas_rf(
                    bed_fixed,
                    pheno,
                    bed_path,
                    train_size_set,
                    estimators,
                    self.gwas_result_name,
                    chrom_mapping,
                    self.add_log,
                    model_nr,
                    nr_jobs,
                    aggregation_method,
                )
            elif self.algorithm == "XGBoost (AI)":
                gwas_df, df_plot = self.gwas.run_gwas_xg(
                    bed_fixed,
                    pheno,
                    bed_path,
                    train_size_set,
                    estimators,
                    self.gwas_result_name,
                    chrom_mapping,
                    self.add_log,
                    model_nr,
                    max_dep_set,
                    nr_jobs,
                    aggregation_method,
                )
            elif self.algorithm == "Ridge Regression":
                gwas_df, df_plot = self.gwas.run_gwas_ridge(
                    bed_fixed,
                    pheno,
                    bed_path,
                    train_size_set,
                    1.0,
                    self.gwas_result_name,
                    chrom_mapping,
                    self.add_log,
                    model_nr,
                    aggregation_method,
                )
            else:
                gwas_df, df_plot = None, None

            if gwas_df is None:
                self.add_log("Error, GWAS Analysis could not be started.", error=True)
                return

            self.add_log("GWAS Analysis done.")
            self.add_log("GWAS Results Plotting...")

            try:
                do_stats = bool(dpg.get_value(self.plot_stats))
            except Exception:
                do_stats = False
            if do_stats:
                self.plot_class.plot_pheno_statistics(pheno_path, self.pheno_stats_name)
                self.plot_class.plot_geno_statistics(bed_fixed, pheno, self.geno_stats_name)

            self.gwas.plot_gwas(
                df_plot,
                snp_limit if (str(snp_limit).strip().isdigit() and int(snp_limit) > 0) else None,
                self.algorithm,
                self.manhatten_plot_name,
                self.qq_plot_name,
                chrom_mapping,
            )

            self._build_top_snps_file()

            try:
                annotate_on = bool(dpg.get_value(self.annotate_enable)) if dpg.does_item_exist(
                    self.annotate_enable) else False
                window_kb = int(dpg.get_value(self.annotate_window_kb)) if dpg.does_item_exist(
                    self.annotate_window_kb) else 50
            except Exception:
                annotate_on, window_kb = False, 50

            if annotate_on:
                gtf_path = None
                try:
                    gtf_path = self._get_appdata_path_safe(getattr(self, "gtf_app_data", None))
                except Exception:
                    gtf_path = None
                if gtf_path and os.path.exists(self.gwas_result_name):
                    try:
                        from annotation_utils import Annotator
                        ann = Annotator()
                        ann.load_gtf_or_gff(gtf_path)
                        ann.build_index()
                        annotated_csv = ann.annotate_and_save(
                            gwas_df,
                            out_csv=self.gwas_result_name,
                            window_kb=window_kb,
                            chr_col="Chr",
                            pos_col="ChrPos",
                        )
                        self.add_log(f"[ANNOT] GWAS annotated → {annotated_csv}")
                    except Exception as e:
                        self.add_log(f"[ANNOT] Annotation failed: {e}", warn=True)

            self.add_log("Done...")
            self.show_results_window(gwas_df, self.algorithm, genomic_predict=False)
            self.bed_app_data = None
            self.pheno_app_data = None
            self.cov_app_data = None

        except TypeError as e:
            self.add_log(f"TypeError in GWAS: {e}", error=True)
            return
        except Exception as e:
            self.add_log(f"Unexpected error in GWAS: {e}", error=True)
            return

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

            self._set_workspace_dir(os.path.dirname(bed_path))

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

                for p in [
                    "gp_validation_plot.png",
                    "correlation_plots.png",
                    os.path.join(os.getcwd(), "gp_validation_plot.png"),
                    os.path.join(os.getcwd(), "correlation_plots.png"),
                ]:
                    if os.path.exists(p):
                        self._show_validation_plot(p)
                        break
                return

            if gp_df is not None:
                self.add_log('Genomic Prediction done.')
                self.add_log('Genomic Prediction Plotting...')
                if self._safe_load_image(self.gp_plot_name, "ba_tag"):
                    pass
                if self._safe_load_image(self.gp_plot_name_scatter, "ba_tag2"):
                    pass
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

    def run_batch_gwas_ui(self, s, data):
        try:
            bed_path = self._get_appdata_path_safe(self.bed_app_data)
            pheno_path = self._get_appdata_path_safe(self.pheno_app_data)
            cov_path = self._get_appdata_path_safe(self.cov_app_data)
            if not bed_path or not pheno_path:
                self.add_log("Please select BED and phenotype file first (Batch GWAS).", error=True)
                return

            self._set_workspace_dir(os.path.dirname(bed_path))

            algorithm = dpg.get_value(self.batch_algo) if self.batch_algo and dpg.does_item_exist(
                self.batch_algo) else "FaST-LMM"

            nr_jobs = int(dpg.get_value(self.nr_jobs)) if self.nr_jobs and dpg.does_item_exist(self.nr_jobs) else -1
            gb_goal = int(dpg.get_value(self.gb_goal)) if self.gb_goal and dpg.does_item_exist(self.gb_goal) else 0
            train_size = float(
                dpg.get_value(self.train_size_set) or 70) / 100.0 if self.train_size_set and dpg.does_item_exist(
                self.train_size_set) else 0.7
            estimators = int(dpg.get_value(self.estim_set) or 200) if self.estim_set and dpg.does_item_exist(
                self.estim_set) else 200
            max_depth = int(dpg.get_value(self.max_dep_set) or 3) if self.max_dep_set and dpg.does_item_exist(
                self.max_dep_set) else 3
            model_nr = int(dpg.get_value(self.model_nr) or 1) if self.model_nr and dpg.does_item_exist(
                self.model_nr) else 1
            aggregation_method = dpg.get_value(
                self.aggregation_method) if self.aggregation_method and dpg.does_item_exist(
                self.aggregation_method) else "sum"
            snp_limit = dpg.get_value(self.snp_limit) if self.snp_limit and dpg.does_item_exist(self.snp_limit) else ""
            snp_limit = int(snp_limit) if str(snp_limit).strip().isdigit() else None

            out_dir = self._workspace_dir

            self.add_log(f"[Batch-GWAS] Starting on: {os.path.basename(pheno_path)} using {algorithm}")
            result = run_batch_gwas_for_all_traits(
                gwas=self.gwas, helper=self.helper,
                bed_path=bed_path, pheno_path=pheno_path, cov_path=cov_path, algorithm=algorithm,
                out_dir=out_dir, log_fn=self.add_log, nr_jobs=nr_jobs, gb_goal=gb_goal,
                train_size=train_size, estimators=estimators, model_nr=model_nr,
                max_depth=max_depth, aggregation_method=aggregation_method, snp_limit=snp_limit
            )
            summary_csv = result["summary_csv"]

            self.ensure_results_window(show=True, title="Batch GWAS – Summary")
            dpg.delete_item(self._results_body, children_only=True)
            dpg.add_button(label="Export Results", parent=self._results_body,
                           callback=lambda: dpg.show_item("select_directory"))
            dpg.add_spacer(height=10, parent=self._results_body)

            import pandas as pd
            try:
                df = pd.read_csv(summary_csv)
                with dpg.table(row_background=True, borders_innerH=True, borders_innerV=True,
                               borders_outerH=True, borders_outerV=True, parent=self._results_body):
                    for c in df.columns:
                        dpg.add_table_column(label=str(c))
                    for _, r in df.iterrows():
                        with dpg.table_row():
                            for c in df.columns:
                                dpg.add_text(str(r[c]))
            except Exception as e:
                self.add_log(f"[Batch-GWAS] Could not render summary table: {e}", warn=True)
                dpg.add_text(f"Summary saved here: {summary_csv}", parent=self._results_body)

            self.add_log(f"[Batch-GWAS] Done. Summary: {summary_csv}")
        except Exception as e:
            self.add_log(f"[Batch-GWAS] Error: {e}", error=True)

    def _show_validation_plot(self, plot_path: str):
        self.ensure_results_window(show=True, title="Validation Plots")
        dpg.delete_item(self._results_body, children_only=True)
        wh = self._safe_load_image(plot_path, "gp_val_plot_tag")
        if not wh:
            self.add_log(f"[VAL] Plot image not found: {plot_path}", warn=True)
            return
        w, h = wh
        dpg.add_text("Validation Correlation Plots", parent=self._results_body)
        dpg.add_spacer(height=6, parent=self._results_body)
        width = min(1200, w)
        dpg.add_image(texture_tag="gp_val_plot_tag", parent=self._results_body,
                      width=width, height=int(h * (width / w)))

    def _save_qc_report(self, vcf_path: str, report):
        try:
            base = os.path.splitext(vcf_path)[0]
            out_path = f"{base}_qc_report.txt"
            text = self.vcf_qc_checker.to_text(report, vcf_path=vcf_path)
            with open(out_path, "w", encoding="utf-8") as f:
                f.write(text)
            self.add_log(f"[OK] QC report saved: {out_path}")
        except Exception as e:
            self.add_log(f"[ERROR] Could not save QC report: {e}", error=True)

    def _add_top_snps_tab(self, parent_tabbar_tag: str, top_n: int = 100):
        import pandas as pd
        path = getattr(self, "gwas_top_snps_name", "gwas_top_snps.csv")
        if not os.path.exists(path):
            self.add_log(f"[TopSNPs] '{path}' not found. Skipping Top SNPs tab.", warn=True)
            return
        try:
            df_top = pd.read_csv(path)
        except Exception as e:
            self.add_log(f"[TopSNPs] Could not read '{path}': {e}", warn=True)
            return

        if isinstance(top_n, int) and top_n > 0:
            df_top = df_top.head(top_n)

        with dpg.tab(label=f"Top SNPs (Top {len(df_top)})", parent=parent_tabbar_tag):
            with dpg.table(row_background=True, borders_innerH=True, borders_outerH=True,
                           borders_innerV=True, borders_outerV=True, resizable=True, sortable=True):
                for col in df_top.columns:
                    dpg.add_table_column(label=str(col))
                for _, r in df_top.iterrows():
                    with dpg.table_row():
                        for col in df_top.columns:
                            v = r[col]
                            dpg.add_text(f"{v:.6g}" if isinstance(v, float) else str(v))

    def show_results_window(self, df, algorithm, genomic_predict):
        self.ensure_results_window(show=True, title="Results")
        dpg.delete_item(self._results_body, children_only=True)

        dpg.add_button(label="Export Results", parent=self._results_body,
                       callback=lambda: dpg.show_item("select_directory"))
        dpg.add_spacer(height=10, parent=self._results_body)

        if genomic_predict:
            _ = self._safe_load_image(self.gp_plot_name, "ba_tag")
            _ = self._safe_load_image(self.gp_plot_name_scatter, "ba_tag2")

            with dpg.tab_bar(label='tabbar_results_gp', parent=self._results_body):
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
                    if dpg.does_item_exist("ba_tag2"):
                        dpg.add_image(texture_tag="ba_tag2", tag="ba_image2", width=750, height=450)
                with dpg.tab(label="Bland-Altman Plot (Model Accuracy)"):
                    if dpg.does_item_exist("ba_tag"):
                        dpg.add_image(texture_tag="ba_tag", tag="ba_image", width=750, height=450)

        else:
            _ = self._safe_load_image(self.manhatten_plot_name, "manhatten_tag")

            with dpg.tab_bar(label='tabbar_results_gwas', parent=self._results_body) as gwas_tabbar:
                with dpg.tab(label="Manhattan Plot"):
                    if dpg.does_item_exist("manhatten_tag"):
                        if algorithm in ("FaST-LMM", "Linear regression"):
                            dpg.add_image(texture_tag="manhatten_tag", tag="manhatten_image", width=950, height=400)
                        else:
                            dpg.add_image(texture_tag="manhatten_tag", tag="manhatten_image", width=900, height=300)
                    else:
                        dpg.add_text("Manhattan plot not available.")

                if algorithm in ("FaST-LMM", "Linear regression"):
                    _ = self._safe_load_image(self.qq_plot_name, "qq_tag")
                    df = df.sort_values(by=['PValue'], ascending=True)
                    with dpg.tab(label="QQ-Plot"):
                        if dpg.does_item_exist("qq_tag"):
                            dpg.add_image(texture_tag="qq_tag", tag="qq_image", height=450, width=450)
                        else:
                            dpg.add_text("QQ plot not available.")
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

                self._add_top_snps_tab(parent_tabbar_tag=gwas_tabbar, top_n=100)

                with dpg.tab(label="QTL / Region tools"):
                    with dpg.group(horizontal=True, horizontal_spacing=8):
                        dpg.add_input_text(tag="qtl_region_inp", label="Region (chr:start-end)", width=260,
                                           hint="e.g., 1:100000-500000")
                        dpg.add_button(label="Export region SNPs (CSV)", callback=self._qtl_export_region)
                        dpg.add_button(label="Plot region Manhattan", callback=self._qtl_plot_region)
                    dpg.add_spacer(height=6)
                    dpg.add_text("", tag="qtl_region_csv_lbl")
                    dpg.add_spacer(height=8)
                    dpg.add_text("Region Manhattan:", color=(200, 200, 220))
                    with dpg.child_window(tag="qtl_region_img_area", width=-1, height=460, border=True):
                        pass

    def _qtl_export_region(self, s=None, a=None):
        try:
            import pandas as pd
            region = (dpg.get_value("qtl_region_inp") or "").strip()
            if not region:
                self.add_log("[QTL] Please enter a region like '1:100000-500000'.", warn=True)
                return
            if not os.path.exists(self.gwas_result_name):
                self.add_log("[QTL] GWAS results CSV not found.", error=True)
                return
            df = pd.read_csv(self.gwas_result_name)
            safe = region.replace(":", "_").replace("-", "_").replace(",", "")
            out_csv = os.path.join(self._workspace_dir, f"gwas_region_{safe}.csv")

            self.gwas.plot_gwas(
                df=df,
                limit=None,
                algorithm=getattr(self, "algorithm", "FaST-LMM"),
                manhatten_plot_name=self.manhatten_plot_name,
                qq_plot_name=self.qq_plot_name,
                chrom_mapping=getattr(self, "_chrom_mapping_last", None),
                region=region,
                region_only_csv=out_csv,
                title_suffix="Region export"
            )
            if dpg.does_item_exist("qtl_region_csv_lbl"):
                dpg.set_value("qtl_region_csv_lbl", f"Saved: {os.path.basename(out_csv)}")
            self.add_log(f"[QTL] Region CSV saved: {out_csv}")
        except Exception as e:
            self.add_log(f"[QTL] Export failed: {e}", error=True)

    def _qtl_plot_region(self, s=None, a=None):
        try:
            import pandas as pd
            region = (dpg.get_value("qtl_region_inp") or "").strip()
            if not region:
                self.add_log("[QTL] Please enter a region like '1:100000-500000'.", warn=True)
                return
            if not os.path.exists(self.gwas_result_name):
                self.add_log("[QTL] GWAS results CSV not found.", error=True)
                return
            df = pd.read_csv(self.gwas_result_name)
            out_png = os.path.join(self._workspace_dir, "manhattan_region.png")

            self.gwas.plot_gwas(
                df=df,
                limit=None,
                algorithm=getattr(self, "algorithm", "FaST-LMM"),
                manhatten_plot_name=out_png,
                qq_plot_name=self.qq_plot_name,
                chrom_mapping=getattr(self, "_chrom_mapping_last", None),
                region=region,
                region_only_csv=None,
                title_suffix="Region view"
            )

            wh = self._safe_load_image(out_png, "manhatten_region_tag")
            if wh and dpg.does_item_exist("qtl_region_img_area"):
                dpg.delete_item("qtl_region_img_area", children_only=True)
                w, h = wh
                width = min(1100, w)
                dpg.add_image(texture_tag="manhatten_region_tag", parent="qtl_region_img_area",
                              width=width, height=int(h * (width / w)))
            self.add_log(f"[QTL] Region Manhattan saved: {out_png}")
        except Exception as e:
            self.add_log(f"[QTL] Plot failed: {e}", error=True)

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
                    "table_gwas", "table_gp", "ba_tag", "ba_tag2", "ba_image", "ba_image2",
                    "gp_val_plot_tag"]:
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
            self.pheno_stats_name, self.geno_stats_name,
            self.gwas_top_snps_name
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
