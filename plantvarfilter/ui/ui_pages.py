# ui/ui_pages.py
from __future__ import annotations
import traceback
from pathlib import Path
import os
import platform as _pf
import subprocess
import dearpygui.dearpygui as dpg
import traceback
import os
import dearpygui.dearpygui as dpg

from ..pangenome_builder import build_pangenome_graph


try:
    from plantvarfilter.preanalysis import (
        ReferenceManager, ReferenceIndexStatus,
        run_fastq_qc, FastqQCReport,
        Aligner, AlignmentResult,
    )
    _HAS_PRE = True
    PLATFORM_DISPLAY = [
        "Illumina (short reads)",
        "Oxford Nanopore (ONT)",
        "PacBio HiFi (CCS)",
        "PacBio CLR",
    ]

    DISPLAY_TO_KEY = {
        "Illumina (short reads)": "illumina",
        "Oxford Nanopore (ONT)": "ont",
        "PacBio HiFi (CCS)": "hifi",
        "PacBio CLR": "pb",
    }

    KEY_TO_MINIMAP2_PRESET = {
        "ont": "map-ont",
        "hifi": "map-hifi",
        "pb": "map-pb",
    }
except Exception:
    _HAS_PRE = False


def _default_start_dir(app=None, prefer_desktop=True) -> str:
    env = os.environ.get("PVF_START_DIR")
    if env and Path(env).exists():
        return env
    if app is not None:
        ws = getattr(app, "workspace_dir", None)
        if ws and Path(ws).exists():
            return str(Path(ws))
    home = Path.home()
    desktop = home / "Desktop"
    if prefer_desktop and desktop.exists():
        return str(desktop)
    return str(home)


def _open_in_os(path: str):
    try:
        if _pf.system() == "Windows":
            os.startfile(path)  # type: ignore
        elif _pf.system() == "Darwin":
            subprocess.Popen(["open", path])
        else:
            subprocess.Popen(["xdg-open", path])
    except Exception:
        pass


def _file_input_row(
    label: str,
    tag_key: str,
    parent,
    file_extensions: tuple[str, ...] = (".*",),
    app=None,
    default_dir: str | None = None,
):
    with dpg.group(parent=parent, horizontal=True, horizontal_spacing=8):
        dpg.add_text(label)
        dpg.add_input_text(tag=f"input_{tag_key}", width=460)
        dpg.add_button(
            label="Browse",
            callback=lambda: dpg.show_item(f"dlg_{tag_key}"),
        )
    base = default_dir or _default_start_dir(app)
    with dpg.file_dialog(
        tag=f"dlg_{tag_key}",
        directory_selector=False,
        show=False,
        default_path=base,
        callback=lambda s, a: dpg.set_value(
            f"input_{tag_key}", a["file_path_name"]
        ),
    ):
        for ext in file_extensions:
            dpg.add_file_extension(ext, color=(150, 150, 150, 255))

def _dir_input_row(
    label: str,
    tag_key: str,
    parent,
    app=None,
    default_dir: str | None = None,
):
    with dpg.group(parent=parent, horizontal=True, horizontal_spacing=8):
        dpg.add_text(label)
        dpg.add_input_text(tag=f"input_{tag_key}", width=460)
        dpg.add_button(
            label="Browse",
            callback=lambda: dpg.show_item(f"dlg_{tag_key}"),
        )
    base = default_dir or _default_start_dir(app)
    with dpg.file_dialog(
        tag=f"dlg_{tag_key}",
        directory_selector=True,
        show=False,
        default_path=base,
        callback=lambda s, a: dpg.set_value(
            f"input_{tag_key}", a["file_path_name"]
        ),
    ):
        pass

def page_reference_manager(app, parent):
    with dpg.group(parent=parent, show=False, tag="page_ref_manager"):
        dpg.add_text("\nReference Manager", indent=10)
        dpg.add_spacer(height=10)

        _file_input_row(
            "Reference FASTA:",
            "ref_fasta",
            parent,
            (".fa", ".fasta", ".fna", ".*"),
            app=app,
        )
        _dir_input_row(
            "Reference out dir (optional):",
            "ref_out_dir",
            parent,
            app=app,
        )

        dpg.add_spacer(height=6)
        with dpg.group(horizontal=True, horizontal_spacing=8):
            btn_build = dpg.add_button(
                tag="btn_ref_build",
                label="Build / Refresh Indexes",
                width=240,
                height=36,
                callback=lambda: _on_build_reference(app),
            )
            with dpg.tooltip(btn_build):
                dpg.add_text(
                    "Build or refresh genome reference indexes (FAI, DICT, minimap2, bowtie2) for the selected FASTA."
                )
            dpg.add_loading_indicator(
                tag="ref_build_spinner",
                radius=6,
                style=2,
                show=False,
            )

    return "page_ref_manager"

def page_pangenome_builder(app, parent):


    with dpg.group(parent=parent, show=False, tag="page_pangenome"):
        dpg.add_text("\nPangenome Builder", indent=10)
        dpg.add_spacer(height=10)

        _file_input_row(
            "Base Reference FASTA:",
            "pan_base_ref",
            parent,
            (".fa", ".fasta", ".fna"),
            app=app,
        )

        _file_input_row(
            "Assemblies / Consensus FASTA (file):",
            "pan_assembly_fasta",
            parent,
            (".fa", ".fasta", ".fna"),
            app=app,
        )

        _dir_input_row(
            "Assemblies / Consensus Folder (many FASTA):",
            "pan_assemblies_dir",
            parent,
            app=app,
        )

        _dir_input_row(
            "Output Folder:",
            "pan_out_dir",
            parent,
            app=app,
        )

        dpg.add_spacer(height=10)
        with dpg.group(parent=parent, horizontal=True, horizontal_spacing=12):
            mode_lbl = dpg.add_text("Mode:")
            mode_combo = dpg.add_combo(
                items=["Fast (subset 25)", "Full (all FASTA in folder)"],
                default_value="Fast (subset 25)",
                tag="pan_mode",
                width=220,
            )

            min_lbl = dpg.add_text("Min contig length:")
            min_inp = dpg.add_input_int(
                default_value=1000,
                min_value=0,
                step=100,
                tag="pan_min_contig_len",
                width=120,
            )

        with dpg.tooltip(mode_lbl):
            dpg.add_text("Fast: builds from a subset for quick testing. Full: uses all assemblies in the folder.")
        with dpg.tooltip(mode_combo):
            dpg.add_text("Choose between a quick subset build or a full build over all assemblies.")
        with dpg.tooltip(min_lbl):
            dpg.add_text("Contigs shorter than this length will be ignored to reduce noise and speed up building.")
        with dpg.tooltip(min_inp):
            dpg.add_text("Minimum contig length threshold (bp).")

        dpg.add_spacer(height=10)

        def _get_vals():
            base_ref = (dpg.get_value("input_pan_base_ref") or "").strip()
            asm_file = (dpg.get_value("input_pan_assembly_fasta") or "").strip()
            asm_dir = (dpg.get_value("input_pan_assemblies_dir") or "").strip()
            out_dir = (dpg.get_value("input_pan_out_dir") or "").strip()
            mode = dpg.get_value("pan_mode")
            min_len = dpg.get_value("pan_min_contig_len")
            return base_ref, asm_file, asm_dir, out_dir, mode, min_len

        def _validate_inputs(base_ref, asm_file, asm_dir, out_dir):
            ok = True
            if not base_ref:
                app.add_log("[PAN] Base Reference FASTA is required.", error=True)
                ok = False
            if not asm_file and not asm_dir:
                app.add_log("[PAN] Provide either a FASTA file OR a folder of FASTA files.", error=True)
                ok = False
            if asm_file and asm_dir:
                app.add_log("[PAN] Choose only one assemblies input: file OR folder.", error=True)
                ok = False
            if not out_dir:
                app.add_log("[PAN] Output folder is required.", error=True)
                ok = False
            return ok

        def _busy_build(is_busy: bool):
            if dpg.does_item_exist("pan_build_btn"):
                dpg.configure_item("pan_build_btn", enabled=not is_busy)
            if dpg.does_item_exist("pan_build_spinner"):
                dpg.configure_item("pan_build_spinner", show=is_busy)

        def _busy_use(is_busy: bool):
            if dpg.does_item_exist("pan_use_ref_btn"):
                dpg.configure_item("pan_use_ref_btn", enabled=not is_busy)
            if dpg.does_item_exist("pan_use_ref_spinner"):
                dpg.configure_item("pan_use_ref_spinner", show=is_busy)

        def on_build_pangenome(sender=None, app_data=None, user_data=None):
            base_ref, asm_file, asm_dir, out_dir, mode, min_len = _get_vals()

            app.ensure_log_window(show=True)
            if not _validate_inputs(base_ref, asm_file, asm_dir, out_dir):
                return

            _busy_build(True)
            app.add_log("[PAN] Build started.")

            def task():
                def logger(msg):
                    app.ui_emit(app.add_log, msg)

                assemblies_input = asm_file if asm_file else asm_dir
                threads = int(os.environ.get("PLANTVARFILTER_THREADS", "8"))
                batch_size = 20 if str(mode).lower().startswith("full") else None

                res = build_pangenome_graph(
                    base_reference_fasta=base_ref,
                    assemblies_input=assemblies_input,
                    output_dir=out_dir,
                    mode=str(mode),
                    subset_n=25,
                    threads=threads,
                    min_contig_len=int(min_len) if min_len is not None else 0,
                    minigraph_preset="ggs",
                    batch_size=batch_size,
                    logger=logger,
                )
                return res

            def on_done(res):
                app._pan_last_gfa = res.pangenome_gfa
                app.add_log(f"[PAN] Output GFA: {res.pangenome_gfa}")
                app.add_log(f"[PAN] Report: {res.report_txt}")
                app.add_log(f"[PAN] Log: {res.log_txt}")
                app.add_log("[PAN] Build finished ✔")
                _busy_build(False)

            def on_error(exc):
                app.add_log(f"[PAN] Build failed: {exc}", error=True)
                _busy_build(False)

            app.run_task("pangenome_build", task, on_done=on_done, on_error=on_error)

        def on_use_as_reference(sender=None, app_data=None, user_data=None):
            gfa = (getattr(app, "_pan_last_gfa", "") or "").strip()
            if not gfa:
                app.ensure_log_window(show=True)
                app.add_log("[PAN] No pangenome GFA found. Build first.", error=True)
                return

            _busy_use(True)
            try:
                app.ensure_log_window(show=True)
                if hasattr(app, "use_reference") and callable(getattr(app, "use_reference")):
                    app.use_reference(reference_path=gfa, reference_type="gfa")
                    app.add_log(f"[PAN] Reference set: {gfa}")
                else:
                    app.add_log("[PAN] Backend not wired yet: implement GWASApp.use_reference(reference_path, reference_type).", warn=True)
            except Exception as exc:
                app.add_log(f"[PAN] Error: {exc}", error=True)
            finally:
                _busy_use(False)

        with dpg.group(parent=parent, horizontal=True, horizontal_spacing=10):
            build_btn = dpg.add_button(
                label="Build Pangenome",
                callback=on_build_pangenome,
                width=180,
                tag="pan_build_btn",
            )
            with dpg.tooltip(build_btn):
                dpg.add_text("Build the pangenome graph (GFA) into the selected output folder.")

            use_btn = dpg.add_button(
                label="Use as Reference",
                callback=on_use_as_reference,
                width=180,
                tag="pan_use_ref_btn",
            )
            with dpg.tooltip(use_btn):
                dpg.add_text("Set the generated pangenome GFA as the active reference for downstream steps.")

            dpg.add_loading_indicator(
                tag="pan_build_spinner",
                radius=6,
                style=2,
                show=False,
            )
            dpg.add_loading_indicator(
                tag="pan_use_ref_spinner",
                radius=6,
                style=2,
                show=False,
            )


def _render_ref_status(app, st: ReferenceIndexStatus):
    lines = [
        f"FASTA: {st.fasta}",
        f"Directory: {st.reference_dir}",
        f"faidx (.fai): {st.faidx or 'missing'}",
        f"dict (.dict): {st.dict or 'missing'}",
        f"minimap2 index (.mmi): {st.mmi or 'missing'}",
        f"bowtie2 prefix: {st.bt2_prefix or 'missing'}",
        "Tools: "
        + ", ".join(
            f"{k}: {'ok' if (v.get('path')) else 'missing'}"
            for k, v in st.tools.items()
        ),
        f"OK: {st.ok}",
    ]
    for ln in lines:
        app.add_log(f"[REF] {ln}")


def _on_build_reference(app):
    dpg.configure_item("ref_build_spinner", show=True)
    dpg.configure_item("btn_ref_build", enabled=False)

    try:
        if not _HAS_PRE:
            app.add_log(
                "[REF] Preanalysis modules not available in environment.",
                error=True,
            )
            return

        fasta = dpg.get_value("input_ref_fasta")
        out_dir_value = dpg.get_value("input_ref_out_dir") or None

        if not fasta or not Path(fasta).exists():
            app.add_log("[REF] Select a valid FASTA file first.", warn=True)
            return

        fasta_path = Path(fasta)

        if out_dir_value is None:
            out_dir = fasta_path.parent / "reference"
        else:
            out_dir = Path(out_dir_value)

        out_dir.mkdir(parents=True, exist_ok=True)

        app.add_log("[REF] Starting reference indexing...")
        app.add_log(f"[REF] FASTA: {fasta_path}")
        app.add_log(f"[REF] Output directory: {out_dir}")

        rm = ReferenceManager(logger=app.add_log, workspace=app.workspace_dir)
        st = rm.build_indices(str(fasta_path), out_dir=str(out_dir))

        _render_ref_status(app, st)

        if st.ok:
            app.add_log("[REF] Reference indexing completed successfully ✔")
        else:
            app.add_log(
                "[REF] Reference indexing finished, but some indexes are missing.",
                warn=True,
            )

    except Exception as exc:
        app.add_log(f"[REF] Error while building reference: {exc}", error=True)
    finally:
        dpg.configure_item("ref_build_spinner", show=False)
        dpg.configure_item("btn_ref_build", enabled=True)

def page_fastq_qc(app, parent):
    with dpg.group(parent=parent, show=False, tag="page_fastq_qc"):
        dpg.add_text("\nFASTQ Quality Control", indent=10)
        dpg.add_spacer(height=10)

        with dpg.group(parent=parent, horizontal=True, horizontal_spacing=12):
            dpg.add_text("Platform:")
            dpg.add_combo(
                items=["illumina", "ont", "hifi", "pb"],
                default_value="illumina",
                width=180,
                tag="fq_platform",
            )

        _file_input_row(
            "Reads #1 (R1 or single):",
            "fq_r1",
            parent,
            (".fastq", ".fq", ".fastq.gz", ".fq.gz", ".*"),
            app=app,
        )
        _file_input_row(
            "Reads #2 (R2, optional):",
            "fq_r2",
            parent,
            (".fastq", ".fq", ".fastq.gz", ".fq.gz", ".*"),
            app=app,
        )
        _dir_input_row("Output dir (optional):", "fq_out", parent, app=app)

        dpg.add_spacer(height=6)
        with dpg.group(parent=parent, horizontal=True, horizontal_spacing=12):
            btn_qc = dpg.add_button(
                tag="btn_fq_qc",
                label="Run QC",
                width=200,
                height=36,
                callback=lambda: _on_run_fastq_qc(app),
            )
            with dpg.tooltip(btn_qc):
                dpg.add_text(
                    "Run basic QC on FASTQ reads.\n"
                    "If successful, a summary and paths to QC reports\n"
                    "will be written to the Log window."
                )

            dpg.add_checkbox(
                label="Use FastQC if available",
                default_value=True,
                tag="fq_use_fastqc",
            )
            dpg.add_loading_indicator(
                tag="fq_qc_spinner",
                radius=6,
                style=2,
                show=False,
            )

    return "page_fastq_qc"


def _render_qc(app, rep: FastqQCReport):
    fields = [
        ("Platform", rep.platform),
        ("Reads (sampled)", rep.n_reads),
        ("Mean length", f"{rep.mean_length:.2f}"),
        ("Median length", f"{rep.median_length:.2f}"),
        ("GC%", f"{rep.gc_percent:.2f}"),
        ("N%", f"{rep.n_percent:.3f}"),
        ("Mean PHRED", "NA" if rep.mean_phred is None else f"{rep.mean_phred:.2f}"),
        ("Verdict", rep.verdict),
        ("Report TXT", rep.report_txt),
    ]

    app.add_log("[FQ-QC] QC summary:")
    for k, v in fields:
        app.add_log(f"[FQ-QC] {k}: {v}")

    if rep.length_hist_png:
        app.add_log(f"[FQ-QC] Length histogram: {rep.length_hist_png}")
    if rep.gc_hist_png:
        app.add_log(f"[FQ-QC] GC% histogram: {rep.gc_hist_png}")
    if rep.per_cycle_q_mean_png:
        app.add_log(
            f"[FQ-QC] Per-cycle mean PHRED plot: {rep.per_cycle_q_mean_png}"
        )

    out_dir = str(Path(rep.report_txt).parent)
    app.add_log(f"[FQ-QC] QC output folder: {out_dir}")


def _on_run_fastq_qc(app):
    if not _HAS_PRE:
        app.add_log(
            "[FQ-QC] Preanalysis modules not available in environment.",
            error=True,
        )
        return

    dpg.configure_item("fq_qc_spinner", show=True)
    dpg.configure_item("btn_fq_qc", enabled=False)

    try:
        r1 = dpg.get_value("input_fq_r1")
        r2 = dpg.get_value("input_fq_r2") or None
        if not r1 or not Path(r1).exists():
            app.add_log("[FQ-QC] Select valid FASTQ file(s).", warn=True)
            return

        out_dir = dpg.get_value("input_fq_out") or None
        platform = (dpg.get_value("fq_platform") or "illumina").lower()
        use_fastqc = bool(dpg.get_value("fq_use_fastqc"))

        app.add_log(
            f"[FQ-QC] Running QC on platform='{platform}' | "
            f"R1='{r1}' | R2='{r2}' | out_dir='{out_dir or 'auto'}'"
        )

        rep = run_fastq_qc(
            r1,
            r2,
            platform=platform,
            out_dir=out_dir,
            use_fastqc_if_available=use_fastqc,
            logger=app.add_log,
        )

        _render_qc(app, rep)
        app.add_log("[FQ-QC] QC run completed successfully.")

    except Exception as exc:
        app.add_log(f"[FQ-QC] Error while running QC: {exc}", error=True)
        app.add_log(traceback.format_exc(), error=True)
    finally:
        dpg.configure_item("fq_qc_spinner", show=False)
        dpg.configure_item("btn_fq_qc", enabled=True)


def _render_qc(rep: FastqQCReport):
    dpg.delete_item("fq_qc_area", children_only=True)
    fields = [
        ("Platform", rep.platform),
        ("Reads (sampled)", rep.n_reads),
        ("Mean length", f"{rep.mean_length:.2f}"),
        ("Median length", f"{rep.median_length:.2f}"),
        ("GC%", f"{rep.gc_percent:.2f}"),
        ("N%", f"{rep.n_percent:.3f}"),
        ("Mean PHRED", "NA" if rep.mean_phred is None else f"{rep.mean_phred:.2f}"),
        ("Verdict", rep.verdict),
        ("Report TXT", rep.report_txt),
    ]
    for k, v in fields:
        dpg.add_text(f"{k}: {v}", parent="fq_qc_area")
    if rep.length_hist_png:
        dpg.add_text(f"Length hist: {rep.length_hist_png}", parent="fq_qc_area")
    if rep.gc_hist_png:
        dpg.add_text(f"GC% hist: {rep.gc_hist_png}", parent="fq_qc_area")
    if rep.per_cycle_q_mean_png:
        dpg.add_text(f"Per-cycle mean PHRED: {rep.per_cycle_q_mean_png}", parent="fq_qc_area")
    dpg.add_spacer(height=6, parent="fq_qc_area")
    out_dir = str(Path(rep.report_txt).parent)
    dpg.add_button(label="Open QC folder", parent="fq_qc_area",
                   callback=lambda: _open_in_os(out_dir))


def _on_run_fastq_qc(app):
    if not _HAS_PRE:
        app.add_log("[FQ-QC] Preanalysis modules not available in environment.", error=True)
        return
    r1 = dpg.get_value("input_fq_r1")
    r2 = dpg.get_value("input_fq_r2") or None
    if not r1 or not Path(r1).exists():
        app.add_log("[FQ-QC] Select valid FASTQ file(s).", warn=True)
        return
    out_dir = dpg.get_value("input_fq_out") or None
    platform = (dpg.get_value("fq_platform") or "illumina").lower()
    use_fastqc = bool(dpg.get_value("fq_use_fastqc"))
    app.add_log(f"[FQ-QC] Running on platform={platform}...")
    rep = run_fastq_qc(r1, r2, platform=platform, out_dir=out_dir, use_fastqc_if_available=use_fastqc, logger=app.add_log)
    _render_qc(rep)
    app.add_log("[FQ-QC] Done.")


def page_alignment(app, parent):
    with dpg.group(parent=parent, show=False, tag="page_alignment"):
        dpg.add_text("\nAlignment", indent=10)
        dpg.add_spacer(height=10)

        # Platform selector
        with dpg.group(parent=parent, horizontal=True, horizontal_spacing=12):
            dpg.add_text("Platform:")
            dpg.add_combo(
                items=["illumina", "ont", "hifi", "pb"],
                default_value="illumina",
                width=180,
                tag="aln_platform",
            )

        # Reference section
        with dpg.collapsing_header(label="Reference", parent=parent, default_open=True):
            dpg.add_text(
                "For ONT/PB: select FASTA or .mmi\n"
                "For Illumina: select bowtie2 prefix base path (no suffix)."
            )
            _file_input_row(
                "Reference (.fa/.mmi or bt2 prefix)",
                "aln_reference",
                parent,
                (".fa", ".fasta", ".fna", ".mmi", ".*"),
                app=app,
            )

        # Reads section
        with dpg.collapsing_header(label="Reads", parent=parent, default_open=True) as reads_header:
            _file_input_row(
                "Reads #1 (R1 or single):",
                "aln_r1",
                reads_header,
                (".fastq", ".fq", ".fastq.gz", ".fq.gz", ".*"),
                app=app,
            )
            _file_input_row(
                "Reads #2 (R2, optional):",
                "aln_r2",
                reads_header,
                (".fastq", ".fq", ".fastq.gz", ".fq.gz", ".*"),
                app=app,
            )

        # Output dir
        _dir_input_row("Output dir (optional):", "aln_out", parent, app=app)

        dpg.add_spacer(height=6)

        # Threads + options
        with dpg.group(parent=parent, horizontal=True, horizontal_spacing=12):
            dpg.add_text("Threads")

            aln_threads = dpg.add_input_int(
                tag="aln_threads",
                default_value=4,          # safer default
                min_value=1,
                max_value=64,
                min_clamped=True,
                max_clamped=True,
                width=100,
            )
            with dpg.tooltip(aln_threads):
                dpg.add_text(
                    "Number of CPU threads used by the aligner.\n"
                    "- Higher = faster but uses more CPU/RAM.\n"
                    "- For ONT on laptops: 2–4 is usually safer.\n"
                    "- For strong servers: 8–16 is often fine."
                )

            save_sam = dpg.add_checkbox(label="Save SAM", tag="aln_save_sam")
            with dpg.tooltip(save_sam):
                dpg.add_text(
                    "If enabled, also write a SAM file in addition to BAM.\n"
                    "This is mainly useful for debugging and will use extra disk."
                )

            markdup = dpg.add_checkbox(
                label="Mark duplicates",
                default_value=True,
                tag="aln_markdup",
            )
            with dpg.tooltip(markdup):
                dpg.add_text(
                    "Mark PCR/optical duplicates in the BAM file.\n"
                    "Recommended for most Illumina data; for ONT it is optional."
                )

        # Read Group
        with dpg.collapsing_header(label="Read Group (optional)", parent=parent, default_open=False):
            for k in ("ID", "SM", "LB", "PL", "PU"):
                with dpg.group(horizontal=True, horizontal_spacing=8):
                    dpg.add_text(k)
                    dpg.add_input_text(tag=f"aln_rg_{k}", width=240)

        dpg.add_spacer(height=8)

        # Run button + spinner
        with dpg.group(parent=parent, horizontal=True, horizontal_spacing=12):
            btn_aln = dpg.add_button(
                tag="btn_aln_run",
                label="Run Alignment",
                width=200,
                height=36,
                callback=lambda: _on_run_alignment(app),
            )
            with dpg.tooltip(btn_aln):
                dpg.add_text(
                    "Align reads to the selected reference.\n"
                    "If successful, BAM/BAI and basic stats\n"
                    "will be written to the output folder and Log."
                )

            dpg.add_loading_indicator(
                tag="aln_spinner",
                radius=6,
                style=2,
                show=False,
            )

    return "page_alignment"

def _bt2_prefix_exists(prefix: str) -> bool:
    p = Path(prefix)
    suff = [".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"]
    return all((p.parent / (p.name + s)).exists() for s in suff)

def _resolve_reference_for_alignment(ref: str, platform: str):
    if not ref:
        return None
    p = Path(ref)
    plat = (platform or "illumina").lower()

    if plat in {"ont", "hifi", "pb"}:
        if p.suffix.lower() in {".fa", ".fasta", ".fna", ".mmi"} and p.exists():
            return str(p)
        return None

    if _bt2_prefix_exists(ref):
        return ref

    if p.suffix.lower() in {".fa", ".fasta", ".fna"} and p.exists():
        cand = [
            str(p.with_suffix("")),
            str((p.parent / (p.stem + "_bt2")).resolve()),
        ]
        for c in cand:
            if _bt2_prefix_exists(c):
                return c
        return None

    return None


def _on_run_alignment(app):
    import traceback
    try:
        if not _HAS_PRE:
            app.add_log(
                "[ALN] Preanalysis modules not available in environment.",
                error=True,
            )
            return

        dpg.configure_item("aln_spinner", show=True)
        dpg.configure_item("btn_aln_run", enabled=False)

        ref_in = dpg.get_value("input_aln_reference")
        r1 = dpg.get_value("input_aln_r1")
        r2 = dpg.get_value("input_aln_r2") or None
        platform = (dpg.get_value("aln_platform") or "illumina").lower()
        threads = dpg.get_value("aln_threads") or 8
        save_sam = bool(dpg.get_value("aln_save_sam"))
        markdup = bool(dpg.get_value("aln_markdup"))

        out_dir_ui = dpg.get_value("input_aln_out") or ""
        if out_dir_ui and not Path(out_dir_ui).exists():
            try:
                Path(out_dir_ui).mkdir(parents=True, exist_ok=True)
            except Exception as e:
                app.add_log(
                    f"[ALN] Could not create output dir '{out_dir_ui}': {e}",
                    warn=True,
                )
                out_dir_ui = ""

        if r1 and Path(r1).exists():
            default_out = str(Path(r1).parent)
        else:
            default_out = app.workspace_dir
        out_dir = out_dir_ui or default_out

        app.add_log(
            f"[ALN] START | ref='{ref_in}' | r1='{r1}' | r2='{r2}' | "
            f"platform={platform} | threads={threads} | out='{out_dir}'"
        )

        ref_resolved = _resolve_reference_for_alignment(ref_in, platform)
        if not ref_resolved:
            app.add_log(
                "[ALN] Select a valid reference (.fa/.mmi or bowtie2 prefix).",
                warn=True,
            )
            return

        r1_ok = bool(r1) and Path(r1).is_file()
        r2_ok = bool(r2) and Path(r2).is_file()
        app.add_log(
            f"[ALN] DEBUG reads: R1 exists={r1_ok} | R2 exists={r2_ok}"
        )

        if not r1_ok:
            app.add_log(
                "[ALN] Select valid FASTQ reads (R1).",
                warn=True,
            )
            return
        if r2 and not r2_ok:
            app.add_log(
                "[ALN] R2 not found or unreadable; continuing as single-end.",
                warn=True,
            )
            r2 = None

        def _sample_prefix_from_r1(p: str) -> str:
            name = Path(p).name
            if name.endswith(".gz"):
                name = name[:-3]
            for ext in (".fastq", ".fq"):
                if name.endswith(ext):
                    name = name[: -len(ext)]
            for tag in ("_R1", ".R1"):
                if name.endswith(tag):
                    name = name[: -len(tag)]
            return name

        base_prefix = _sample_prefix_from_r1(r1)
        out_prefix = (
            f"{base_prefix}.minimap2"
            if platform in {"ont", "hifi", "pb"}
            else f"{base_prefix}.bowtie2"
        )

        rg = {}
        for k in ("ID", "SM", "LB", "PL", "PU"):
            v = dpg.get_value(f"aln_rg_{k}")
            if v:
                rg[k] = v

        aln = Aligner(logger=app.add_log, workspace=app.workspace_dir)

        if platform in {"ont", "hifi", "pb"}:
            preset = (
                "map-ont"
                if platform == "ont"
                else ("map-hifi" if platform == "hifi" else "map-pb")
            )
            reads = [r1] if not r2 else [r1, r2]
            res = aln.minimap2(
                ref_resolved,
                reads,
                preset=preset,
                threads=threads,
                read_group=rg or None,
                save_sam=save_sam,
                mark_duplicates=markdup,
                out_dir=out_dir,
                out_prefix=out_prefix,
            )
        else:
            res = aln.bowtie2(
                ref_resolved,
                r1,
                r2,
                threads=threads,
                read_group=rg or None,
                save_sam=save_sam,
                mark_duplicates=markdup,
                out_dir=out_dir,
                out_prefix=out_prefix,
            )

        # Log summary instead of using a Result window
        app.add_log(f"[ALN] Tool: {res.tool}")
        if res.sam:
            app.add_log(f"[ALN] SAM: {res.sam}")
        app.add_log(f"[ALN] BAM: {res.bam}")
        app.add_log(f"[ALN] BAI: {res.bai}")
        app.add_log(f"[ALN] flagstat: {res.flagstat}")
        app.add_log(f"[ALN] Elapsed: {res.elapsed_sec:.1f} sec")
        app.add_log(
            f"[ALN] Output folder: {Path(res.bam).parent}"
        )
        app.add_log("[ALN] Alignment finished.")

    except Exception as e:
        app.add_log(f"[ALN] ERROR: {e}", error=True)
        app.add_log(traceback.format_exc(), error=True)
    finally:
        dpg.configure_item("aln_spinner", show=False)
        dpg.configure_item("btn_aln_run", enabled=True)


def page_preprocess_samtools(app, parent):
    with dpg.group(parent=parent, show=False, tag="page_pre_sam"):
        dpg.add_text(
            "\nClean BAM from SAM/BAM: sort / fixmate / markdup / index + QC reports",
            indent=10,
        )
        dpg.add_spacer(height=10)

        with dpg.group(horizontal=True, horizontal_spacing=60):
            with dpg.group():
                bam_btn = dpg.add_button(
                    label="Choose SAM/BAM",
                    callback=lambda: dpg.show_item("file_dialog_bam"),
                    width=220,
                    tag="tooltip_bam_sam",
                )
                app._secondary_buttons.append(bam_btn)
                dpg.add_text("No file", tag="sam_bam_path_lbl", wrap=500)

            with dpg.group():
                app.sam_threads = dpg.add_input_int(
                    label="Threads",
                    width=220,
                    default_value=4,
                    min_value=1,
                    min_clamped=True,
                )
                app._inputs.append(app.sam_threads)

                app.sam_remove_dups = dpg.add_checkbox(
                    label="Remove duplicates (instead of marking)",
                    default_value=False,
                )
                app._inputs.append(app.sam_remove_dups)

                app.sam_compute_stats = dpg.add_checkbox(
                    label="Compute QC reports (flagstat/stats/idxstats/depth)",
                    default_value=True,
                )
                app._inputs.append(app.sam_compute_stats)

                dpg.add_spacer(height=6)
                app.sam_out_prefix = dpg.add_input_text(
                    label="Output prefix (optional)",
                    hint="Leave empty to auto-generate next to the input file",
                    width=320,
                )
                app._inputs.append(app.sam_out_prefix)

                dpg.add_spacer(height=12)
                run_sam = dpg.add_button(
                    tag="btn_sam_preprocess",
                    label="Run samtools preprocess",
                    callback=lambda: _on_run_sam_preprocess(app),
                    width=240,
                    height=38,
                )
                app._primary_buttons.append(run_sam)

                with dpg.tooltip(run_sam):
                    dpg.add_text(
                        "Run samtools pipeline to clean the BAM: sort, fixmate, "
                        "mark or remove duplicates, index and generate QC reports."
                    )

                dpg.add_loading_indicator(
                    tag="sam_pre_spinner",
                    radius=6,
                    style=2,
                    show=False,
                )

    return "page_pre_sam"


def _on_run_sam_preprocess(app):
    dpg.configure_item("sam_pre_spinner", show=True)
    dpg.configure_item("btn_sam_preprocess", enabled=False)

    try:
        app.run_samtools_preprocess()
    except Exception as exc:
        app.add_log(f"[SAMTOOLS] Error during preprocess: {exc}", error=True)
        app.add_log(traceback.format_exc(), error=True)
    finally:
        dpg.configure_item("sam_pre_spinner", show=False)
        dpg.configure_item("btn_sam_preprocess", enabled=True)


def page_variant_calling(app, parent):
    with dpg.group(parent=parent, show=False, tag="page_vc"):
        dpg.add_text("\nCall variants with bcftools mpileup + call", indent=10)
        dpg.add_spacer(height=10)

        with dpg.group(horizontal=True, horizontal_spacing=60):
            with dpg.group():
                b1 = dpg.add_button(
                    label="Choose BAM (single)",
                    callback=lambda: dpg.show_item("file_dialog_bam_vc"),
                    width=220,
                    tag="tooltip_bam_vc",
                )
                app._secondary_buttons.append(b1)
                dpg.add_text("", tag="vc_bam_path_lbl", wrap=500)

                dpg.add_spacer(height=6)
                b2 = dpg.add_button(
                    label="Choose BAM-list (.list)",
                    callback=lambda: dpg.show_item("file_dialog_bamlist"),
                    width=220,
                    tag="tooltip_bamlist_vc",
                )
                app._secondary_buttons.append(b2)
                dpg.add_text("", tag="vc_bamlist_path_lbl", wrap=500)

                dpg.add_spacer(height=6)
                fa_btn = dpg.add_button(
                    label="Choose reference FASTA",
                    callback=lambda: dpg.show_item("file_dialog_fasta"),
                    width=220,
                    tag="tooltip_fa_vc",
                )
                app._secondary_buttons.append(fa_btn)
                dpg.add_text("", tag="vc_ref_path_lbl", wrap=500)

                dpg.add_spacer(height=6)
                reg_btn2 = dpg.add_button(
                    label="Choose regions BED (optional)",
                    callback=lambda: dpg.show_item("file_dialog_blacklist"),
                    width=220,
                    tag="tooltip_reg_vc",
                )
                app._secondary_buttons.append(reg_btn2)
                dpg.add_text("", tag="vc_regions_path_lbl", wrap=500)

            with dpg.group():
                app.vc_threads = dpg.add_input_int(
                    label="Threads",
                    width=220,
                    default_value=4,
                    min_value=1,
                    min_clamped=True,
                )
                app._inputs.append(app.vc_threads)

                app.vc_ploidy = dpg.add_input_int(
                    label="Ploidy",
                    width=220,
                    default_value=2,
                    min_value=1,
                    min_clamped=True,
                    tag="tooltip_ploidy",
                )
                app._inputs.append(app.vc_ploidy)

                app.vc_min_bq = dpg.add_input_int(
                    label="Min BaseQ",
                    width=220,
                    default_value=20,
                    min_value=0,
                    min_clamped=True,
                    tag="tooltip_bq",
                )
                app._inputs.append(app.vc_min_bq)

                app.vc_min_mq = dpg.add_input_int(
                    label="Min MapQ",
                    width=220,
                    default_value=20,
                    min_value=0,
                    min_clamped=True,
                    tag="tooltip_mq",
                )
                app._inputs.append(app.vc_min_mq)

                dpg.add_spacer(height=6)
                app.vc_out_prefix = dpg.add_input_text(
                    label="Output prefix (optional)",
                    hint="Leave empty to auto-generate next to the BAM",
                    width=320,
                )
                app._inputs.append(app.vc_out_prefix)

                dpg.add_spacer(height=6)
                app.vc_split_after = dpg.add_checkbox(
                    label="Split VCF by variant type (SNPs / INDELs)",
                    default_value=False,
                    tag="vc_split_after_calling",
                )
                app._inputs.append(app.vc_split_after)

                dpg.add_spacer(height=12)
                run_vc = dpg.add_button(
                    tag="btn_vc_call",
                    label="Call variants (bcftools)",
                    callback=lambda: _on_run_variant_calling(app),
                    width=240,
                    height=38,
                )
                app._primary_buttons.append(run_vc)

                with dpg.tooltip(run_vc):
                    dpg.add_text(
                        "Run bcftools mpileup + call on the selected BAM/BAM-list "
                        "to generate a VCF file (optionally split into SNPs / INDELs)."
                    )

                dpg.add_loading_indicator(
                    tag="vc_call_spinner",
                    radius=6,
                    style=2,
                    show=False,
                )

    return "page_vc"


def _on_run_variant_calling(app):
    dpg.configure_item("vc_call_spinner", show=True)
    dpg.configure_item("btn_vc_call", enabled=False)

    try:
        app.run_variant_calling()
    except Exception as exc:
        app.add_log(f"[VC] Error during variant calling: {exc}", error=True)
        app.add_log(traceback.format_exc(), error=True)
    finally:
        dpg.configure_item("vc_call_spinner", show=False)
        dpg.configure_item("btn_vc_call", enabled=True)


def page_preprocess_bcftools(app, parent):
    with dpg.group(parent=parent, show=False, tag="page_pre_bcf"):
        dpg.add_text("\nNormalize / split multiallelic / sort / filter / set IDs (bcftools)", indent=10)
        dpg.add_spacer(height=10)

        with dpg.group(horizontal=True, horizontal_spacing=60):
            with dpg.group():
                vcf_btn_bcf = dpg.add_button(
                    label="Choose a VCF file",
                    callback=lambda: dpg.show_item("file_dialog_vcf"),
                    width=220,
                    tag="tooltip_vcf_bcf",
                )
                app._secondary_buttons.append(vcf_btn_bcf)
                dpg.add_text("", tag="bcf_vcf_path_lbl", wrap=500)

                dpg.add_spacer(height=6)
                fasta_btn = dpg.add_button(
                    label="Choose reference FASTA (for left-align)",
                    callback=lambda: dpg.show_item("file_dialog_fasta"),
                    width=220,
                    tag="tooltip_fa_bcf",
                )
                app._secondary_buttons.append(fasta_btn)
                dpg.add_text("", tag="bcf_ref_path_lbl", wrap=500)

                dpg.add_spacer(height=6)
                reg_btn = dpg.add_button(
                    label="Choose regions BED (optional)",
                    callback=lambda: dpg.show_item("file_dialog_blacklist"),
                    width=220,
                    tag="tooltip_reg_bcf",
                )
                app._secondary_buttons.append(reg_btn)
                dpg.add_text("", tag="bcf_regions_path_lbl", wrap=500)

            with dpg.group():
                app.bcf_split = dpg.add_checkbox(label="Split multiallelic", default_value=True)
                app._inputs.append(app.bcf_split)

                app.bcf_left = dpg.add_checkbox(label="Left-align indels (needs FASTA)", default_value=True)
                app._inputs.append(app.bcf_left)

                app.bcf_sort = dpg.add_checkbox(label="Sort", default_value=True)
                app._inputs.append(app.bcf_sort)

                app.bcf_setid = dpg.add_checkbox(label="Set ID to CHR:POS:REF:ALT", default_value=True)
                app._inputs.append(app.bcf_setid)

                app.bcf_compr = dpg.add_checkbox(label="Compress output (.vcf.gz)", default_value=True)
                app._inputs.append(app.bcf_compr)

                app.bcf_index = dpg.add_checkbox(label="Index output (tabix)", default_value=True)
                app._inputs.append(app.bcf_index)

                app.bcf_rmflt = dpg.add_checkbox(label="Keep only PASS (remove filtered)", default_value=False)
                app._inputs.append(app.bcf_rmflt)

                app.bcf_filltags = dpg.add_checkbox(
                    label="Fill tags (AC, AN, AF, MAF, HWE) before filtering",
                    default_value=True,
                )
                app._inputs.append(app.bcf_filltags)

                dpg.add_spacer(height=6)
                app.bcf_filter_expr = dpg.add_input_text(
                    label="bcftools filter expression (optional)",
                    hint="Example: QUAL>=30 && INFO/DP>=10",
                    width=320,
                )
                app._inputs.append(app.bcf_filter_expr)

                dpg.add_spacer(height=6)
                app.bcf_out_prefix = dpg.add_input_text(
                    label="Output prefix (optional)",
                    hint="Leave empty to auto-generate next to the input file",
                    width=320,
                )
                app._inputs.append(app.bcf_out_prefix)

                dpg.add_spacer(height=6)
                app.bcf_make_snps = dpg.add_checkbox(label="Produce SNP-only VCF", default_value=False)
                app._inputs.append(app.bcf_make_snps)

                app.bcf_make_svs = dpg.add_checkbox(label="Produce SV-only VCF", default_value=False)
                app._inputs.append(app.bcf_make_svs)

                dpg.add_spacer(height=12)
                run_bcf = dpg.add_button(
                    tag="btn_bcf_pre",
                    label="Run bcftools preprocess",
                    callback=lambda: _on_run_bcftools_preprocess(app),
                    width=240,
                    height=38,
                )
                app._primary_buttons.append(run_bcf)

                with dpg.tooltip(run_bcf):
                    dpg.add_text(
                        "Normalize and clean the selected VCF using bcftools "
                        "(split multiallelic, left-align, sort, filter, compress, index)."
                    )

                dpg.add_loading_indicator(
                    tag="bcf_pre_spinner",
                    radius=6,
                    style=2,
                    show=False,
                )

    return "page_pre_bcf"


def _on_run_bcftools_preprocess(app):
    dpg.configure_item("bcf_pre_spinner", show=True)
    dpg.configure_item("btn_bcf_pre", enabled=False)

    try:
        app.run_bcftools_preprocess()
    except Exception as exc:
        app.add_log(f"[BCF] Error during bcftools preprocess: {exc}", error=True)
        app.add_log(traceback.format_exc(), error=True)
    finally:
        dpg.configure_item("bcf_pre_spinner", show=False)
        dpg.configure_item("btn_bcf_pre", enabled=True)


def page_check_vcf(app, parent):
    with dpg.group(parent=parent, show=False, tag="page_check_vcf"):
        dpg.add_text("\nCheck VCF quality before conversion/analysis", indent=10)
        dpg.add_spacer(height=10)

        with dpg.group(horizontal=True, horizontal_spacing=60):
            # -------- LEFT: file selectors --------
            with dpg.group():
                vcf_btn_qc = dpg.add_button(
                    label="Choose a VCF file",
                    callback=lambda: dpg.show_item("file_dialog_vcf"),
                    width=220,
                    tag="btn_qc_vcf",
                )
                app._secondary_buttons.append(vcf_btn_qc)
                dpg.add_text("", tag="qc_vcf_path_lbl", wrap=500)

                with dpg.tooltip(vcf_btn_qc):
                    dpg.add_text("Main VCF file to be checked for basic quality issues.")

                dpg.add_spacer(height=6)
                vcf_btn_qc2 = dpg.add_button(
                    label="Choose another VCF (optional)",
                    callback=lambda: dpg.show_item("file_dialog_vcf2"),
                    width=220,
                    tag="btn_qc_vcf2",
                )
                app._secondary_buttons.append(vcf_btn_qc2)
                dpg.add_text("", tag="qc_vcf2_path_lbl", wrap=500)

                with dpg.tooltip(vcf_btn_qc2):
                    dpg.add_text("Optional second VCF to compare against the first one.")

                dpg.add_spacer(height=6)
                bl_btn = dpg.add_button(
                    label="Choose blacklist BED (optional)",
                    callback=lambda: dpg.show_item("file_dialog_blacklist"),
                    width=220,
                    tag="btn_qc_bl",
                )
                app._secondary_buttons.append(bl_btn)
                dpg.add_text("", tag="qc_bl_path_lbl", wrap=500)

                with dpg.tooltip(bl_btn):
                    dpg.add_text("Optional BED file with regions to be treated as blacklist.")

            # -------- RIGHT: options + run --------
            with dpg.group():
                app.deep_scan = dpg.add_checkbox(
                    label="Deep scan",
                    default_value=False,
                    tag="qc_deep_scan",
                )
                app._inputs.append(app.deep_scan)

                with dpg.tooltip(app.deep_scan):
                    dpg.add_text(
                        "If enabled, performs more detailed checks (slower but more thorough)."
                    )

                dpg.add_spacer(height=12)

                run_qc_btn = dpg.add_button(
                    label="Run Quality Check",
                    tag="btn_run_vcf_qc",
                    callback=lambda: _on_run_vcf_qc(app),
                    width=200,
                    height=36,
                )
                app._primary_buttons.append(run_qc_btn)

                with dpg.tooltip(run_qc_btn):
                    dpg.add_text(
                        "Run bcftools-based quality checks on the selected VCF file(s)."
                    )

                dpg.add_loading_indicator(
                    tag="vcf_qc_spinner",
                    radius=6,
                    style=2,
                    show=False,
                )

    return "page_check_vcf"


def _on_run_vcf_qc(app):
    """Wrapper to show loader while VCF QC is running."""
    dpg.configure_item("vcf_qc_spinner", show=True)
    dpg.configure_item("btn_run_vcf_qc", enabled=False)

    try:
        app.run_vcf_qc()
    except Exception as exc:
        app.add_log(f"[VCF-QC] Error during quality check: {exc}", error=True)
        app.add_log(traceback.format_exc(), error=True)
    finally:
        dpg.configure_item("vcf_qc_spinner", show=False)
        dpg.configure_item("btn_run_vcf_qc", enabled=True)


def page_convert_plink(app, parent):
    if app is None:
        raise ValueError("app must not be None in page_convert_plink")

    with dpg.group(parent=parent, show=False, tag="page_plink"):
        dpg.add_text(
            "\nConvert a VCF file into PLINK BED and apply MAF/missing genotype filters.",
            indent=10,
        )
        dpg.add_spacer(height=10)

        with dpg.group(horizontal=True, horizontal_spacing=60):
            with dpg.group():
                dpg.add_text("Select files:", indent=0)

                vcf = dpg.add_button(
                    label="Choose a VCF file",
                    callback=lambda: dpg.show_item("file_dialog_vcf"),
                    width=220,
                    tag="btn_plink_vcf",
                )
                app._secondary_buttons.append(vcf)
                dpg.add_text("", tag="conv_vcf_path_lbl", wrap=500)

                with dpg.tooltip(vcf):
                    dpg.add_text("Main VCF file that will be converted to PLINK BED.")

                dpg.add_spacer(height=6)
                variant_ids = dpg.add_button(
                    label="Choose IDs file (optional)",
                    callback=lambda: dpg.show_item("file_dialog_variants"),
                    width=220,
                    tag="btn_plink_ids",
                )
                app._secondary_buttons.append(variant_ids)
                dpg.add_text("", tag="ids_path_lbl", wrap=500)

                with dpg.tooltip(variant_ids):
                    dpg.add_text(
                        "Optional list of variant IDs to keep during conversion."
                    )

            with dpg.group():
                dpg.add_text("Apply filters:", indent=0)

                maf_input = dpg.add_input_float(
                    label="Minor allele frequency (MAF)",
                    width=220,
                    default_value=0.05,
                    step=0.005,
                    tag="plink_maf_input",
                )
                app._inputs.append(maf_input)

                with dpg.tooltip(maf_input):
                    dpg.add_text(
                        "Variants with MAF below this threshold will be removed."
                    )

                dpg.add_spacer(height=6)
                geno_input = dpg.add_input_float(
                    label="Missing genotype rate",
                    width=220,
                    default_value=0.10,
                    step=0.005,
                    tag="plink_missing_input",
                )
                app._inputs.append(geno_input)

                with dpg.tooltip(geno_input):
                    dpg.add_text(
                        "Variants with missing genotype rate above this value are filtered out."
                    )

                dpg.add_spacer(height=14)
                convert_btn = dpg.add_button(
                    tag="convert_vcf_btn",
                    label="Convert VCF",
                    callback=app.convert_vcf,
                    width=160,
                    height=36,
                    enabled=True,
                )

                app._primary_buttons.append(convert_btn)

                with dpg.tooltip(convert_btn):
                    dpg.add_text(
                        "Convert the selected VCF to PLINK BED applying the filters above."
                    )

                dpg.add_loading_indicator(
                    tag="plink_convert_spinner",
                    radius=6,
                    style=2,
                    show=False,
                )

    return "page_plink"


def _on_convert_vcf(app, sender, app_data, user_data):
    """Wrapper to run app.convert_vcf with UI feedback (spinner + button disable)."""
    try:
        dpg.configure_item("plink_convert_spinner", show=True)
        dpg.configure_item("convert_vcf_btn", enabled=False)
    except Exception:
        pass

    try:
        if app is None:
            raise RuntimeError("app is None in _on_convert_vcf")

        if not hasattr(app, "convert_vcf"):
            raise AttributeError("App object has no method 'convert_vcf'")

        app.convert_vcf(sender, app_data, user_data)

    except Exception as exc:
        msg = f"[PLINK] Error while converting VCF: {exc}"

        if app is not None and hasattr(app, "add_log"):
            try:
                app.add_log(msg, error=True)
                app.add_log(traceback.format_exc(), error=True)
            except Exception:
                print(msg)
                print(traceback.format_exc())
        else:
            print(msg)
            print(traceback.format_exc())

    finally:
        try:
            dpg.configure_item("plink_convert_spinner", show=False)
            dpg.configure_item("convert_vcf_btn", enabled=True)
        except Exception:
            pass

def page_ld_analysis(app, parent):
    with dpg.group(parent=parent, show=False, tag="page_ld"):
        dpg.add_text("\nLD analysis: LD decay, LD heatmap, and diversity metrics", indent=10)
        dpg.add_spacer(height=10)

        with dpg.group(horizontal=True, horizontal_spacing=60):
            # LEFT: input files
            with dpg.group():
                dpg.add_text("Input files", color=(200, 220, 200))
                dpg.add_spacer(height=6)

                btn_bed = dpg.add_button(
                    label="Choose PLINK BED",
                    callback=lambda: dpg.show_item("file_dialog_bed"),
                    width=220,
                    tag="ld_btn_bed",
                )
                app._secondary_buttons.append(btn_bed)
                dpg.add_text("", tag="ld_bed_path_lbl", wrap=520)

                with dpg.tooltip(btn_bed):
                    dpg.add_text("PLINK BED/BIM/FAM set that will be used for LD analysis.")

                dpg.add_spacer(height=6)
                btn_vcf = dpg.add_button(
                    label="Choose VCF (optional)",
                    callback=lambda: dpg.show_item("file_dialog_vcf"),
                    width=220,
                    tag="ld_btn_vcf",
                )
                app._secondary_buttons.append(btn_vcf)
                dpg.add_text("", tag="ld_vcf_path_lbl", wrap=520)

                with dpg.tooltip(btn_vcf):
                    dpg.add_text("Optional VCF file; used only for some plots if provided.")

                dpg.add_spacer(height=6)
                app.ld_region = dpg.add_input_text(
                    label="Region (optional)",
                    hint="chr:start-end (e.g., 1:1000000-2000000)",
                    width=320,
                    tag="ld_region_input",
                )
                app._inputs.append(app.ld_region)

                with dpg.tooltip(app.ld_region):
                    dpg.add_text("Restrict LD analysis to this genomic interval if desired.")

            # RIGHT: options
            with dpg.group():
                dpg.add_text("Options", color=(200, 220, 200))
                dpg.add_spacer(height=6)

                app.ld_window_kb = dpg.add_input_int(
                    label="LD window (kb)",
                    default_value=500,
                    min_value=1,
                    min_clamped=True,
                    width=220,
                    tag="ld_window_kb_input",
                )
                app._inputs.append(app.ld_window_kb)

                with dpg.tooltip(app.ld_window_kb):
                    dpg.add_text("Physical window (in kb) around each SNP for LD calculation.")

                app.ld_window_snp = dpg.add_input_int(
                    label="LD window size (SNPs)",
                    default_value=5000,
                    min_value=10,
                    min_clamped=True,
                    width=220,
                    tag="ld_window_snp_input",
                )
                app._inputs.append(app.ld_window_snp)

                with dpg.tooltip(app.ld_window_snp):
                    dpg.add_text("Maximum number of SNPs considered in the LD window.")

                app.ld_max_kb = dpg.add_input_int(
                    label="Max distance (kb)",
                    default_value=1000,
                    min_value=1,
                    min_clamped=True,
                    width=220,
                    tag="ld_max_kb_input",
                )
                app._inputs.append(app.ld_max_kb)

                with dpg.tooltip(app.ld_max_kb):
                    dpg.add_text("Pairs of SNPs farther than this distance are ignored.")

                app.ld_min_r2 = dpg.add_input_float(
                    label="Min r²",
                    default_value=0.1,
                    min_value=0.0,
                    max_value=1.0,
                    min_clamped=True,
                    max_clamped=True,
                    step=0.05,
                    width=220,
                    tag="ld_min_r2_input",
                )
                app._inputs.append(app.ld_min_r2)

                with dpg.tooltip(app.ld_min_r2):
                    dpg.add_text("Only SNP pairs with r² at or above this threshold are kept.")

                dpg.add_spacer(height=6)
                app.ld_do_decay = dpg.add_checkbox(
                    label="Compute LD decay",
                    default_value=True,
                    tag="ld_decay_checkbox",
                )
                app._inputs.append(app.ld_do_decay)

                with dpg.tooltip(app.ld_do_decay):
                    dpg.add_text("Generate LD decay curves across distance.")

                app.ld_do_heatmap = dpg.add_checkbox(
                    label="LD heatmap",
                    default_value=True,
                    tag="ld_heatmap_checkbox",
                )
                app._inputs.append(app.ld_do_heatmap)

                with dpg.tooltip(app.ld_do_heatmap):
                    dpg.add_text("Produce LD heatmap figures for the selected region.")

                app.ld_do_div = dpg.add_checkbox(
                    label="Diversity metrics",
                    default_value=True,
                    tag="ld_div_checkbox",
                )
                app._inputs.append(app.ld_do_div)

                with dpg.tooltip(app.ld_do_div):
                    dpg.add_text("Compute nucleotide diversity and related metrics.")

                dpg.add_spacer(height=12)
                run_btn = dpg.add_button(
                    label="Run LD analysis",
                    tag="ld_run_btn",
                    callback=_on_run_ld_analysis,
                    user_data={"app": app},
                    width=240,
                    height=38,
                )

                app._primary_buttons.append(run_btn)

                with dpg.tooltip(run_btn):
                    dpg.add_text(
                        "Start LD analysis. If everything is OK, summary files and plots "
                        "will appear in the chosen output directory."
                    )

                dpg.add_loading_indicator(
                    tag="ld_run_spinner",
                    radius=6,
                    style=2,
                    show=False,
                )

    return "page_ld"


def _on_run_ld_analysis(sender, app_data, user_data):
    """
    Callback for 'Run LD analysis' button.
    user_data: {"app": app}
    """
    app = None
    if isinstance(user_data, dict):
        app = user_data.get("app")

    if app is None:
        print("[LD] ERROR: app is None in _on_run_ld_analysis")
        return

    try:
        dpg.configure_item("ld_run_spinner", show=True)
        dpg.configure_item("ld_run_btn", enabled=False)
    except Exception:
        pass

    try:
        app.run_ld_analysis(sender, app_data)
    except Exception as exc:
        try:
            app.add_log(f"[LD] Error while running LD analysis: {exc}", error=True)
            app.add_log(traceback.format_exc(), error=True)
        except Exception:
            print(f"[LD] Error while running LD analysis: {exc}")
            print(traceback.format_exc())
    finally:
        try:
            dpg.configure_item("ld_run_spinner", show=False)
            dpg.configure_item("ld_run_btn", enabled=True)
        except Exception:
            pass


def page_gwas(app, parent):
    with dpg.group(parent=parent, show=False, tag="page_gwas"):
        dpg.add_text("\nStart GWAS Analysis", indent=10)
        dpg.add_spacer(height=10)

        with dpg.group(horizontal=True, horizontal_spacing=60):
            # LEFT: file inputs
            with dpg.group():
                geno = dpg.add_button(
                    label="Choose a BED file",
                    callback=lambda: dpg.show_item("file_dialog_bed"),
                    width=220,
                    tag="gwas_btn_bed",
                )
                app._secondary_buttons.append(geno)
                dpg.add_text("", tag="gwas_bed_path_lbl", wrap=500)

                with dpg.tooltip(geno):
                    dpg.add_text("Genotype PLINK files (BED/BIM/FAM) used as GWAS input.")

                dpg.add_spacer(height=6)
                pheno = dpg.add_button(
                    label="Choose a phenotype file",
                    callback=lambda: dpg.show_item("file_dialog_pheno"),
                    width=220,
                    tag="gwas_btn_pheno",
                )
                app._secondary_buttons.append(pheno)
                dpg.add_text("", tag="gwas_pheno_path_lbl", wrap=500)

                with dpg.tooltip(pheno):
                    dpg.add_text("Phenotype file (one or more traits) for the samples.")

                dpg.add_spacer(height=6)
                cov_file = dpg.add_button(
                    label="Choose covariate file (optional)",
                    callback=lambda: dpg.show_item("file_dialog_cov"),
                    width=220,
                    tag="gwas_btn_cov",
                )
                app._secondary_buttons.append(cov_file)
                dpg.add_text("", tag="gwas_cov_path_lbl", wrap=500)

                with dpg.tooltip(cov_file):
                    dpg.add_text("Optional covariates (e.g. PCs, batch, environment).")

                dpg.add_spacer(height=6)
                kin_file = dpg.add_button(
                    label="Choose kinship file (optional)",
                    callback=lambda: dpg.show_item("file_dialog_kinship"),
                    width=220,
                    tag="gwas_btn_kin",
                )
                app._secondary_buttons.append(kin_file)
                dpg.add_text("", tag="gwas_kinship_path_lbl", wrap=500)

                with dpg.tooltip(kin_file):
                    dpg.add_text("Pre-computed kinship/GRM matrix for mixed models.")

                dpg.add_spacer(height=6)
                kin_compute = dpg.add_button(
                    label="Compute kinship from BED",
                    callback=lambda: app.compute_kinship_from_bed(),
                    width=220,
                    tag="gwas_btn_compute_kin",
                )
                app._secondary_buttons.append(kin_compute)

                with dpg.tooltip(kin_compute):
                    dpg.add_text("Estimate kinship matrix directly from the selected BED.")

                dpg.add_spacer(height=8)
                anno_btn = dpg.add_button(
                    label="Choose GTF/GFF annotation",
                    callback=lambda: dpg.show_item("file_dialog_gtf"),
                    width=220,
                    tag="gwas_btn_gtf",
                )
                app._secondary_buttons.append(anno_btn)
                dpg.add_text("", tag="gwas_gtf_path_lbl", wrap=500)

                with dpg.tooltip(anno_btn):
                    dpg.add_text("Genome annotation file for mapping SNPs to genes.")

            # RIGHT: options and run
            with dpg.group():
                app.gwas_combo = dpg.add_combo(
                    label="Analysis Algorithms",
                    items=[
                        "FaST-LMM",
                        "Linear regression",
                        "Ridge Regression",
                        "Random Forest (AI)",
                        "XGBoost (AI)",
                        "GLM (PLINK2)",
                        "SAIGE (mixed model)",
                    ],
                    width=260,
                    default_value="FaST-LMM",
                    tag="gwas_algorithm_combo",
                )
                app._inputs.append(app.gwas_combo)

                with dpg.tooltip(app.gwas_combo):
                    dpg.add_text("Choose the GWAS method to run on the selected data.")

                dpg.add_spacer(height=8)
                app.snp_limit = dpg.add_input_int(
                    label="Limit SNPs in plots (optional)",
                    width=260,
                    min_value=0,
                    step=1000,
                    default_value=0,
                    tag="gwas_snp_limit",
                )
                app._inputs.append(app.snp_limit)

                with dpg.tooltip(app.snp_limit):
                    dpg.add_text("Optional cap on number of SNPs shown in plots (0 = no limit).")

                dpg.add_spacer(height=6)
                app.plot_stats = dpg.add_checkbox(
                    label="Produce pheno/geno statistics (PDF)",
                    default_value=False,
                    tag="gwas_plot_stats",
                )

                with dpg.tooltip(app.plot_stats):
                    dpg.add_text("Generate summary statistics for phenotypes and genotypes.")

                dpg.add_spacer(height=6)
                app.annotate_enable = dpg.add_checkbox(
                    label="Annotate GWAS with GTF/GFF",
                    default_value=False,
                    tag="gwas_annotate_checkbox",
                )

                with dpg.tooltip(app.annotate_enable):
                    dpg.add_text("If enabled, top hits will be annotated with nearby genes.")

                app.annotate_window_kb = dpg.add_input_int(
                    label="Window around TSS (kb)",
                    width=260,
                    min_value=0,
                    step=10,
                    default_value=50,
                    tag="gwas_annot_window",
                )

                with dpg.tooltip(app.annotate_window_kb):
                    dpg.add_text("Distance around gene TSS used when assigning SNPs to genes.")

                dpg.add_spacer(height=10)
                dpg.add_text("Optional Region Filter", color=(150, 150, 255))

                app.region_chr = dpg.add_input_text(
                    label="Chromosome",
                    width=260,
                    default_value="",
                    tag="gwas_region_chr",
                )
                with dpg.tooltip(app.region_chr):
                    dpg.add_text("Restrict analysis to a single chromosome (optional).")

                app.region_start = dpg.add_input_int(
                    label="Start position",
                    width=260,
                    min_value=0,
                    step=1000,
                    default_value=0,
                    tag="gwas_region_start",
                )
                with dpg.tooltip(app.region_start):
                    dpg.add_text("Start genomic coordinate for region-based GWAS (optional).")

                app.region_end = dpg.add_input_int(
                    label="End position",
                    width=260,
                    min_value=0,
                    step=1000,
                    default_value=0,
                    tag="gwas_region_end",
                )
                with dpg.tooltip(app.region_end):
                    dpg.add_text("End genomic coordinate for region-based GWAS (optional).")

                dpg.add_spacer(height=14)
                gwas_btn = dpg.add_button(
                    label="Run GWAS",
                    callback=_on_run_gwas,
                    user_data={"app": app, "controls": [geno, pheno]},
                    width=200,
                    height=36,
                    tag="gwas_run_btn",
                )

                app._primary_buttons.append(gwas_btn)

                with dpg.tooltip(gwas_btn):
                    dpg.add_text(
                        "Start the GWAS pipeline. If everything is configured correctly, "
                        "results and plots will be written to the output directory."
                    )

                dpg.add_loading_indicator(
                    tag="gwas_run_spinner",
                    radius=6,
                    style=2,
                    show=False,
                )

        dpg.add_spacer(height=12)
        dpg.add_separator()
        dpg.add_spacer(height=8)
        dpg.add_text(
            "Results will appear in the Results window (tables, Manhattan/QQ plots).",
            wrap=900,
        )

    return "page_gwas"


def _on_run_gwas(sender, app_data, user_data):
    app = None
    controls = None

    if isinstance(user_data, dict):
        app = user_data.get("app")
        controls = user_data.get("controls")

    if app is None:
        print("[GWAS] ERROR: app is None inside _on_run_gwas. Check button user_data wiring.")
        return

    dpg.configure_item("gwas_run_spinner", show=True)
    dpg.configure_item("gwas_run_btn", enabled=False)

    try:
        app.run_gwas(sender, app_data, controls)
    except Exception as exc:
        app.add_log(f"[GWAS] Error while running GWAS: {exc}", error=True)
        app.add_log(traceback.format_exc(), error=True)
    finally:
        dpg.configure_item("gwas_run_spinner", show=False)
        dpg.configure_item("gwas_run_btn", enabled=True)


def page_pca(app, parent):
    with dpg.group(parent=parent, show=False, tag="page_pca"):
        dpg.add_text("Population Structure (PCA) & Kinship", indent=10)
        dpg.add_spacer(height=10)

        with dpg.group(horizontal=True, horizontal_spacing=60):
            # Left side: inputs
            with dpg.group():
                btn_bed = dpg.add_button(
                    label="Choose PLINK BED",
                    callback=lambda: dpg.show_item("file_dialog_bed"),
                    width=220,
                    tag="pca_btn_bed",
                )
                app._secondary_buttons.append(btn_bed)
                dpg.add_text("", tag="pca_bed_path_lbl", wrap=500)

                with dpg.tooltip(btn_bed):
                    dpg.add_text("Select the PLINK BED/BIM/FAM set used for PCA/kinship.")

                dpg.add_spacer(height=8)
                app.pca_out_prefix = dpg.add_input_text(
                    label="Output prefix (optional)",
                    hint="Default: next to BED prefix",
                    width=320,
                    tag="pca_out_prefix",
                )
                app._inputs.append(app.pca_out_prefix)

                with dpg.tooltip(app.pca_out_prefix):
                    dpg.add_text("Optional prefix for PCA and kinship output files.")

            # Right side: options + run
            with dpg.group():
                app.pca_npcs = dpg.add_input_int(
                    label="Number of PCs",
                    default_value=10,
                    min_value=2,
                    max_value=50,
                    width=160,
                    tag="pca_npcs",
                )
                app._inputs.append(app.pca_npcs)

                with dpg.tooltip(app.pca_npcs):
                    dpg.add_text("How many principal components to compute (2–50).")

                app.pca_kinship = dpg.add_checkbox(
                    label="Compute kinship matrix",
                    default_value=True,
                    tag="pca_kinship",
                )
                app._inputs.append(app.pca_kinship)

                with dpg.tooltip(app.pca_kinship):
                    dpg.add_text("If enabled, a kinship/GRM matrix will also be generated.")

                dpg.add_spacer(height=12)
                run_btn = dpg.add_button(
                    label="Run PCA",
                    tag="pca_run_btn",
                    callback=_on_run_pca,
                    user_data={"app": app},
                    width=200,
                    height=36,
                )

                app._primary_buttons.append(run_btn)

                with dpg.tooltip(run_btn):
                    dpg.add_text(
                        "Start PCA on the selected BED file. Results will be saved and "
                        "previewed in the Results window."
                    )

                dpg.add_loading_indicator(
                    tag="pca_run_spinner",
                    radius=6,
                    style=2,
                    show=False,
                )

        dpg.add_spacer(height=12)
        dpg.add_separator()
        dpg.add_spacer(height=6)
        # dpg.add_text(
        #     "Results preview will appear in the Results window (PC1 vs PC2 plot, files list)."
        # )

    return "page_pca"


def _on_run_pca(sender, app_data, user_data):
    app = user_data.get("app")
    if app is None:
        print("[PCA] app is None in _on_run_pca")
        return

    dpg.configure_item("pca_run_spinner", show=True)
    dpg.configure_item("pca_run_btn", enabled=False)
    try:
        app.run_pca_module(sender, app_data)
    except Exception as exc:
        app.add_log(f"[PCA] Error while running PCA: {exc}", error=True)
        app.add_log(traceback.format_exc(), error=True)
    finally:
        dpg.configure_item("pca_run_spinner", show=False)
        dpg.configure_item("pca_run_btn", enabled=True)


def page_genomic_prediction(app, parent):
    with dpg.group(parent=parent, show=False, tag="page_gp"):
        dpg.add_text("\nStart Genomic Prediction", indent=10)
        dpg.add_spacer(height=10)

        with dpg.group(horizontal=True, horizontal_spacing=60):
            with dpg.group():
                geno = dpg.add_button(
                    label="Choose a BED file",
                    callback=lambda: dpg.show_item("file_dialog_bed"),
                    width=220,
                    tag="gp_bed_btn",
                )
                app._secondary_buttons.append(geno)
                dpg.add_text("", tag="gp_bed_path_lbl", wrap=500)

                with dpg.tooltip(geno):
                    dpg.add_text("Select the PLINK BED/BIM/FAM set for genomic prediction.")

                dpg.add_spacer(height=6)
                pheno = dpg.add_button(
                    label="Choose a phenotype file",
                    callback=lambda: dpg.show_item("file_dialog_pheno"),
                    width=220,
                    tag="gp_pheno_btn",
                )
                app._secondary_buttons.append(pheno)
                dpg.add_text("", tag="gp_pheno_path_lbl", wrap=500)

                with dpg.tooltip(pheno):
                    dpg.add_text("Select the phenotype file matching the individuals in the BED set.")

            with dpg.group():
                app.gwas_gp = dpg.add_combo(
                    label="Analysis Algorithms",
                    items=[
                        "XGBoost (AI)",
                        "Random Forest (AI)",
                        "Ridge Regression",
                        "GP_LMM",
                        "val",
                    ],
                    width=240,
                    default_value="XGBoost (AI)",
                    tag="gp_algorithm_combo",
                )
                app._inputs.append(app.gwas_gp)

                with dpg.tooltip(app.gwas_gp):
                    dpg.add_text("Choose the model used to predict phenotypes from genotypes.")

                dpg.add_spacer(height=14)
                gp_btn = dpg.add_button(
                    label="Run Genomic Prediction",
                    tag="gp_run_btn",
                    callback=_on_run_genomic_prediction,
                    user_data={"app": app, "deps": [geno, pheno]},
                    width=200,
                    height=36,
                )

                app._primary_buttons.append(gp_btn)

                with dpg.tooltip(gp_btn):
                    dpg.add_text("Train and evaluate the selected prediction model on the provided data.")

                dpg.add_loading_indicator(
                    tag="gp_run_spinner",
                    radius=6,
                    style=2,
                    thickness=2,
                    show=False,
                )

    return "page_gp"


def _on_run_genomic_prediction(sender, app_data, user_data):
    """
    DearPyGui callback for the 'Run Genomic Prediction' button.
    user_data should be a dict containing:
        - "app":  main application instance
        - "deps": optional dependencies (e.g. [geno_button, pheno_button])
    """
    app = None
    deps = None

    if isinstance(user_data, dict):
        app = user_data.get("app")
        deps = user_data.get("deps")

    if app is None:
        print("[GenomicPrediction] ERROR: app is None in _on_run_genomic_prediction")
        return

    # show spinner / disable button if tags exist
    try:
        dpg.configure_item("gp_run_spinner", show=True)
        dpg.configure_item("gp_run_btn", enabled=False)
    except Exception:
        pass

    try:
        # call the real implementation
        app.run_genomic_prediction(sender, app_data, deps)
    except Exception as exc:
        try:
            app.add_log(f"[GenomicPrediction] Error: {exc}", error=True)
            app.add_log(traceback.format_exc(), error=True)
        except Exception:
            print(f"[GenomicPrediction] Error: {exc}")
            print(traceback.format_exc())
    finally:
        try:
            dpg.configure_item("gp_run_spinner", show=False)
            dpg.configure_item("gp_run_btn", enabled=True)
        except Exception:
            pass


def page_batch_gwas(app, parent):
    with dpg.group(parent=parent, show=False, tag="page_batch"):
        dpg.add_text(
            "\nBatch GWAS for all traits in a phenotype file (FID IID + multiple traits).",
            indent=10,
        )
        dpg.add_spacer(height=10)

        with dpg.group(horizontal=True, horizontal_spacing=60):
            with dpg.group():
                geno_btn = dpg.add_button(
                    label="Choose a BED file (SNP or SV)",
                    callback=lambda: dpg.show_item("file_dialog_bed"),
                    width=240,
                    tag="batch_bed_btn",
                )
                app._secondary_buttons.append(geno_btn)
                dpg.add_text("", tag="batch_bed_path_lbl", wrap=520)

                with dpg.tooltip(geno_btn):
                    dpg.add_text("Select the PLINK BED/BIM/FAM set for batch GWAS.")

                dpg.add_spacer(height=6)
                pheno_btn = dpg.add_button(
                    label="Choose a multi-trait phenotype file",
                    callback=lambda: dpg.show_item("file_dialog_pheno"),
                    width=240,
                    tag="batch_pheno_btn",
                )
                app._secondary_buttons.append(pheno_btn)
                dpg.add_text("", tag="batch_pheno_path_lbl", wrap=520)

                with dpg.tooltip(pheno_btn):
                    dpg.add_text("Phenotype file with FID/IID and multiple trait columns.")

                dpg.add_spacer(height=6)
                cov_btn = dpg.add_button(
                    label="Choose covariates (optional)",
                    callback=lambda: dpg.show_item("file_dialog_cov"),
                    width=240,
                    tag="batch_cov_btn",
                )
                app._secondary_buttons.append(cov_btn)
                dpg.add_text("", tag="batch_cov_path_lbl", wrap=520)

                with dpg.tooltip(cov_btn):
                    dpg.add_text("Optional covariates to include in all batch GWAS runs.")

            with dpg.group():
                app.batch_algo = dpg.add_combo(
                    label="Algorithm",
                    items=[
                        "FaST-LMM",
                        "Linear regression",
                        "Ridge Regression",
                        "Random Forest (AI)",
                        "XGBoost (AI)",
                    ],
                    width=260,
                    default_value="FaST-LMM",
                    tag="batch_algo_combo",
                )
                app._inputs.append(app.batch_algo)

                with dpg.tooltip(app.batch_algo):
                    dpg.add_text("Choose the GWAS method used for each trait in the batch.")

                dpg.add_spacer(height=8)
                dpg.add_text(
                    "Uses settings from 'Settings' page (trees, depth, train %, jobs ...)."
                )

                dpg.add_spacer(height=12)
                run_btn = dpg.add_button(
                    label="Run Batch GWAS",
                    callback=lambda s, a, _app=app: _on_run_batch_gwas(_app, s, a),
                    width=220,
                    height=38,
                    tag="batch_run_btn",
                )
                app._primary_buttons.append(run_btn)

                with dpg.tooltip(run_btn):
                    dpg.add_text("Launch GWAS for all traits in the phenotype file.")

                dpg.add_loading_indicator(
                    tag="batch_run_spinner",
                    radius=6,
                    style=2,
                    thickness=2,
                    show=False,
                )

    return "page_batch"


def _on_run_batch_gwas(app, sender, app_data):
    dpg.configure_item("batch_run_spinner", show=True)
    dpg.configure_item("batch_run_btn", enabled=False)
    try:
        app.run_batch_gwas_ui(sender, app_data)
    except Exception as exc:
        app.add_log(f"[BatchGWAS] Error: {exc}", error=True)
        app.add_log(traceback.format_exc(), error=True)
    finally:
        dpg.configure_item("batch_run_spinner", show=False)
        dpg.configure_item("batch_run_btn", enabled=True)


def page_settings(app, parent):
    with dpg.group(parent=parent, show=False, tag="page_settings"):
        dpg.add_spacer(height=10)

        with dpg.table(
            header_row=False,
            borders_innerH=False,
            borders_outerH=False,
            borders_innerV=False,
            borders_outerV=False,
            resizable=False,
        ):
            dpg.add_table_column()
            dpg.add_table_column()

            with dpg.table_row():
                # General settings
                with dpg.group():
                    dpg.add_text("General Settings", color=(200, 180, 90))
                    dpg.add_spacer(height=8)

                    dark_toggle = dpg.add_checkbox(
                        label="Night Mode (Dark)",
                        default_value=True,
                        tag="settings_dark_toggle",
                    )
                    with dpg.tooltip(dark_toggle):
                        dpg.add_text(
                            "Switch between dark/light UI themes.\n"
                            "No impact on performance, only on appearance and eye comfort."
                        )

                    dpg.add_spacer(height=6)
                    app.nr_jobs = dpg.add_input_int(
                        label="Number of jobs to run",
                        width=220,
                        default_value=-1,
                        step=1,
                        min_value=-1,
                        max_value=50,
                        min_clamped=True,
                        max_clamped=True,
                        tag="tooltip_nr_jobs",
                    )
                    app._inputs.append(app.nr_jobs)
                    with dpg.tooltip(app.nr_jobs):
                        dpg.add_text(
                            "Maximum number of analysis jobs to run in parallel.\n"
                            "-1 = auto (uses available CPU cores).\n"
                            "Higher values can speed up workflows but increase CPU usage."
                        )

                    app.gb_goal = dpg.add_input_int(
                        label="Gigabytes of memory per run",
                        width=220,
                        default_value=0,
                        step=4,
                        min_value=0,
                        max_value=512,
                        min_clamped=True,
                        max_clamped=True,
                        tag="tooltip_gb_goal",
                    )
                    app._inputs.append(app.gb_goal)
                    with dpg.tooltip(app.gb_goal):
                        dpg.add_text(
                            "Approximate RAM budget per run in GB.\n"
                            "Use this to avoid overloading machines with limited memory.\n"
                            "0 = no explicit limit (fastest but less safe on small systems)."
                        )

                    app.plot_stats = dpg.add_checkbox(
                        label="Advanced Plotting",
                        default_value=False,
                        tag="tooltip_stats",
                    )
                    app._inputs.append(app.plot_stats)
                    with dpg.tooltip(app.plot_stats):
                        dpg.add_text(
                            "Enable richer plots and extra summary statistics.\n"
                            "This may increase runtime and memory for very large datasets."
                        )

                    app.snp_limit = dpg.add_input_text(
                        label="SNP limit",
                        width=220,
                        default_value="",
                        tag="tooltip_limit",
                    )
                    app._inputs.append(app.snp_limit)
                    with dpg.tooltip(app.snp_limit):
                        dpg.add_text(
                            "Optional hard cap on the number of SNPs used in plots/AI.\n"
                            "Leave empty for all SNPs.\n"
                            "Lower values reduce runtime and memory for extremely dense data."
                        )

                # ML settings
                with dpg.group():
                    dpg.add_text("Machine Learning Settings", color=(200, 180, 90))
                    dpg.add_spacer(height=8)

                    app.train_size_set = dpg.add_input_int(
                        label="Training size %",
                        width=220,
                        default_value=70,
                        step=10,
                        min_value=0,
                        max_value=100,
                        min_clamped=True,
                        max_clamped=True,
                        tag="tooltip_training",
                    )
                    app._inputs.append(app.train_size_set)
                    with dpg.tooltip(app.train_size_set):
                        dpg.add_text(
                            "Percentage of samples used for training models.\n"
                            "Higher training % can improve accuracy but leaves fewer samples for testing."
                        )

                    app.estim_set = dpg.add_input_int(
                        label="Number of trees",
                        width=220,
                        default_value=200,
                        step=10,
                        min_value=1,
                        min_clamped=True,
                        tag="tooltip_trees",
                    )
                    app._inputs.append(app.estim_set)
                    with dpg.tooltip(app.estim_set):
                        dpg.add_text(
                            "Number of trees/estimators for tree-based models (RF/XGBoost).\n"
                            "More trees usually improve stability and accuracy but increase runtime and memory."
                        )

                    app.max_dep_set = dpg.add_input_int(
                        label="Max depth",
                        width=220,
                        default_value=3,
                        step=10,
                        min_value=0,
                        max_value=100,
                        min_clamped=True,
                        max_clamped=True,
                        tag="tooltip_depth",
                    )
                    app._inputs.append(app.max_dep_set)
                    with dpg.tooltip(app.max_dep_set):
                        dpg.add_text(
                            "Maximum depth of each tree.\n"
                            "Shallow trees are faster and less prone to overfitting.\n"
                            "Very deep trees capture complex patterns but are slower and heavier."
                        )

                    app.model_nr = dpg.add_input_int(
                        label="Nr. of models",
                        width=220,
                        default_value=1,
                        step=1,
                        min_value=1,
                        max_value=50,
                        min_clamped=True,
                        tag="tooltip_model",
                    )
                    app._inputs.append(app.model_nr)
                    with dpg.tooltip(app.model_nr):
                        dpg.add_text(
                            "Number of independent models to train and ensemble.\n"
                            "Higher values smooth results but multiply runtime."
                        )

                    app.aggregation_method = dpg.add_combo(
                        ("sum", "median", "mean"),
                        label="Aggregation Method",
                        width=220,
                        default_value="sum",
                        tag="tooltip_aggr",
                    )
                    app._inputs.append(app.aggregation_method)
                    with dpg.tooltip(app.aggregation_method):
                        dpg.add_text(
                            "How predictions from multiple models are combined.\n"
                            "Median/mean can make ensembles more robust to outliers."
                        )

        dpg.add_spacer(height=18)
        dpg.add_text("Large-file handling", color=(200, 180, 90))
        dpg.add_spacer(height=8)

        with dpg.group(horizontal=True, horizontal_spacing=60):
            with dpg.group():
                app.large_enable = dpg.add_checkbox(
                    label="Enable Large-file mode",
                    default_value=True,
                    tag="large_enable",
                )
                app._inputs.append(app.large_enable)
                with dpg.tooltip(app.large_enable):
                    dpg.add_text(
                        "Use chunked processing for very large VCF/BCF files.\n"
                        "Reduces peak memory usage at the cost of extra I/O and slightly longer runtimes."
                    )

                dpg.add_spacer(height=6)
                app.large_chunk_lines = dpg.add_input_int(
                    label="Chunk size (VCF lines per part)",
                    width=260,
                    default_value=500_000,
                    step=50_000,
                    min_value=10_000,
                    min_clamped=True,
                    tag="large_chunk_lines",
                )
                app._inputs.append(app.large_chunk_lines)
                with dpg.tooltip(app.large_chunk_lines):
                    dpg.add_text(
                        "Number of variant lines per chunk when splitting large VCFs.\n"
                        "Smaller chunks lower memory usage but create more temporary files and overhead."
                    )

                dpg.add_spacer(height=6)
                app.large_max_workers = dpg.add_input_int(
                    label="Max workers",
                    width=260,
                    default_value=2,
                    min_value=1,
                    max_value=64,
                    min_clamped=True,
                    max_clamped=True,
                    tag="large_max_workers",
                )
                app._inputs.append(app.large_max_workers)
                with dpg.tooltip(app.large_max_workers):
                    dpg.add_text(
                        "Maximum number of worker processes used for chunked jobs.\n"
                        "Higher values speed up processing but increase CPU and disk contention."
                    )

            with dpg.group():
                app.large_merge_strategy = dpg.add_combo(
                    items=["bcftools", "cat"],
                    label="Merge strategy",
                    width=260,
                    default_value="bcftools",
                    tag="large_merge_strategy",
                )
                app._inputs.append(app.large_merge_strategy)
                with dpg.tooltip(app.large_merge_strategy):
                    dpg.add_text(
                        "Tool/strategy for merging processed chunks.\n"
                        "bcftools is safer and feature-rich; cat is simpler and faster but less flexible."
                    )

                dpg.add_spacer(height=6)
                app.large_resume = dpg.add_checkbox(
                    label="Resume interrupted runs",
                    default_value=True,
                    tag="large_resume",
                )
                app._inputs.append(app.large_resume)
                with dpg.tooltip(app.large_resume):
                    dpg.add_text(
                        "Reuse existing intermediate results when rerunning a job.\n"
                        "Saves time after crashes or manual interruption."
                    )

                dpg.add_spacer(height=6)
                app.large_temp_dir = dpg.add_input_text(
                    label="Temp folder (optional)",
                    width=260,
                    default_value="",
                    tag="large_temp_dir",
                )
                app._inputs.append(app.large_temp_dir)
                with dpg.tooltip(app.large_temp_dir):
                    dpg.add_text(
                        "Custom directory for temporary files.\n"
                        "Use a fast local disk/SSD if possible. Leave empty for system default."
                    )

        dpg.add_spacer(height=18)
        dpg.add_text("Appearance", color=(200, 180, 90))
        dpg.add_spacer(height=8)

        font_scale = dpg.add_slider_float(
            label="Font scale",
            min_value=0.85,
            max_value=1.35,
            default_value=1.10,
            width=520,
            tag="settings_font_scale",
        )
        with dpg.tooltip(font_scale):
            dpg.add_text(
                "Global font scaling factor.\n"
                "Larger values improve readability on high-resolution screens."
            )

        dpg.add_spacer(height=8)

        accent_combo = dpg.add_combo(
            items=[
                "Evergreen (Green)",
                "Teal",
                "Blue",
                "Amber",
                "Purple",
            ],
            default_value="Evergreen (Green)",
            width=260,
            label="Accent color",
            tag="settings_accent_combo",
        )
        with dpg.tooltip(accent_combo):
            dpg.add_text(
                "Accent color used for buttons, highlights and progress indicators.\n"
                "No impact on performance."
            )

    return "page_settings"


def build_pages(app, parent):
    pages = {}

    def _mount(key, builder_fn):
        container = f"view_{key}"
        if dpg.does_item_exist(container):
            dpg.delete_item(container)
        with dpg.child_window(tag=container, parent=parent, show=False, border=False):
            inner = builder_fn(app, container)
            if isinstance(inner, str) and dpg.does_item_exist(inner):
                dpg.configure_item(inner, show=True)
        pages[key] = container

    _mount("ref_manager", page_reference_manager)
    _mount("pangenome", page_pangenome_builder)
    _mount("fastq_qc", page_fastq_qc)
    _mount("alignment", page_alignment)
    _mount("pre_sam", page_preprocess_samtools)
    _mount("vc", page_variant_calling)
    _mount("pre_bcf", page_preprocess_bcftools)
    _mount("check_vcf", page_check_vcf)
    _mount("plink", page_convert_plink)
    _mount("ld", page_ld_analysis)
    _mount("gwas", page_gwas)
    _mount("pca", page_pca)
    _mount("gp", page_genomic_prediction)
    _mount("batch", page_batch_gwas)
    _mount("settings", page_settings)

    return pages
