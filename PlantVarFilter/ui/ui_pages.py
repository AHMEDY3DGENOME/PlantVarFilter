# ui/ui_pages.py
from __future__ import annotations

from pathlib import Path
import os
import platform as _pf
import subprocess

import dearpygui.dearpygui as dpg

try:
    from PlantVarFilter.preanalysis import (
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


def _file_input_row(label: str, tag_key: str, parent, file_extensions: tuple[str, ...] = (".*",), app=None):
    with dpg.group(parent=parent, horizontal=True, horizontal_spacing=8):
        dpg.add_text(label)
        dpg.add_input_text(tag=f"input_{tag_key}", width=460)
        dpg.add_button(label="Browse", callback=lambda: dpg.show_item(f"dlg_{tag_key}"))
    with dpg.file_dialog(tag=f"dlg_{tag_key}", directory_selector=False, show=False,
                         callback=lambda s, a: dpg.set_value(f"input_{tag_key}", a["file_path_name"])):
        for ext in file_extensions:
            dpg.add_file_extension(ext, color=(150, 150, 150, 255))


def _dir_input_row(label: str, tag_key: str, parent):
    with dpg.group(parent=parent, horizontal=True, horizontal_spacing=8):
        dpg.add_text(label)
        dpg.add_input_text(tag=f"input_{tag_key}", width=460)
        dpg.add_button(label="Select", callback=lambda: dpg.show_item(f"dlg_{tag_key}"))
    with dpg.file_dialog(tag=f"dlg_{tag_key}", directory_selector=True, show=False,
                         callback=lambda s, a: dpg.set_value(f"input_{tag_key}", a["file_path_name"])):
        pass


def page_reference_manager(app, parent):
    with dpg.group(parent=parent, show=False, tag="page_ref_manager"):
        dpg.add_text("\nReference Manager", indent=10)
        dpg.add_spacer(height=10)

        _file_input_row("Reference FASTA:", "ref_fasta", parent, (".fa", ".fasta", ".fna", ".*"), app=app)
        _dir_input_row("Reference out dir (optional):", "ref_out_dir", parent)

        dpg.add_spacer(height=6)
        dpg.add_button(label="Build / Refresh Indexes", width=240, height=36,
                       callback=lambda: _on_build_reference(app))

        dpg.add_spacer(height=8)
        dpg.add_text("Status:", parent=parent)
        dpg.add_child_window(tag="ref_status_area", parent=parent, width=-1, height=240, border=True)
    return "page_ref_manager"


def _render_ref_status(st: ReferenceIndexStatus):
    dpg.delete_item("ref_status_area", children_only=True)
    lines = [
        f"FASTA: {st.fasta}",
        f"Directory: {st.reference_dir}",
        f"faidx (.fai): {st.faidx or 'missing'}",
        f"dict (.dict): {st.dict or 'missing'}",
        f"minimap2 index (.mmi): {st.mmi or 'missing'}",
        f"bowtie2 prefix: {st.bt2_prefix or 'missing'}",
        "Tools: " + ", ".join(f"{k}: {'ok' if (v.get('path')) else 'missing'}" for k, v in st.tools.items()),
        f"OK: {st.ok}",
    ]
    for ln in lines:
        dpg.add_text(ln, parent="ref_status_area")
    if st.reference_dir:
        dpg.add_spacer(height=6, parent="ref_status_area")
        dpg.add_button(label="Open reference folder", parent="ref_status_area",
                       callback=lambda: _open_in_os(st.reference_dir))


def _on_build_reference(app):
    if not _HAS_PRE:
        app.add_log("[REF] Preanalysis modules not available in environment.", error=True)
        return
    fasta = dpg.get_value("input_ref_fasta")
    out_dir = dpg.get_value("input_ref_out_dir") or None
    if not fasta or not Path(fasta).exists():
        app.add_log("[REF] Select a valid FASTA file first.", warn=True)
        return
    rm = ReferenceManager(logger=app.add_log, workspace=app.workspace_dir)
    st = rm.build_indices(fasta, out_dir=out_dir)
    _render_ref_status(st)
    app.add_log("[REF] Reference indexing finished.")


def page_fastq_qc(app, parent):
    with dpg.group(parent=parent, show=False, tag="page_fastq_qc"):
        dpg.add_text("\nFASTQ Quality Control", indent=10)
        dpg.add_spacer(height=10)

        with dpg.group(parent=parent, horizontal=True, horizontal_spacing=12):
            dpg.add_text("Platform:")
            dpg.add_combo(items=["illumina", "ont", "hifi", "pb"], default_value="illumina", width=180,
                          tag="fq_platform")

        _file_input_row("Reads #1 (R1 or single):", "fq_r1", parent, (".fastq", ".fq", ".fastq.gz", ".fq.gz", ".*"), app=app)
        _file_input_row("Reads #2 (R2, optional):", "fq_r2", parent, (".fastq", ".fq", ".fastq.gz", ".fq.gz", ".*"), app=app)
        _dir_input_row("Output dir (optional):", "fq_out", parent)

        with dpg.group(parent=parent, horizontal=True, horizontal_spacing=12):
            dpg.add_button(label="Run QC", width=200, height=36, callback=lambda: _on_run_fastq_qc(app))
            dpg.add_checkbox(label="Use FastQC if available", default_value=True, tag="fq_use_fastqc")

        dpg.add_spacer(height=8)
        dpg.add_text("QC Summary:", parent=parent)
        dpg.add_child_window(tag="fq_qc_area", parent=parent, width=-1, height=260, border=True)
    return "page_fastq_qc"


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

        with dpg.group(parent=parent, horizontal=True, horizontal_spacing=12):
            dpg.add_text("Platform:")
            dpg.add_combo(items=["illumina", "ont", "hifi", "pb"], default_value="illumina", width=180,
                          tag="aln_platform")

        with dpg.collapsing_header(label="Reference", parent=parent, default_open=True):
            dpg.add_text("For ONT/PB: select FASTA or .mmi\nFor Illumina: select bowtie2 prefix base path (no suffix).")
            _file_input_row("Reference (.fa/.mmi or bt2 prefix)", "aln_reference", parent,
                            (".fa", ".fasta", ".fna", ".mmi", ".*"))

        with dpg.collapsing_header(label="Reads", parent=parent, default_open=True):
            _file_input_row("Reads #1 (R1 or single):", "aln_r1", parent,
                            (".fastq", ".fq", ".fastq.gz", ".fq.gz", ".*"))
            _file_input_row("Reads #2 (R2, optional):", "aln_r2", parent,
                            (".fastq", ".fq", ".fastq.gz", ".fq.gz", ".*"))

        with dpg.group(parent=parent, horizontal=True, horizontal_spacing=12):
            dpg.add_text("Threads")
            dpg.add_input_int(tag="aln_threads", default_value=8, min_value=1, width=80)
            dpg.add_checkbox(label="Save SAM", tag="aln_save_sam")
            dpg.add_checkbox(label="Mark duplicates", default_value=True, tag="aln_markdup")

        with dpg.collapsing_header(label="Read Group (optional)", parent=parent, default_open=False):
            for k in ("ID", "SM", "LB", "PL", "PU"):
                with dpg.group(horizontal=True, horizontal_spacing=8):
                    dpg.add_text(k)
                    dpg.add_input_text(tag=f"aln_rg_{k}", width=240)

        dpg.add_spacer(height=6)
        dpg.add_button(label="Run Alignment", width=200, height=36, callback=lambda: _on_run_alignment(app))

        dpg.add_spacer(height=8)
        dpg.add_text("Result:", parent=parent)
        dpg.add_child_window(tag="aln_result_area", parent=parent, width=-1, height=240, border=True)
    return "page_alignment"


def _on_run_alignment(app):
    if not _HAS_PRE:
        app.add_log("[ALN] Preanalysis modules not available in environment.", error=True)
        return
    ref = dpg.get_value("input_aln_reference")
    r1 = dpg.get_value("input_aln_r1")
    r2 = dpg.get_value("input_aln_r2") or None
    platform = (dpg.get_value("aln_platform") or "illumina").lower()
    threads = dpg.get_value("aln_threads") or 8
    save_sam = bool(dpg.get_value("aln_save_sam"))
    markdup = bool(dpg.get_value("aln_markdup"))
    if not ref or not Path(ref).exists():
        app.add_log("[ALN] Select a valid reference (.fa/.mmi or bowtie2 prefix path).", warn=True)
        return
    if not r1 or not Path(r1).exists():
        app.add_log("[ALN] Select valid FASTQ reads.", warn=True)
        return
    rg = {}
    for k in ("ID", "SM", "LB", "PL", "PU"):
        v = dpg.get_value(f"aln_rg_{k}")
        if v:
            rg[k] = v
    aln = Aligner(logger=app.add_log, workspace=app.workspace_dir)
    if platform in {"ont", "hifi", "pb"}:
        preset = "map-ont" if platform == "ont" else ("map-hifi" if platform == "hifi" else "map-pb")
        reads = [r1] if not r2 else [r1, r2]
        res = aln.minimap2(ref, reads, preset=preset, threads=threads, read_group=rg or None,
                           save_sam=save_sam, mark_duplicates=markdup, out_dir=app.workspace_dir)
    else:
        res = aln.bowtie2(ref, r1, r2, threads=threads, read_group=rg or None,
                          save_sam=save_sam, mark_duplicates=markdup, out_dir=app.workspace_dir)
    dpg.delete_item("aln_result_area", children_only=True)
    dpg.add_text(f"Tool: {res.tool}", parent="aln_result_area")
    if res.sam:
        dpg.add_text(f"SAM: {res.sam}", parent="aln_result_area")
    dpg.add_text(f"BAM: {res.bam}", parent="aln_result_area")
    dpg.add_text(f"BAI: {res.bai}", parent="aln_result_area")
    dpg.add_text(f"flagstat: {res.flagstat}", parent="aln_result_area")
    dpg.add_text(f"Elapsed: {res.elapsed_sec:.1f} sec", parent="aln_result_area")
    dpg.add_spacer(height=6, parent="aln_result_area")
    dpg.add_button(label="Open alignment folder", parent="aln_result_area",
                   callback=lambda: _open_in_os(str(Path(res.bam).parent)))
    app.add_log("[ALN] Alignment finished.")


def page_preprocess_samtools(app, parent):
    with dpg.group(parent=parent, show=False, tag="page_pre_sam"):
        dpg.add_text("\nClean BAM from SAM/BAM: sort / fixmate / markdup / index + QC reports", indent=10)
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
                    label="Threads", width=220, default_value=4, min_value=1, min_clamped=True
                )
                app._inputs.append(app.sam_threads)

                app.sam_remove_dups = dpg.add_checkbox(
                    label="Remove duplicates (instead of marking)", default_value=False
                )
                app._inputs.append(app.sam_remove_dups)

                app.sam_compute_stats = dpg.add_checkbox(
                    label="Compute QC reports (flagstat/stats/idxstats/depth)", default_value=True
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
                    label="Run samtools preprocess",
                    callback=app.run_samtools_preprocess,
                    width=240,
                    height=38,
                )
                app._primary_buttons.append(run_sam)
    return "page_pre_sam"


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
                    label="Threads", width=220, default_value=4, min_value=1, min_clamped=True
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
                    label="Call variants (bcftools)",
                    callback=app.run_variant_calling,
                    width=240,
                    height=38,
                )
                app._primary_buttons.append(run_vc)

    return "page_vc"


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
                    label="Fill tags (AC, AN, AF, MAF, HWE) before filtering", default_value=True
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
                    label="Run bcftools preprocess",
                    callback=app.run_bcftools_preprocess,
                    width=240,
                    height=38,
                )
                app._primary_buttons.append(run_bcf)
    return "page_pre_bcf"


def page_check_vcf(app, parent):
    with dpg.group(parent=parent, show=False, tag="page_check_vcf"):
        dpg.add_text("\nCheck VCF quality before conversion/analysis", indent=10)
        dpg.add_spacer(height=10)
        with dpg.group(horizontal=True, horizontal_spacing=60):
            with dpg.group():
                vcf_btn_qc = dpg.add_button(
                    label="Choose a VCF file",
                    callback=lambda: dpg.show_item("file_dialog_vcf"),
                    width=220,
                    tag="tooltip_vcf_qc",
                )
                app._secondary_buttons.append(vcf_btn_qc)
                dpg.add_text("", tag="qc_vcf_path_lbl", wrap=500)

                dpg.add_spacer(height=6)
                vcf_btn_qc2 = dpg.add_button(
                    label="Choose another VCF (optional)",
                    callback=lambda: dpg.show_item("file_dialog_vcf2"),
                    width=220,
                )
                app._secondary_buttons.append(vcf_btn_qc2)
                dpg.add_text("", tag="qc_vcf2_path_lbl", wrap=500)

                dpg.add_spacer(height=6)
                bl_btn = dpg.add_button(
                    label="Choose blacklist BED (optional)",
                    callback=lambda: dpg.show_item("file_dialog_blacklist"),
                    width=220,
                    tag="tooltip_bl_qc",
                )
                app._secondary_buttons.append(bl_btn)
                dpg.add_text("", tag="qc_bl_path_lbl", wrap=500)

            with dpg.group():
                app.deep_scan = dpg.add_checkbox(label="Deep scan", default_value=False)
                app._inputs.append(app.deep_scan)

                dpg.add_spacer(height=12)
                run_qc_btn = dpg.add_button(
                    label="Run Quality Check",
                    callback=app.run_vcf_qc,
                    width=200,
                    height=36,
                )
                app._primary_buttons.append(run_qc_btn)
    return "page_check_vcf"


def page_convert_plink(app, parent):
    with dpg.group(parent=parent, show=False, tag="page_plink"):
        dpg.add_text("\nConvert a VCF file into PLINK BED and apply MAF/missing genotype filters.", indent=10)
        dpg.add_spacer(height=10)
        with dpg.group(horizontal=True, horizontal_spacing=60):
            with dpg.group():
                dpg.add_text("Select files:", indent=0)
                vcf = dpg.add_button(
                    label="Choose a VCF file",
                    callback=lambda: dpg.show_item("file_dialog_vcf"),
                    width=220,
                    tag="tooltip_vcf",
                )
                app._secondary_buttons.append(vcf)
                dpg.add_text("", tag="conv_vcf_path_lbl", wrap=500)

                dpg.add_spacer(height=6)
                variant_ids = dpg.add_button(
                    label="Choose IDs file (optional)",
                    callback=lambda: dpg.show_item("file_dialog_variants"),
                    width=220,
                    tag="tooltip_variant",
                )
                app._secondary_buttons.append(variant_ids)
                dpg.add_text("", tag="ids_path_lbl", wrap=500)

            with dpg.group():
                dpg.add_text("Apply filters:", indent=0)
                maf_input = dpg.add_input_float(
                    label="Minor allele frequency (MAF)",
                    width=220,
                    default_value=0.05,
                    step=0.005,
                    tag="tooltip_maf",
                )
                app._inputs.append(maf_input)

                dpg.add_spacer(height=6)
                geno_input = dpg.add_input_float(
                    label="Missing genotype rate",
                    width=220,
                    default_value=0.10,
                    step=0.005,
                    tag="tooltip_missing",
                )
                app._inputs.append(geno_input)

                dpg.add_spacer(height=14)
                convert_btn = dpg.add_button(
                    tag="convert_vcf_btn",
                    label="Convert VCF",
                    callback=app.convert_vcf,
                    user_data={"maf": maf_input, "geno": geno_input},
                    width=160,
                    height=36,
                    enabled=True,
                )
                app._primary_buttons.append(convert_btn)
    return "page_plink"


def page_ld_analysis(app, parent):
    with dpg.group(parent=parent, show=False, tag="page_ld"):
        dpg.add_text("\nLD analysis: LD decay, LD heatmap, and diversity metrics", indent=10)
        dpg.add_spacer(height=10)

        with dpg.group(horizontal=True, horizontal_spacing=60):
            with dpg.group():
                dpg.add_text("Input files", color=(200, 220, 200))
                dpg.add_spacer(height=6)

                btn_bed = dpg.add_button(
                    label="Choose PLINK BED",
                    callback=lambda: dpg.show_item("file_dialog_bed"),
                    width=220,
                )
                app._secondary_buttons.append(btn_bed)
                dpg.add_text("", tag="ld_bed_path_lbl", wrap=520)

                dpg.add_spacer(height=6)
                btn_vcf = dpg.add_button(
                    label="Choose VCF (optional)",
                    callback=lambda: dpg.show_item("file_dialog_vcf"),
                    width=220,
                )
                app._secondary_buttons.append(btn_vcf)
                dpg.add_text("", tag="ld_vcf_path_lbl", wrap=520)

                dpg.add_spacer(height=6)
                app.ld_region = dpg.add_input_text(
                    label="Region (optional)",
                    hint="chr:start-end (e.g., 1:1000000-2000000)",
                    width=320,
                )
                app._inputs.append(app.ld_region)

            with dpg.group():
                dpg.add_text("Options", color=(200, 220, 200))
                dpg.add_spacer(height=6)

                app.ld_window_kb = dpg.add_input_int(
                    label="LD window (kb)",
                    default_value=500,
                    min_value=1,
                    min_clamped=True,
                    width=220,
                )
                app._inputs.append(app.ld_window_kb)

                app.ld_window_snp = dpg.add_input_int(
                    label="LD window size (SNPs)",
                    default_value=5000,
                    min_value=10,
                    min_clamped=True,
                    width=220,
                )
                app._inputs.append(app.ld_window_snp)

                app.ld_max_kb = dpg.add_input_int(
                    label="Max distance (kb)",
                    default_value=1000,
                    min_value=1,
                    min_clamped=True,
                    width=220,
                )
                app._inputs.append(app.ld_max_kb)

                app.ld_min_r2 = dpg.add_input_float(
                    label="Min r²",
                    default_value=0.1,
                    min_value=0.0,
                    max_value=1.0,
                    min_clamped=True,
                    max_clamped=True,
                    step=0.05,
                    width=220,
                )
                app._inputs.append(app.ld_min_r2)

                dpg.add_spacer(height=6)
                app.ld_do_decay = dpg.add_checkbox(label="Compute LD decay", default_value=True)
                app._inputs.append(app.ld_do_decay)

                app.ld_do_heatmap = dpg.add_checkbox(label="LD heatmap", default_value=True)
                app._inputs.append(app.ld_do_heatmap)

                app.ld_do_div = dpg.add_checkbox(label="Diversity metrics", default_value=True)
                app._inputs.append(app.ld_do_div)

                dpg.add_spacer(height=12)
                run_btn = dpg.add_button(
                    label="Run LD analysis",
                    callback=app.run_ld_analysis,
                    width=240,
                    height=38,
                )
                app._primary_buttons.append(run_btn)

    return "page_ld"


def page_gwas(app, parent):
    with dpg.group(parent=parent, show=False, tag="page_gwas"):
        dpg.add_text("\nStart GWAS Analysis", indent=10)
        dpg.add_spacer(height=10)

        with dpg.group(horizontal=True, horizontal_spacing=60):
            with dpg.group():
                geno = dpg.add_button(
                    label="Choose a BED file",
                    callback=lambda: dpg.show_item("file_dialog_bed"),
                    width=220,
                    tag="tooltip_bed",
                )
                app._secondary_buttons.append(geno)
                dpg.add_text("", tag="gwas_bed_path_lbl", wrap=500)

                dpg.add_spacer(height=6)
                pheno = dpg.add_button(
                    label="Choose a phenotype file",
                    callback=lambda: dpg.show_item("file_dialog_pheno"),
                    width=220,
                    tag="tooltip_pheno",
                )
                app._secondary_buttons.append(pheno)
                dpg.add_text("", tag="gwas_pheno_path_lbl", wrap=500)

                dpg.add_spacer(height=6)
                cov_file = dpg.add_button(
                    label="Choose covariate file (optional)",
                    callback=lambda: dpg.show_item("file_dialog_cov"),
                    width=220,
                    tag="tooltip_cov",
                )
                app._secondary_buttons.append(cov_file)
                dpg.add_text("", tag="gwas_cov_path_lbl", wrap=500)

                dpg.add_spacer(height=6)
                kin_file = dpg.add_button(
                    label="Choose kinship file (optional)",
                    callback=lambda: dpg.show_item("file_dialog_kinship"),
                    width=220,
                    tag="tooltip_kinship",
                )
                app._secondary_buttons.append(kin_file)
                dpg.add_text("", tag="gwas_kinship_path_lbl", wrap=500)

                dpg.add_spacer(height=6)
                kin_compute = dpg.add_button(
                    label="Compute kinship from BED",
                    callback=lambda: app.compute_kinship_from_bed(),
                    width=220,
                )
                app._secondary_buttons.append(kin_compute)

                dpg.add_spacer(height=8)
                anno_btn = dpg.add_button(
                    label="Choose GTF/GFF annotation",
                    callback=lambda: dpg.show_item("file_dialog_gtf"),
                    width=220,
                )
                app._secondary_buttons.append(anno_btn)
                dpg.add_text("", tag="gwas_gtf_path_lbl", wrap=500)

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
                        "SAIGE (mixed model)"
                    ],
                    width=260,
                    default_value="FaST-LMM",
                    tag="tooltip_algorithm",
                )
                app._inputs.append(app.gwas_combo)

                dpg.add_spacer(height=8)
                app.snp_limit = dpg.add_input_int(
                    label="Limit SNPs in plots (optional)",
                    width=260,
                    min_value=0,
                    step=1000,
                    default_value=0
                )
                app._inputs.append(app.snp_limit)

                dpg.add_spacer(height=6)
                app.plot_stats = dpg.add_checkbox(
                    label="Produce pheno/geno statistics (PDF)",
                    default_value=False
                )

                dpg.add_spacer(height=6)
                app.annotate_enable = dpg.add_checkbox(
                    label="Annotate GWAS with GTF/GFF",
                    default_value=False
                )
                app.annotate_window_kb = dpg.add_input_int(
                    label="Window around TSS (kb)",
                    width=260,
                    min_value=0,
                    step=10,
                    default_value=50
                )

                dpg.add_spacer(height=10)
                dpg.add_text("Optional Region Filter", color=(150, 150, 255))
                app.region_chr = dpg.add_input_text(
                    label="Chromosome",
                    width=260,
                    default_value=""
                )
                app.region_start = dpg.add_input_int(
                    label="Start position",
                    width=260,
                    min_value=0,
                    step=1000,
                    default_value=0
                )
                app.region_end = dpg.add_input_int(
                    label="End position",
                    width=260,
                    min_value=0,
                    step=1000,
                    default_value=0
                )

                dpg.add_spacer(height=14)
                gwas_btn = dpg.add_button(
                    label="Run GWAS",
                    callback=lambda s, a: app.run_gwas(s, a, [geno, pheno]),
                    width=200,
                    height=36,
                )
                app._primary_buttons.append(gwas_btn)

        dpg.add_spacer(height=12)
        dpg.add_separator()
        dpg.add_spacer(height=8)
        dpg.add_text(
            "Results will appear in the Results window (tables, Manhattan/QQ plots).",
            wrap=900
        )

    return "page_gwas"


def page_pca(app, parent):
    with dpg.group(parent=parent, show=False, tag="page_pca"):
        dpg.add_text("Population Structure (PCA) & Kinship", indent=10)
        dpg.add_spacer(height=10)

        with dpg.group(horizontal=True, horizontal_spacing=60):
            with dpg.group():
                btn_bed = dpg.add_button(
                    label="Choose PLINK BED",
                    callback=lambda: dpg.show_item("file_dialog_bed"),
                    width=220,
                    tag="tooltip_bed_pca"
                )
                app._secondary_buttons.append(btn_bed)
                dpg.add_text("", tag="pca_bed_path_lbl", wrap=500)

                dpg.add_spacer(height=8)
                app.pca_out_prefix = dpg.add_input_text(
                    label="Output prefix (optional)",
                    hint="Default: next to BED prefix",
                    width=320
                )
                app._inputs.append(app.pca_out_prefix)

            with dpg.group():
                app.pca_npcs = dpg.add_input_int(label="Number of PCs", default_value=10, min_value=2, max_value=50, width=160)
                app._inputs.append(app.pca_npcs)

                app.pca_kinship = dpg.add_checkbox(label="Compute kinship matrix", default_value=True)
                app._inputs.append(app.pca_kinship)

                dpg.add_spacer(height=12)
                run_btn = dpg.add_button(
                    label="Run PCA",
                    callback=app.run_pca_module,
                    width=200, height=36
                )
                app._primary_buttons.append(run_btn)

        dpg.add_spacer(height=12)
        dpg.add_separator()
        dpg.add_spacer(height=6)
        dpg.add_text("Results preview will appear in the Results window (PC1 vs PC2 plot, files list).")

    return "page_pca"


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
                    tag="tooltip_bed_gp",
                )
                app._secondary_buttons.append(geno)
                dpg.add_text("", tag="gp_bed_path_lbl", wrap=500)

                dpg.add_spacer(height=6)
                pheno = dpg.add_button(
                    label="Choose a phenotype file",
                    callback=lambda: dpg.show_item("file_dialog_pheno"),
                    width=220,
                    tag="tooltip_pheno_gp",
                )
                app._secondary_buttons.append(pheno)
                dpg.add_text("", tag="gp_pheno_path_lbl", wrap=500)

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
                    tag="tooltip_algorithm_gp",
                )
                app._inputs.append(app.gwas_gp)

                dpg.add_spacer(height=14)
                gp_btn = dpg.add_button(
                    label="Run Genomic Prediction",
                    callback=lambda s, a: app.run_genomic_prediction(s, a, [geno, pheno]),
                    width=220,
                    height=36,
                )
                app._primary_buttons.append(gp_btn)
    return "page_gp"


def page_batch_gwas(app, parent):
    with dpg.group(parent=parent, show=False, tag="page_batch"):
        dpg.add_text("\nBatch GWAS for all traits in a phenotype file (FID IID + multiple traits).", indent=10)
        dpg.add_spacer(height=10)
        with dpg.group(horizontal=True, horizontal_spacing=60):
            with dpg.group():
                geno_btn = dpg.add_button(
                    label="Choose a BED file (SNP or SV)",
                    callback=lambda: dpg.show_item("file_dialog_bed"),
                    width=240,
                )
                app._secondary_buttons.append(geno_btn)
                dpg.add_text("", tag="batch_bed_path_lbl", wrap=520)

                dpg.add_spacer(height=6)
                pheno_btn = dpg.add_button(
                    label="Choose a multi-trait phenotype file",
                    callback=lambda: dpg.show_item("file_dialog_pheno"),
                    width=240,
                )
                app._secondary_buttons.append(pheno_btn)
                dpg.add_text("", tag="batch_pheno_path_lbl", wrap=520)

                dpg.add_spacer(height=6)
                cov_btn = dpg.add_button(
                    label="Choose covariates (optional)",
                    callback=lambda: dpg.show_item("file_dialog_cov"),
                    width=240,
                )
                app._secondary_buttons.append(cov_btn)
                dpg.add_text("", tag="batch_cov_path_lbl", wrap=520)

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
                )
                app._inputs.append(app.batch_algo)

                dpg.add_spacer(height=8)
                dpg.add_text("Uses settings from 'Settings' page (trees, depth, train %, jobs ...).")

                dpg.add_spacer(height=12)
                run_btn = dpg.add_button(
                    label="Run Batch GWAS",
                    callback=app.run_batch_gwas_ui,
                    width=220,
                    height=38,
                )
                app._primary_buttons.append(run_btn)

    return "page_batch"


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
                with dpg.group():
                    dpg.add_text("General Settings", color=(200, 180, 90))
                    dpg.add_spacer(height=8)

                    dpg.add_checkbox(
                        label="Night Mode (Dark)",
                        default_value=True,
                        tag="settings_dark_toggle",
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

                    app.plot_stats = dpg.add_checkbox(
                        label="Advanced Plotting", default_value=False, tag="tooltip_stats"
                    )
                    app._inputs.append(app.plot_stats)

                    app.snp_limit = dpg.add_input_text(
                        label="SNP limit", width=220, default_value="", tag="tooltip_limit"
                    )
                    app._inputs.append(app.snp_limit)

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

                    app.aggregation_method = dpg.add_combo(
                        ("sum", "median", "mean"),
                        label="Aggregation Method",
                        width=220,
                        default_value="sum",
                        tag="tooltip_aggr",
                    )
                    app._inputs.append(app.aggregation_method)

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

            with dpg.group():
                app.large_merge_strategy = dpg.add_combo(
                    items=["bcftools", "cat"],
                    label="Merge strategy",
                    width=260,
                    default_value="bcftools",
                    tag="large_merge_strategy",
                )
                app._inputs.append(app.large_merge_strategy)

                dpg.add_spacer(height=6)
                app.large_resume = dpg.add_checkbox(
                    label="Resume interrupted runs",
                    default_value=True,
                    tag="large_resume",
                )
                app._inputs.append(app.large_resume)

                dpg.add_spacer(height=6)
                app.large_temp_dir = dpg.add_input_text(
                    label="Temp folder (optional)",
                    width=260,
                    default_value="",
                    tag="large_temp_dir",
                )
                app._inputs.append(app.large_temp_dir)

        dpg.add_spacer(height=18)
        dpg.add_text("Appearance", color=(200, 180, 90))
        dpg.add_spacer(height=8)

        dpg.add_slider_float(
            label="Font scale",
            min_value=0.85,
            max_value=1.35,
            default_value=1.10,
            width=520,
            tag="settings_font_scale",
        )

        dpg.add_spacer(height=8)

        dpg.add_combo(
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
