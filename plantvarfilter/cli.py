"""
plantvarfilter CLI
==================
Headless command-line interface for PlantOmicsGWAS / plantvarfilter.

Wraps the existing analytical engine (gwas_pipeline, vcf_quality, samtools_utils,
bcftools_utils, etc.) so the full GUI workflow runs on Linux servers via SSH
without a display.

Entry point (registered in pyproject.toml as 'plantvarfilter'):
    plantvarfilter <subcommand> [options]
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import sys
import traceback
from dataclasses import asdict, is_dataclass
from pathlib import Path
from typing import List, Optional


# ============================================================================
# Logging
# ============================================================================

class CLILogger:
    """Drop-in replacement for the GUI's add_log callback.

    The engine modules call `add_log(msg, error=True/False, warn=True/False)`.
    This logger forwards everything to stdout/stderr so it shows up in the
    user's terminal (or Slurm/nohup logfile).
    """

    def __init__(self, verbose: bool = True):
        self.verbose = verbose

    def __call__(self, msg: str, *, error: bool = False, warn: bool = False, **_):
        if error:
            tag, stream = "[ERROR]", sys.stderr
        elif warn:
            tag, stream = "[WARN]",  sys.stderr
        else:
            tag, stream = "[INFO]",  sys.stdout
            if not self.verbose:
                return
        print(f"{tag} {msg}", file=stream, flush=True)


def make_logger(args: argparse.Namespace) -> CLILogger:
    return CLILogger(verbose=not getattr(args, "quiet", False))


# ============================================================================
# Common helpers
# ============================================================================

def _abspath(p: Optional[str]) -> Optional[str]:
    return str(Path(p).expanduser().resolve()) if p else None


def _ensure_dir(p: Optional[str]) -> Optional[str]:
    if not p:
        return None
    Path(p).mkdir(parents=True, exist_ok=True)
    return _abspath(p)


def _bed_prefix(p: str) -> str:
    """Accept either a PLINK prefix or any of .bed/.bim/.fam; return the prefix."""
    s = str(p)
    for ext in (".bed", ".bim", ".fam"):
        if s.endswith(ext):
            return s[: -len(ext)]
    return s


def _add_common_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("--quiet", action="store_true", help="Suppress INFO logs.")
    p.add_argument("--debug", action="store_true",
                   help="Print full Python traceback on error.")


def _dump_json(obj) -> None:
    """Pretty-print a dataclass / dict / mixed object as JSON."""
    if is_dataclass(obj):
        obj = asdict(obj)
    print(json.dumps(obj, default=str, indent=2))


def _exit_ok(msg: str = "Done.") -> int:
    print(f"[OK] {msg}", flush=True)
    return 0


def _exit_err(exc: BaseException, debug: bool = False) -> int:
    sys.stderr.write(f"[FAIL] {type(exc).__name__}: {exc}\n")
    if debug:
        traceback.print_exc()
    return 2


def _build_chrom_mapping(bed_prefix: str):
    """Use the engine's HELPERS to map non-numeric chromosome labels to integers."""
    from plantvarfilter.helpers import HELPERS
    return HELPERS().replace_with_integers(f"{bed_prefix}.bim")


def _load_bed_pheno(bed_prefix: str, pheno_path: str):
    """Open BED + Pheno via pysnptools, apply filter_out_missing.

    Returns (gwas_obj, bed, bed_fixed, pheno, chrom_mapping).
    """
    from pysnptools.snpreader import Bed, Pheno
    from plantvarfilter.gwas_pipeline import GWAS

    chrom_mapping = _build_chrom_mapping(bed_prefix)
    bed = Bed(bed_prefix, count_A1=False, chrom_map=chrom_mapping)
    pheno = Pheno(pheno_path)
    gwas = GWAS()
    bed_fixed = gwas.filter_out_missing(bed)
    return gwas, bed, bed_fixed, pheno, chrom_mapping


# ============================================================================
# Subcommand: ref-index
# ============================================================================

def cmd_ref_index(args: argparse.Namespace) -> int:
    from plantvarfilter.preanalysis.reference_manager import ReferenceManager

    log = make_logger(args)
    rm = ReferenceManager(logger=log, workspace=_ensure_dir(args.out) or os.getcwd())
    status = rm.build_indices(
        fasta=_abspath(args.fasta),
        out_dir=_abspath(args.out),
        build_mmi=not args.no_mmi,
        build_bt2=not args.no_bt2,
        build_dict=not args.no_dict,
    )
    _dump_json(status)
    return _exit_ok("Reference indexes ready.")


def _parser_ref_index(sp: argparse._SubParsersAction) -> None:
    p = sp.add_parser("ref-index", help="Build samtools/bowtie2/minimap2 indexes.")
    p.add_argument("--fasta", required=True, help="Reference FASTA.")
    p.add_argument("--out", default=None, help="Output directory.")
    p.add_argument("--no-mmi",  action="store_true", help="Skip minimap2 index.")
    p.add_argument("--no-bt2",  action="store_true", help="Skip bowtie2 index.")
    p.add_argument("--no-dict", action="store_true", help="Skip sequence dictionary.")
    _add_common_args(p)
    p.set_defaults(func=cmd_ref_index)


# ============================================================================
# Subcommand: fastq-qc
# ============================================================================

def cmd_fastq_qc(args: argparse.Namespace) -> int:
    from plantvarfilter.preanalysis.fastq_qc import run_fastq_qc

    log = make_logger(args)
    reads1 = _abspath(args.reads[0])
    reads2 = _abspath(args.reads[1]) if len(args.reads) > 1 else None

    report = run_fastq_qc(
        reads1=reads1,
        reads2=reads2,
        platform=args.platform,
        out_dir=_abspath(args.out),
        sample_max=args.sample_max,
        use_fastqc_if_available=not args.no_fastqc,
        logger=log,
    )
    _dump_json(report)
    return _exit_ok("FASTQ QC complete.")


def _parser_fastq_qc(sp: argparse._SubParsersAction) -> None:
    p = sp.add_parser("fastq-qc", help="FASTQ quality report (length, GC, PHRED).")
    p.add_argument("--reads", nargs="+", required=True,
                   help="One (single-end) or two (paired-end) FASTQ files.")
    p.add_argument("--platform", default="illumina",
                   choices=["illumina", "ont", "nanopore", "pb", "pacbio", "hifi"],
                   help="Sequencing platform (default: illumina).")
    p.add_argument("--out", default=None, help="Output directory for QC artifacts.")
    p.add_argument("--sample-max", type=int, default=1_000_000,
                   help="Cap reads sampled for stats (default: 1,000,000).")
    p.add_argument("--no-fastqc", action="store_true",
                   help="Skip running FastQC even if available on PATH.")
    _add_common_args(p)
    p.set_defaults(func=cmd_fastq_qc)


# ============================================================================
# Subcommand: align
# ============================================================================

def cmd_align(args: argparse.Namespace) -> int:
    from plantvarfilter.preanalysis.aligner import Aligner

    log = make_logger(args)
    aln = Aligner(logger=log, workspace=_ensure_dir(args.out) or os.getcwd())

    rg = None
    if args.read_group:
        # Format: KEY=VAL,KEY=VAL,...
        rg = dict(item.split("=", 1) for item in args.read_group.split(",") if "=" in item)

    reads1 = _abspath(args.reads[0])
    reads2 = _abspath(args.reads[1]) if len(args.reads) > 1 else None

    result = aln.align(
        platform=args.platform,
        reference=_abspath(args.ref),
        reads1=reads1,
        reads2=reads2,
        threads=args.threads,
        read_group=rg,
        save_sam=args.keep_sam,
        out_dir=_abspath(args.out),
        out_prefix=args.prefix,
    )
    _dump_json(result)
    return _exit_ok("Alignment complete.")


def _parser_align(sp: argparse._SubParsersAction) -> None:
    p = sp.add_parser("align", help="Align reads (Bowtie2 short / Minimap2 long).")
    p.add_argument("--ref", required=True, help="Reference FASTA or Minimap2 .mmi.")
    p.add_argument("--reads", nargs="+", required=True,
                   help="Single-end (R1) or paired-end (R1 R2).")
    p.add_argument("--platform", default="illumina",
                   choices=["illumina", "short", "ont", "nanopore", "pb", "pacbio",
                            "hifi", "long"],
                   help="Read type. 'illumina'/'short' -> Bowtie2; others -> Minimap2.")
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--read-group", default=None,
                   help="Comma-separated key=value pairs, e.g. ID=S1,SM=S1,PL=ILLUMINA")
    p.add_argument("--keep-sam", action="store_true", help="Keep intermediate SAM file.")
    p.add_argument("--out", default=None, help="Output directory.")
    p.add_argument("--prefix", default=None, help="Output filename prefix.")
    _add_common_args(p)
    p.set_defaults(func=cmd_align)


# ============================================================================
# Subcommand: preprocess-bam
# ============================================================================

def cmd_preprocess_bam(args: argparse.Namespace) -> int:
    from plantvarfilter.samtools_utils import Samtools

    log = make_logger(args)
    st = Samtools()
    out = st.preprocess(
        input_path=_abspath(args.bam),
        out_prefix=args.prefix,
        threads=args.threads,
        remove_dups=args.remove_dups,
        compute_stats=not args.no_stats,
        log=log,
        keep_temps=args.keep_temps,
    )
    _dump_json(out)
    return _exit_ok("BAM preprocessing complete.")


def _parser_preprocess_bam(sp: argparse._SubParsersAction) -> None:
    p = sp.add_parser("preprocess-bam",
                      help="samtools sort/fixmate/markdup/index pipeline.")
    p.add_argument("--bam", required=True, help="Input BAM/SAM/CRAM.")
    p.add_argument("--prefix", default=None, help="Output prefix (default: <bam>.sm).")
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--remove-dups", action="store_true",
                   help="Remove duplicates instead of marking them.")
    p.add_argument("--no-stats", action="store_true",
                   help="Skip flagstat/stats/idxstats/depth.")
    p.add_argument("--keep-temps", action="store_true")
    _add_common_args(p)
    p.set_defaults(func=cmd_preprocess_bam)


# ============================================================================
# Subcommand: call
# ============================================================================

def cmd_call(args: argparse.Namespace) -> int:
    from plantvarfilter.variant_caller_utils import VariantCaller

    log = make_logger(args)
    vc = VariantCaller()
    vcf_gz, stats = vc.call_bcftools(
        bams=[_abspath(b) for b in args.bams],
        ref_fasta=_abspath(args.ref),
        out_prefix=_abspath(args.out),
        regions_bed=_abspath(args.regions),
        threads=args.threads,
        min_baseq=args.min_baseq,
        min_mapq=args.min_mapq,
        ploidy=args.ploidy,
        log=log,
        split_after_calling=args.split,
    )
    _dump_json({"vcf_gz": vcf_gz, "stats": stats})
    return _exit_ok("Variant calling complete.")


def _parser_call(sp: argparse._SubParsersAction) -> None:
    p = sp.add_parser("call", help="Variant calling with bcftools mpileup | call.")
    p.add_argument("--ref", required=True, help="Reference FASTA (must be indexed).")
    p.add_argument("--bams", nargs="+", required=True, help="One or more sorted BAMs.")
    p.add_argument("--out", default=None,
                   help="Output prefix (default: <bam_dir>/calls).")
    p.add_argument("--regions", default=None, help="BED of regions to restrict calling.")
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--min-baseq", type=int, default=13)
    p.add_argument("--min-mapq", type=int, default=20)
    p.add_argument("--ploidy", type=int, default=2)
    p.add_argument("--split", action="store_true",
                   help="Split multi-sample VCF after calling.")
    _add_common_args(p)
    p.set_defaults(func=cmd_call)


# ============================================================================
# Subcommand: vcf-prep
# ============================================================================

def cmd_vcf_prep(args: argparse.Namespace) -> int:
    from plantvarfilter.bcftools_utils import BCFtools

    log = make_logger(args)
    bt = BCFtools()
    final_vcf, stats = bt.preprocess(
        input_vcf=_abspath(args.vcf),
        out_prefix=_abspath(args.out),
        log=log,
        ref_fasta=_abspath(args.ref),
        regions_bed=_abspath(args.regions),
        split_multiallelic=not args.no_split_multi,
        left_align=not args.no_left_align,
        do_sort=not args.no_sort,
        set_id_from_fields=not args.no_set_id,
        filter_expr=args.filter_expr,
        remove_filtered=args.remove_filtered,
        compress_output=True,
        index_output=not args.no_index,
        keep_temps=args.keep_temps,
        fill_tags=args.fill_tags,
    )
    _dump_json({"final_vcf": final_vcf, "stats": stats})
    return _exit_ok("VCF preprocessing complete.")


def _parser_vcf_prep(sp: argparse._SubParsersAction) -> None:
    p = sp.add_parser("vcf-prep",
                      help="bcftools-based normalization, sort, ID, filter.")
    p.add_argument("--vcf", required=True, help="Input VCF/VCF.gz.")
    p.add_argument("--out", default=None, help="Output prefix.")
    p.add_argument("--ref", default=None, help="Reference FASTA (for left-align).")
    p.add_argument("--regions", default=None)
    p.add_argument("--no-split-multi", action="store_true")
    p.add_argument("--no-left-align",  action="store_true")
    p.add_argument("--no-sort",        action="store_true")
    p.add_argument("--no-set-id",      action="store_true")
    p.add_argument("--no-index",       action="store_true")
    p.add_argument("--filter-expr", default=None,
                   help="bcftools view filter expression, e.g. 'QUAL>=30 && DP>=10'")
    p.add_argument("--remove-filtered", action="store_true")
    p.add_argument("--keep-temps",      action="store_true")
    p.add_argument("--fill-tags",       action="store_true")
    _add_common_args(p)
    p.set_defaults(func=cmd_vcf_prep)


# ============================================================================
# Subcommand: vcf-qc
# ============================================================================

def cmd_vcf_qc(args: argparse.Namespace) -> int:
    from plantvarfilter.vcf_quality import VCFQualityChecker

    log = make_logger(args)
    qc = VCFQualityChecker(
        max_sites_scan=args.max_sites,
        min_sites_required=args.min_sites,
    )
    report = qc.evaluate(_abspath(args.vcf), log_fn=log)

    text = qc.to_text(report, vcf_path=_abspath(args.vcf))
    print(text)

    if args.out_text:
        Path(_abspath(args.out_text)).write_text(text, encoding="utf-8")
    if args.out_json:
        Path(_abspath(args.out_json)).write_text(
            json.dumps(asdict(report), default=str, indent=2),
            encoding="utf-8",
        )
    print(f"\n[VCF-QAScore] {report.score:.1f}  Verdict: {report.verdict}")
    return _exit_ok("VCF QC complete.")


def _parser_vcf_qc(sp: argparse._SubParsersAction) -> None:
    p = sp.add_parser("vcf-qc", help="VCF-QAScore quality assessment (0-100).")
    p.add_argument("--vcf", required=True)
    p.add_argument("--out-text", default=None, help="Write human-readable report.")
    p.add_argument("--out-json", default=None, help="Write JSON report.")
    p.add_argument("--max-sites", type=int, default=200_000,
                   help="Maximum sites to scan.")
    p.add_argument("--min-sites", type=int, default=5_000,
                   help="Minimum sites required for full score.")
    _add_common_args(p)
    p.set_defaults(func=cmd_vcf_qc)


# ============================================================================
# Subcommand: vcf2bed
# ============================================================================

def cmd_vcf2bed(args: argparse.Namespace) -> int:
    from plantvarfilter.gwas_pipeline import GWAS

    gwas = GWAS()
    out_prefix = _abspath(args.out)
    _ensure_dir(os.path.dirname(out_prefix))
    msg = gwas.vcf_to_bed(
        vcf_file=_abspath(args.vcf),
        id_file=_abspath(args.keep),
        file_out=out_prefix,
        maf=args.maf,
        geno=args.geno,
    )
    print(msg)
    return _exit_ok(f"PLINK BED prefix: {out_prefix}")


def _parser_vcf2bed(sp: argparse._SubParsersAction) -> None:
    p = sp.add_parser("vcf2bed",
                      help="Convert VCF to PLINK BED with MAF/geno filtering.")
    p.add_argument("--vcf", required=True)
    p.add_argument("--out", required=True, help="Output PLINK prefix.")
    p.add_argument("--maf",  type=float, default=0.05)
    p.add_argument("--geno", type=float, default=0.10)
    p.add_argument("--keep", default=None, help="--keep ID file (FID IID).")
    _add_common_args(p)
    p.set_defaults(func=cmd_vcf2bed)


# ============================================================================
# Subcommand: gwas
# ============================================================================

GWAS_ALGORITHMS = ["FaST-LMM", "Linear regression", "RandomForest", "XGBoost",
                   "Ridge", "GLM-PLINK2", "SAIGE"]


def cmd_gwas(args: argparse.Namespace) -> int:
    log = make_logger(args)
    bed_prefix = _bed_prefix(_abspath(args.bed))
    bed_path   = f"{bed_prefix}.bed"
    out_dir    = _ensure_dir(args.out) or os.getcwd()
    out_csv    = os.path.join(out_dir, "gwas_results.csv")
    manh_png   = os.path.join(out_dir, "manhatten_plot.png")
    qq_png     = os.path.join(out_dir, "qq_plot.png")

    algo = args.algo

    # GLM-PLINK2 and SAIGE bypass pysnptools/fastlmm and operate on file paths.
    if algo == "GLM-PLINK2":
        from plantvarfilter.gwas_pipeline import GWAS
        gwas = GWAS()
        result = gwas.run_gwas_glm_plink2(
            bed_file=bed_path,
            pheno_file=_abspath(args.pheno),
            cov_file=_abspath(args.cov),
            out_csv=out_csv,
            add_log=log,
        )
        _dump_json({"results_csv": out_csv, "details": str(result)})
        return _exit_ok("GLM-PLINK2 GWAS complete.")

    if algo == "SAIGE":
        from plantvarfilter.gwas_pipeline import GWAS
        gwas = GWAS()
        result = gwas.run_gwas_saige(
            bed_file=bed_path,
            pheno_file=_abspath(args.pheno),
            cov_file=_abspath(args.cov),
            out_csv=out_csv,
            add_log=log,
        )
        _dump_json({"results_csv": out_csv, "details": str(result)})
        return _exit_ok("SAIGE GWAS complete.")

    # All other algorithms share the pysnptools-based loading path.
    gwas, bed, bed_fixed, pheno, chrom_mapping = _load_bed_pheno(
        bed_prefix, _abspath(args.pheno))

    df_result = None
    df_plot   = None

    if algo in ("FaST-LMM", "Linear regression"):
        df_result, df_plot = gwas.run_gwas_lmm(
            bed_fixed=bed_fixed,
            pheno=pheno,
            chrom_mapping=chrom_mapping,
            add_log=log,
            gwas_result_name=out_csv,
            algorithm=algo,
            bed_file=bed_path,
            cov_file=_abspath(args.cov),
            gb_goal=args.gb_goal,
            kinship_path=_abspath(args.kinship),
        )

    elif algo == "XGBoost":
        df_result, df_plot = gwas.run_gwas_xg(
            bed_fixed=bed_fixed, pheno=pheno, bed_file=bed_path,
            test_size=args.test_size, estimators=args.estimators,
            gwas_result_name=out_csv, chrom_mapping=chrom_mapping,
            add_log=log, model_nr=args.models,
            max_dep_set=args.max_depth, nr_jobs=args.jobs,
            method=args.aggregation,
        )

    elif algo == "RandomForest":
        df_result, df_plot = gwas.run_gwas_rf(
            bed_fixed=bed_fixed, pheno=pheno, bed_file=bed_path,
            test_size=args.test_size, estimators=args.estimators,
            gwas_result_name=out_csv, chrom_mapping=chrom_mapping,
            add_log=log, model_nr=args.models, nr_jobs=args.jobs,
            method=args.aggregation,
        )

    elif algo == "Ridge":
        df_result, df_plot = gwas.run_gwas_ridge(
            bed_fixed=bed_fixed, pheno=pheno, bed_file=bed_path,
            test_size=args.test_size, alpha=args.alpha,
            gwas_result_name=out_csv, chrom_mapping=chrom_mapping,
            add_log=log, model_nr=args.models,
            method=args.aggregation,
        )

    else:
        raise ValueError(f"Unknown algorithm: {algo}")

    # Plot Manhattan + QQ
    if not args.no_plots and df_plot is not None:
        gwas.plot_gwas(
            df=df_plot,
            limit=args.plot_limit,
            algorithm=algo,
            manhatten_plot_name=manh_png,
            qq_plot_name=qq_png,
            chrom_mapping=chrom_mapping,
            region=args.region,
            region_only_csv=os.path.join(out_dir, "gwas_region.csv") if args.region else None,
            title_suffix=args.title_suffix,
        )

    _dump_json({
        "results_csv": out_csv,
        "manhattan":   manh_png if not args.no_plots else None,
        "qq":          qq_png   if not args.no_plots else None,
        "rows":        int(len(df_result)) if df_result is not None else 0,
    })
    return _exit_ok(f"GWAS ({algo}) complete.")


def _parser_gwas(sp: argparse._SubParsersAction) -> None:
    p = sp.add_parser("gwas", help="Single-trait GWAS.")
    p.add_argument("--bed", required=True, help="PLINK prefix or .bed/.bim/.fam path.")
    p.add_argument("--pheno", required=True, help="Phenotype file (FID IID Value).")
    p.add_argument("--algo", required=True, choices=GWAS_ALGORITHMS,
                   help="Statistical model to use.")
    p.add_argument("--out", required=True, help="Output directory for results & plots.")
    p.add_argument("--cov", default=None, help="Covariates file (optional).")

    # FaST-LMM / linear options
    p.add_argument("--kinship", default=None, help="Pre-computed kinship .npy.")
    p.add_argument("--gb-goal", type=int, default=0,
                   help="FaST-LMM RAM target in GB (0 = auto).")

    # ML options (RF/XGBoost/Ridge)
    p.add_argument("--test-size",   type=float, default=0.3)
    p.add_argument("--estimators",  type=int,   default=200)
    p.add_argument("--models",      type=int,   default=1, dest="models",
                   help="Number of model iterations (cross-validation folds).")
    p.add_argument("--max-depth",   type=int,   default=3)
    p.add_argument("--alpha",       type=float, default=1.0, help="Ridge alpha.")
    p.add_argument("--jobs",        type=int,   default=-1, help="Parallel jobs (-1 = all).")
    p.add_argument("--aggregation", default="sum",
                   choices=["sum", "mean", "median", "max", "min"],
                   help="How to aggregate per-iteration scores.")

    # Plotting options
    p.add_argument("--no-plots",   action="store_true")
    p.add_argument("--plot-limit", type=float, default=5e-8,
                   help="Genome-wide significance threshold line.")
    p.add_argument("--region", default=None,
                   help="Plot region only, e.g. '1:1000000-2000000' or 'chr5'.")
    p.add_argument("--title-suffix", default=None)
    _add_common_args(p)
    p.set_defaults(func=cmd_gwas)


# ============================================================================
# Subcommand: batch-gwas
# ============================================================================

def cmd_batch_gwas(args: argparse.Namespace) -> int:
    from plantvarfilter.batch_gwas import run_batch_gwas_for_all_traits
    from plantvarfilter.gwas_pipeline import GWAS
    from plantvarfilter.helpers import HELPERS

    log = make_logger(args)
    bed_prefix = _bed_prefix(_abspath(args.bed))
    out_dir = _ensure_dir(args.out)

    summary = run_batch_gwas_for_all_traits(
        gwas=GWAS(),
        helper=HELPERS(),
        bed_path=f"{bed_prefix}.bed",
        pheno_path=_abspath(args.pheno),
        cov_path=_abspath(args.cov),
        algorithm=args.algo,
        out_dir=out_dir,
        log_fn=log,
        nr_jobs=args.jobs,
        gb_goal=args.gb_goal,
        train_size=args.test_size,
        estimators=args.estimators,
        model_nr=args.models,
        max_depth=args.max_depth,
        aggregation_method=args.aggregation,
        snp_limit=args.snp_limit,
    )
    _dump_json(summary)
    return _exit_ok("Batch GWAS complete.")


def _parser_batch_gwas(sp: argparse._SubParsersAction) -> None:
    p = sp.add_parser("batch-gwas", help="GWAS over every trait column.")
    p.add_argument("--bed",   required=True)
    p.add_argument("--pheno", required=True,
                   help="Multi-trait phenotype file: FID,IID,Trait1,Trait2,...")
    p.add_argument("--algo",  required=True, choices=GWAS_ALGORITHMS)
    p.add_argument("--out",   required=True)
    p.add_argument("--cov",   default=None)
    p.add_argument("--jobs",  type=int, default=-1)
    p.add_argument("--gb-goal", type=int, default=0)
    p.add_argument("--test-size",  type=float, default=0.3)
    p.add_argument("--estimators", type=int,   default=200)
    p.add_argument("--models",     type=int,   default=1)
    p.add_argument("--max-depth",  type=int,   default=3)
    p.add_argument("--aggregation", default="sum",
                   choices=["sum", "mean", "median", "max", "min"])
    p.add_argument("--snp-limit", type=int, default=None,
                   help="Optional cap on number of SNPs (debug/speed).")
    _add_common_args(p)
    p.set_defaults(func=cmd_batch_gwas)


# ============================================================================
# Subcommand: gwas-plot
# ============================================================================

def cmd_gwas_plot(args: argparse.Namespace) -> int:
    import pandas as pd
    from plantvarfilter.gwas_pipeline import GWAS

    log = make_logger(args)
    df = pd.read_csv(_abspath(args.gwas_csv))
    out_dir  = _ensure_dir(args.out) or os.getcwd()
    manh_png = os.path.join(out_dir, "manhatten_plot.png")
    qq_png   = os.path.join(out_dir, "qq_plot.png")

    gwas = GWAS()
    gwas.plot_gwas(
        df=df,
        limit=args.limit,
        algorithm=args.algo,
        manhatten_plot_name=manh_png,
        qq_plot_name=qq_png,
        chrom_mapping={},
        region=args.region,
        region_only_csv=os.path.join(out_dir, "gwas_region.csv") if args.region else None,
        title_suffix=args.title_suffix,
    )
    _dump_json({"manhattan": manh_png, "qq": qq_png})
    return _exit_ok("Plots written.")


def _parser_gwas_plot(sp: argparse._SubParsersAction) -> None:
    p = sp.add_parser("gwas-plot", help="Manhattan + QQ from a GWAS results CSV.")
    p.add_argument("--gwas-csv", required=True)
    p.add_argument("--out",      required=True, help="Output directory.")
    p.add_argument("--algo",     default="FaST-LMM",
                   help="Algorithm label (used for axis/title).")
    p.add_argument("--limit", type=float, default=5e-8,
                   help="Genome-wide significance threshold.")
    p.add_argument("--region", default=None)
    p.add_argument("--title-suffix", default=None)
    _add_common_args(p)
    p.set_defaults(func=cmd_gwas_plot)


# ============================================================================
# Subcommand: annotate
# ============================================================================

def cmd_annotate(args: argparse.Namespace) -> int:
    import pandas as pd
    from plantvarfilter.annotation_utils import Annotator

    df = pd.read_csv(_abspath(args.gwas_csv))

    ann = Annotator()
    ann.load_gtf_or_gff(_abspath(args.gtf))
    ann.build_index()

    out_csv = _abspath(args.out) or os.path.splitext(_abspath(args.gwas_csv))[0] + "_annotated.csv"
    written = ann.annotate_and_save(
        gwas_df=df,
        out_csv=out_csv,
        window_kb=args.window_kb,
        chr_col=args.chr_col,
        pos_col=args.pos_col,
    )
    _dump_json({"annotated_csv": written})
    return _exit_ok("Annotation written.")


def _parser_annotate(sp: argparse._SubParsersAction) -> None:
    p = sp.add_parser("annotate", help="Annotate GWAS CSV against a GTF/GFF.")
    p.add_argument("--gwas-csv", required=True)
    p.add_argument("--gtf",      required=True, help="GTF or GFF (gz ok).")
    p.add_argument("--out",      default=None,  help="Output CSV path.")
    p.add_argument("--window-kb", type=int, default=50)
    p.add_argument("--chr-col",  default="Chr")
    p.add_argument("--pos-col",  default="ChrPos")
    _add_common_args(p)
    p.set_defaults(func=cmd_annotate)


# ============================================================================
# Subcommand: predict
# ============================================================================

PREDICT_MODELS = ["LMM", "RandomForest", "XGBoost", "Ridge"]


def cmd_predict(args: argparse.Namespace) -> int:
    from plantvarfilter.genomic_prediction_pipeline import GenomicPrediction

    log = make_logger(args)
    bed_prefix = _bed_prefix(_abspath(args.bed))
    bed_path   = f"{bed_prefix}.bed"
    out_dir    = _ensure_dir(args.out) or os.getcwd()
    out_csv    = os.path.join(out_dir, args.out_name or "gp_predictions.csv")

    _, _, bed_fixed, pheno, chrom_mapping = _load_bed_pheno(
        bed_prefix, _abspath(args.pheno))

    gp = GenomicPrediction()

    if args.model == "LMM":
        df = gp.run_lmm_gp(
            bed_fixed=bed_fixed, pheno=pheno,
            genomic_predict_name=out_csv,
            model_nr=args.models, add_log=log,
            bed_file=bed_path, chrom_mapping=chrom_mapping,
        )
    elif args.model == "RandomForest":
        df = gp.run_gp_rf(
            bed_fixed=bed_fixed, pheno=pheno, bed_file=bed_path,
            test_size=args.test_size, estimators=args.estimators,
            genomic_predict_name=out_csv, chrom_mapping=chrom_mapping,
            add_log=log, model_nr=args.models, nr_jobs=args.jobs,
        )
    elif args.model == "XGBoost":
        df = gp.run_gp_xg(
            bed_fixed=bed_fixed, pheno=pheno, bed_file=bed_path,
            test_size=args.test_size, estimators=args.estimators,
            genomic_predict_name=out_csv, chrom_mapping=chrom_mapping,
            add_log=log, model_nr=args.models, nr_jobs=args.jobs,
        )
    elif args.model == "Ridge":
        df = gp.run_gp_ridge(
            bed_fixed=bed_fixed, pheno=pheno, bed_file=bed_path,
            test_size=args.test_size, alpha=args.alpha,
            genomic_predict_name=out_csv, chrom_mapping=chrom_mapping,
            add_log=log, model_nr=args.models,
        )
    else:
        raise ValueError(f"Unknown model: {args.model}")

    # Optional plots
    if not args.no_plots:
        ba_png = os.path.join(out_dir, "Bland_Altman_plot.png")
        sc_png = os.path.join(out_dir, "GP_scatter_plot.png")
        try:
            gp.plot_gp(df, ba_png, args.model)
            gp.plot_gp_scatter(df, sc_png, args.model)
        except Exception as e:
            log(f"Plotting skipped: {e}", warn=True)

    _dump_json({
        "predictions_csv": out_csv,
        "rows":           int(len(df)) if df is not None else 0,
    })
    return _exit_ok(f"Genomic prediction ({args.model}) complete.")


def _parser_predict(sp: argparse._SubParsersAction) -> None:
    p = sp.add_parser("predict", help="Genomic prediction (GEBVs).")
    p.add_argument("--bed",   required=True)
    p.add_argument("--pheno", required=True)
    p.add_argument("--model", required=True, choices=PREDICT_MODELS)
    p.add_argument("--out",   required=True, help="Output directory.")
    p.add_argument("--out-name", default=None,
                   help="Output CSV filename (default: gp_predictions.csv).")
    p.add_argument("--test-size",  type=float, default=0.3)
    p.add_argument("--estimators", type=int,   default=200)
    p.add_argument("--models",     type=int,   default=5,
                   help="Cross-validation folds / iterations.")
    p.add_argument("--alpha",      type=float, default=1.0)
    p.add_argument("--jobs",       type=int,   default=-1)
    p.add_argument("--no-plots",   action="store_true")
    _add_common_args(p)
    p.set_defaults(func=cmd_predict)


# ============================================================================
# Subcommand: ld-decay
# ============================================================================

def cmd_ld_decay(args: argparse.Namespace) -> int:
    from plantvarfilter.ld_utils import LDAnalyzer

    out_prefix = _abspath(args.out)
    _ensure_dir(os.path.dirname(out_prefix))

    ana = LDAnalyzer()
    res = ana.ld_decay(
        out_prefix=out_prefix,
        bfile_prefix=_bed_prefix(_abspath(args.bed)) if args.bed else None,
        vcf_path=_abspath(args.vcf),
        window_kb=args.window_kb,
        window_snps=args.window_snps,
        max_dist_kb=args.max_dist_kb,
        min_r2=args.min_r2,
        region=args.region,
        keep_samples=_abspath(args.keep_samples),
    )
    _dump_json(res)
    return _exit_ok("LD decay computed.")


def _parser_ld_decay(sp: argparse._SubParsersAction) -> None:
    p = sp.add_parser("ld-decay", help="LD decay (mean r^2 vs distance).")
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--bed", help="PLINK prefix.")
    src.add_argument("--vcf", help="VCF as input (PLINK will read it directly).")
    p.add_argument("--out", required=True, help="Output prefix.")
    p.add_argument("--window-kb",   type=int,   default=500)
    p.add_argument("--window-snps", type=int,   default=5000)
    p.add_argument("--max-dist-kb", type=int,   default=1000)
    p.add_argument("--min-r2",      type=float, default=0.1)
    p.add_argument("--region", default=None)
    p.add_argument("--keep-samples", default=None)
    _add_common_args(p)
    p.set_defaults(func=cmd_ld_decay)


# ============================================================================
# Subcommand: ld-heatmap
# ============================================================================

def cmd_ld_heatmap(args: argparse.Namespace) -> int:
    from plantvarfilter.ld_utils import LDAnalyzer

    out_prefix = _abspath(args.out)
    _ensure_dir(os.path.dirname(out_prefix))

    ana = LDAnalyzer()
    res = ana.ld_heatmap(
        out_prefix=out_prefix,
        bfile_prefix=_bed_prefix(_abspath(args.bed)) if args.bed else None,
        vcf_path=_abspath(args.vcf),
        region=args.region,
        window_snps=args.window_snps,
        min_r2=args.min_r2,
        keep_samples=_abspath(args.keep_samples),
        max_snps_for_plot=args.max_snps,
    )
    _dump_json(res)
    return _exit_ok("LD heatmap written.")


def _parser_ld_heatmap(sp: argparse._SubParsersAction) -> None:
    p = sp.add_parser("ld-heatmap", help="LD heatmap (square r^2 matrix).")
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--bed", help="PLINK prefix.")
    src.add_argument("--vcf", help="VCF input.")
    p.add_argument("--out", required=True, help="Output prefix.")
    p.add_argument("--region", default=None,
                   help="Region to plot, e.g. '1:1000000-2000000'.")
    p.add_argument("--window-snps", type=int,   default=1000)
    p.add_argument("--min-r2",      type=float, default=0.0)
    p.add_argument("--max-snps",    type=int,   default=500,
                   help="Down-sample to this many SNPs for plotting.")
    p.add_argument("--keep-samples", default=None)
    _add_common_args(p)
    p.set_defaults(func=cmd_ld_heatmap)


# ============================================================================
# Subcommand: pangenome
# ============================================================================

def cmd_pangenome(args: argparse.Namespace) -> int:
    from plantvarfilter.pangenome_builder import build_pangenome_graph

    log = make_logger(args)
    res = build_pangenome_graph(
        assemblies_input=_abspath(args.assemblies),
        output_dir=_ensure_dir(args.out),
        mode=args.mode,
        subset_n=args.subset_n,
        threads=args.threads,
        minigraph_preset=args.preset,
        backbone_strategy=args.backbone,
        logger=log,
        vg_gzip_only=args.vg_gzip_only,
    )
    _dump_json(res)
    return _exit_ok("Pangenome built.")


def _parser_pangenome(sp: argparse._SubParsersAction) -> None:
    p = sp.add_parser("pangenome",
                      help="Build a graph-based pangenome (minigraph + vg).")
    p.add_argument("--assemblies", required=True,
                   help="Multi-FASTA file OR a directory containing FASTAs.")
    p.add_argument("--out",        required=True, help="Output directory.")
    p.add_argument("--mode", choices=["full", "fast"], default="full")
    p.add_argument("--subset-n", type=int, default=25,
                   help="Top-N largest inputs in 'fast' mode.")
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--preset", default="ggs", help="minigraph preset (-x).")
    p.add_argument("--backbone", default="largest",
                   choices=["largest", "first"], help="Backbone selection strategy.")
    p.add_argument("--vg-gzip-only", action="store_true")
    _add_common_args(p)
    p.set_defaults(func=cmd_pangenome)


# ============================================================================
# Top-level dispatcher
# ============================================================================

def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="plantomicsgwas",
        description="Headless CLI for the PlantOmicsGWAS / plantvarfilter pipeline.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  plantomicsgwas ref-index --fasta ref.fa --out ref/\n"
            "  plantomicsgwas vcf-qc    --vcf calls.vcf.gz --out-json qc.json\n"
            "  plantomicsgwas vcf2bed   --vcf calls.vcf.gz --out data --maf 0.05\n"
            "  plantomicsgwas gwas      --bed data --pheno pheno.txt "
            "--algo FaST-LMM --out results/\n"
        ),
    )
    parser.add_argument("--version", action="store_true",
                        help="Print package version and exit.")

    sp = parser.add_subparsers(dest="subcommand", metavar="<subcommand>")

    _parser_ref_index(sp)
    _parser_fastq_qc(sp)
    _parser_align(sp)
    _parser_preprocess_bam(sp)
    _parser_call(sp)
    _parser_vcf_prep(sp)
    _parser_vcf_qc(sp)
    _parser_vcf2bed(sp)
    _parser_gwas(sp)
    _parser_batch_gwas(sp)
    _parser_gwas_plot(sp)
    _parser_annotate(sp)
    _parser_predict(sp)
    _parser_ld_decay(sp)
    _parser_ld_heatmap(sp)
    _parser_pangenome(sp)

    return parser

def main(argv: Optional[List[str]] = None) -> int:
    # Make the engine's logging.info() messages visible by default.
    logging.basicConfig(
        level=logging.INFO,
        format="[%(levelname)s] %(message)s",
        stream=sys.stdout,
    )

    parser = _build_parser()
    args = parser.parse_args(argv)

    if getattr(args, "version", False):
        try:
            from plantvarfilter import __version__
            print(__version__)
        except Exception:
            print("unknown")
        return 0

    if not getattr(args, "subcommand", None):
        parser.print_help(sys.stderr)
        return 1

    try:
        return args.func(args)
    except KeyboardInterrupt:
        sys.stderr.write("\n[ABORT] Interrupted by user.\n")
        return 130
    except Exception as e:
        return _exit_err(e, debug=getattr(args, "debug", False))


if __name__ == "__main__":
    sys.exit(main())