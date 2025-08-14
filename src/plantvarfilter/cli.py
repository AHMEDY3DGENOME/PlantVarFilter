# src/plantvarfilter/cli.py

import argparse, gzip, sys, os, json, time, logging, re
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

from scipy import stats
import pyarrow.feather as feather

from plantvarfilter.filter import improved_filter_variants
from plantvarfilter.annotator import build_gene_db, annotate_variants_with_genes
from plantvarfilter.regression_gwas import run_snp_gwas_matrix
from plantvarfilter.parser import read_gene_traits


# ---------------------------------------------------------------------
# logging
# ---------------------------------------------------------------------
log_file = None
logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler(sys.stdout)]
)


# ---------------------------------------------------------------------
# helpers: config / IO
# ---------------------------------------------------------------------
def load_config_file(path):
    if not path:
        logging.error("No config file path provided.")
        sys.exit(1)
    if not os.path.exists(path):
        logging.error(f"Config file not found: {path}")
        sys.exit(1)
    try:
        with open(path, "r") as f:
            return json.load(f)
    except Exception as e:
        logging.error(f"Failed to load config file: {e}")
        sys.exit(1)


def _to_path(config, key, required=True):
    val = config.get(key)
    if not val or str(val).strip() == "":
        if required:
            logging.error(f"Missing '{key}' in config.")
            sys.exit(1)
        return None
    return Path(val).expanduser().resolve()


def initialize_user_data(path: str):
    base = Path(path).expanduser().resolve()
    (base / "input").mkdir(parents=True, exist_ok=True)
    (base / "output").mkdir(parents=True, exist_ok=True)
    cfg = {
        "vcf": str(base / "input/your_data.vcf.gz"),
        "gff": "",
        "traits": str(base / "input/your_traits.csv"),
        "include_intergenic": True,
        "consequence_types": None,
        "output_format": "csv",
        "output": str(base / "output/filtered_variants.csv"),
        "output_dir": str(base / "output"),
        "plot": True,
        "gwas": True,
        "force_snp_gwas": True,

        "trait_missing_values": [-9, "-9", "NA", "na", ""],
        "sample_id_column": "Sample",
        "maf_min": 0.01,
        "callrate_min": 0.9
    }
    with open(base / "config.json", "w") as f:
        json.dump(cfg, f, indent=2)
    print(f"Project initialized at {base}")


# ---------------------------------------------------------------------
# plotting helpers
# ---------------------------------------------------------------------
def generate_plots(df: pd.DataFrame, out_dir: Path):
    plot_dir = out_dir / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)

    if "Consequence" in df.columns and not df["Consequence"].empty:
        plt.figure(figsize=(10, 6))
        sns.countplot(y=df["Consequence"], order=df["Consequence"].value_counts().index)
        plt.title("Distribution of Variant Consequences")
        plt.tight_layout()
        plt.savefig(plot_dir / "consequence_distribution.png")
        plt.close()

    if "Variant_Type" in df.columns and not df["Variant_Type"].empty:
        plt.figure(figsize=(6, 6))
        df["Variant_Type"].value_counts().plot.pie(autopct='%1.1f%%')
        plt.title("Variant Type Proportions")
        plt.ylabel("")
        plt.tight_layout()
        plt.savefig(plot_dir / "variant_type_pie.png")
        plt.close()


def plot_manhattan_coord(df: pd.DataFrame, out_png: Path, title="Manhattan Plot (SNP-level)"):
    if df.empty or "P" not in df.columns:
        return
    d = df.copy()
    d["CHR"] = d["CHROM"].astype(str)
    d["CHR"] = d["CHR"].str.replace("^chr", "", regex=True)

    def chr_key(c):
        try:
            return (0, int(re.sub(r"\D", "", c)))
        except Exception:
            return (1, c)

    chrs = sorted(d["CHR"].unique(), key=chr_key)
    offsets = {}
    cum = 0
    tick_pos, tick_lbl = [], []
    for c in chrs:
        sub = d[d["CHR"] == c]
        if sub.empty:
            continue
        offsets[c] = cum
        tick_pos.append(cum + sub["POS"].median())
        tick_lbl.append(c)
        cum += sub["POS"].max() + 1

    d["BPcum"] = d.apply(lambda r: r["POS"] + offsets.get(r["CHR"], 0), axis=1)
    d.sort_values("BPcum", inplace=True)

    plt.figure(figsize=(12, 6))
    y = -np.log10(d["P"].clip(lower=np.finfo(float).tiny))
    plt.scatter(d["BPcum"], y, s=6)
    plt.xticks(tick_pos, tick_lbl, rotation=0)
    plt.xlabel("Chromosome")
    plt.ylabel("-log10(P-value)")
    plt.title(title)
    plt.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_png)
    plt.close()


# ---------------------------------------------------------------------
# VCF → dosage + meta
# ---------------------------------------------------------------------
def _open_vcf_text(p: Path):
    return gzip.open(p, "rt", encoding="utf-8") if str(p).endswith(".gz") else open(p, "rt", encoding="utf-8")


def vcf_to_dosage_matrix(vcf_path: Path, max_variants=None):
    samples, rows, meta = [], [], []
    header_found = False
    with _open_vcf_text(vcf_path) as f:
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header_found = True
                hdr = line.strip().lstrip("#").split("\t")
                samples = hdr[9:]
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 10:
                continue
            chrom, pos, vid, ref, alt, qual, flt, info, fmt = parts[:9]
            smp = parts[9:]
            fmt_keys = fmt.split(":")
            try:
                gt_idx = fmt_keys.index("GT")
            except ValueError:
                continue

            dos = []
            for sf in smp:
                fields = sf.split(":")
                gt = fields[gt_idx] if gt_idx < len(fields) else "."
                if gt is None or gt == "." or gt.startswith("."):
                    dos.append(np.nan)
                    continue
                a = gt.replace("|", "/").split("/")
                if len(a) != 2:
                    dos.append(np.nan)
                    continue
                try:
                    x = (0 if a[0] == "." else int(a[0])) + (0 if a[1] == "." else int(a[1]))
                    if "." in a:
                        dos.append(np.nan)
                    else:
                        dos.append(float(x))
                except Exception:
                    dos.append(np.nan)

            key = f"{chrom}:{pos}:{ref}:{alt}"
            rows.append(dos)
            meta.append((chrom, int(pos), ref, alt, vid if vid != "." else key))
            if max_variants and len(rows) >= max_variants:
                break

    if not rows:
        return pd.DataFrame(), pd.DataFrame(), []

    if not header_found or len(samples) < len(rows[0]):
        samples = [f"S{i + 1}" for i in range(len(rows[0]))]

    X = pd.DataFrame(rows, columns=samples)
    idx = [f"{m[0]}:{m[1]}:{m[2]}:{m[3]}" for m in meta]
    X.index = idx
    meta_df = pd.DataFrame(meta, columns=["CHROM", "POS", "REF", "ALT", "ID"], index=idx)
    return X, meta_df, samples


# ---------------------------------------------------------------------
# traits cleaning / mapping + QC helpers
# ---------------------------------------------------------------------
def clean_traits(df: pd.DataFrame, sample_col: str, missing_values):
    t = df.copy()
    if sample_col not in t.columns:
        for cand in ["Sample", "ID", "Name"]:
            if cand in t.columns:
                t.rename(columns={cand: sample_col}, inplace=True)
                break
    if "Trait_Score" in t.columns:
        t["Trait_Score"] = pd.to_numeric(t["Trait_Score"], errors="coerce")
        numeric_missing = {
            float(x) for x in missing_values
            if str(x).replace('.', '', 1).lstrip('-').isdigit()
        }
        t.loc[t["Trait_Score"].isin(numeric_missing), "Trait_Score"] = np.nan
        t.loc[t["Trait_Score"] == -9.0, "Trait_Score"] = np.nan
        t = t.dropna(subset=["Trait_Score"])
    t[sample_col] = t[sample_col].astype(str).str.strip().str.upper()
    return t


def compute_maf_callrate(X: pd.DataFrame):
    n_called = X.notna().sum(axis=1)
    callrate = n_called / X.shape[1]
    denom = 2 * n_called.replace(0, np.nan)
    p = X.sum(axis=1) / denom
    maf = np.minimum(p, 1 - p)
    return maf, callrate


# ---------------------------------------------------------------------
# SNP-level GWAS with QC (delegates to run_snp_gwas_matrix)
# ---------------------------------------------------------------------
def run_snp_gwas_with_qc(vcf_path: Path,
                         traits_df: pd.DataFrame,
                         sample_col: str,
                         missing_values,
                         maf_min: float,
                         callrate_min: float,
                         out_dir: Path):
    X, meta, vcf_samples = vcf_to_dosage_matrix(vcf_path)
    if X.empty:
        logging.error("No variants parsed from VCF for SNP-GWAS.")
        sys.exit(1)

    t = clean_traits(traits_df, sample_col=sample_col, missing_values=missing_values)
    if sample_col not in t.columns or "Trait_Score" not in t.columns:
        logging.error("Traits file must contain columns for sample IDs and Trait_Score after parsing.")
        sys.exit(1)

    vcf_cols_upper = [str(s).upper() for s in X.columns]
    t_samples_set = set(t[sample_col].astype(str).str.upper())
    name_map = {s: su for s, su in zip(X.columns, vcf_cols_upper)}
    inter = [s for s in X.columns if name_map[s] in t_samples_set]
    if not inter:
        logging.error("No overlapping samples between VCF header and traits. Check sample naming (e.g. V001..).")
        sys.exit(1)

    y = t.set_index(sample_col).loc[[name_map[s] for s in inter], "Trait_Score"].astype(float)
    X = X[inter]

    maf, callrate = compute_maf_callrate(X)
    keep = (callrate >= float(callrate_min)) & (maf >= float(maf_min))
    Xq = X.loc[keep]
    meta_q = meta.loc[Xq.index]
    logging.info(f"SNP QC: kept {Xq.shape[0]}/{X.shape[0]} variants (MAF≥{maf_min}, callrate≥{callrate_min}).")
    pd.DataFrame({
        "SNP": Xq.index,
        "MAF": maf.loc[keep].values,
        "CallRate": callrate.loc[keep].values,
        "CHROM": meta_q["CHROM"].values,
        "POS": meta_q["POS"].values
    }).to_csv(out_dir / "gwas_snp_qc_variants.csv", index=False)

    # run GWAS via regression module (writes results + QQ internally)
    res = run_snp_gwas_matrix(
        X=Xq,          # rows = SNP, cols = samples (already intersected)
        meta=meta_q,
        y=y,           # pandas Series indexed by sample
        covariates=None,
        out_dir=str(out_dir)
    )

    # Manhattan (from returned results joined with meta)
    jj = res.dropna(subset=["P"])
    if not jj.empty and {"CHROM", "POS", "P"}.issubset(jj.columns):
        plot_manhattan_coord(
            jj[["CHROM", "POS", "P"]],
            out_png=(out_dir / "plots" / "manhattan_plot_snp.png"),
            title="Manhattan Plot (SNP-level)"
        )
        logging.info("Manhattan plot written.")


# ---------------------------------------------------------------------
# pipeline
# ---------------------------------------------------------------------
def run_pipeline(config):
    get = lambda k, d=None: config.get(k, d)

    try:
        vcf_path = _to_path(config, "vcf", required=True)
        gff_path = _to_path(config, "gff", required=False)
        traits_path = _to_path(config, "traits", required=False)
    except Exception as e:
        logging.error(f"Invalid path in config: {e}")
        sys.exit(1)

    for label, path, req in [("VCF", vcf_path, True), ("GFF", gff_path, False), ("Traits", traits_path, False)]:
        if req and not path.exists():
            logging.error(f"{label} file not found at: {path}")
            sys.exit(1)
        if not req and path and not path.exists():
            logging.warning(f"{label} file not found at: {path}")
            if label == "GFF":
                gff_path = None
            if label == "Traits":
                traits_path = None

    output_path = Path(get("output"))
    out_dir = Path(get("output_dir")) if get("output_dir") else output_path.parent
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "plots").mkdir(parents=True, exist_ok=True)

    global log_file
    log_file = out_dir / "run.log"
    logging.getLogger().handlers.clear()
    logging.basicConfig(
        level=logging.INFO,
        format='[%(asctime)s] %(levelname)s - %(message)s',
        handlers=[logging.StreamHandler(sys.stdout), logging.FileHandler(log_file)]
    )

    logging.info("Reading VCF...")
    open_vcf = (lambda p: gzip.open(p, 'rt')) if str(vcf_path).endswith(".gz") else (lambda p: open(p, "r", encoding="utf-8"))
    with open_vcf(vcf_path) as vcf_stream:
        feather_path = improved_filter_variants(
            vcf_stream,
            include_intergenic=get("include_intergenic", True),
            store_as_feather=True,
            consequence_types=get("consequence_types")
        )
    variants_df = pd.read_feather(feather_path)
    if variants_df.empty:
        logging.warning("No variants found after filtering.")
        sys.exit(1)

    annotated_df = variants_df
    if gff_path:
        logging.info("Building gene database...")
        open_gff = (lambda p: gzip.open(p, 'rt')) if str(gff_path).endswith(".gz") else (lambda p: open(p, "r", encoding="utf-8"))
        with open_gff(gff_path) as gff_stream:
            gene_db = build_gene_db(gff_stream)
        logging.info("Annotating variants with genes...")
        annotated_df = annotate_variants_with_genes(
            variants_df, gene_db, include_intergenic=get("include_intergenic", True)
        )
    else:
        logging.info("No GFF provided — skipping gene annotation.")

    traits_df = None
    if traits_path:
        logging.info("Reading trait data...")
        traits_df = read_gene_traits(traits_path)
        try:
            logging.info(f"Traits columns after parsing: {list(traits_df.columns)}")
            logging.info("Traits preview:\n" + str(traits_df.head().to_string(index=False)))
        except Exception:
            pass

    fmt = (get("output_format") or "csv").lower()
    if fmt == "csv":
        annotated_df.to_csv(output_path, index=False)
    elif fmt == "tsv":
        annotated_df.to_csv(output_path, sep="\t", index=False)
    elif fmt == "json":
        annotated_df.to_json(output_path, orient="records", lines=True)
    elif fmt == "feather":
        feather.write_feather(annotated_df, output_path)
    elif fmt == "xlsx":
        annotated_df.to_excel(output_path, index=False)
    else:
        logging.warning(f"Unknown output_format '{fmt}', defaulting to CSV")
        annotated_df.to_csv(output_path, index=False)

    if get("plot", True):
        logging.info("Generating plots...")
        generate_plots(annotated_df, out_dir)

    if get("gwas", False) and traits_df is not None and not traits_df.empty:
        logging.info("Running SNP-level GWAS with QC, Manhattan, and QQ...")
        run_snp_gwas_with_qc(
            vcf_path=vcf_path,
            traits_df=traits_df,
            sample_col=config.get("sample_id_column", "Sample"),
            missing_values=config.get("trait_missing_values", [-9, "-9", "NA", "na", ""]),
            maf_min=float(config.get("maf_min", 0.01)),
            callrate_min=float(config.get("callrate_min", 0.9)),
            out_dir=out_dir
        )


# ---------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser(description="PlantVarFilter - Variant Filtering for Plant Genomics")
    sp = p.add_subparsers(dest="command", required=True)

    p_init = sp.add_parser("init", help="Initialize a project in a custom directory")
    p_init.add_argument("path", type=str, help="Target directory")

    p_run = sp.add_parser("run", help="Run the full analysis pipeline")
    p_run.add_argument("--config", type=str, help="Path to config.json")

    p_plot = sp.add_parser("plot-only", help="Plot from an existing GWAS CSV (P column required)")
    p_plot.add_argument("--config", type=str, help="Path to config.json")
    p_plot.add_argument("--file", type=str, help="CSV path (overrides config)")

    args = p.parse_args()

    if args.command == "init":
        initialize_user_data(args.path)
    elif args.command == "run":
        cfg_path = args.config or (Path.home() / ".plantvarfilter_data" / "config.json")
        cfg = load_config_file(cfg_path)
        cfg["start_time"] = time.time()
        run_pipeline(cfg)
    elif args.command == "plot-only":
        cfg_path = args.config or (Path.home() / ".plantvarfilter_data" / "config.json")
        cfg = load_config_file(cfg_path)
        out_dir = Path(cfg.get("output_dir", "."))
        f = args.file or cfg.get("gwas_results")
        if not f or not os.path.exists(f):
            logging.error(f"GWAS results file not found: {f}")
            sys.exit(1)
        df = pd.read_csv(f)
        if "P" not in df.columns:
            logging.error("CSV must contain a 'P' column.")
            sys.exit(1)
        out_dir.mkdir(exist_ok=True, parents=True)
        plot_manhattan_coord(
            df[["CHROM", "POS", "P"]].dropna() if {"CHROM", "POS", "P"}.issubset(df.columns) else df,
            out_dir / "plots" / "manhattan_plot_from_file.png",
            "Manhattan Plot (from existing file)"
        )
        plt.close()
        logging.info("Plots written.")
    else:
        p.print_help()


if __name__ == "__main__":
    main()
