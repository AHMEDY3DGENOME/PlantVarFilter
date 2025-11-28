# ld_utils.py
import os
import re
import math
import csv
import subprocess
from typing import Optional, Tuple, List, Dict

try:
    from plantvarfilter.linux import resolve_tool
except Exception:
    from shutil import which as _which
    def resolve_tool(name: str) -> str:
        return _which(name) or name


class LDAnalysisError(Exception):
    pass


def _run(cmd: List[str]) -> str:
    p = subprocess.run(cmd, text=True, capture_output=True)
    if p.returncode != 0:
        raise LDAnalysisError(f"Command failed ({p.returncode}): {' '.join(cmd)}\n{p.stdout}\n{p.stderr}")
    return p.stdout or ""


def _parse_region(region: Optional[str]) -> Tuple[Optional[str], Optional[int], Optional[int]]:
    if not region or not str(region).strip():
        return None, None, None
    m = re.match(r"^(\w+):(\d+)[-:](\d+)$", region.replace(",", ""))
    if not m:
        raise LDAnalysisError("Region must look like: chr1:100000-200000")
    chrom, start, end = m.group(1), int(m.group(2)), int(m.group(3))
    if end < start:
        start, end = end, start
    return chrom, start, end


def _ensure_dir(path: str):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)


class LDAnalyzer:
    def __init__(self, plink_bin: Optional[str] = None):
        self.plink = plink_bin or resolve_tool("plink")

    def _plink_base(
        self,
        bfile_prefix: Optional[str] = None,
        vcf_path: Optional[str] = None,
        out_prefix: Optional[str] = None,
        keep_samples: Optional[str] = None,
        region: Optional[str] = None,
    ) -> List[str]:
        if not out_prefix:
            raise LDAnalysisError("out_prefix is required")
        _ensure_dir(out_prefix + ".x")
        cmd = [self.plink, "--allow-no-sex", "--silent", "--out", out_prefix]
        if bfile_prefix:
            cmd += ["--bfile", bfile_prefix]
        elif vcf_path:
            cmd += ["--vcf", vcf_path]
        else:
            raise LDAnalysisError("Either bfile_prefix or vcf_path must be provided")
        if keep_samples:
            cmd += ["--keep", keep_samples]
        if region:
            chrom, start, end = _parse_region(region)
            if chrom:
                cmd += ["--chr", str(chrom), "--from-bp", str(start), "--to-bp", str(end)]
        return cmd

    # ---------- LD decay ----------
    def ld_decay(
        self,
        out_prefix: str,
        bfile_prefix: Optional[str] = None,
        vcf_path: Optional[str] = None,
        window_kb: int = 500,
        window_snps: int = 5000,
        max_dist_kb: int = 1000,
        min_r2: float = 0.1,
        region: Optional[str] = None,
        keep_samples: Optional[str] = None,
    ) -> Dict[str, str]:
        """
        Produces:
          - <out_prefix>.ld.gz        (pairwise LD list)
          - <out_prefix>.decay.csv    (binned mean r^2 vs distance)
          - <out_prefix>.decay.png    (plot)
        """
        # 1) pairwise LD with distances (bp)
        cmd = self._plink_base(bfile_prefix, vcf_path, out_prefix, keep_samples, region)
        cmd += [
            "--r2", "yes-really", "gz",
            "--ld-window-kb", str(window_kb),
            "--ld-window", str(window_snps),
            "--ld-window-r2", str(min_r2),
            # "--ld-window-bp", str(max_dist_kb * 1000),
        ]
        _run(cmd)

        # 2) Summarize decay curve from .ld.gz
        ld_gz = out_prefix + ".ld.gz"
        if not os.path.exists(ld_gz):
            raise LDAnalysisError("plink did not create .ld.gz")

        # Use zcat or python gzip to read; prefer python for portability
        import gzip
        dists, r2s = [], []
        with gzip.open(ld_gz, "rt") as fh:
            header = fh.readline().strip().split()
            # Expect columns: CHR_A BP_A SNP_A CHR_B BP_B SNP_B R2
            # plink v1.9 has: CHR_A BP_A SNP_A CHR_B BP_B SNP_B R^2
            # normalize key:
            r2_key = "R2" if "R2" in header else "R^2"
            idx_bp_a = header.index("BP_A")
            idx_bp_b = header.index("BP_B")
            idx_r2 = header.index(r2_key)
            for ln in fh:
                parts = ln.strip().split()
                try:
                    d = abs(int(parts[idx_bp_b]) - int(parts[idx_bp_a]))
                    r = float(parts[idx_r2])
                except Exception:
                    continue
                dists.append(d)
                r2s.append(r)

        if not dists:
            raise LDAnalysisError("No LD pairs returned by plink (check filters/region)")

        # Bin by distance (kb) and compute mean r^2
        max_bp = max(dists)
        bin_size_bp = max(1, int(max(1000, (max_bp // 50))))  # ~50 bins
        bins: Dict[int, List[float]] = {}
        for d, r in zip(dists, r2s):
            b = (d // bin_size_bp) * bin_size_bp
            bins.setdefault(b, []).append(r)

        rows = []
        for b in sorted(bins.keys()):
            arr = bins[b]
            mean_r2 = sum(arr) / len(arr)
            kb = b / 1000.0
            rows.append((kb, mean_r2))

        decay_csv = out_prefix + ".decay.csv"
        with open(decay_csv, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["distance_kb", "mean_r2"])
            for kb, r in rows:
                w.writerow([f"{kb:.6f}", f"{r:.6f}"])

        # 3) Plot
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            xs = [r[0] for r in rows]
            ys = [r[1] for r in rows]
            plt.figure(figsize=(6, 4))
            plt.plot(xs, ys)
            plt.xlabel("Distance (kb)")
            plt.ylabel("Mean r²")
            plt.title("LD decay")
            plt.tight_layout()
            decay_png = out_prefix + ".decay.png"
            plt.savefig(decay_png, dpi=150)
            plt.close()
        except Exception:
            decay_png = ""

        return {"pairs_gz": ld_gz, "decay_csv": decay_csv, "decay_png": decay_png}

    # ---------- LD heatmap ----------
    def ld_heatmap(
        self,
        out_prefix: str,
        bfile_prefix: Optional[str] = None,
        vcf_path: Optional[str] = None,
        region: Optional[str] = None,
        window_snps: int = 1000,
        min_r2: float = 0.0,
        keep_samples: Optional[str] = None,
        max_snps_for_plot: int = 500,
    ) -> Dict[str, str]:
        """
        Produces a square r^2 matrix heatmap for the given region.
        If region is None, plink will use all SNPs up to window_snps.
        """
        # 1) get pairwise LD list
        cmd = self._plink_base(bfile_prefix, vcf_path, out_prefix, keep_samples, region)
        cmd += ["--r2", "yes-really", "gz", "--ld-window", str(window_snps), "--ld-window-r2", str(min_r2)]
        _run(cmd)
        ld_gz = out_prefix + ".ld.gz"

        # 2) build matrix from list
        import gzip
        snp_index: Dict[str, int] = {}
        pairs: List[Tuple[int, int, float]] = []

        with gzip.open(ld_gz, "rt") as fh:
            header = fh.readline().strip().split()
            r2_key = "R2" if "R2" in header else "R^2"
            i_snp_a = header.index("SNP_A")
            i_snp_b = header.index("SNP_B")
            i_r2 = header.index(r2_key)
            for ln in fh:
                p = ln.strip().split()
                a, b = p[i_snp_a], p[i_snp_b]
                r2 = float(p[i_r2])
                if a not in snp_index:
                    snp_index[a] = len(snp_index)
                if b not in snp_index:
                    snp_index[b] = len(snp_index)
                pairs.append((snp_index[a], snp_index[b], r2))

        n = len(snp_index)
        # downsample if too big for plotting
        if n > max_snps_for_plot:
            # keep first N by their order of appearance
            keep = set(i for i, _ in zip(range(n), range(max_snps_for_plot)))
            pairs = [(i, j, r) for (i, j, r) in pairs if i in keep and j in keep]
            n = max_snps_for_plot

        import numpy as np
        mat = np.ones((n, n), dtype=float)
        mat[:] = np.nan
        for i, j, r in pairs:
            if i < n and j < n:
                mat[i, j] = r
                mat[j, i] = r
        for i in range(n):
            mat[i, i] = 1.0

        # 3) plot
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            plt.figure(figsize=(6, 5))
            im = plt.imshow(mat, interpolation="nearest", origin="lower", vmin=0.0, vmax=1.0)
            plt.colorbar(im, fraction=0.046, pad=0.04, label="r²")
            plt.title("LD heatmap")
            plt.tight_layout()
            heat_png = out_prefix + ".heatmap.png"
            plt.savefig(heat_png, dpi=150)
            plt.close()
        except Exception:
            heat_png = ""

        return {"pairs_gz": ld_gz, "heat_png": heat_png}

    # ---------- Diversity metrics ----------
    def diversity(
        self,
        out_prefix: str,
        bfile_prefix: Optional[str] = None,
        vcf_path: Optional[str] = None,
        keep_samples: Optional[str] = None,
        region: Optional[str] = None,
    ) -> Dict[str, str]:
        """
        Uses plink to compute common per-variant and per-sample metrics:
          --freq  (MAF)
          --hardy (HWE)
          --het   (observed/expected het per individual)
        Produces:
          <out_prefix>.frq, <out_prefix>.hwe, <out_prefix>.het and a summary CSV.
        """
        base = self._plink_base(bfile_prefix, vcf_path, out_prefix, keep_samples, region)
        # freq
        _run(base + ["--freq"])
        # hardy
        _run(base + ["--hardy"])
        # sample het
        _run(base + ["--het"])

        # Summarize
        frq, hwe, het = out_prefix + ".frq", out_prefix + ".hwe", out_prefix + ".het"
        summary_csv = out_prefix + ".diversity.csv"

        # quick aggregation: mean MAF, fraction of SNPs with MAF<0.05, HWE fails, mean F
        try:
            import pandas as pd
            df_frq = pd.read_csv(frq, delim_whitespace=True) if os.path.exists(frq) else None
            df_hwe = pd.read_csv(hwe, delim_whitespace=True) if os.path.exists(hwe) else None
            df_het = pd.read_csv(het, delim_whitespace=True) if os.path.exists(het) else None

            rows = []
            if df_frq is not None and "MAF" in df_frq.columns:
                mean_maf = float(df_frq["MAF"].replace(".", "0").astype(float).mean())
                low05 = float((df_frq["MAF"].replace(".", "0").astype(float) < 0.05).mean())
                rows.append(("mean_maf", mean_maf))
                rows.append(("fraction_maf_lt_0.05", low05))
            if df_hwe is not None and "P" in df_hwe.columns:
                fails = float((df_hwe["P"].replace(".", "1").astype(float) < 1e-6).mean())
                rows.append(("hwe_fail_fraction_p_lt_1e-6", fails))
            if df_het is not None and "F" in df_het.columns:
                mean_F = float(df_het["F"].astype(float).mean())
                rows.append(("mean_inbreeding_coeff_F", mean_F))

            with open(summary_csv, "w", newline="") as f:
                w = csv.writer(f)
                w.writerow(["metric", "value"])
                for k, v in rows:
                    w.writerow([k, f"{v:.6g}"])
        except Exception:
            # fallback: just create an empty summary
            with open(summary_csv, "w", newline="") as f:
                csv.writer(f).writerow(["metric", "value"])

        return {"frq": frq, "hwe": hwe, "het": het, "summary_csv": summary_csv}
