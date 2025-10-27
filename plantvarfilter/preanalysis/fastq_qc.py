# plantvarfilter/preanalysis/fastq_qc.py

from __future__ import annotations

from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Tuple
import gzip
import io
import os
import shutil
import subprocess
import statistics
import json

import matplotlib.pyplot as plt


@dataclass
class FastqQCReport:
    """Lightweight summary for FASTQ quality control."""
    platform: str                     # "illumina", "ont", "pb"
    n_reads: int
    mean_length: float
    median_length: float
    gc_percent: float                 # overall (per-base)
    n_percent: float                  # overall (per-base)
    mean_phred: Optional[float]       # None when qualities are missing
    per_cycle_q_mean_png: Optional[str]
    length_hist_png: Optional[str]
    gc_hist_png: Optional[str]
    fastqc_html: Optional[str]
    report_txt: str
    verdict: str                      # "Pass" | "Warning" | "Fail"
    details_json: str                 # machine-readable dump of metrics

    def to_dict(self) -> Dict:
        d = asdict(self)
        return d


def _open_text_fastq(path: str) -> io.TextIOBase:
    p = Path(path)
    if p.suffix == ".gz":
        return io.TextIOWrapper(gzip.open(p, "rb"))
    return open(p, "r")


def _iter_fastq_reads(path: str) -> Iterator[Tuple[str, str]]:
    """Yield (seq, qual) tuples; qual may be '' for missing qualities."""
    with _open_text_fastq(path) as fh:
        while True:
            h = fh.readline()
            if not h:
                break
            seq = fh.readline().strip()
            plus = fh.readline()
            qual = fh.readline().strip()
            # If the record is malformed, stop gracefully
            if not plus:
                break
            # Some long-read files may have empty qualities
            if len(qual) != len(seq):
                qual = ""
            yield seq, qual


def _consume_fastq(paths: List[str],
                   sample_max: int,
                   compute_per_cycle_q: bool = True
                   ) -> Tuple[int, List[int], int, int, int, List[float], Optional[List[float]]]:
    """
    Stream FASTQ and collect metrics.

    Returns:
        n_reads
        lengths (list)
        gc_bases
        n_bases
        total_bases
        q_values_sample (subset of per-read mean phred for histogram)
        per_cycle_q_mean (list) or None
    """
    n_reads = 0
    lengths: List[int] = []
    gc_bases = 0
    n_bases = 0
    total_bases = 0
    q_values_sample: List[float] = []

    # per-cycle accumulators
    cycle_q_sum: List[int] = []
    cycle_q_cnt: List[int] = []

    for p in paths:
        for seq, qual in _iter_fastq_reads(p):
            n_reads += 1
            L = len(seq)
            lengths.append(L)
            total_bases += L
            gc_bases += sum(1 for c in seq if c in "GCgc")
            n_bases += seq.count("N") + seq.count("n")

            if qual:
                # mean phred for this read (phred+33)
                q_sum = sum((ord(ch) - 33) for ch in qual)
                q_values_sample.append(q_sum / L)
                if compute_per_cycle_q:
                    # grow accumulators if needed
                    if len(cycle_q_sum) < L:
                        extra = L - len(cycle_q_sum)
                        cycle_q_sum.extend([0] * extra)
                        cycle_q_cnt.extend([0] * extra)
                    for i, ch in enumerate(qual):
                        cycle_q_sum[i] += (ord(ch) - 33)
                        cycle_q_cnt[i] += 1

            if n_reads >= sample_max:
                break
        if n_reads >= sample_max:
            break

    per_cycle_q_mean: Optional[List[float]] = None
    if cycle_q_cnt:
        per_cycle_q_mean = [
            (cycle_q_sum[i] / cycle_q_cnt[i]) if cycle_q_cnt[i] else 0.0
            for i in range(len(cycle_q_cnt))
        ]

    return (
        n_reads,
        lengths,
        gc_bases,
        n_bases,
        total_bases,
        q_values_sample,
        per_cycle_q_mean,
    )


def _save_histogram(data: List[float], out_png: Path, title: str, xlabel: str) -> None:
    if not data:
        return
    plt.figure()
    plt.hist(data, bins=60)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel("Count")
    plt.tight_layout()
    plt.savefig(out_png)
    plt.close()


def _save_line(values: List[float], out_png: Path, title: str, ylabel: str, xlabel: str = "Cycle") -> None:
    if not values:
        return
    plt.figure()
    plt.plot(range(1, len(values) + 1), values)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.tight_layout()
    plt.savefig(out_png)
    plt.close()


def _run_fastqc_if_available(paths: List[str], out_dir: Path, logger=print) -> Optional[str]:
    fastqc = shutil.which("fastqc")
    if not fastqc:
        return None
    try:
        cmd = [fastqc, "-o", str(out_dir), "-f", "fastq", "-q", *paths]
        logger(f"[FQ-QC] $ {' '.join(cmd)}")
        subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        # Return the first HTML we find
        htmls = sorted(out_dir.glob("*.html"))
        return str(htmls[0]) if htmls else None
    except Exception as e:
        logger(f"[FQ-QC] fastqc failed: {e}")
        return None


def _decide_verdict(platform: str,
                    n_reads: int,
                    mean_len: float,
                    gc_p: float,
                    n_p: float,
                    mean_phred: Optional[float]) -> str:
    """
    Simple rule-based verdict:
      - Illumina: mean_phred >= 20 and N% < 1 and 50bp <= mean_len <= 300 -> Pass
      - ONT/PB: N% < 5 and mean_len >= 1000 -> Pass
      - Otherwise Warning; Fail if reads too few (<1000) or clear red flags
    """
    platform = platform.lower()
    if platform in ("illumina", "short", "sr"):
        if mean_phred is not None and mean_phred >= 20.0 and n_p < 1.0 and 50 <= mean_len <= 300:
            return "Pass"
        if n_reads < 1000 or (mean_phred is not None and mean_phred < 15) or n_p > 5.0:
            return "Fail"
        return "Warning"
    else:  # ont / pb / long
        if n_p < 5.0 and mean_len >= 1000:
            return "Pass"
        if n_reads < 500 or n_p > 20.0:
            return "Fail"
        return "Warning"


def run_fastq_qc(
    reads1: str,
    reads2: Optional[str] = None,
    *,
    platform: str = "illumina",
    out_dir: Optional[str] = None,
    sample_max: int = 1_000_000,
    use_fastqc_if_available: bool = True,
    logger=print,
) -> FastqQCReport:
    """
    Perform lightweight FASTQ QC and optionally call FastQC if available.
    Returns a FastqQCReport with paths to saved artifacts.
    """
    paths = [reads1] + ([reads2] if reads2 else [])
    if not all(Path(p).exists() for p in paths):
        missing = [p for p in paths if not Path(p).exists()]
        raise FileNotFoundError(f"Missing FASTQ files: {missing}")

    out = Path(out_dir or Path(reads1).parent / "qc").resolve()
    out.mkdir(parents=True, exist_ok=True)

    n_reads, lengths, gc_bases, n_bases, total_bases, q_values_sample, per_cycle_q_mean = _consume_fastq(
        paths, sample_max=sample_max, compute_per_cycle_q=True
    )

    mean_len = statistics.mean(lengths) if lengths else 0.0
    median_len = statistics.median(lengths) if lengths else 0.0
    gc_percent = (100.0 * gc_bases / total_bases) if total_bases else 0.0
    n_percent = (100.0 * n_bases / total_bases) if total_bases else 0.0
    mean_phred = (statistics.mean(q_values_sample) if q_values_sample else None)

    # Plots
    len_png: Optional[str] = None
    gc_hist_png: Optional[str] = None
    per_cycle_png: Optional[str] = None

    if lengths:
        len_png_path = out / "length_hist.png"
        _save_histogram([float(x) for x in lengths], len_png_path, "Read length distribution", "Length (bp)")
        len_png = str(len_png_path)

        # per-read GC% distribution
        gc_per_read: List[float] = []
        # limit for speed
        cap = min(len(lengths), 50000)
        i = 0
        for p in paths:
            for seq, _ in _iter_fastq_reads(p):
                if i >= cap:
                    break
                if seq:
                    gc_per_read.append(100.0 * (sum(ch in "GCgc" for ch in seq) / len(seq)))
                i += 1
            if i >= cap:
                break
        gc_png_path = out / "gc_percent_hist.png"
        _save_histogram(gc_per_read, gc_png_path, "Per-read GC% distribution", "GC%")
        gc_hist_png = str(gc_png_path)

    if per_cycle_q_mean:
        per_cycle_png_path = out / "per_cycle_quality_mean.png"
        _save_line(per_cycle_q_mean, per_cycle_png_path, "Per-cycle mean PHRED", "Mean PHRED")
        per_cycle_png = str(per_cycle_png_path)

    # Optional FastQC
    fastqc_html: Optional[str] = None
    if use_fastqc_if_available:
        fastqc_html = _run_fastqc_if_available(paths, out, logger=logger)

    verdict = _decide_verdict(platform, n_reads, mean_len, gc_percent, n_percent, mean_phred)

    # Text report
    txt = out / "fastq_qc_report.txt"
    with open(txt, "w") as fw:
        fw.write(f"Platform: {platform}\n")
        fw.write(f"Reads (sampled): {n_reads}\n")
        fw.write(f"Mean length: {mean_len:.2f}\n")
        fw.write(f"Median length: {median_len:.2f}\n")
        fw.write(f"GC% (overall): {gc_percent:.2f}\n")
        fw.write(f"N% (overall): {n_percent:.3f}\n")
        if mean_phred is not None:
            fw.write(f"Mean PHRED (approx.): {mean_phred:.2f}\n")
        fw.write(f"Verdict: {verdict}\n")
        if fastqc_html:
            fw.write(f"FastQC: {fastqc_html}\n")

    # Machine-readable details
    details = {
        "platform": platform,
        "n_reads": n_reads,
        "mean_length": mean_len,
        "median_length": median_len,
        "gc_percent_overall": gc_percent,
        "n_percent_overall": n_percent,
        "mean_phred": mean_phred,
        "per_cycle_quality_mean": per_cycle_q_mean or [],
        "plots": {
            "length_hist": len_png,
            "gc_hist": gc_hist_png,
            "per_cycle_quality_mean": per_cycle_png,
        },
        "fastqc_html": fastqc_html,
    }
    details_json_path = out / "fastq_qc_details.json"
    with open(details_json_path, "w") as jf:
        json.dump(details, jf, indent=2)

    return FastqQCReport(
        platform=platform,
        n_reads=n_reads,
        mean_length=mean_len,
        median_length=median_len,
        gc_percent=gc_percent,
        n_percent=n_percent,
        mean_phred=mean_phred,
        per_cycle_q_mean_png=per_cycle_png,
        length_hist_png=len_png,
        gc_hist_png=gc_hist_png,
        fastqc_html=fastqc_html,
        report_txt=str(txt),
        verdict=verdict,
        details_json=str(details_json_path),
    )


__all__ = ["FastqQCReport", "run_fastq_qc"]
