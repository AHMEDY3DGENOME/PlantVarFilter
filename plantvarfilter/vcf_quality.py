# vcf_quality.py
# Lightweight VCF Quality Checker (+ distributions for interactive plots)

from __future__ import annotations
import os, gzip, math, statistics as stats
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional

DNA = set("ACGT")

@dataclass
class VCFQualityReport:
    score: float
    verdict: str                     # "Pass" / "Caution" / "Fail"
    metrics: Dict[str, float]        # key -> value
    recommendations: List[str]
    hard_fail_reasons: List[str]
    dists: Optional[Dict[str, List[float]]] = None   # raw distributions (for plots)
    data_type: Optional[str] = None                  # "SNPs" or "SVs"

class VCFQualityChecker:
    """
    Header sanity + streaming sampling of sites to compute core metrics:
      - counts: total_sites, snps, indels, multiallelic
      - Ti/Tv for SNPs
      - site/sample missingness (approx)
      - DP / GQ summaries (median, p10, p90) when available
      - Heterozygote allele balance (from AD or AB when available)
      - Data type detection: "SNPs" vs "SVs" (SNP-only vs otherwise)
    Returns a VCF-QAScore (0..100) + verdict + recommendations + raw dists for plotting.
    """
    def __init__(self, max_sites_scan: int = 200_000, min_sites_required: int = 5_000):
        self.max_sites_scan = max_sites_scan
        self.min_sites_required = min_sites_required

    # ---------- Public API ----------
    def evaluate(self, vcf_path: str, log_fn=print) -> VCFQualityReport:
        reasons = self._quick_header_checks(vcf_path, log_fn)
        if reasons:
            return VCFQualityReport(0.0, "Fail", {}, [], reasons, dists=None, data_type=None)

        feats, dists, data_type = self._scan_features(vcf_path, log_fn)

        sites = int(feats.get("sites_scanned", 0))
        small_penalty = 0.0
        recs: List[str] = []

        if sites < 150:
            return VCFQualityReport(
                0.0, "Fail", feats,
                ["VCF is extremely small (<150 sites). Add more variants then re-check."],
                [], dists, data_type=data_type
            )
        if sites < self.min_sites_required:
            small_penalty = min(40.0, (self.min_sites_required - sites) / self.min_sites_required * 30.0)
            recs.append(f"Small VCF (n={sites}) — score penalized; treat results with caution.")

        base_score, _base_verdict, recs_base = self._score_and_recommend(feats)
        score = max(0.0, base_score - small_penalty)
        verdict = "Pass" if score >= 80 else ("Caution" if score >= 60 else "Fail")
        return VCFQualityReport(score, verdict, feats, recs + recs_base, [], dists, data_type=data_type)

    def to_text(self, report: VCFQualityReport, vcf_path: Optional[str] = None) -> str:
        """Render a human-readable QC report string."""
        lines: List[str] = []
        if vcf_path:
            lines.append(f"VCF: {vcf_path}")
        lines.append(f"VCF-QAScore: {report.score:.1f}  |  Verdict: {report.verdict}")
        if "samples" in report.metrics:
            lines.append(f"Samples: {int(report.metrics['samples'])}")
        if report.data_type:
            lines.append(f"Data type: {report.data_type}")
        lines.append("")
        lines.append("Recommendations:")
        if report.recommendations:
            for r in report.recommendations:
                lines.append(f"- {r}")
        else:
            lines.append("- None")
        lines.append("")
        lines.append("Metrics:")
        headline = ["samples", "sites_scanned", "snps", "indels", "multiallelic",
                    "multiallelic_ratio", "snp_indel_ratio", "titv",
                    "site_missing_mean", "site_missing_p90",
                    "sample_missing_mean", "sample_missing_p90",
                    "dp_median", "dp_p10", "dp_p90",
                    "gq_median", "gq_p10", "gq_p90",
                    "ab_dev_mean", "ab_dev_p90"]
        shown = set()
        for k in headline:
            if k in report.metrics:
                lines.append(f"{k:>22}: {report.metrics[k]}")
                shown.add(k)
        for k, v in sorted(report.metrics.items()):
            if k not in shown:
                lines.append(f"{k:>22}: {v}")
        if report.hard_fail_reasons:
            lines.append("")
            lines.append("Hard Fail Reasons:")
            for r in report.hard_fail_reasons:
                lines.append(f"- {r}")
        return "\n".join(lines)

    # ---------- Internals ----------
    def _is_gz(self, path: str) -> bool:
        return path.endswith(".gz")

    def _opener(self, path: str):
        return gzip.open(path, "rt") if self._is_gz(path) else open(path, "r")

    def _quick_header_checks(self, vcf_path: str, log_fn) -> List[str]:
        reasons = []
        if not os.path.exists(vcf_path):
            return ["File not found"]

        samples_found = False
        has_format = False
        contigs = []
        try:
            with self._opener(vcf_path) as fh:
                for line in fh:
                    if not line.startswith("#"):
                        break
                    if line.startswith("##contig="):
                        contigs.append(line.strip())
                    if line.startswith("##FORMAT="):
                        has_format = True
                    if line.startswith("#CHROM"):
                        parts = line.strip().split("\t")
                        if len(parts) < 8:
                            reasons.append("Header columns incomplete (#CHROM line).")
                        if len(parts) >= 10:
                            samples_found = True
        except Exception as e:
            reasons.append(f"Error reading header: {e}")

        if not contigs:
            reasons.append("No contigs declared in header.")
        if not has_format:
            reasons.append("No FORMAT fields declared.")
        if not samples_found:
            reasons.append("No samples detected (need at least one).")

        return reasons

    def _scan_features(self, vcf_path: str, log_fn) -> Tuple[Dict[str, float], Dict[str, List[float]], str]:
        total = snps = indels = multiallelic = 0
        ti = tv = 0
        site_missing_rates: List[float] = []
        sample_missing_counts: Optional[List[int]] = None
        dp_vals: List[float] = []
        gq_vals: List[float] = []
        ab_deviations: List[float] = []

        # Counters for data-type
        non_snv_sites = 0
        snv_sites = 0

        samples_n = 0
        format_keys: List[str] = []

        with self._opener(vcf_path) as fh:
            for line in fh:
                if line.startswith("##"):
                    continue
                if line.startswith("#CHROM"):
                    cols = line.strip().split("\t")
                    samples_n = max(0, len(cols) - 9)  # 9 = up to FORMAT
                    if samples_n > 0:
                        sample_missing_counts = [0] * samples_n
                    continue

                if total >= self.max_sites_scan:
                    break

                parts = line.rstrip("\n").split("\t")
                if len(parts) < 8:
                    continue

                _chrom, _pos, _id, REF, ALT, _QUAL, _FILTER, _INFO = parts[:8]
                sample_fields = parts[9:] if len(parts) > 9 else []
                FORMAT = parts[8] if len(parts) > 8 else ""

                alts = ALT.split(",")
                if len(alts) > 1:
                    multiallelic += 1

                ref_len = len(REF)
                alt_is_snv = all(len(a) == 1 and a in DNA for a in alts if a != ".")
                ref_is_snv = (ref_len == 1 and REF in DNA)
                is_snv_site = (ref_is_snv and alt_is_snv)

                if is_snv_site:
                    snv_sites += 1
                else:
                    non_snv_sites += 1

                if is_snv_site:
                    snps += 1
                    if len(alts) == 1 and alts[0] in DNA:
                        if self._is_transition(REF, alts[0]):
                            ti += 1
                        else:
                            tv += 1
                else:
                    indels += 1

                site_missing = 0.0
                if samples_n and FORMAT:
                    if not format_keys:
                        format_keys = FORMAT.split(":")
                    key_idx = {k: i for i, k in enumerate(format_keys)}
                    has_GT = "GT" in key_idx
                    idx_AD = key_idx.get("AD")
                    idx_AB = key_idx.get("AB")
                    idx_DP = key_idx.get("DP")
                    idx_GQ = key_idx.get("GQ")

                    for si, cell in enumerate(sample_fields):
                        if not cell or cell == ".":
                            if sample_missing_counts is not None:
                                sample_missing_counts[si] += 1
                            site_missing += 1
                            continue

                        toks = cell.split(":")
                        gt = toks[key_idx["GT"]] if has_GT and key_idx["GT"] < len(toks) else None

                        if gt is None or gt == "." or gt == "./.":
                            if sample_missing_counts is not None:
                                sample_missing_counts[si] += 1
                            site_missing += 1
                        else:
                            if ("0/1" in gt) or ("1/0" in gt) or ("0|1" in gt) or ("1|0" in gt):
                                if idx_AB is not None and idx_AB < len(toks):
                                    try:
                                        ab = float(toks[idx_AB])
                                        if 0 <= ab <= 1:
                                            ab_deviations.append(abs(ab - 0.5))
                                    except Exception:
                                        pass
                                elif idx_AD is not None and idx_AD < len(toks):
                                    try:
                                        ad = toks[idx_AD]
                                        rs = ad.split(",")
                                        if len(rs) >= 2:
                                            rc = float(rs[0]) if rs[0].isdigit() else None
                                            ac = float(rs[1]) if rs[1].isdigit() else None
                                            if rc is not None and ac is not None and (rc + ac) > 0:
                                                ab = ac / (rc + ac)
                                                ab_deviations.append(abs(ab - 0.5))
                                    except Exception:
                                        pass

                            if idx_DP is not None and idx_DP < len(toks):
                                try:
                                    dp_val = float(toks[idx_DP])
                                    if dp_val >= 0:
                                        dp_vals.append(dp_val)
                                except Exception:
                                    pass
                            if idx_GQ is not None and idx_GQ < len(toks):
                                try:
                                    gq_val = float(toks[idx_GQ])
                                    if gq_val >= 0:
                                        gq_vals.append(gq_val)
                                except Exception:
                                    pass

                if samples_n > 0:
                    site_missing_rates.append(site_missing / samples_n)

                total += 1

        metrics: Dict[str, float] = {}
        metrics["samples"] = float(samples_n)
        metrics["sites_scanned"] = float(total)
        metrics["snps"] = float(snv_sites)  # single-nucleotide sites
        metrics["indels"] = float(indels)   # includes non-SNVs seen in counting branch
        metrics["multiallelic"] = float(multiallelic)
        metrics["multiallelic_ratio"] = (multiallelic / total) if total else 0.0
        metrics["snp_indel_ratio"] = (snv_sites / indels) if indels else float("inf")
        metrics["titv"] = (ti / tv) if tv else float("inf")

        if site_missing_rates:
            metrics["site_missing_mean"] = float(stats.mean(site_missing_rates))
            metrics["site_missing_p90"] = float(self._percentile(site_missing_rates, 90))
        else:
            metrics["site_missing_mean"] = 0.0
            metrics["site_missing_p90"] = 0.0

        if sample_missing_counts is not None and total > 0:
            per_sample_missing = [c / total for c in sample_missing_counts]
            metrics["sample_missing_mean"] = float(stats.mean(per_sample_missing))
            metrics["sample_missing_p90"] = float(self._percentile(per_sample_missing, 90))
        else:
            metrics["sample_missing_mean"] = 0.0
            metrics["sample_missing_p90"] = 0.0

        for name, arr in (("dp", dp_vals), ("gq", gq_vals)):
            if arr:
                metrics[f"{name}_median"] = float(stats.median(arr))
                metrics[f"{name}_p10"] = float(self._percentile(arr, 10))
                metrics[f"{name}_p90"] = float(self._percentile(arr, 90))
            else:
                metrics[f"{name}_median"] = 0.0
                metrics[f"{name}_p10"] = 0.0
                metrics[f"{name}_p90"] = 0.0

        metrics["ab_dev_mean"] = float(stats.mean(ab_deviations)) if ab_deviations else 0.0
        metrics["ab_dev_p90"] = float(self._percentile(ab_deviations, 90)) if ab_deviations else 0.0

        dists = {
            "dp": dp_vals,
            "gq": gq_vals,
            "ab_dev": ab_deviations,
            "site_missing": site_missing_rates,
        }

        # Data type decision (two-label policy):
        # - SNPs: all sites are single-nucleotide (no non-SNVs encountered)
        # - SVs : any site is non-SNV (indel/MNV/symbolic etc.)
        data_type = "SNPs" if non_snv_sites == 0 else "SVs"

        return metrics, dists, data_type

    def _score_and_recommend(self, m: Dict[str, float]) -> Tuple[float, str, List[str]]:
        penalties = []
        recs: List[str] = []

        titv = m["titv"]
        if math.isinf(titv) or titv == 0:
            p = 1.0
            recs.append("Ti/Tv not computable → check SNP calling / ALT parsing.")
        else:
            p = min(abs(titv - 2.0) / 1.0, 1.0)
        penalties.append(("TiTv", 10, p))

        p_site = min(m["site_missing_mean"] * 5.0, 1.0)
        p_samp = min(m["sample_missing_mean"] * 5.0, 1.0)
        penalties.append(("SiteMissing", 10, p_site))
        penalties.append(("SampleMissing", 10, p_samp))
        if m["site_missing_mean"] > 0.05:
            recs.append("High site missingness → consider --geno 0.05 or higher.")
        if m["sample_missing_p90"] > 0.1:
            recs.append("Some samples have high missingness → consider --mind 0.1 and remove worst samples.")

        snps, indels = m["snps"], m["indels"]
        snp_indel_ratio = m["snp_indel_ratio"]
        if indels == 0:
            p = 0.3
        else:
            p = min(abs(snp_indel_ratio - 3.0) / 2.0, 1.0)
        penalties.append(("SNP_INDEL", 6, p))

        mar = m["multiallelic_ratio"]
        p = 0.0 if mar <= 0.1 else min((mar - 0.1) / 0.2, 1.0)
        penalties.append(("Multiallelic", 6, p))

        dp_med = m["dp_median"]
        gq_med = m["gq_median"]
        p_dp = 0.0 if 10 <= dp_med <= 200 else 0.6
        p_gq = 0.0 if gq_med >= 20 else 0.8
        penalties.append(("DepthMedian", 8, p_dp))
        penalties.append(("GQMedian", 8, p_gq))
        if dp_med and dp_med < 10:
            recs.append("Low depth median → consider stricter QUAL/DP filters or re-calling.")
        if gq_med and gq_med < 20:
            recs.append("Low GQ median → consider filtering low GQ genotypes.")

        ab_dev = m["ab_dev_mean"]
        p = min(ab_dev / 0.15, 1.0)
        penalties.append(("ABDeviation", 10, p))
        if ab_dev > 0.2:
            recs.append("Heterozygote allele balance is skewed → could indicate contamination or mapping bias.")

        max_score = 100.0
        total_penalty = sum(w * p for _, w, p in penalties)
        score = max(0.0, max_score - total_penalty)

        verdict = "Pass" if score >= 80 else ("Caution" if score >= 60 else "Fail")

        if snps > 0 and verdict != "Pass":
            recs.append("Try MAF filter around 0.03–0.05 and re-evaluate.")
        return score, verdict, recs

    # ---------- Helpers ----------
    def _is_transition(self, a: str, b: str) -> bool:
        pair = {a, b}
        return pair == {"A", "G"} or pair == {"C", "T"}

    def _percentile(self, arr: List[float], p: float) -> float:
        if not arr:
            return 0.0
        arr_sorted = sorted(arr)
        k = (len(arr_sorted) - 1) * (p / 100.0)
        f = math.floor(k); c = math.ceil(k)
        if f == c:
            return arr_sorted[int(k)]
        d0 = arr_sorted[f] * (c - k)
        d1 = arr_sorted[c] * (k - f)
        return d0 + d1
