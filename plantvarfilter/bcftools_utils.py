# bcftools_utils.py
# Lightweight wrapper around bcftools/bgzip/tabix for VCF preprocessing

import os
import shutil
import subprocess
from typing import Optional, Tuple, List, Callable, Dict

from plantvarfilter.linux import resolve_tool

LogFn = Callable[[str], None]


class BCFtoolsError(Exception):
    pass


class BCFtools:
    def __init__(
        self,
        bcftools_bin: Optional[str] = None,
        bgzip_bin: Optional[str] = None,
        tabix_bin: Optional[str] = None,
    ):
        self.bcftools = bcftools_bin or resolve_tool("bcftools")
        self.bgzip = bgzip_bin or resolve_tool("bgzip")
        self.tabix = tabix_bin or resolve_tool("tabix")

    def _emit(self, log: LogFn, msg: str):
        try:
            log(msg)
        except TypeError:
            log(msg)

    def ensure_bins(self, log: LogFn = print):
        missing = []
        if not self.bcftools:
            missing.append("bcftools")
        if not self.bgzip:
            missing.append("bgzip")
        if not self.tabix:
            missing.append("tabix")
        if missing:
            raise BCFtoolsError(
                "Missing binaries: "
                + ", ".join(missing)
                + ". Bundle them under plantvarfilter/linux or install them in PATH."
            )
        self._emit(log, f"bcftools: {self.bcftools}")
        self._emit(log, f"bgzip   : {self.bgzip}")
        self._emit(log, f"tabix   : {self.tabix}")

    def _run(self, cmd: List[str], log: LogFn, workdir: Optional[str] = None):
        self._emit(log, " ".join(cmd))
        res = subprocess.run(
            cmd,
            cwd=workdir,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
        if res.returncode != 0:
            if res.stdout:
                self._emit(log, res.stdout)
            raise BCFtoolsError(f"Command failed ({res.returncode})")
        if res.stdout and res.stdout.strip():
            self._emit(log, res.stdout.strip())

    def _index_vcf(self, path: str, log: LogFn):
        if not os.path.exists(path):
            raise BCFtoolsError(f"VCF not found to index: {path}")
        cmd = [self.tabix, "-f", "-p", "vcf", path]
        self._run(cmd, log)

    def preprocess(
        self,
        input_vcf: str,
        out_prefix: Optional[str],
        log: LogFn = print,
        ref_fasta: Optional[str] = None,
        regions_bed: Optional[str] = None,
        split_multiallelic: bool = True,
        left_align: bool = True,
        do_sort: bool = True,
        set_id_from_fields: bool = True,
        filter_expr: Optional[str] = None,
        remove_filtered: bool = False,
        compress_output: bool = True,
        index_output: bool = True,
        keep_temps: bool = False,
        fill_tags: bool = False,
    ) -> Tuple[str, Optional[str]]:
        """
        Runs a bcftools-based preprocessing pipeline and returns:
          (final_vcf_path, stats_path)
        stats_path may be None if indexing/statistics are skipped.
        """
        self.ensure_bins(log)
        if not os.path.exists(input_vcf):
            raise BCFtoolsError(f"Input VCF not found: {input_vcf}")

        workdir = os.path.dirname(os.path.abspath(input_vcf))
        base = os.path.basename(input_vcf)
        if out_prefix is None or not out_prefix.strip():
            stem = os.path.splitext(base)[0]
            out_prefix = os.path.join(workdir, f"{stem}.prep")

        def _vcf_ext(gz: bool) -> str:
            return ".vcf.gz" if gz else ".vcf"

        temps: List[str] = []
        current = input_vcf

        if left_align or split_multiallelic:
            out1 = f"{out_prefix}.norm{_vcf_ext(True)}"
            cmd = [self.bcftools, "norm"]
            if left_align:
                if not ref_fasta:
                    self._emit(
                        log,
                        "WARN: left-align requested but no reference FASTA provided → skipping left-align",
                    )
                else:
                    cmd += ["-f", ref_fasta]
            if split_multiallelic:
                cmd += ["-m", "-both"]
            cmd += ["-Oz", "-o", out1, current]
            self._run(cmd, log)
            current = out1
            temps.append(out1)

        if do_sort:
            out2 = f"{out_prefix}.sort{_vcf_ext(True)}"
            cmd = [self.bcftools, "sort", "-Oz", "-o", out2, current]
            self._run(cmd, log)
            current = out2
            temps.append(out2)

        if set_id_from_fields:
            out3 = f"{out_prefix}.id{_vcf_ext(True)}"
            cmd = [
                self.bcftools,
                "annotate",
                "-x",
                "ID",
                "-I",
                "+%CHROM:%POS:%REF:%ALT",
                "-Oz",
                "-o",
                out3,
                current,
            ]
            self._run(cmd, log)
            current = out3
            temps.append(out3)

        if fill_tags:
            out_tags = f"{out_prefix}.tags{_vcf_ext(True)}"
            cmd = [
                self.bcftools,
                "+fill-tags",
                current,
                "-Oz",
                "-o",
                out_tags,
                "--",
                "-t",
                "AC,AN,AF,MAF,HWE",
            ]
            self._run(cmd, log)
            current = out_tags
            temps.append(out_tags)

        if regions_bed and os.path.exists(regions_bed):
            out4 = f"{out_prefix}.region{_vcf_ext(True)}"
            cmd = [self.bcftools, "view", "-R", regions_bed, "-Oz", "-o", out4, current]
            self._run(cmd, log)
            current = out4
            temps.append(out4)

        if filter_expr and filter_expr.strip():
            out5 = f"{out_prefix}.filt{_vcf_ext(True)}"
            cmd = [self.bcftools, "filter", "-i", filter_expr, "-Oz", "-o", out5, current]
            self._run(cmd, log)
            current = out5
            temps.append(out5)

        if remove_filtered:
            out6 = f"{out_prefix}.pass{_vcf_ext(True)}"
            cmd = [self.bcftools, "view", "-f", ".,PASS", "-Oz", "-o", out6, current]
            self._run(cmd, log)
            current = out6
            temps.append(out6)

        final_path = f"{out_prefix}{_vcf_ext(compress_output)}"
        if os.path.abspath(current) != os.path.abspath(final_path):
            if compress_output:
                if current.endswith(".vcf.gz"):
                    shutil.copyfile(current, final_path)
                else:
                    cmd = [self.bgzip, "-c", current]
                    with open(final_path, "wb") as fout:
                        self._emit(log, " ".join(cmd) + f" > {final_path}")
                        proc = subprocess.run(cmd, stdout=fout)
                        if proc.returncode != 0:
                            raise BCFtoolsError("bgzip failed")
            else:
                if current.endswith(".vcf"):
                    shutil.copyfile(current, final_path)
                else:
                    cmd = [self.bcftools, "view", "-Ov", "-o", final_path, current]
                    self._run(cmd, log)

        stats_path: Optional[str] = None
        if compress_output and index_output:
            self._index_vcf(final_path, log)

            stats_path = f"{out_prefix}.stats.txt"
            cmd = [self.bcftools, "stats", "-s", "-", final_path]
            self._emit(log, " ".join(cmd) + f" > {stats_path}")
            with open(stats_path, "w") as fout:
                proc = subprocess.run(cmd, stdout=fout, text=True)
                if proc.returncode != 0:
                    self._emit(log, "WARN: bcftools stats failed")

        if not keep_temps:
            for t in temps:
                try:
                    os.remove(t)
                except Exception:
                    pass
                for ext in (".tbi", ".csi"):
                    try:
                        os.remove(t + ext)
                    except Exception:
                        pass

        self._emit(log, f"Final VCF: {final_path}")
        if stats_path and os.path.exists(stats_path):
            self._emit(log, f"Stats    : {stats_path}")
        return final_path, stats_path

    def split_vcf(
        self,
        input_vcf: str,
        out_prefix: Optional[str] = None,
        log: LogFn = print,
        regions_bed: Optional[str] = None,
    ) -> Dict[str, str]:
        """
        Split a VCF into SNPs and INDEL/other variant types.
        Returns a dict with keys: {'snps': path, 'indels': path}
        """
        self.ensure_bins(log)
        if not os.path.exists(input_vcf):
            raise BCFtoolsError(f"Input VCF not found: {input_vcf}")

        workdir = os.path.dirname(os.path.abspath(input_vcf))
        base = os.path.basename(input_vcf)
        if out_prefix is None or not out_prefix.strip():
            stem = base
            if stem.endswith(".vcf.gz"):
                stem = stem[:-7]
            elif stem.endswith(".vcf"):
                stem = stem[:-4]
            out_prefix = os.path.join(workdir, stem)

        snps_out = f"{out_prefix}.snps.vcf.gz"
        indels_out = f"{out_prefix}.indels.vcf.gz"

        cmd_snps = [self.bcftools, "view", "-v", "snps"]
        cmd_indels = [self.bcftools, "view", "-v", "indels,other"]

        if regions_bed and os.path.exists(regions_bed):
            cmd_snps += ["-R", regions_bed]
            cmd_indels += ["-R", regions_bed]

        cmd_snps += ["-Oz", "-o", snps_out, input_vcf]
        cmd_indels += ["-Oz", "-o", indels_out, input_vcf]

        self._run(cmd_snps, log)
        self._index_vcf(snps_out, log)

        self._run(cmd_indels, log)
        self._index_vcf(indels_out, log)

        self._emit(log, f"SNPs VCF   : {snps_out}")
        self._emit(log, f"INDELs VCF : {indels_out}")
        return {"snps": snps_out, "indels": indels_out}
