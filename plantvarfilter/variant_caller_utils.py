# variant_caller_utils.py
import os
import shutil
import subprocess
import tempfile
from typing import List, Optional, Tuple, Callable, Union, Dict

try:
    from .bcftools_utils import BCFtools
except ImportError:
    from bcftools_utils import BCFtools

try:
    from .linux import resolve_tool
except Exception:
    try:
        from plantvarfilter.linux import resolve_tool
    except Exception:
        resolve_tool = None


LogFn = Callable[[str], None]


class VariantCallerError(RuntimeError):
    pass


def _resolve_exe(name: str) -> str:
    if resolve_tool is not None:
        p = resolve_tool(name)
        if p:
            return p
    p = shutil.which(name)
    if p:
        return p
    raise VariantCallerError(f"Executable not found: {name}. Bundle it under PlantVarFilter/linux or add to PATH.")


class VariantCaller:
    def __init__(
        self,
        bcftools_bin: Optional[str] = None,
        bgzip_bin: Optional[str] = None,
        tabix_bin: Optional[str] = None,
    ):
        self.bcftools = bcftools_bin or _resolve_exe("bcftools")
        self.bgzip = bgzip_bin or _resolve_exe("bgzip")
        self.tabix = tabix_bin or _resolve_exe("tabix")
        self.last_split: Optional[Dict[str, str]] = None

    def _run(self, cmd: List[str], log: LogFn):
        log(" ".join(cmd))
        res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        out = res.stdout.decode("utf-8", errors="replace")
        if res.returncode != 0:
            log(out)
            raise VariantCallerError(f"Command failed: {' '.join(cmd)}")
        if out.strip():
            for ln in out.strip().splitlines():
                log(ln)

    def get_last_split(self) -> Optional[Dict[str, str]]:
        return self.last_split

    def call_bcftools(
        self,
        bams: Union[List[str], str],
        ref_fasta: str,
        out_prefix: Optional[str] = None,
        regions_bed: Optional[str] = None,
        threads: int = 8,
        min_baseq: int = 13,
        min_mapq: int = 20,
        ploidy: int = 2,
        log: LogFn = print,
        split_after_calling: bool = False,
    ) -> Tuple[str, str]:
        if not os.path.exists(ref_fasta):
            raise VariantCallerError(f"Reference FASTA not found: {ref_fasta}")

        workdir = os.getcwd()
        if isinstance(bams, list) and len(bams) > 0:
            workdir = os.path.dirname(os.path.abspath(bams[0])) or workdir
        elif isinstance(bams, str) and os.path.exists(bams):
            workdir = os.path.dirname(os.path.abspath(bams)) or workdir

        if not out_prefix:
            out_prefix = os.path.join(workdir, "calls")

        vcf_gz = f"{out_prefix}.vcf.gz"
        out_dir = os.path.dirname(vcf_gz) or "."
        os.makedirs(out_dir, exist_ok=True)

        bam_list_path = None
        single_bam: Optional[str] = None
        is_temp_bamlist = False

        if isinstance(bams, list):
            if len(bams) == 0:
                raise VariantCallerError("No BAMs provided.")
            if len(bams) == 1:
                single_bam = bams[0]
            else:
                fd, tmp = tempfile.mkstemp(prefix="bamlist_", suffix=".list", dir=workdir)
                os.close(fd)
                with open(tmp, "w") as fh:
                    for b in bams:
                        fh.write(b + "\n")
                bam_list_path = tmp
                single_bam = None
                is_temp_bamlist = True
        else:
            if not os.path.exists(bams):
                raise VariantCallerError(f"BAM list not found: {bams}")
            bam_list_path = bams
            single_bam = None
            is_temp_bamlist = False

        try:
            mp_cmd = [
                self.bcftools,
                "mpileup",
                "-Ou",
                "-f",
                ref_fasta,
                "-q",
                str(min_mapq),
                "-Q",
                str(min_baseq),
                "-a",
                "FORMAT/AD,FORMAT/DP",
            ]
            if threads and threads > 1:
                mp_cmd += ["--threads", str(threads)]
            if regions_bed and os.path.exists(regions_bed):
                mp_cmd += ["-R", regions_bed]
            if bam_list_path:
                mp_cmd += ["-b", bam_list_path]
            elif single_bam:
                mp_cmd += [single_bam]

            call_cmd = [
                self.bcftools,
                "call",
                "-mv",
                "--ploidy",
                str(ploidy),
                "-Oz",
                "-o",
                vcf_gz,
            ]
            if threads and threads > 1:
                call_cmd = ["taskset", "-c", f"0-{max(0, threads - 1)}"] + call_cmd

            log(" | ".join([" ".join(mp_cmd), " ".join(call_cmd)]))
            p1 = subprocess.Popen(mp_cmd, stdout=subprocess.PIPE)
            p2 = subprocess.Popen(call_cmd, stdin=p1.stdout, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            p1.stdout.close()
            out = p2.communicate()[0]
            rc = p2.returncode
            if rc != 0:
                out_txt = (out or b"").decode("utf-8", errors="replace")
                raise VariantCallerError(f"bcftools call failed:\n{out_txt}")

            self._run([self.tabix, "-f", "-p", "vcf", vcf_gz], log)

            self.last_split = None
            if split_after_calling:
                try:
                    bcf = BCFtools()
                    paths = bcf.split_vcf(
                        input_vcf=vcf_gz,
                        out_prefix=out_prefix,
                        log=log,
                        regions_bed=regions_bed,
                    )
                    self.last_split = paths
                    if paths.get("snps"):
                        log(f"SNPs VCF: {paths['snps']}")
                    if paths.get("indels"):
                        log(f"INDELs VCF: {paths['indels']}")
                except Exception as e:
                    log(f"WARN: split_vcf failed: {e}")

            return vcf_gz, vcf_gz + ".tbi"
        finally:
            if is_temp_bamlist and bam_list_path and os.path.exists(bam_list_path):
                try:
                    os.remove(bam_list_path)
                except Exception:
                    pass
