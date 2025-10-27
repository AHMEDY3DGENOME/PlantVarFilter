# samtools_utils.py
import os
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from typing import Callable, Dict, Optional

try:
    from plantvarfilter.linux import resolve_tool
except Exception:
    resolve_tool = None

LogFn = Callable[[str], None]


class SamtoolsError(RuntimeError):
    """Exception raised for samtools-related failures."""
    pass


def _resolve_exe(name: str) -> str:
    """Resolve executable path from bundled tools first, then PATH."""
    if resolve_tool is not None:
        exe = resolve_tool(name)
        if exe:
            return exe
    path = shutil.which(name)
    if path:
        return path
    raise SamtoolsError(
        f"Executable not found: {name}. Bundle it under PlantVarFilter/linux or add to PATH."
    )


def _run(cmd: list, log: Optional[LogFn] = None, cwd: Optional[str] = None) -> str:
    """Run a command, return stdout (decoded safely)."""
    if log:
        log(" ".join(cmd))
    proc = subprocess.run(cmd, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    out = proc.stdout.decode("utf-8", errors="replace")
    if proc.returncode != 0:
        raise SamtoolsError(f"Command failed ({proc.returncode}): {' '.join(cmd)}\n{out}")
    if out.strip() and log:
        for line in out.strip().splitlines():
            log(line)
    return out


def _run_to_file(cmd: list, out_path: str, log: Optional[LogFn] = None, cwd: Optional[str] = None):
    """Run a command and write its stdout directly to a file."""
    if log:
        log(" ".join(cmd) + f" > {out_path}")
    with open(out_path, "w") as fh:
        proc = subprocess.run(cmd, cwd=cwd, stdout=fh, stderr=subprocess.PIPE, text=True)
    if proc.returncode != 0:
        raise SamtoolsError(
            f"Command failed ({proc.returncode}): {' '.join(cmd)}\nSTDERR:\n{proc.stderr}"
        )


@dataclass
class PreprocessOutputs:
    final_bam: str
    bai: Optional[str]
    stats_files: Dict[str, str]  # keys: flagstat, stats, idxstats, depth, markdup


class Samtools:
    def __init__(self, exe: str = "samtools"):
        self.exe = _resolve_exe(exe)

    # -------- Basic ops --------
    def version(self) -> str:
        try:
            proc = subprocess.run([self.exe, "--version"],
                                  stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            out = proc.stdout.decode("utf-8", errors="replace")
            return out.splitlines()[0].strip() if out else "samtools (unknown)"
        except Exception as e:
            return f"samtools (unknown): {e}"

    def sort_name(self, in_bam: str, out_bam: str, threads: int = 4, log: Optional[LogFn] = None):
        _run([self.exe, "sort", "-n", "-@", str(threads), "-o", out_bam, in_bam], log=log)

    def sort_coord(self, in_bam: str, out_bam: str, threads: int = 4,
                   tmpdir: Optional[str] = None, log: Optional[LogFn] = None):
        args = [self.exe, "sort", "-@", str(threads)]
        if tmpdir:
            args += ["-T", os.path.join(tmpdir, "samtools_sort_tmp")]
        args += ["-o", out_bam, in_bam]
        _run(args, log=log)

    def fixmate(self, in_bam_namesorted: str, out_bam: str, add_ms_tag: bool = True,
                log: Optional[LogFn] = None):
        args = [self.exe, "fixmate"]
        if add_ms_tag:
            args.append("-m")
        args += [in_bam_namesorted, out_bam]
        _run(args, log=log)

    def markdup(self, in_bam_coordsorted: str, out_bam: str, threads: int = 4,
                remove_dups: bool = False, stats_file: Optional[str] = None,
                log: Optional[LogFn] = None):
        args = [self.exe, "markdup", "-@", str(threads)]
        if remove_dups:
            args.append("-r")
        if stats_file:
            args += ["-f", stats_file]
        args += [in_bam_coordsorted, out_bam]
        _run(args, log=log)

    def index(self, bam: str, threads: int = 2, log: Optional[LogFn] = None) -> str:
        _run([self.exe, "index", "-@", str(threads), bam], log=log)
        # Typical outputs are bam.bai or basename.bai
        candidates = [bam + ".bai", os.path.splitext(bam)[0] + ".bai"]
        for p in candidates:
            if os.path.exists(p):
                return p
        return ""

    # -------- QC helpers --------
    def flagstat(self, bam: str, out_txt: str, log: Optional[LogFn] = None):
        _run_to_file([self.exe, "flagstat", bam], out_txt, log=log)

    def stats(self, bam: str, out_txt: str, log: Optional[LogFn] = None):
        _run_to_file([self.exe, "stats", bam], out_txt, log=log)

    def idxstats(self, bam: str, out_txt: str, log: Optional[LogFn] = None):
        _run_to_file([self.exe, "idxstats", bam], out_txt, log=log)

    def depth_all(self, bam: str, out_txt: str, max_depth: Optional[int] = None,
                  log: Optional[LogFn] = None):
        args = [self.exe, "depth", "-a"]
        if max_depth is not None:
            args += ["-d", str(max_depth)]
        args.append(bam)
        _run_to_file(args, out_txt, log=log)

    # -------- High-level pipeline --------
    def preprocess(self,
                   input_path: str,
                   out_prefix: Optional[str] = None,
                   threads: int = 4,
                   remove_dups: bool = False,
                   compute_stats: bool = True,
                   log: Optional[LogFn] = None,
                   keep_temps: bool = False) -> PreprocessOutputs:
        """
        Pipeline:
          0) If input is SAM/CRAM -> convert to BAM
          1) sort -n
          2) fixmate -m
          3) sort (coordinate)
          4) markdup [-r]
          5) index
          6) optional: flagstat, stats, idxstats, depth
        """
        if not os.path.exists(input_path):
            raise SamtoolsError(f"Input file not found: {input_path}")

        workdir = os.path.dirname(os.path.abspath(input_path))
        stem = os.path.splitext(os.path.basename(input_path))[0]
        base = out_prefix or f"{stem}.sm"
        base = os.path.join(workdir, base)

        tmpdir = tempfile.mkdtemp(prefix="samtools_", dir=workdir)
        try:
            # Step 0: ensure BAM
            ext = os.path.splitext(input_path)[1].lower()
            start_bam = base + ".start.bam"
            if ext == ".bam":
                if os.path.abspath(input_path) != os.path.abspath(start_bam):
                    shutil.copyfile(input_path, start_bam)
            elif ext == ".sam":
                _run([self.exe, "view", "-bS", "-o", start_bam, input_path], log=log)
            elif ext == ".cram":
                _run([self.exe, "view", "-b", "-o", start_bam, input_path], log=log)
            else:
                raise SamtoolsError(f"Unsupported input extension: {ext}")

            ns_bam = base + ".namesort.bam"
            fx_bam = base + ".fixmate.bam"
            cs_bam = base + ".coordsort.bam"
            dd_bam = base + ".dedup.bam"
            md_stats = base + ".markdup.txt"

            # 1) name sort
            self.sort_name(start_bam, ns_bam, threads=threads, log=log)
            # 2) fixmate
            self.fixmate(ns_bam, fx_bam, add_ms_tag=True, log=log)
            # 3) coordinate sort
            self.sort_coord(fx_bam, cs_bam, threads=threads, tmpdir=tmpdir, log=log)
            # 4) markdup
            self.markdup(cs_bam, dd_bam, threads=threads, remove_dups=remove_dups,
                         stats_file=md_stats, log=log)
            # 5) index
            bai = self.index(dd_bam, threads=max(1, threads // 2), log=log)

            stats_files: Dict[str, str] = {"markdup": md_stats}
            if compute_stats:
                flagstat_p = base + ".flagstat.txt"
                stats_p = base + ".stats.txt"
                idxstats_p = base + ".idxstats.txt"
                depth_p = base + ".depth.txt"
                self.flagstat(dd_bam, flagstat_p, log=log)
                self.stats(dd_bam, stats_p, log=log)
                self.idxstats(dd_bam, idxstats_p, log=log)
                self.depth_all(dd_bam, depth_p, log=log)
                stats_files.update({
                    "flagstat": flagstat_p,
                    "stats": stats_p,
                    "idxstats": idxstats_p,
                    "depth": depth_p
                })

            if not keep_temps:
                for p in (start_bam, ns_bam, fx_bam, cs_bam):
                    try:
                        if os.path.exists(p):
                            os.remove(p)
                    except Exception:
                        pass

            return PreprocessOutputs(final_bam=dd_bam, bai=bai, stats_files=stats_files)
        finally:
            try:
                os.rmdir(tmpdir)
            except Exception:
                pass
