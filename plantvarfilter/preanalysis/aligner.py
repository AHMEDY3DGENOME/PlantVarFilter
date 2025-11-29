# plantvarfilter/preanalysis/aligner.py
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional, List
import os
import shutil
import subprocess
import time

try:
    from plantvarfilter.variant_caller_utils import resolve_tool  # type: ignore
except Exception:
    resolve_tool = shutil.which  # fallback

try:
    from plantvarfilter.samtools_utils import Samtools  # type: ignore
except Exception:
    Samtools = None  # type: ignore

try:
    from plantvarfilter.preanalysis.reference_manager import ReferenceManager  # type: ignore
except Exception:
    ReferenceManager = None  # type: ignore


@dataclass
class AlignmentResult:
    tool: str
    sam: Optional[str]
    bam: str
    bai: str
    flagstat: str
    elapsed_sec: float
    cmdline: str


class Aligner:
    def __init__(self, logger=print, workspace: Optional[str] = None):
        self.log = logger
        self.workspace = Path(workspace or os.getcwd())
        self.workspace.mkdir(parents=True, exist_ok=True)

        self._paths: Dict[str, Optional[str]] = {
            "samtools": self._tool_path("samtools"),
            "minimap2": self._tool_path("minimap2"),
            "bowtie2": self._tool_path("bowtie2"),
            "bowtie2-build": self._tool_path("bowtie2-build"),
        }

        if Samtools:
            try:
                self.st = Samtools(logger=self.log)  # type: ignore
            except TypeError:
                self.st = Samtools()  # type: ignore
            self._samtools_bin = (
                getattr(self.st, "binary", None)
                or getattr(self.st, "exe", None)
                or self._paths["samtools"]
            )
        else:
            self.st = None
            self._samtools_bin = self._paths["samtools"]

        if not self._samtools_bin:
            raise RuntimeError("samtools not found in PATH")

    def minimap2(
        self,
        reference_mmi_or_fasta: str,
        reads: List[str],
        *,
        preset: str = "map-ont",
        threads: int = 8,
        read_group: Optional[Dict[str, str]] = None,
        save_sam: bool = False,
        mark_duplicates: bool = False,
        out_dir: Optional[str] = None,
        out_prefix: Optional[str] = None,
    ) -> AlignmentResult:
        if not self._paths["minimap2"]:
            raise RuntimeError("minimap2 not found in PATH")

        t0 = time.time()
        out = Path(out_dir or (self.workspace / "alignment")).resolve()
        out.mkdir(parents=True, exist_ok=True)

        ref = Path(reference_mmi_or_fasta)
        if ref.suffix != ".mmi":
            if ReferenceManager is None:
                raise RuntimeError("FASTA provided but .mmi is missing and ReferenceManager is unavailable.")
            rm = ReferenceManager(logger=self.log, workspace=str(self.workspace))  # type: ignore
            status = rm.build_indices(str(ref), out_dir=out.as_posix())
            if not status.mmi:
                raise RuntimeError("Failed to build minimap2 index (.mmi).")
            ref = Path(status.mmi)

        prefix = out_prefix or "aln.minimap2"
        sam_path = out / f"{prefix}.sam"
        bam_path = out / f"{prefix}.sorted.bam"
        bai_path = out / f"{prefix}.sorted.bam.bai"
        flagstat_path = out / f"{prefix}.flagstat.txt"

        rg_args: List[str] = []
        if read_group:
            rg_str = "@RG" + "".join([f"\\t{k}:{v}" for k, v in read_group.items()])
            rg_args = ["-R", rg_str]

        align_cmd = [
            self._paths["minimap2"], "-t", str(threads), "-ax", preset, str(ref), *reads, *rg_args  # type: ignore
        ]
        self.log(f"[ALN] $ {' '.join(align_cmd)}")

        # Always write SAM to disk, then convert to sorted BAM
        self._run(align_cmd + ["-o", str(sam_path)])
        self._sam_to_sorted_bam(str(sam_path), str(bam_path), threads=threads)

        if not save_sam:
            try:
                sam_path.unlink()
            except OSError:
                pass

        if mark_duplicates:
            self._markdup_inplace(str(bam_path))

        self._index_bam(str(bam_path))
        self._flagstat(str(bam_path), str(flagstat_path))

        return AlignmentResult(
            tool="minimap2",
            sam=str(sam_path) if save_sam else None,
            bam=str(bam_path),
            bai=str(bai_path),
            flagstat=str(flagstat_path),
            elapsed_sec=time.time() - t0,
            cmdline=" ".join(align_cmd),
        )

    def bowtie2(
        self,
        bt2_prefix: str,
        reads1: str,
        reads2: Optional[str] = None,
        *,
        threads: int = 8,
        read_group: Optional[Dict[str, str]] = None,
        very_sensitive: bool = True,
        save_sam: bool = False,
        mark_duplicates: bool = True,
        out_dir: Optional[str] = None,
        out_prefix: Optional[str] = None,
        extra_args: Optional[List[str]] = None,
    ) -> AlignmentResult:
        if not self._paths["bowtie2"]:
            raise RuntimeError("bowtie2 not found in PATH")

        t0 = time.time()
        out = Path(out_dir or (self.workspace / "alignment")).resolve()
        out.mkdir(parents=True, exist_ok=True)

        prefix = out_prefix or "aln.bowtie2"
        sam_path = out / f"{prefix}.sam"
        bam_path = out / f"{prefix}.sorted.bam"
        bai_path = out / f"{prefix}.sorted.bam.bai"
        flagstat_path = out / f"{prefix}.flagstat.txt"

        args = [self._paths["bowtie2"], "-p", str(threads), "-x", bt2_prefix]  # type: ignore
        if very_sensitive:
            args += ["--very-sensitive"]

        if read_group:
            if "ID" in read_group:
                args += ["--rg-id", str(read_group["ID"])]
            for k in ("SM", "LB", "PL", "PU"):
                if k in read_group:
                    args += ["--rg", f"{k}:{read_group[k]}"]

        if reads2:
            args += ["-1", reads1, "-2", reads2]
        else:
            args += ["-U", reads1]

        if extra_args:
            args += extra_args

        self.log(f"[ALN] $ {' '.join(args)}")

        # Always write SAM to disk, then convert to sorted BAM
        self._run(args + ["-S", str(sam_path)])
        self._sam_to_sorted_bam(str(sam_path), str(bam_path), threads=threads)

        if not save_sam:
            try:
                sam_path.unlink()
            except OSError:
                pass

        if mark_duplicates:
            self._markdup_inplace(str(bam_path))

        self._index_bam(str(bam_path))
        self._flagstat(str(bam_path), str(flagstat_path))

        return AlignmentResult(
            tool="bowtie2",
            sam=str(sam_path) if save_sam else None,
            bam=str(bam_path),
            bai=str(bai_path),
            flagstat=str(flagstat_path),
            elapsed_sec=time.time() - t0,
            cmdline=" ".join(args),
        )

    def align(
        self,
        platform: str,
        reference: str,
        reads1: str,
        reads2: Optional[str] = None,
        *,
        threads: int = 8,
        read_group: Optional[Dict[str, str]] = None,
        save_sam: bool = False,
        out_dir: Optional[str] = None,
        out_prefix: Optional[str] = None,
    ) -> AlignmentResult:
        p = platform.lower().strip()
        if p in {"ont", "nanopore", "pb", "pacbio", "hifi", "long"}:
            preset = "map-ont" if p in {"ont", "nanopore", "long"} else ("map-hifi" if p in {"hifi"} else "map-pb")
            return self.minimap2(
                reference,
                [reads1] if not reads2 else [reads1, reads2],
                preset=preset,
                threads=threads,
                read_group=read_group,
                save_sam=save_sam,
                out_dir=out_dir,
                out_prefix=out_prefix,
            )
        return self.bowtie2(
            reference,
            reads1,
            reads2,
            threads=threads,
            read_group=read_group,
            save_sam=save_sam,
            out_dir=out_dir,
            out_prefix=out_prefix,
        )

    def _tool_path(self, name: str) -> Optional[str]:
        p = resolve_tool(name)
        return str(p) if p else None

    def _run(self, cmd: List[str]) -> None:
        # Log the command being executed
        self.log(f"[ALN] RUN: {' '.join(cmd)}")

        proc = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )

        if proc.stdout:
            self.log(proc.stdout.rstrip())

        if proc.stderr:
            self.log(proc.stderr.rstrip())

        if proc.returncode != 0:
            output = (proc.stdout or "").strip()
            err = (proc.stderr or "").strip()

            msg_parts = [f"Command failed (rc={proc.returncode}): {' '.join(cmd)}"]
            if output:
                msg_parts.append("STDOUT:\n" + output[-4000:])
            if err and err not in output:
                msg_parts.append("STDERR:\n" + err[-4000:])

            raise RuntimeError("\n".join(msg_parts))

    def _pipe_to_sorted_bam(self, align_cmd: List[str], bam_out: str, threads: int = 8) -> None:
        assert self._samtools_bin

        import shlex

        align_str = " ".join(shlex.quote(x) for x in align_cmd)
        view_str = f"{shlex.quote(self._samtools_bin)} view -b -"
        sort_str = f"{shlex.quote(self._samtools_bin)} sort -@ {int(threads)} -o {shlex.quote(bam_out)} -"

        pipeline = f"{align_str} | {view_str} | {sort_str}"

        self.log(f"[ALN] pipe (shell): {pipeline}")

        proc = subprocess.run(
            pipeline,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )

        if proc.stdout:
            self.log("[ALN] pipe stdout (truncated):\n" + proc.stdout[:1000])

        if proc.stderr:
            self.log("[ALN] pipe stderr:\n" + proc.stderr)

        if proc.returncode != 0:
            raise RuntimeError(
                f"Failed to produce sorted BAM (pipeline rc={proc.returncode})"
            )

    def _sam_to_sorted_bam(self, sam_path: str, bam_out: str, threads: int = 8) -> None:
        assert self._samtools_bin
        self._run([self._samtools_bin, "view", "-b", sam_path, "-o", bam_out])  # type: ignore
        tmp_sorted = bam_out + ".tmp"
        self._run([
            self._samtools_bin,
            "sort",
            "-@",
            str(threads),
            "-m",
            "1G",
            "-o",
            tmp_sorted,
            bam_out,
        ])  # type: ignore
        os.replace(tmp_sorted, bam_out)

    def _markdup_inplace(self, bam_path: str) -> None:
        assert self._samtools_bin
        tmp = bam_path + ".mkdup.tmp.bam"
        try:
            self._run([self._samtools_bin, "markdup", "-r", bam_path, tmp])  # type: ignore
            os.replace(tmp, bam_path)
        except Exception:
            if os.path.exists(tmp):
                try:
                    os.remove(tmp)
                except Exception:
                    pass

    def _index_bam(self, bam_path: str) -> None:
        assert self._samtools_bin
        self._run([self._samtools_bin, "index", bam_path])  # type: ignore

    def _flagstat(self, bam_path: str, out_txt: str) -> None:
        assert self._samtools_bin
        with open(out_txt, "w") as fw:
            proc = subprocess.run(
                [self._samtools_bin, "flagstat", bam_path],
                stdout=fw,
                stderr=subprocess.STDOUT,
                text=True,
            )  # type: ignore
            if proc.returncode != 0:
                raise RuntimeError("samtools flagstat failed")


__all__ = ["Aligner", "AlignmentResult"]
