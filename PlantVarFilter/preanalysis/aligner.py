# PlantVarFilter/preanalysis/aligner.py

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional, List
import os
import shutil
import subprocess
import time

# Prefer project helpers when available
try:
    from PlantVarFilter.variant_caller_utils import resolve_tool  # type: ignore
except Exception:
    resolve_tool = shutil.which  # fallback

try:
    # Optional: use Samtools wrapper if the project provides it
    from PlantVarFilter.samtools_utils import Samtools  # type: ignore
except Exception:
    Samtools = None  # type: ignore

# Optional: to autogenerate .mmi if only FASTA is provided
try:
    from PlantVarFilter.preanalysis.reference_manager import ReferenceManager  # type: ignore
except Exception:
    ReferenceManager = None  # type: ignore


@dataclass
class AlignmentResult:
    """Outputs of an alignment run."""
    tool: str                  # "minimap2" or "bowtie2"
    sam: Optional[str]         # may be None when streamed directly to BAM
    bam: str
    bai: str
    flagstat: str
    elapsed_sec: float
    cmdline: str               # primary aligner command used


class Aligner:
    """
    Run read alignment with minimap2 (long reads: ONT/PacBio) or bowtie2 (short reads: Illumina).
    Produces a coordinate-sorted, indexed BAM suitable for downstream analysis.
    """

    def __init__(self, logger=print, workspace: Optional[str] = None):
        self.log = logger
        self.workspace = Path(workspace or os.getcwd())
        self.workspace.mkdir(parents=True, exist_ok=True)

        # Resolve tools
        self._paths: Dict[str, Optional[str]] = {
            "samtools": self._tool_path("samtools"),
            "minimap2": self._tool_path("minimap2"),
            "bowtie2": self._tool_path("bowtie2"),
            "bowtie2-build": self._tool_path("bowtie2-build"),
        }

        # Prefer project Samtools wrapper if available
        if Samtools:
            self.st = Samtools(logger=self.log)
            # Expose binary path for piping helpers
            self._samtools_bin = getattr(self.st, "binary", self._paths["samtools"])
        else:
            self.st = None
            self._samtools_bin = self._paths["samtools"]

        if not self._samtools_bin:
            raise RuntimeError("samtools not found in PATH")

    # ------------------------------- public API ------------------------------- #

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
        """
        Long-read alignment. 'reference_mmi_or_fasta' may be a .mmi index or a FASTA path.
        """
        if not self._paths["minimap2"]:
            raise RuntimeError("minimap2 not found in PATH")

        t0 = time.time()
        out = Path(out_dir or (self.workspace / "alignment")).resolve()
        out.mkdir(parents=True, exist_ok=True)

        # Ensure .mmi exists if a FASTA was provided
        ref = Path(reference_mmi_or_fasta)
        if ref.suffix != ".mmi":
            if ReferenceManager is None:
                raise RuntimeError("FASTA provided but .mmi is missing and ReferenceManager is unavailable.")
            rm = ReferenceManager(logger=self.log, workspace=str(self.workspace))
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

        if save_sam:
            self._run(align_cmd + ["-o", str(sam_path)])
            self._sam_to_sorted_bam(str(sam_path), str(bam_path), threads=threads)
        else:
            self._pipe_to_sorted_bam(align_cmd, str(bam_path), threads=threads)

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
        """
        Short-read alignment for Illumina paired-end or single-end reads.
        """
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

        # Read group formatting
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

        if save_sam:
            self._run(args + ["-S", str(sam_path)])
            self._sam_to_sorted_bam(str(sam_path), str(bam_path), threads=threads)
        else:
            self._pipe_to_sorted_bam(args + ["-S", "-"], str(bam_path), threads=threads)

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
        """
        Convenience dispatcher based on platform:
          - "ont", "pb", "hifi", "long" -> minimap2 (preset auto-chosen)
          - otherwise -> bowtie2
        'reference' is a .mmi or FASTA for minimap2, or bowtie2 prefix for bowtie2.
        """
        p = platform.lower().strip()
        if p in {"ont", "nanopore", "pb", "pacbio", "hifi", "long"}:
            preset = "map-ont" if p in {"ont", "nanopore", "long"} else ("map-hifi" if p in {"hifi"} else "map-pb")
            return self.minimap2(
                reference, [reads1] if not reads2 else [reads1, reads2],
                preset=preset, threads=threads, read_group=read_group,
                save_sam=save_sam, out_dir=out_dir, out_prefix=out_prefix
            )
        # bowtie2 for short reads
        return self.bowtie2(
            reference, reads1, reads2,
            threads=threads, read_group=read_group, save_sam=save_sam,
            out_dir=out_dir, out_prefix=out_prefix
        )

    # ------------------------------- helpers ------------------------------- #

    def _tool_path(self, name: str) -> Optional[str]:
        p = resolve_tool(name)
        return str(p) if p else None

    def _run(self, cmd: List[str]) -> None:
        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        if proc.stdout:
            self.log(proc.stdout.rstrip())
        if proc.returncode != 0:
            raise RuntimeError(f"Command failed: {' '.join(cmd)}")

    def _pipe_to_sorted_bam(self, align_cmd: List[str], bam_out: str, threads: int = 8) -> None:
        """aligner | samtools view -bS - | samtools sort -@ t -o bam_out"""
        assert self._samtools_bin
        p1 = subprocess.Popen(align_cmd, stdout=subprocess.PIPE)
        p2 = subprocess.Popen([self._samtools_bin, "view", "-bS", "-"], stdin=p1.stdout, stdout=subprocess.PIPE)  # type: ignore
        p3 = subprocess.Popen([self._samtools_bin, "sort", "-@", str(threads), "-o", bam_out], stdin=p2.stdout)  # type: ignore
        p3.communicate()
        if p3.returncode != 0:
            raise RuntimeError("Failed to produce sorted BAM")

    def _sam_to_sorted_bam(self, sam_path: str, bam_out: str, threads: int = 8) -> None:
        """samtools view -> sort"""
        assert self._samtools_bin
        self._run([self._samtools_bin, "view", "-bS", sam_path, "-o", bam_out])  # type: ignore
        tmp_sorted = bam_out + ".tmp"
        self._run([self._samtools_bin, "sort", "-@", str(threads), "-o", tmp_sorted, bam_out])  # type: ignore
        os.replace(tmp_sorted, bam_out)

    def _markdup_inplace(self, bam_path: str) -> None:
        """samtools markdup -r"""
        assert self._samtools_bin
        tmp = bam_path + ".mkdup.tmp.bam"
        try:
            self._run([self._samtools_bin, "markdup", "-r", bam_path, tmp])  # type: ignore
            os.replace(tmp, bam_path)
        except Exception:
            # If markdup is unavailable or fails, keep original BAM
            if os.path.exists(tmp):
                try: os.remove(tmp)
                except Exception: pass

    def _index_bam(self, bam_path: str) -> None:
        assert self._samtools_bin
        self._run([self._samtools_bin, "index", bam_path])  # type: ignore

    def _flagstat(self, bam_path: str, out_txt: str) -> None:
        assert self._samtools_bin
        with open(out_txt, "w") as fw:
            proc = subprocess.run([self._samtools_bin, "flagstat", bam_path], stdout=fw, stderr=subprocess.STDOUT, text=True)  # type: ignore
            if proc.returncode != 0:
                raise RuntimeError("samtools flagstat failed")


__all__ = ["Aligner", "AlignmentResult"]
