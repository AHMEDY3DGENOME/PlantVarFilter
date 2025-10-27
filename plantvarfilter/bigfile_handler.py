# bigfile_handler.py
import os
import io
import gzip
import json
import math
import shutil
import hashlib
import tempfile
import subprocess
from dataclasses import dataclass
from typing import Callable, Iterator, List, Optional, Tuple

# -----------------------------
# Settings and small utilities
# -----------------------------

@dataclass
class LargeFileSettings:
    enabled: bool = True
    # Split by number of variant lines per chunk (header lines are copied to each chunk)
    chunk_lines: int = 500_000
    # Max number of parallel workers when mapping chunks
    max_workers: int = 2
    # Optional temp directory (falls back to system tmp if None)
    temp_dir: Optional[str] = None
    # Resume from a previous run (uses checkpoint file under temp_dir)
    resume: bool = True
    # Merge method for VCFs: "bcftools" or "cat" (bcftools is safer for VCFs)
    merge_strategy: str = "bcftools"
    # Gzip-level for chunk files created by Python writer (0-9)
    compress_level: int = 5
    # Buffer size in bytes for streaming
    io_buffer_bytes: int = 2 * 1024 * 1024

class BigFileError(RuntimeError):
    pass

class TempManager:
    def __init__(self, base: Optional[str] = None, prefix: str = "pvf_"):
        self._root = base or tempfile.gettempdir()
        os.makedirs(self._root, exist_ok=True)
        self._dir = tempfile.mkdtemp(prefix=prefix, dir=self._root)
        self._alive = True

    @property
    def path(self) -> str:
        return self._dir

    def cleanup(self):
        if self._alive:
            try:
                shutil.rmtree(self._dir, ignore_errors=True)
            finally:
                self._alive = False

    def __del__(self):
        try:
            self.cleanup()
        except Exception:
            pass

class Checkpoint:
    def __init__(self, path: str):
        self.path = path
        self.state = {"done": [], "meta": {}}
        self._load()

    def _load(self):
        try:
            if os.path.exists(self.path):
                with open(self.path, "r", encoding="utf-8") as f:
                    self.state = json.load(f)
        except Exception:
            self.state = {"done": [], "meta": {}}

    def save(self):
        try:
            with open(self.path, "w", encoding="utf-8") as f:
                json.dump(self.state, f, ensure_ascii=False, indent=2)
        except Exception:
            pass

    def mark_done(self, key: str):
        if key not in self.state["done"]:
            self.state["done"].append(key)
            self.save()

    def is_done(self, key: str) -> bool:
        return key in self.state["done"]

class ExternalRunner:
    @staticmethod
    def run(cmd: List[str], log: Callable[[str], None], check: bool = True) -> int:
        log(" ".join(cmd))
        p = subprocess.run(cmd, text=True, capture_output=True)
        if p.stdout:
            for ln in p.stdout.splitlines():
                log(ln)
        if p.stderr:
            for ln in p.stderr.splitlines():
                log(ln)
        if check and p.returncode != 0:
            raise BigFileError(f"Command failed: {' '.join(cmd)} (rc={p.returncode})")
        return p.returncode

# --------------------------------
# Core: BigFileProcessor (generic)
# --------------------------------

class BigFileProcessor:
    """
    Generic large-file helper for:
      - Streaming and chunking text files (VCF/VCF.GZ)
      - Mapping a function over chunks (in sequence or with limited parallelism)
      - Merging chunk outputs (prefer bcftools for VCF)
      - Checkpoint/Resume
    """
    def __init__(self, settings: LargeFileSettings, logger: Callable[[str], None] = print):
        self.s = settings
        self.log = logger
        self.tmp = TempManager(self.s.temp_dir, prefix="pvf_big_")

    # ---------- Basics ----------

    @staticmethod
    def _is_gzip(path: str) -> bool:
        return path.lower().endswith(".gz")

    @staticmethod
    def _hash_path(path: str) -> str:
        h = hashlib.sha1(path.encode("utf-8", errors="ignore")).hexdigest()[:10]
        return h

    def _open_text(self, path: str):
        # Returns a file-like object in text mode for .gz or plain
        if self._is_gzip(path):
            return io.TextIOWrapper(
                gzip.open(path, "rb"),
                encoding="utf-8",
                newline=""
            )
        return open(path, "r", encoding="utf-8", newline="")

    # ---------- VCF helpers ----------

    def _read_vcf_header(self, path: str) -> List[str]:
        header = []
        with self._open_text(path) as f:
            for line in f:
                if line.startswith("#"):
                    header.append(line)
                else:
                    break
        if not header or not header[-1].startswith("#CHROM"):
            # Fallback: try to re-scan to guarantee #CHROM exists
            with self._open_text(path) as f:
                for line in f:
                    if line.startswith("#CHROM"):
                        header.append(line)
                        break
        if not header:
            raise BigFileError("Invalid VCF: missing header lines.")
        return header

    def _iter_vcf_records(self, path: str) -> Iterator[str]:
        with self._open_text(path) as f:
            for line in f:
                if line and (line[0] != "#"):
                    yield line

    # ---------- Chunking ----------

    def chunk_vcf_by_lines(self, input_vcf: str, chunk_lines: Optional[int] = None) -> List[str]:
        """
        Split a VCF/VCF.GZ into N chunk files (gzipped), copying header to each chunk.
        Returns a list of chunk paths.
        """
        clines = int(chunk_lines or self.s.chunk_lines)
        header = self._read_vcf_header(input_vcf)
        base = os.path.basename(input_vcf)
        stem = base.replace(".vcf.gz", "").replace(".vcf", "")
        out_paths: List[str] = []
        idx = 0
        written = 0

        def _new_chunk_path(i: int) -> str:
            return os.path.join(self.tmp.path, f"{stem}.part{i:05d}.vcf.gz")

        writer = None
        gz = None

        try:
            for rec in self._iter_vcf_records(input_vcf):
                if writer is None:
                    part_path = _new_chunk_path(idx)
                    gz = gzip.open(part_path, "wb", compresslevel=self.s.compress_level)
                    writer = io.TextIOWrapper(gz, encoding="utf-8")
                    for h in header:
                        writer.write(h)
                    out_paths.append(part_path)
                    written = 0
                writer.write(rec)
                written += 1
                if written >= clines:
                    writer.flush()
                    writer.detach()
                    gz.close()
                    writer = None
                    gz = None
                    idx += 1
            # Close last
            if writer is not None:
                writer.flush()
                writer.detach()
                gz.close()
        finally:
            try:
                if writer is not None:
                    writer.detach()
                if gz is not None:
                    gz.close()
            except Exception:
                pass

        self.log(f"[BigFile] Chunked VCF: {len(out_paths)} part(s)")
        return out_paths

    # ---------- Merge ----------

    def merge_vcfs(self, vcfs: List[str], out_vcf_gz: str, bcftools_path: Optional[str] = "bcftools") -> str:
        if not vcfs:
            raise BigFileError("No VCFs to merge.")
        os.makedirs(os.path.dirname(out_vcf_gz) or ".", exist_ok=True)

        if self.s.merge_strategy.lower() == "bcftools":
            # Safe merge using bcftools concat -a
            cmd = [bcftools_path or "bcftools", "concat", "-a", "-Oz", "-o", out_vcf_gz] + vcfs
            ExternalRunner.run(cmd, self.log, check=True)
            ExternalRunner.run([bcftools_path or "bcftools", "index", "-f", "-p", "vcf", out_vcf_gz], self.log, check=False)
            return out_vcf_gz

        # Fallback simple concatenation (assumes compatible headers)
        # Extract header from first file + all body lines from each
        header = self._read_vcf_header(vcfs[0])
        with gzip.open(out_vcf_gz, "wb", compresslevel=self.s.compress_level) as gz:
            w = io.TextIOWrapper(gz, encoding="utf-8")
            for h in header:
                w.write(h)
            for path in vcfs:
                with self._open_text(path) as f:
                    for line in f:
                        if line and (line[0] != "#"):
                            w.write(line)
            w.flush()
        return out_vcf_gz

    # ---------- Map over chunks ----------

    def map_chunks(
        self,
        chunk_paths: List[str],
        fn: Callable[[str], str],
        checkpoint_key: str,
        checkpoint_dir: Optional[str] = None
    ) -> List[str]:
        """
        Apply fn(chunk_path) -> output_path on each chunk sequentially (worker pool can be added later).
        Respects resume when enabled.
        """
        ck_dir = checkpoint_dir or self.tmp.path
        os.makedirs(ck_dir, exist_ok=True)
        ck_path = os.path.join(ck_dir, f"{checkpoint_key}.json")
        ck = Checkpoint(ck_path)

        outputs: List[str] = []
        total = len(chunk_paths)
        for i, part in enumerate(chunk_paths, 1):
            key = os.path.basename(part)
            if self.s.resume and ck.is_done(key):
                self.log(f"[BigFile] Skip (resume): {key} ({i}/{total})")
                # Assume output present next to chunk under naming rule:
                out_guess = self._guess_output_path(part)
                if out_guess and os.path.exists(out_guess):
                    outputs.append(out_guess)
                continue
            self.log(f"[BigFile] Processing: {key} ({i}/{total})")
            out = fn(part)
            outputs.append(out)
            ck.mark_done(key)
        return outputs

    @staticmethod
    def _guess_output_path(chunk_path: str) -> Optional[str]:
        # Heuristic: replace ".partXXXXX.vcf.gz" with ".partXXXXX.proc.vcf.gz"
        if chunk_path.endswith(".vcf.gz"):
            return chunk_path[:-7] + ".proc.vcf.gz"
        return None

    # ---------- High-level pipelines ----------

    def preprocess_vcf_in_chunks(
        self,
        input_vcf: str,
        out_prefix: str,
        build_chunk_cmd: Callable[[str, str], List[str]],
        bcftools_path: Optional[str] = "bcftools",
    ) -> Tuple[str, List[str]]:
        """
        Generic chunked preprocessing:
          - split input into chunks
          - run external command per chunk (provided by build_chunk_cmd)
          - merge outputs with bcftools concat
        build_chunk_cmd(part_vcf, out_vcf_gz) -> command list
        """
        chunks = self.chunk_vcf_by_lines(input_vcf, self.s.chunk_lines)

        def _apply(part: str) -> str:
            out_part = part[:-7] + ".proc.vcf.gz" if part.endswith(".vcf.gz") else part + ".proc.vcf.gz"
            cmd = build_chunk_cmd(part, out_part)
            ExternalRunner.run(cmd, self.log, check=True)
            # index per-part to speed up concat
            ExternalRunner.run([bcftools_path or "bcftools", "index", "-f", "-p", "vcf", out_part], self.log, check=False)
            return out_part

        outputs = self.map_chunks(chunks, _apply, checkpoint_key=f"ck_{self._hash_path(input_vcf)}")
        final_out = out_prefix + ".vcf.gz" if not out_prefix.endswith(".vcf.gz") else out_prefix
        merged = self.merge_vcfs(outputs, final_out, bcftools_path=bcftools_path)
        self.log(f"[BigFile] Final merged: {merged}")
        return merged, outputs
