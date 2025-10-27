# plantvarfilter/preanalysis/reference_manager.py

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import os
import shutil
import subprocess
from typing import Dict, Optional


# Prefer the project's resolver if available
try:
    from plantvarfilter.variant_caller_utils import resolve_tool  # type: ignore
except Exception:
    resolve_tool = shutil.which  # fallback


@dataclass
class ReferenceIndexStatus:
    """Summary of reference indexing artifacts and tool availability."""
    fasta: str
    reference_dir: str
    faidx: Optional[str]
    dict: Optional[str]
    mmi: Optional[str]
    bt2_prefix: Optional[str]
    tools: Dict[str, Dict[str, Optional[str]]]
    ok: bool


class ReferenceManager:
    """
    Build and validate reference indexes used by downstream steps.

    Artifacts (co-located in a dedicated 'reference' directory):
      - FASTA copy/link
      - <fasta>.fai                  (samtools faidx)
      - <fasta_stem>.dict            (samtools dict or Picard)
      - <fasta_stem>.mmi             (minimap2 index)
      - <fasta_stem>_bt2.*.bt2[l]?   (bowtie2-build index)
    """

    def __init__(self, logger=print, workspace: Optional[str] = None):
        self.log = logger
        self.workspace = Path(workspace) if workspace else None

    # ---------------------------- public API ---------------------------- #

    def build_indices(
        self,
        fasta: str,
        out_dir: Optional[str] = None,
        build_mmi: bool = True,
        build_bt2: bool = True,
        build_dict: bool = True,
    ) -> ReferenceIndexStatus:
        """
        Build any missing indexes for the given FASTA.
        Returns a ReferenceIndexStatus describing what exists.
        """
        fa_src = Path(fasta).expanduser().resolve()
        if not fa_src.exists():
            raise FileNotFoundError(f"FASTA not found: {fa_src}")

        ref_dir = Path(out_dir or self._default_reference_dir_for(fa_src)).resolve()
        ref_dir.mkdir(parents=True, exist_ok=True)

        # Place a copy/hardlink of the FASTA inside reference dir so sidecar files sit together.
        fa_dst = ref_dir / fa_src.name
        if not fa_dst.exists():
            try:
                os.link(str(fa_src), str(fa_dst))  # hard link when possible
                self.log(f"[REF] Linked FASTA → {fa_dst}")
            except Exception:
                shutil.copy2(str(fa_src), str(fa_dst))
                self.log(f"[REF] Copied FASTA → {fa_dst}")

        # Build faidx
        samtools = self._tool_path("samtools")
        if samtools:
            fai = fa_dst.with_suffix(fa_dst.suffix + ".fai")
            if not fai.exists():
                self._run([samtools, "faidx", str(fa_dst)], tag="samtools")
        else:
            self.log("[REF] samtools not found in PATH (faidx will be missing)")

        # Build dict (samtools dict preferred; Picard fallback if available)
        if build_dict:
            dict_path = ref_dir / (fa_dst.stem + ".dict")
            if not dict_path.exists():
                if samtools and self._supports_samtools_dict(samtools):
                    try:
                        self._run([samtools, "dict", str(fa_dst), "-o", str(dict_path)], tag="samtools")
                    except Exception as e:
                        self.log(f"[REF] samtools dict failed ({e}). Will try Picard if present.")
                        self._try_picard_dict(fa_dst, dict_path)
                else:
                    self._try_picard_dict(fa_dst, dict_path)

        # Build minimap2 .mmi
        if build_mmi:
            minimap2 = self._tool_path("minimap2")
            mmi = ref_dir / (fa_dst.stem + ".mmi")
            if minimap2:
                if not mmi.exists():
                    self._run([minimap2, "-d", str(mmi), str(fa_dst)], tag="minimap2")
            else:
                self.log("[REF] minimap2 not found in PATH (mmi will be missing)")

        # Build bowtie2 index
        if build_bt2:
            bt2_build = self._tool_path("bowtie2-build")
            bt2_prefix = ref_dir / (fa_dst.stem + "_bt2")
            if bt2_build:
                if not self._bowtie2_index_exists(bt2_prefix):
                    self._run([bt2_build, str(fa_dst), str(bt2_prefix)], tag="bowtie2-build")
            else:
                self.log("[REF] bowtie2-build not found in PATH (bt2 index will be missing)")

        return self.check_status(str(fa_dst))

    def check_status(self, fasta_in_reference_dir: str) -> ReferenceIndexStatus:
        """Return a status snapshot of all reference artifacts and tools."""
        fa = Path(fasta_in_reference_dir).expanduser().resolve()
        ref_dir = fa.parent

        tools = {
            "samtools": self._tool_info("samtools"),
            "minimap2": self._tool_info("minimap2"),
            "bowtie2-build": self._tool_info("bowtie2-build"),
            "picard": self._tool_info("picard"),  # optional
        }

        faidx = fa.with_suffix(fa.suffix + ".fai")
        dict_path = ref_dir / (fa.stem + ".dict")
        mmi = ref_dir / (fa.stem + ".mmi")
        bt2_prefix = ref_dir / (fa.stem + "_bt2")

        status = ReferenceIndexStatus(
            fasta=str(fa),
            reference_dir=str(ref_dir),
            faidx=str(faidx) if faidx.exists() else None,
            dict=str(dict_path) if dict_path.exists() else None,
            mmi=str(mmi) if mmi.exists() else None,
            bt2_prefix=str(bt2_prefix) if self._bowtie2_index_exists(bt2_prefix) else None,
            tools=tools,
            ok=False,
        )
        status.ok = bool(status.faidx and (status.mmi or status.bt2_prefix))
        return status

    # ---------------------------- helpers ---------------------------- #

    def _default_reference_dir_for(self, fasta: Path) -> str:
        # Prefer project workspace if provided; otherwise sibling "reference" next to the FASTA.
        if self.workspace:
            return str(self.workspace / "reference")
        return str(fasta.parent / "reference")

    def _run(self, cmd: list[str], tag: str = "cmd") -> None:
        self.log(f"[REF] $ {' '.join(cmd)}")
        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        if proc.stdout:
            self.log(proc.stdout.strip())
        if proc.returncode != 0:
            raise RuntimeError(f"[REF] {tag} failed with exit code {proc.returncode}")

    def _tool_path(self, name: str) -> Optional[str]:
        path = resolve_tool(name)
        return str(path) if path else None

    def _tool_info(self, name: str) -> Dict[str, Optional[str]]:
        path = self._tool_path(name)
        version = None
        if path:
            version = self._probe_version(name, path)
        return {"path": path, "version": version}

    def _probe_version(self, name: str, path: str) -> Optional[str]:
        try:
            if name == "samtools":
                out = subprocess.run([path, "--version"], stdout=subprocess.PIPE, text=True).stdout.strip()
                return out.splitlines()[0] if out else None
            if name == "minimap2":
                out = subprocess.run([path, "--version"], stdout=subprocess.PIPE, text=True).stdout.strip()
                return out or None
            if name == "bowtie2-build":
                out = subprocess.run([path, "--version"], stdout=subprocess.PIPE, text=True).stdout.strip()
                return out or None
            if name == "picard":
                out = subprocess.run([path], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True).stdout.strip()
                # Picard typically prints usage with version on first line
                return out.splitlines()[0] if out else None
        except Exception:
            return None
        return None

    def _supports_samtools_dict(self, samtools_path: str) -> bool:
        try:
            out = subprocess.run([samtools_path, "dict", "-h"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
            return out.returncode == 0
        except Exception:
            return False

    def _try_picard_dict(self, fa_dst: Path, dict_path: Path) -> None:
        picard = self._tool_path("picard")
        if not picard:
            self.log("[REF] Picard not found; skipping .dict creation")
            return
        try:
            self._run([picard, "CreateSequenceDictionary", f"R={fa_dst}", f"O={dict_path}"], tag="picard")
        except Exception as e:
            self.log(f"[REF] Picard dict failed: {e}")

    def _bowtie2_index_exists(self, prefix: Path) -> bool:
        # Bowtie2 produces either .bt2 or .bt2l (large index) suffixes
        base = prefix.parent / prefix.name
        patterns = [
            f"{base}.1.bt2", f"{base}.2.bt2", f"{base}.3.bt2", f"{base}.4.bt2",
            f"{base}.rev.1.bt2", f"{base}.rev.2.bt2",
            f"{base}.1.bt2l", f"{base}.2.bt2l", f"{base}.3.bt2l", f"{base}.4.bt2l",
            f"{base}.rev.1.bt2l", f"{base}.rev.2.bt2l",
        ]
        return all(Path(p).exists() for p in patterns[:6]) or all(Path(p).exists() for p in patterns[6:])
