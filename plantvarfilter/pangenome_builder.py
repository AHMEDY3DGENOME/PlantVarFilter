from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import List, Tuple, Optional, Dict
import os
import time
import subprocess


FASTA_EXTS = (".fa", ".fasta", ".fna")


@dataclass
class PangenomeBuildResult:
    pangenome_fasta: str
    report_txt: str
    included_files: List[str]
    skipped_files: List[str]
    total_sequences_written: int
    total_bases_written: int


def _iter_fasta_records(path: Path):
    header = None
    seq_chunks = []
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_chunks)
                header = line[1:].strip()
                seq_chunks = []
            else:
                seq_chunks.append(line)
        if header is not None:
            yield header, "".join(seq_chunks)


def _sanitize_sample_name(p: Path) -> str:
    name = p.stem
    safe = []
    for ch in name:
        if ch.isalnum() or ch in ("_", "-", "."):
            safe.append(ch)
        else:
            safe.append("_")
    return "".join(safe)[:80] if safe else "sample"


def _list_assemblies(dir_path: Path) -> List[Path]:
    files = []
    for p in dir_path.iterdir():
        if p.is_file() and p.suffix.lower() in FASTA_EXTS:
            files.append(p)
    return sorted(files, key=lambda x: x.name.lower())


def _merge_intervals(intervals: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    if not intervals:
        return []
    intervals.sort()
    merged = [intervals[0]]
    for s, e in intervals[1:]:
        ps, pe = merged[-1]
        if s <= pe:
            merged[-1] = (ps, max(pe, e))
        else:
            merged.append((s, e))
    return merged


def _paf_query_coverage(paf_path: Path) -> Dict[str, float]:
    intervals_by_q: Dict[str, List[Tuple[int, int]]] = {}
    qlen_by_q: Dict[str, int] = {}

    with paf_path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 12:
                continue

            qname = parts[0]
            qlen = int(parts[1])
            qstart = int(parts[2])
            qend = int(parts[3])

            qlen_by_q[qname] = qlen
            intervals_by_q.setdefault(qname, []).append((qstart, qend))

    cov: Dict[str, float] = {}
    for qname, intervals in intervals_by_q.items():
        qlen = qlen_by_q.get(qname, 0)
        if qlen <= 0:
            cov[qname] = 0.0
            continue
        merged = _merge_intervals(intervals)
        covered = sum(e - s for s, e in merged)
        cov[qname] = covered / float(qlen)

    return cov


def _run_minimap2_paf(ref_fa: Path, qry_fa: Path, paf_out: Path, threads: int, preset: str, logger=None) -> None:
    cmd = [
        "minimap2",
        "-x", preset,
        "-t", str(max(1, int(threads))),
        str(ref_fa),
        str(qry_fa),
    ]
    if logger:
        logger(f"[PAN] minimap2: {' '.join(cmd)}")

    with paf_out.open("w", encoding="utf-8") as out:
        p = subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, text=True)

    if p.returncode != 0:
        raise RuntimeError(f"minimap2 failed for {qry_fa.name}: {p.stderr.strip()}")


def build_pangenome_mvp(
    base_reference_fasta: str,
    assemblies_dir: str,
    output_dir: str,
    mode: str = "Fast (subset 25)",
    min_contig_len: int = 1000,
    logger=None,
) -> PangenomeBuildResult:
    t0 = time.time()

    base_ref = Path(base_reference_fasta).expanduser().resolve()
    asm_dir = Path(assemblies_dir).expanduser().resolve()
    out_dir = Path(output_dir).expanduser().resolve()

    if not base_ref.exists() or not base_ref.is_file():
        raise FileNotFoundError(f"Base reference FASTA not found: {base_ref}")
    if not asm_dir.exists() or not asm_dir.is_dir():
        raise FileNotFoundError(f"Assemblies folder not found: {asm_dir}")

    out_dir.mkdir(parents=True, exist_ok=True)

    def log(msg: str):
        if logger:
            logger(msg)
        else:
            print(msg)

    assemblies = _list_assemblies(asm_dir)
    if not assemblies:
        raise ValueError(f"No FASTA assemblies found in: {asm_dir}")

    chosen = assemblies
    if str(mode).lower().startswith("fast"):
        scored = sorted(
            assemblies,
            key=lambda p: (p.stat().st_size if p.exists() else 0),
            reverse=True,
        )
        chosen = scored[:25]

    pan_fa = out_dir / "pangenome.fa"
    report = out_dir / "pangenome_build_report.txt"
    work_dir = out_dir / "_pan_work"
    work_dir.mkdir(parents=True, exist_ok=True)

    threads = int(os.environ.get("PLANTVARFILTER_THREADS", "8"))
    novelty_cov_threshold = float(os.environ.get("PLANTVARFILTER_PAN_COV_THR", "0.70"))
    preset = os.environ.get("PLANTVARFILTER_MINIMAP2_PRESET", "asm5")

    included: List[str] = []
    skipped: List[str] = []
    total_seqs = 0
    total_bases = 0
    novel_total = 0

    log("[PAN] Reference-guided pangenome build (minimap2 novelty filter + write novel contigs).")
    log(f"[PAN] Base reference: {base_ref}")
    log(f"[PAN] Assemblies dir:  {asm_dir}")
    log(f"[PAN] Output dir:     {out_dir}")
    log(f"[PAN] Mode:           {mode}")
    log(f"[PAN] Min contig len: {min_contig_len}")
    log(f"[PAN] Novelty cov thr: {novelty_cov_threshold}")
    log(f"[PAN] Threads:        {threads}")
    log(f"[PAN] minimap2 preset: {preset}")
    log(f"[PAN] Assemblies found: {len(assemblies)} | Selected: {len(chosen)}")

    with pan_fa.open("w", encoding="utf-8") as out:
        for h, s in _iter_fasta_records(base_ref):
            if min_contig_len and len(s) < int(min_contig_len):
                continue
            out.write(f">{h}\n{s}\n")
            total_seqs += 1
            total_bases += len(s)

        for asm in chosen:
            sample = _sanitize_sample_name(asm)
            paf_path = work_dir / f"{sample}.paf"

            try:
                _run_minimap2_paf(
                    ref_fa=base_ref,
                    qry_fa=asm,
                    paf_out=paf_path,
                    threads=threads,
                    preset=preset,
                    logger=log,
                )

                cov = _paf_query_coverage(paf_path)

                wrote_any = False
                kept_this = 0
                seen = 0

                for h, s in _iter_fasta_records(asm):
                    if min_contig_len and len(s) < int(min_contig_len):
                        continue

                    seen += 1
                    c = cov.get(h, 0.0)
                    is_novel = (h not in cov) or (c < novelty_cov_threshold)
                    if not is_novel:
                        continue

                    nh = f"{sample}|{h}"
                    out.write(f">{nh}\n{s}\n")
                    total_seqs += 1
                    total_bases += len(s)
                    wrote_any = True
                    kept_this += 1
                    novel_total += 1

                    if novel_total % 500 == 0:
                        log(f"[PAN] Novel contigs written: {novel_total}")

                if wrote_any:
                    included.append(str(asm))
                    log(f"[PAN] {asm.name}: scanned={seen} kept_novel={kept_this}")
                else:
                    skipped.append(str(asm))
                    log(f"[PAN] {asm.name}: scanned={seen} kept_novel=0 (coverage >= {novelty_cov_threshold})")

            except Exception as exc:
                skipped.append(str(asm))
                log(f"[PAN] {asm.name}: skipped due to error: {exc}")

    dt = time.time() - t0

    with report.open("w", encoding="utf-8") as rep:
        rep.write("PlantVarFilter - Pangenome Builder (Reference-guided)\n")
        rep.write(f"Base reference: {base_ref}\n")
        rep.write(f"Assemblies dir:  {asm_dir}\n")
        rep.write(f"Output dir:      {out_dir}\n")
        rep.write(f"Mode:            {mode}\n")
        rep.write(f"Min contig len:  {min_contig_len}\n")
        rep.write(f"Novelty cov thr: {novelty_cov_threshold}\n")
        rep.write(f"Threads:         {threads}\n")
        rep.write(f"minimap2 preset: {preset}\n")
        rep.write(f"Assemblies found: {len(assemblies)}\n")
        rep.write(f"Assemblies selected: {len(chosen)}\n")
        rep.write(f"Included files:  {len(included)}\n")
        rep.write(f"Skipped files:   {len(skipped)}\n")
        rep.write(f"Total sequences written: {total_seqs}\n")
        rep.write(f"Total bases written:     {total_bases}\n")
        rep.write(f"Novel contigs written:   {novel_total}\n")
        rep.write(f"Elapsed seconds:         {dt:.2f}\n")
        rep.write("\nIncluded:\n")
        for p in included:
            rep.write(f"- {p}\n")
        rep.write("\nSkipped:\n")
        for p in skipped:
            rep.write(f"- {p}\n")

    log(f"[PAN] Done ✔ Output: {pan_fa}")
    log(f"[PAN] Report: {report}")

    return PangenomeBuildResult(
        pangenome_fasta=str(pan_fa),
        report_txt=str(report),
        included_files=included,
        skipped_files=skipped,
        total_sequences_written=total_seqs,
        total_bases_written=total_bases,
    )
