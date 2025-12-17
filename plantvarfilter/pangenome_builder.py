from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import List, Tuple, Optional, Dict, Union, Callable
import os
import time
import subprocess


FASTA_EXTS = (".fa", ".fasta", ".fna")


# -----------------------------
# Results
# -----------------------------
@dataclass
class PangenomeBuildResult:
    # output artifacts (either gfa or fasta depending on engine)
    pangenome_gfa: Optional[str] = None
    pangenome_fasta: Optional[str] = None

    report_txt: str = ""
    log_txt: str = ""

    inputs_used: List[str] = None
    renamed_inputs: List[str] = None
    included_files: List[str] = None
    skipped_files: List[str] = None

    total_sequences_written: int = 0
    total_bases_written: int = 0
    elapsed_seconds: float = 0.0


# -----------------------------
# Helpers
# -----------------------------
def _log(logger: Optional[Callable[[str], None]], msg: str) -> None:
    if logger:
        logger(msg)
    else:
        print(msg)


def _require_tool(tool: str) -> None:
    p = subprocess.run(["bash", "-lc", f"command -v {tool} >/dev/null 2>&1"], text=True)
    if p.returncode != 0:
        raise RuntimeError(f"Required tool not found in PATH: {tool}")


def _run_cmd(cmd: List[str], log_file: Path, cwd: Optional[Path] = None,
             logger: Optional[Callable[[str], None]] = None) -> None:
    _log(logger, f"[PAN] RUN: {' '.join(cmd)}")
    with log_file.open("a", encoding="utf-8") as lf:
        lf.write("\n" + "=" * 90 + "\n")
        lf.write("CMD: " + " ".join(cmd) + "\n")
        p = subprocess.run(
            cmd,
            cwd=str(cwd) if cwd else None,
            stdout=lf,
            stderr=lf,
            text=True
        )
    if p.returncode != 0:
        raise RuntimeError(f"Command failed: {' '.join(cmd)} (see log: {log_file})")


def _is_fasta(p: Path) -> bool:
    return p.is_file() and p.suffix.lower() in FASTA_EXTS


def _list_fastas_in_dir(dir_path: Path) -> List[Path]:
    files = [p for p in dir_path.iterdir() if _is_fasta(p)]
    return sorted(files, key=lambda x: x.name.lower())


def _resolve_assemblies_input(assemblies_input: Union[str, List[str]]) -> List[Path]:
    """
    Supports:
      - folder containing FASTA files
      - single FASTA file (multi-fasta allowed)
      - list of FASTA file paths
    """
    if isinstance(assemblies_input, list):
        paths = [Path(x).expanduser().resolve() for x in assemblies_input]
        for p in paths:
            if not p.exists() or not p.is_file():
                raise FileNotFoundError(f"FASTA file not found: {p}")
            if p.suffix.lower() not in FASTA_EXTS:
                raise ValueError(f"Not a FASTA file: {p}")
        return paths

    p = Path(assemblies_input).expanduser().resolve()
    if not p.exists():
        raise FileNotFoundError(f"Assemblies input not found: {p}")

    if p.is_dir():
        files = _list_fastas_in_dir(p)
        if not files:
            raise ValueError(f"No FASTA files found in folder: {p}")
        return files

    if p.is_file():
        if p.suffix.lower() not in FASTA_EXTS:
            raise ValueError(f"Assemblies file must be FASTA: {p}")
        return [p]

    raise ValueError(f"Unsupported assemblies_input: {assemblies_input}")


def _iter_fasta_records(path: Path):
    header = None
    seq_chunks: List[str] = []
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
    s = "".join(safe)
    return s[:80] if s else "sample"


def _samtools_faidx(fa: Path, log_file: Path, logger: Optional[Callable[[str], None]] = None) -> None:
    _run_cmd(["samtools", "faidx", str(fa)], log_file=log_file, logger=logger)


def _write_renamed_fasta(
    in_fa: Path,
    out_fa: Path,
    sample_prefix: str,
    min_contig_len: int = 0,
) -> Tuple[int, int]:
    """
    Rename headers to avoid collisions:
      >{sample_prefix}__{original_header}
    """
    seqs = 0
    bases = 0
    with out_fa.open("w", encoding="utf-8") as out:
        for h, s in _iter_fasta_records(in_fa):
            if min_contig_len and len(s) < int(min_contig_len):
                continue
            nh = f"{sample_prefix}__{h}"
            out.write(f">{nh}\n{s}\n")
            seqs += 1
            bases += len(s)
    return seqs, bases


# -----------------------------
# Engine 1: Graph (GFA) via minigraph
# -----------------------------
def build_pangenome_graph(
    base_reference_fasta: str,
    assemblies_input: Union[str, List[str]],
    output_dir: str,
    mode: str = "full",               # "full" or "fast"
    subset_n: int = 25,
    threads: int = 8,
    min_contig_len: int = 0,
    minigraph_preset: str = "ggs",    # uses -xggs
    batch_size: Optional[int] = None, # helpful for 100+ samples
    logger=None,
) -> PangenomeBuildResult:
    """
    Builds a TRUE pangenome graph reference (GFA) using minigraph:
      minigraph -c -xggs -t THREADS ref.fa sample1.renamed.fa ... > pangenome.gfa
    """

    t0 = time.time()
    _require_tool("minigraph")
    _require_tool("samtools")

    base_ref = Path(base_reference_fasta).expanduser().resolve()
    if not base_ref.exists() or not base_ref.is_file():
        raise FileNotFoundError(f"Base reference FASTA not found: {base_ref}")

    out_dir = Path(output_dir).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    log_txt = out_dir / "pangenome_graph.log"
    report_txt = out_dir / "pangenome_graph_report.txt"
    log_txt.write_text("", encoding="utf-8")

    work_dir = out_dir / "_pan_work"
    renamed_dir = out_dir / "renamed_inputs"
    work_dir.mkdir(parents=True, exist_ok=True)
    renamed_dir.mkdir(parents=True, exist_ok=True)

    _log(logger, f"[PAN] Graph build (minigraph) started")
    _log(logger, f"[PAN] Reference: {base_ref}")
    _samtools_faidx(base_ref, log_file=log_txt, logger=logger)

    inputs = _resolve_assemblies_input(assemblies_input)

    selected = inputs
    if str(mode).lower().startswith("fast") and len(inputs) > subset_n:
        selected = sorted(inputs, key=lambda p: p.stat().st_size, reverse=True)[:subset_n]

    renamed_inputs: List[str] = []
    inputs_used: List[str] = []

    total_seqs = 0
    total_bases = 0
    skipped: List[str] = []

    # prepare renamed sample fastas
    for fa in selected:
        sample = _sanitize_sample_name(fa)
        out_fa = renamed_dir / f"{sample}.renamed.fa"
        try:
            seqs, bases = _write_renamed_fasta(fa, out_fa, sample_prefix=sample, min_contig_len=min_contig_len)
            if seqs == 0:
                skipped.append(str(fa))
                _log(logger, f"[PAN] Skipped (empty after filters): {fa.name}")
                continue
            _samtools_faidx(out_fa, log_file=log_txt, logger=logger)
            renamed_inputs.append(str(out_fa))
            inputs_used.append(str(fa))
            total_seqs += seqs
            total_bases += bases
            _log(logger, f"[PAN] Prepared: {fa.name} -> {out_fa.name} | seqs={seqs} bases={bases}")
        except Exception as exc:
            skipped.append(str(fa))
            _log(logger, f"[PAN] Skipped due to error: {fa.name} | {exc}")

    if not renamed_inputs:
        raise ValueError("No valid input sequences after preparation (check min_contig_len and FASTA validity).")

    gfa_out = out_dir / "pangenome.gfa"
    xflag = f"-x{minigraph_preset}"
    if not xflag.startswith("-x"):
        xflag = "-xggs"

    def _run_minigraph(sample_fastas: List[str], out_gfa: Path) -> None:
        cmd = ["minigraph", "-c", xflag, "-t", str(max(1, int(threads))), str(base_ref)] + sample_fastas
        _log(logger, f"[PAN] RUN: {' '.join(cmd)} > {out_gfa}")
        with out_gfa.open("w", encoding="utf-8") as out_handle:
            with log_txt.open("a", encoding="utf-8") as lf:
                lf.write("\n" + "=" * 90 + "\n")
                lf.write("CMD: " + " ".join(cmd) + f" > {out_gfa}\n")
                p = subprocess.run(cmd, stdout=out_handle, stderr=lf, text=True)
        if p.returncode != 0:
            raise RuntimeError(f"minigraph failed (see log): {log_txt}")

    # batching (optional, helps with memory)
    if batch_size and batch_size > 0 and len(renamed_inputs) > batch_size:
        _log(logger, f"[PAN] Batch mode ON (batch_size={batch_size})")
        current: List[str] = []
        batch_idx = 0
        for i in range(0, len(renamed_inputs), batch_size):
            batch_idx += 1
            current.extend(renamed_inputs[i:i + batch_size])
            tmp = out_dir / f"pangenome.batch_{batch_idx}.gfa"
            _run_minigraph(current, tmp)
        # last batch becomes final
        tmp.replace(gfa_out)
    else:
        _run_minigraph(renamed_inputs, gfa_out)

    if not gfa_out.exists() or gfa_out.stat().st_size < 1024:
        raise RuntimeError("GFA output missing/too small. Check log.")

    dt = time.time() - t0

    with report_txt.open("w", encoding="utf-8") as rep:
        rep.write("PlantVarFilter - Pangenome Graph Builder (minigraph)\n")
        rep.write(f"Reference FASTA: {base_ref}\n")
        rep.write(f"Assemblies input: {assemblies_input}\n")
        rep.write(f"Output dir: {out_dir}\n")
        rep.write(f"Mode: {mode}\n")
        rep.write(f"Subset N: {subset_n}\n")
        rep.write(f"Threads: {threads}\n")
        rep.write(f"Min contig len: {min_contig_len}\n")
        rep.write(f"minigraph preset: {minigraph_preset} ({xflag})\n")
        rep.write(f"Batch size: {batch_size}\n")
        rep.write(f"Inputs found: {len(inputs)}\n")
        rep.write(f"Inputs selected: {len(selected)}\n")
        rep.write(f"Inputs used: {len(inputs_used)}\n")
        rep.write(f"Inputs skipped: {len(skipped)}\n")
        rep.write(f"Sample sequences kept: {total_seqs}\n")
        rep.write(f"Sample bases kept: {total_bases}\n")
        rep.write(f"Pangenome GFA: {gfa_out}\n")
        rep.write(f"Log: {log_txt}\n")
        rep.write(f"Elapsed seconds: {dt:.2f}\n\n")
        rep.write("Inputs used:\n")
        for p in inputs_used:
            rep.write(f"- {p}\n")
        if skipped:
            rep.write("\nSkipped:\n")
            for p in skipped:
                rep.write(f"- {p}\n")
        rep.write("\nRenamed inputs:\n")
        for p in renamed_inputs:
            rep.write(f"- {p}\n")

    _log(logger, f"[PAN] Done ✔ GFA: {gfa_out}")
    _log(logger, f"[PAN] Report: {report_txt}")
    _log(logger, f"[PAN] Log: {log_txt}")

    return PangenomeBuildResult(
        pangenome_gfa=str(gfa_out),
        pangenome_fasta=None,
        report_txt=str(report_txt),
        log_txt=str(log_txt),
        inputs_used=inputs_used,
        renamed_inputs=renamed_inputs,
        included_files=inputs_used,
        skipped_files=skipped,
        total_sequences_written=total_seqs,
        total_bases_written=total_bases,
        elapsed_seconds=dt,
    )


# -----------------------------
# Engine 2: Pan-FASTA (novel contigs) - optional legacy/quick
# (You can keep your older idea but name it correctly.)
# -----------------------------
def build_pan_fasta_novel_contigs(
    base_reference_fasta: str,
    assemblies_dir: str,
    output_dir: str,
    mode: str = "Fast (subset 25)",
    min_contig_len: int = 1000,
    novelty_cov_threshold: float = 0.70,
    threads: int = 8,
    minimap2_preset: str = "asm5",
    logger=None,
) -> PangenomeBuildResult:
    """
    Builds a LINEAR pan-FASTA by appending contigs considered "novel" (low coverage vs reference).
    This is NOT a graph pangenome.
    """

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

    def _run_minimap2_paf(ref_fa: Path, qry_fa: Path, paf_out: Path) -> None:
        cmd = ["minimap2", "-x", minimap2_preset, "-t", str(max(1, int(threads))),
               "--secondary=no", str(ref_fa), str(qry_fa)]
        _run_cmd(cmd, log_file=log_txt, logger=logger)
        # minimap2 PAF to stdout; re-run capturing stdout
        with paf_out.open("w", encoding="utf-8") as out:
            p = subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, text=True)
        if p.returncode != 0:
            raise RuntimeError(p.stderr.strip())

    t0 = time.time()
    _require_tool("minimap2")

    base_ref = Path(base_reference_fasta).expanduser().resolve()
    asm_dir = Path(assemblies_dir).expanduser().resolve()
    out_dir = Path(output_dir).expanduser().resolve()

    if not base_ref.exists() or not base_ref.is_file():
        raise FileNotFoundError(f"Base reference FASTA not found: {base_ref}")
    if not asm_dir.exists() or not asm_dir.is_dir():
        raise FileNotFoundError(f"Assemblies folder not found: {asm_dir}")

    out_dir.mkdir(parents=True, exist_ok=True)
    log_txt = out_dir / "pan_fasta.log"
    report = out_dir / "pan_fasta_report.txt"
    log_txt.write_text("", encoding="utf-8")

    assemblies = _list_fastas_in_dir(asm_dir)
    if not assemblies:
        raise ValueError(f"No FASTA assemblies found in: {asm_dir}")

    chosen = assemblies
    if str(mode).lower().startswith("fast"):
        chosen = sorted(assemblies, key=lambda p: p.stat().st_size, reverse=True)[:25]

    pan_fa = out_dir / "pangenome.fa"
    work_dir = out_dir / "_pan_work"
    work_dir.mkdir(parents=True, exist_ok=True)

    included: List[str] = []
    skipped: List[str] = []
    total_seqs = 0
    total_bases = 0
    novel_total = 0

    with pan_fa.open("w", encoding="utf-8") as out:
        # write reference
        for h, s in _iter_fasta_records(base_ref):
            if min_contig_len and len(s) < int(min_contig_len):
                continue
            out.write(f">{h}\n{s}\n")
            total_seqs += 1
            total_bases += len(s)

        # append "novel" contigs
        for asm in chosen:
            sample = _sanitize_sample_name(asm)
            paf_path = work_dir / f"{sample}.paf"
            try:
                _run_minimap2_paf(base_ref, asm, paf_path)
                cov = _paf_query_coverage(paf_path)

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
                    kept_this += 1
                    novel_total += 1

                if kept_this > 0:
                    included.append(str(asm))
                else:
                    skipped.append(str(asm))

                _log(logger, f"[PAN-FASTA] {asm.name}: scanned={seen} kept_novel={kept_this}")

            except Exception as exc:
                skipped.append(str(asm))
                _log(logger, f"[PAN-FASTA] {asm.name}: skipped due to error: {exc}")

    dt = time.time() - t0

    with report.open("w", encoding="utf-8") as rep:
        rep.write("PlantVarFilter - Pan-FASTA Builder (Novel contigs)\n")
        rep.write(f"Base reference: {base_ref}\n")
        rep.write(f"Assemblies dir: {asm_dir}\n")
        rep.write(f"Output dir: {out_dir}\n")
        rep.write(f"Mode: {mode}\n")
        rep.write(f"Min contig len: {min_contig_len}\n")
        rep.write(f"Novelty cov thr: {novelty_cov_threshold}\n")
        rep.write(f"Threads: {threads}\n")
        rep.write(f"minimap2 preset: {minimap2_preset}\n")
        rep.write(f"Assemblies found: {len(assemblies)}\n")
        rep.write(f"Assemblies selected: {len(chosen)}\n")
        rep.write(f"Included files: {len(included)}\n")
        rep.write(f"Skipped files: {len(skipped)}\n")
        rep.write(f"Total sequences written: {total_seqs}\n")
        rep.write(f"Total bases written: {total_bases}\n")
        rep.write(f"Novel contigs written: {novel_total}\n")
        rep.write(f"Elapsed seconds: {dt:.2f}\n")

    return PangenomeBuildResult(
        pangenome_gfa=None,
        pangenome_fasta=str(pan_fa),
        report_txt=str(report),
        log_txt=str(log_txt),
        inputs_used=included + skipped,
        renamed_inputs=[],
        included_files=included,
        skipped_files=skipped,
        total_sequences_written=total_seqs,
        total_bases_written=total_bases,
        elapsed_seconds=dt,
    )
