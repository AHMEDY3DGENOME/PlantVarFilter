from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Union, Callable
import time
import subprocess

FASTA_EXTS = (".fa", ".fasta", ".fna")


@dataclass
class PangenomeBuildResult:
    pangenome_gfa: Optional[str] = None
    pangenome_vg: Optional[str] = None
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


def _log(logger: Optional[Callable[[str], None]], msg: str) -> None:
    if logger:
        logger(msg)
    else:
        print(msg)


def _require_tool(tool: str) -> str:
    p = subprocess.run(["bash", "-lc", f"command -v {tool}"], text=True, capture_output=True)
    if p.returncode != 0:
        raise RuntimeError(f"Required tool not found in PATH: {tool}")
    return p.stdout.strip()


def _is_fasta(p: Path) -> bool:
    return p.is_file() and p.suffix.lower() in FASTA_EXTS


def _list_fastas_in_dir(dir_path: Path) -> List[Path]:
    files = [p for p in dir_path.iterdir() if _is_fasta(p)]
    return sorted(files, key=lambda x: x.name.lower())


def _resolve_assemblies_input(assemblies_input: Union[str, List[str]]) -> List[Path]:
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


def _pick_backbone(paths: List[Path], strategy: str = "largest") -> Path:
    if not paths:
        raise ValueError("No FASTA inputs provided.")
    st = (strategy or "").strip().lower()
    if st in ("first", "input_order"):
        return paths[0]
    if st in ("largest", "biggest", "max_size"):
        return max(paths, key=lambda p: p.stat().st_size)
    raise ValueError(f"Unsupported backbone_strategy: {strategy}")


def _rc_human(rc: int) -> str:
    if rc == 0:
        return "0 (success)"
    if rc > 0:
        return f"{rc} (non-zero exit code)"
    # negative => killed by signal
    return f"{rc} (killed by signal {abs(rc)})"


def _likely_oom(rc: int, log_tail: str) -> bool:
    # Typical SIGKILL patterns and hints in logs
    if rc in (-9, 137):
        return True
    lt = (log_tail or "").lower()
    if "out of memory" in lt or "killed process" in lt or "oom" in lt:
        return True
    return False


def _tail_text(path: Path, n_lines: int = 80) -> str:
    if not path.exists():
        return ""
    try:
        lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()
        return "\n".join(lines[-n_lines:])
    except Exception:
        return ""


def _run_cmd_capture(
    cmd: List[str],
    log_file: Path,
    stdout_to_file: Optional[Path] = None,
    cwd: Optional[Path] = None,
    logger: Optional[Callable[[str], None]] = None,
    use_time_v: bool = True,
) -> int:
    """
    Runs a command and appends verbose logs.
    - Writes CMD + timestamps + EXIT_CODE into log_file
    - If stdout_to_file is provided, stdout is redirected to that file, stderr goes to log_file
    - If /usr/bin/time exists and use_time_v=True, wraps the command with `/usr/bin/time -v`
      to capture Max RSS/CPU/elapsed in the same log file (very useful to debug OOM).
    """
    log_file.parent.mkdir(parents=True, exist_ok=True)
    if stdout_to_file is not None:
        stdout_to_file.parent.mkdir(parents=True, exist_ok=True)

    start_ts = time.strftime("%Y-%m-%d %H:%M:%S")
    _log(logger, f"[PAN] RUN: {' '.join(cmd)}")
    _log(logger, f"[PAN] Log file: {log_file}")

    # Build a bash command string so we can do redirects + optional /usr/bin/time -v cleanly
    # NOTE: we keep stderr in log_file; stdout either goes to file or log.
    safe_cmd = " ".join([subprocess.list2cmdline([c]) if " " in c else c for c in cmd])
    # The above is conservative but not perfect for edge chars; your inputs are file paths so OK.

    time_prefix = ""
    if use_time_v:
        # use absolute /usr/bin/time if available; fallback silently if missing
        time_prefix = 'command -v /usr/bin/time >/dev/null 2>&1 && TIMEP="/usr/bin/time -v" || TIMEP="" ; '

    if stdout_to_file is None:
        bash_line = f'{time_prefix} $TIMEP {safe_cmd}'
    else:
        bash_line = f'{time_prefix} $TIMEP {safe_cmd} 1> "{stdout_to_file}"'

    with log_file.open("a", encoding="utf-8") as lf:
        lf.write("\n" + "=" * 110 + "\n")
        lf.write(f"START: {start_ts}\n")
        lf.write("CWD:   " + (str(cwd) if cwd else "<none>") + "\n")
        lf.write("CMD:   " + " ".join(cmd) + "\n")
        lf.write("-" * 110 + "\n")
        lf.flush()

        p = subprocess.run(
            ["bash", "-lc", bash_line],
            cwd=str(cwd) if cwd else None,
            stdout=lf if stdout_to_file is None else subprocess.DEVNULL,
            stderr=lf,
            text=True,
        )

        end_ts = time.strftime("%Y-%m-%d %H:%M:%S")
        lf.write("-" * 110 + "\n")
        lf.write(f"END:       {end_ts}\n")
        lf.write(f"EXIT_CODE: {p.returncode}\n")
        lf.write(f"EXIT_HUMAN:{_rc_human(p.returncode)}\n")
        lf.flush()

    # Push a short, visible status into the UI log immediately
    _log(logger, f"[PAN] EXIT: {_rc_human(p.returncode)}")

    if stdout_to_file is not None:
        sz = stdout_to_file.stat().st_size if stdout_to_file.exists() else 0
        _log(logger, f"[PAN] STDOUT file: {stdout_to_file} (size={sz})")

    # If failed, surface a short tail of the log to UI (very important for “not fake logs”)
    if p.returncode != 0:
        tail = _tail_text(log_file, n_lines=60)
        if tail.strip():
            _log(logger, "[PAN] --- log tail (last 60 lines) ---")
            for line in tail.splitlines()[-60:]:
                _log(logger, line)
            _log(logger, "[PAN] --- end log tail ---")

        if _likely_oom(p.returncode, tail):
            _log(
                logger,
                "[PAN][OOM] The process was likely killed by the kernel due to Out-Of-Memory (OOM). "
                "Try: reduce threads (PLANTVARFILTER_THREADS=1 or 2), close heavy apps, or add swap.",
            )

    return p.returncode


def build_pangenome_graph(
    assemblies_input: Union[str, List[str]],
    output_dir: str,
    mode: str = "full",
    subset_n: int = 25,
    threads: int = 8,
    minigraph_preset: str = "ggs",
    backbone_strategy: str = "largest",
    logger=None,
    vg_gzip_only: bool = False,
) -> PangenomeBuildResult:
    t0 = time.time()

    minigraph_path = _require_tool("minigraph")
    vg_path = _require_tool("vg")

    out_dir = Path(output_dir).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    log_txt = out_dir / "pangenome_graph.log"
    report_txt = out_dir / "pangenome_graph_report.txt"
    log_txt.write_text("", encoding="utf-8")

    _log(logger, "[PAN] Direct build started (minigraph -> GFA, vg convert -> VG)")
    _log(logger, f"[PAN] which minigraph: {minigraph_path}")
    _log(logger, f"[PAN] which vg:       {vg_path}")
    _log(logger, f"[PAN] output_dir:     {out_dir}")
    _log(logger, f"[PAN] threads:        {threads}")
    _log(logger, f"[PAN] vg_gzip_only:   {vg_gzip_only}")

    inputs = _resolve_assemblies_input(assemblies_input)

    selected = inputs
    if str(mode).lower().startswith("fast") and len(inputs) > subset_n:
        selected = sorted(inputs, key=lambda p: p.stat().st_size, reverse=True)[:subset_n]
        _log(logger, f"[PAN] Mode=FAST: selected top {len(selected)}/{len(inputs)} by size")
    else:
        _log(logger, f"[PAN] Mode=FULL: selected all {len(selected)} inputs")

    if len(selected) < 2:
        raise ValueError("Need at least 2 FASTA inputs to build a pangenome graph with minigraph.")

    backbone = _pick_backbone(selected, strategy=backbone_strategy)
    samples = [p for p in selected if p != backbone]

    _log(logger, f"[PAN] Backbone: {backbone}")
    _log(logger, f"[PAN] Samples:  {len(samples)} file(s)")

    gfa_out = out_dir / "pangenome.gfa"
    vg_out = out_dir / "pangenome.vg"
    vgz_out = out_dir / "pangenome.vg.gz"

    xflag = f"-x{minigraph_preset}".strip()
    if not xflag.startswith("-x"):
        xflag = "-xggs"

    cmd_mg = ["minigraph", "-c", xflag, "-t", str(max(1, int(threads))), str(backbone)] + [str(p) for p in samples]

    _log(logger, f"[PAN] Step 1/2: minigraph -> {gfa_out}")
    rc = _run_cmd_capture(cmd_mg, log_file=log_txt, stdout_to_file=gfa_out, logger=logger, use_time_v=True)
    if rc != 0:
        raise RuntimeError(f"minigraph failed: {_rc_human(rc)}. See log: {log_txt}")

    gfa_size = gfa_out.stat().st_size if gfa_out.exists() else 0
    _log(logger, f"[PAN] GFA written: {gfa_out} (size={gfa_size})")
    if gfa_size < 1024:
        tail = _tail_text(log_txt, n_lines=80)
        raise RuntimeError(
            f"GFA output is missing/too small (size={gfa_size}). "
            f"See log: {log_txt}\n"
            f"--- log tail ---\n{tail}\n--- end ---"
        )

    _log(logger, f"[PAN] Step 2/2: vg convert -> {vg_out}")
    cmd_vg = ["vg", "convert", "-g", "-p", str(gfa_out)]
    rc2 = _run_cmd_capture(cmd_vg, log_file=log_txt, stdout_to_file=vg_out, logger=logger, use_time_v=True)
    if rc2 != 0:
        raise RuntimeError(f"vg convert failed: {_rc_human(rc2)}. See log: {log_txt}")

    vg_size = vg_out.stat().st_size if vg_out.exists() else 0
    _log(logger, f"[PAN] VG written: {vg_out} (size={vg_size})")
    if vg_size < 1024:
        tail = _tail_text(log_txt, n_lines=80)
        raise RuntimeError(
            f"VG output is missing/too small (size={vg_size}). See log: {log_txt}\n"
            f"--- log tail ---\n{tail}\n--- end ---"
        )

    final_vg_path = vg_out

    if vg_gzip_only:
        _log(logger, f"[PAN] Step 3/3: gzip -> {vgz_out}")
        cmd_gz = ["gzip", "-c", str(vg_out)]
        rc3 = _run_cmd_capture(cmd_gz, log_file=log_txt, stdout_to_file=vgz_out, logger=logger, use_time_v=True)
        if rc3 != 0:
            raise RuntimeError(f"gzip failed: {_rc_human(rc3)}. See log: {log_txt}")

        vgz_size = vgz_out.stat().st_size if vgz_out.exists() else 0
        _log(logger, f"[PAN] VG.GZ written: {vgz_out} (size={vgz_size})")
        if vgz_size < 1024:
            tail = _tail_text(log_txt, n_lines=80)
            raise RuntimeError(
                f"VG.GZ output is missing/too small (size={vgz_size}). See log: {log_txt}\n"
                f"--- log tail ---\n{tail}\n--- end ---"
            )

        try:
            vg_out.unlink(missing_ok=True)
        except Exception:
            pass

        final_vg_path = vgz_out

    dt = time.time() - t0

    with report_txt.open("w", encoding="utf-8") as rep:
        rep.write("PlantOmicsGwas - Pangenome Graph Builder (Direct)\n")
        rep.write(f"Assemblies input: {assemblies_input}\n")
        rep.write(f"Output dir: {out_dir}\n")
        rep.write(f"Mode: {mode}\n")
        rep.write(f"Subset N: {subset_n}\n")
        rep.write(f"Threads: {threads}\n")
        rep.write(f"minigraph preset: {minigraph_preset} ({xflag})\n")
        rep.write(f"Backbone strategy: {backbone_strategy}\n")
        rep.write(f"Backbone: {backbone}\n")
        rep.write(f"Inputs found: {len(inputs)}\n")
        rep.write(f"Inputs selected: {len(selected)}\n")
        rep.write(f"Pangenome GFA: {gfa_out} (size={gfa_size})\n")
        if vg_gzip_only:
            rep.write(f"Pangenome VG.GZ: {vgz_out} (size={vgz_out.stat().st_size if vgz_out.exists() else 0})\n")
        else:
            rep.write(f"Pangenome VG: {vg_out} (size={vg_size})\n")
        rep.write(f"Log: {log_txt}\n")
        rep.write(f"Elapsed seconds: {dt:.2f}\n")

    _log(logger, f"[PAN] Done ✔ elapsed={dt:.2f}s")
    _log(logger, f"[PAN] Outputs: {gfa_out} , {final_vg_path}")

    return PangenomeBuildResult(
        pangenome_gfa=str(gfa_out),
        pangenome_vg=str(final_vg_path),
        pangenome_fasta=None,
        report_txt=str(report_txt),
        log_txt=str(log_txt),
        inputs_used=[str(p) for p in selected],
        renamed_inputs=[],
        included_files=[str(p) for p in selected],
        skipped_files=[],
        total_sequences_written=0,
        total_bases_written=0,
        elapsed_seconds=dt,
    )

