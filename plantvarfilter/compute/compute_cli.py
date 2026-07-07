"""
PlantOmicsGWAS Compute CLI

Command-line entry point for the headless Compute/HPC engine.
"""

import argparse
import platform
import sys

from plantvarfilter.compute.config_loader import load_config
from plantvarfilter.compute.workflow_engine import run_workflow
from plantvarfilter.compute.workflow_registry import list_steps
from plantvarfilter.hpc import HPCManager


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="plantomicsgwas-compute",
        description="PlantOmicsGWAS Compute Engine for CLI, server, and HPC execution.",
    )

    subparsers = parser.add_subparsers(dest="command", required=True)

    run_parser = subparsers.add_parser(
        "run",
        help="Run a PlantOmicsGWAS workflow from a YAML config file.",
    )
    run_parser.add_argument("--config", "-c", required=True)
    run_parser.add_argument("--dry-run", action="store_true")

    write_job_parser = subparsers.add_parser(
        "write-job",
        help="Generate an HPC job script from a YAML config file.",
    )
    write_job_parser.add_argument("--config", "-c", required=True)
    write_job_parser.add_argument("--out", "-o", default=None)

    submit_parser = subparsers.add_parser(
        "submit",
        help="Generate and submit an HPC job script.",
    )
    submit_parser.add_argument("--config", "-c", required=True)
    submit_parser.add_argument("--out", "-o", default=None)
    submit_parser.add_argument("--no-submit", action="store_true")

    cluster_parser = subparsers.add_parser(
        "cluster-info",
        help="Show detected HPC scheduler information.",
    )
    cluster_parser.add_argument("--config", "-c", required=True)

    status_parser = subparsers.add_parser(
        "status",
        help="Show scheduler job status.",
    )
    status_parser.add_argument("--config", "-c", required=True)
    status_parser.add_argument("--job-id", "-j", default=None)

    cancel_parser = subparsers.add_parser(
        "cancel",
        help="Cancel a scheduler job.",
    )
    cancel_parser.add_argument("--config", "-c", required=True)
    cancel_parser.add_argument("--job-id", "-j", required=True)

    subparsers.add_parser(
        "list-steps",
        help="List all available Compute Engine workflow steps.",
    )

    return parser


def _manager(config_path: str) -> HPCManager:
    config = load_config(config_path)
    return HPCManager(config=config, config_path=config_path)


def _is_windows() -> bool:
    return platform.system().lower() == "windows"


def command_run(args: argparse.Namespace) -> int:
    try:
        run_workflow(config_path=args.config, dry_run=args.dry_run)
        return 0
    except Exception as exc:
        print(f"\nERROR: {exc}", file=sys.stderr)
        return 1


def command_write_job(args: argparse.Namespace) -> int:
    try:
        manager = _manager(args.config)
        job_path = manager.write_job(output_path=args.out)
        print(f"HPC job script written: {job_path}")
        return 0
    except Exception as exc:
        print(f"\nERROR: {exc}", file=sys.stderr)
        return 1


def command_submit(args: argparse.Namespace) -> int:
    try:
        if _is_windows() and not args.no_submit:
            print()
            print("HPC scheduler submission is not supported directly on Windows.")
            print()
            print("Use this command to generate the job script:")
            print("  plantomicsgwas-compute write-job --config <config.yaml>")
            print()
            print("Then copy the generated .slurm/.pbs/.lsf file to your HPC cluster")
            print("and submit it there using sbatch, qsub, or bsub.")
            print()
            return 1

        manager = _manager(args.config)
        result = manager.submit(
            output_path=args.out,
            no_submit=args.no_submit,
        )

        print(f"HPC job script written: {result.job_script}")

        if result.command:
            print("Command:")
            print(" ".join(result.command))

        if result.stdout:
            print(result.stdout)

        if result.job_id:
            print(f"Job ID: {result.job_id}")

        return result.return_code

    except Exception as exc:
        print(f"\nERROR: {exc}", file=sys.stderr)
        return 1


def command_cluster_info(args: argparse.Namespace) -> int:
    try:
        manager = _manager(args.config)
        print(manager.cluster_info())
        return 0
    except Exception as exc:
        print(f"\nERROR: {exc}", file=sys.stderr)
        return 1


def command_status(args: argparse.Namespace) -> int:
    try:
        manager = _manager(args.config)
        print(manager.status(job_id=args.job_id))
        return 0
    except Exception as exc:
        print(f"\nERROR: {exc}", file=sys.stderr)
        return 1


def command_cancel(args: argparse.Namespace) -> int:
    try:
        manager = _manager(args.config)
        print(manager.cancel(job_id=args.job_id))
        return 0
    except Exception as exc:
        print(f"\nERROR: {exc}", file=sys.stderr)
        return 1


def command_list_steps() -> int:
    print("\nAvailable PlantOmicsGWAS Compute Engine steps")
    print("=" * 55)

    for step in list_steps():
        print(f"{step.step_id:22s} | {step.category:20s} | {step.name}")

    print("=" * 55)
    return 0


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    if args.command == "run":
        return command_run(args)

    if args.command == "write-job":
        return command_write_job(args)

    if args.command == "submit":
        return command_submit(args)

    if args.command == "cluster-info":
        return command_cluster_info(args)

    if args.command == "status":
        return command_status(args)

    if args.command == "cancel":
        return command_cancel(args)

    if args.command == "list-steps":
        return command_list_steps()

    parser.print_help()
    return 1


if __name__ == "__main__":
    raise SystemExit(main())