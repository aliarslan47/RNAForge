"""Komut satırı girişi."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from rnaforge import __version__
from rnaforge.config import ConfigError, load_config
from rnaforge.metadata import MetadataError
from rnaforge.modules.m01_validate import run_validation
from rnaforge.platform import UnsupportedPlatformError
from rnaforge.state import new_run_dir


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="rnaforge", description="Bulk RNA-seq pipeline")
    parser.add_argument("--version", action="version", version=__version__)
    sub = parser.add_subparsers(dest="command")

    validate = sub.add_parser("validate", help="validate config, metadata and FASTQ inputs")
    validate.add_argument("--config", required=True, type=Path)
    validate.add_argument("--metadata", required=True, type=Path)
    validate.add_argument("--runs-dir", type=Path, default=Path("runs"))
    validate.add_argument("--run-id", default="run")
    return parser


def _cmd_validate(args) -> int:
    config = load_config(args.config)
    run_dir = new_run_dir(args.runs_dir, args.run_id)
    summary = run_validation(config, args.metadata, run_dir)
    print(
        f"validation OK: {summary['n_samples']} sample(s), "
        f"platform={summary['platform']}, "
        f"organism_type={summary['organism_type']}, "
        f"conditions={summary['conditions']}"
    )
    print(f"run directory: {run_dir}")
    return 0


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command is None:
        print("error: no command given (try: rnaforge validate --help)")
        return 2
    try:
        return _cmd_validate(args)
    except (ConfigError, MetadataError, UnsupportedPlatformError,
            FileNotFoundError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
