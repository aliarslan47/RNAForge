"""Komut satırı girişi."""
from __future__ import annotations

import argparse

from rnaforge import __version__


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="rnaforge", description="Bulk RNA-seq pipeline")
    parser.add_argument("--version", action="version", version=__version__)
    parser.add_argument("command", nargs="?", help="validate")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command is None:
        print("error: no command given (try: rnaforge validate)")
        return 2
    return 0
