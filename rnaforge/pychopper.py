"""Pychopper (ONT full-length cDNA orient/trim) runner + stats parser.

Pychopper 2.7.10 writes the oriented FASTQ and the -S stats TSV, then crashes
in its end-of-run PDF report under pandas 3 (_plot_stats -> float(Series)).
The real work is complete before that crash, so run_pychopper tolerates exactly
that signature and raises on anything else (feedback_gurultulu_hata)."""
from __future__ import annotations

import gzip
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path


class PychopperParseError(ValueError):
    """pychopper stats TSV could not be parsed."""


class PychopperRunError(RuntimeError):
    """pychopper failed for a reason other than the known plotting crash."""


@dataclass(frozen=True)
class PychopperStats:
    pass_reads: int
    primers_found: int
    rescue: int
    unusable: int
    len_fail: int


def parse_pychopper_stats(tsv_text: str) -> PychopperStats:
    table: dict[tuple[str, str], float] = {}
    for line in tsv_text.splitlines():
        parts = line.split("\t")
        if len(parts) != 3 or parts[0] == "Category":
            continue
        try:
            table[(parts[0], parts[1])] = float(parts[2])
        except ValueError:
            continue

    def get(cat: str, name: str) -> int:
        if (cat, name) not in table:
            raise PychopperParseError(f"pychopper stats missing {cat}/{name}")
        return int(table[(cat, name)])

    return PychopperStats(
        pass_reads=get("ReadStats", "PassReads"),
        primers_found=get("Classification", "Primers_found"),
        rescue=get("Classification", "Rescue"),
        unusable=get("Classification", "Unusable"),
        len_fail=get("ReadStats", "LenFail"),
    )


def run_pychopper(in_fastq: Path, out_fastq: Path, stats_tsv: Path,
                  env: str = "rnaforge-longread", kit: str | None = None,
                  threads: int = 4) -> PychopperStats:
    in_fastq = Path(in_fastq)
    out_fastq = Path(out_fastq)
    stats_tsv = Path(stats_tsv)
    out_fastq.parent.mkdir(parents=True, exist_ok=True)

    tmp = None
    src = in_fastq
    if in_fastq.name.endswith(".gz"):
        tmp = Path(tempfile.mkstemp(suffix=".fastq")[1])
        with gzip.open(in_fastq, "rt") as fh, tmp.open("w") as out:
            shutil.copyfileobj(fh, out)
        src = tmp
    try:
        cmd = ["conda", "run", "-n", env, "pychopper",
               "-t", str(threads), "-S", str(stats_tsv)]
        if kit:
            cmd += ["-k", kit]
        cmd += [str(src), str(out_fastq)]
        proc = subprocess.run(cmd, capture_output=True, text=True)
    finally:
        if tmp is not None:
            tmp.unlink(missing_ok=True)

    outputs_ok = out_fastq.exists() and stats_tsv.exists()
    known_plot_crash = "_plot_stats" in proc.stderr
    if proc.returncode != 0 and not (outputs_ok and known_plot_crash):
        raise PychopperRunError(
            f"pychopper failed (exit {proc.returncode}) on {in_fastq}\n"
            f"stderr: {proc.stderr[-800:]}"
        )
    if proc.returncode != 0:
        sys.stderr.write(
            f"WARNING: pychopper exited {proc.returncode} on its PDF report "
            f"(known pandas-3 _plot_stats bug); outputs are complete, continuing.\n"
        )
    return parse_pychopper_stats(stats_tsv.read_text())
