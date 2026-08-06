"""NanoPlot (ONT/long-read QC) output parsing + runner.

Mirrors rnaforge/fastqc.py: a pure parser over NanoPlot's --tsv_stats
NanoStats.txt, plus a real runner in the rnaforge-longread env. Long-read QC
is diagnostic (like FastQC) — it never stops the run."""
from __future__ import annotations

import re
import subprocess
from dataclasses import dataclass
from pathlib import Path


class NanoStatsParseError(ValueError):
    """NanoStats.txt could not be parsed."""


class NanoPlotRunError(RuntimeError):
    """NanoPlot failed to run."""


@dataclass(frozen=True)
class NanoStats:
    number_of_reads: int
    number_of_bases: int
    mean_read_length: float
    median_read_length: float
    read_length_stdev: float
    n50: float
    mean_qual: float
    median_qual: float
    reads_above_q10_pct: float | None = None


_CORE = {
    "number_of_reads": int,
    "number_of_bases": lambda v: int(float(v)),
    "mean_read_length": float,
    "median_read_length": float,
    "read_length_stdev": float,
    "n50": float,
    "mean_qual": float,
    "median_qual": float,
}


def parse_nanostats(text: str) -> NanoStats:
    raw: dict[str, str] = {}
    for line in text.splitlines():
        if "\t" not in line:
            continue
        key, value = line.split("\t", 1)
        raw[key.strip()] = value.strip()

    parsed: dict[str, object] = {}
    for key, caster in _CORE.items():
        if key not in raw:
            raise NanoStatsParseError(f"NanoStats.txt missing required metric: {key!r}")
        try:
            parsed[key] = caster(raw[key])
        except (TypeError, ValueError) as exc:
            raise NanoStatsParseError(f"bad value for {key!r}: {raw[key]!r}") from exc

    q10 = None
    m = re.search(r"\(([\d.]+)%\)", raw.get("Reads >Q10:", ""))
    if m:
        q10 = float(m.group(1))

    return NanoStats(reads_above_q10_pct=q10, **parsed)  # type: ignore[arg-type]


def run_nanoplot(fastq: Path, out_dir: Path, env: str = "rnaforge-longread") -> Path:
    """Run NanoPlot on a single long-read FASTQ; return the NanoStats.txt path.

    --tsv_stats gives the parseable stats file; --no_static skips the
    kaleido/orca static-image dependency (HTML plots are still written)."""
    fastq = Path(fastq)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        "conda", "run", "-n", env, "NanoPlot",
        "--fastq", str(fastq),
        "--outdir", str(out_dir),
        "--tsv_stats", "--no_static",
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    stats_path = out_dir / "NanoStats.txt"
    if proc.returncode != 0 or not stats_path.exists():
        raise NanoPlotRunError(
            f"NanoPlot failed (exit {proc.returncode}) on {fastq}\n"
            f"stdout: {proc.stdout[-500:]}\nstderr: {proc.stderr[-500:]}"
        )
    return stats_path
