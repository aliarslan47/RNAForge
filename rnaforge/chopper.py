"""chopper (ONT length/quality filter) runner. Reads stdin, writes stdout."""
from __future__ import annotations

import re
import subprocess
from pathlib import Path


class ChopperRunError(RuntimeError):
    """chopper failed to run."""


def parse_kept(stderr: str) -> int | None:
    m = re.search(r"Kept (\d+) reads out of", stderr)
    return int(m.group(1)) if m else None


def run_chopper(in_fastq: Path, out_fastq: Path, env: str = "rnaforge-longread",
                min_qual: int = 7, min_len: int = 50, threads: int = 4) -> int:
    in_fastq = Path(in_fastq)
    out_fastq = Path(out_fastq)
    out_fastq.parent.mkdir(parents=True, exist_ok=True)
    decomp = "zcat" if in_fastq.name.endswith(".gz") else "cat"
    pipe = (
        f"{decomp} {in_fastq!s} | "
        f"chopper -q {min_qual} -l {min_len} --threads {threads} > {out_fastq!s}"
    )
    cmd = ["conda", "run", "-n", env, "bash", "-lc", pipe]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0 or not out_fastq.exists():
        raise ChopperRunError(
            f"chopper failed (exit {proc.returncode}) on {in_fastq}\n"
            f"stderr: {proc.stderr[-500:]}"
        )
    kept = parse_kept(proc.stderr)
    if kept is None:
        kept = sum(1 for _ in out_fastq.open()) // 4
    return kept
