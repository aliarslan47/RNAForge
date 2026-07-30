"""bowtie2 hizalama: çalıştırır ve özetini parse eder. Parser saftır."""
from __future__ import annotations

import re
import subprocess
from dataclasses import dataclass
from pathlib import Path

_RATE_RE = re.compile(r"([0-9]+(?:\.[0-9]+)?)%\s+overall alignment rate")


class Bowtie2ParseError(ValueError):
    """bowtie2 özeti beklenen biçimde değil."""


class Bowtie2RunError(RuntimeError):
    """bowtie2/samtools çalıştırılamadı ya da beklenen çıktıyı üretmedi."""


@dataclass(frozen=True)
class AlignmentResult:
    bam: Path
    alignment_rate: float


def parse_bowtie2_summary(stderr_text: str) -> float:
    match = None
    for m in _RATE_RE.finditer(stderr_text):
        match = m
    if match is None:
        raise Bowtie2ParseError(
            "bowtie2 stderr has no 'overall alignment rate' line — run may have failed"
        )
    return float(match.group(1)) / 100.0


def _run(cmd: list[str]) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, capture_output=True, text=True)


def build_index(genome_fasta: Path, index_dir: Path,
                env: str = "rnaforge-quant-prok") -> Path:
    genome_fasta = Path(genome_fasta)
    if not genome_fasta.exists():
        raise Bowtie2RunError(f"genome FASTA does not exist: {genome_fasta}")
    index_dir = Path(index_dir)
    index_dir.mkdir(parents=True, exist_ok=True)
    prefix = index_dir / "genome"
    r = _run(["conda", "run", "-n", env, "bowtie2-build", "-q",
              str(genome_fasta), str(prefix)])
    if r.returncode != 0:
        raise Bowtie2RunError(
            f"bowtie2-build failed (exit {r.returncode}) for {genome_fasta}\n"
            f"stderr: {r.stderr.strip()}"
        )
    return prefix


def run_bowtie2(index_prefix: Path, out_dir: Path, fastq_1: Path,
                fastq_2: Path | None = None, threads: int = 4,
                env: str = "rnaforge-quant-prok") -> AlignmentResult:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    sam = out_dir / "aligned.sam"
    bam = out_dir / "aligned.sorted.bam"
    log = out_dir / "bowtie2.log"

    bt = ["conda", "run", "-n", env, "bowtie2", "-x", str(index_prefix),
          "-p", str(threads), "-S", str(sam)]
    if fastq_2 is not None:
        bt += ["-1", str(fastq_1), "-2", str(fastq_2)]
    else:
        bt += ["-U", str(fastq_1)]
    r = _run(bt)
    log.write_text(r.stderr)
    if r.returncode != 0:
        raise Bowtie2RunError(
            f"bowtie2 failed (exit {r.returncode})\ncmd: {' '.join(bt)}\n"
            f"stderr: {r.stderr.strip()}"
        )
    rate = parse_bowtie2_summary(r.stderr)

    sort = _run(["conda", "run", "-n", env, "samtools", "sort", "-o", str(bam), str(sam)])
    if sort.returncode != 0:
        raise Bowtie2RunError(f"samtools sort failed: {sort.stderr.strip()}")
    idx = _run(["conda", "run", "-n", env, "samtools", "index", str(bam)])
    if idx.returncode != 0:
        raise Bowtie2RunError(f"samtools index failed: {idx.stderr.strip()}")
    sam.unlink(missing_ok=True)
    return AlignmentResult(bam=bam, alignment_rate=rate)
