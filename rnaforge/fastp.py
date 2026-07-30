"""fastp çıktısını parse eder ve çalıştırır. Parser saftır: string girer,
FastpResult çıkar — I/O yok, hızlı ve deterministik test edilir."""
from __future__ import annotations

import json
import subprocess
from dataclasses import dataclass
from pathlib import Path


class FastpParseError(ValueError):
    """fastp JSON çıktısı beklenen biçimde değil."""


class FastpRunError(RuntimeError):
    """fastp çalıştırılamadı ya da beklenen çıktıyı üretmedi."""


@dataclass(frozen=True)
class FastpResult:
    reads_before: int
    reads_after: int
    survival_rate: float
    out1: Path | None = None
    out2: Path | None = None


def parse_fastp_json(json_text: str, out1: Path | None = None,
                     out2: Path | None = None) -> FastpResult:
    try:
        data = json.loads(json_text)
    except json.JSONDecodeError as exc:
        raise FastpParseError(f"fastp JSON is not valid JSON: {exc}") from None
    summary = data.get("summary")
    if not isinstance(summary, dict):
        raise FastpParseError("fastp JSON has no 'summary' object — output is malformed")
    try:
        before = int(summary["before_filtering"]["total_reads"])
        after = int(summary["after_filtering"]["total_reads"])
    except (KeyError, TypeError, ValueError) as exc:
        raise FastpParseError(
            f"fastp JSON summary missing before/after total_reads: {exc}"
        ) from None
    survival = (after / before) if before > 0 else 0.0
    return FastpResult(reads_before=before, reads_after=after,
                       survival_rate=survival, out1=out1, out2=out2)


def trimmed_name(fastq: Path) -> str:
    stem = fastq.name
    for ext in (".gz", ".fastq", ".fq"):
        if stem.endswith(ext):
            stem = stem[: -len(ext)]
    return f"{stem}.trimmed.fastq"


def run_fastp(fastq_1: Path, out_dir: Path, min_length: int,
              fastq_2: Path | None = None, aggressive_quality: bool = False,
              env: str = "rnaforge-qc") -> FastpResult:
    fastq_1 = Path(fastq_1)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out1 = out_dir / trimmed_name(fastq_1)
    json_path = out_dir / "fastp.json"
    html_path = out_dir / "fastp.html"
    cmd = ["conda", "run", "-n", env, "fastp",
           "-i", str(fastq_1), "-o", str(out1),
           "-l", str(min_length),
           "-j", str(json_path), "-h", str(html_path)]
    out2 = None
    if fastq_2 is not None:
        fastq_2 = Path(fastq_2)
        out2 = out_dir / trimmed_name(fastq_2)
        cmd += ["-I", str(fastq_2), "-O", str(out2)]
    if aggressive_quality:
        cmd += ["-r"]  # BİLİNÇLİ agresif (sliding-window); nazik varsayılandan sapma
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise FastpRunError(
            f"fastp failed (exit {result.returncode}) for {fastq_1}\n"
            f"cmd: {' '.join(cmd)}\nstderr: {result.stderr.strip()}"
        )
    if not json_path.exists():
        raise FastpRunError(
            f"fastp reported success but produced no JSON at {json_path}"
        )
    return parse_fastp_json(json_path.read_text(), out1=out1, out2=out2)
