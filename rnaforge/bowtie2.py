"""bowtie2 hizalama: çalıştırır ve özetini parse eder. Parser saftır."""
from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

_RATE_RE = re.compile(r"([0-9]+(?:\.[0-9]+)?)%\s+overall alignment rate")


class Bowtie2ParseError(ValueError):
    """bowtie2 özeti beklenen biçimde değil."""


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
