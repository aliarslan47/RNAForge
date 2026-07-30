"""fastp çıktısını parse eder ve çalıştırır. Parser saftır: string girer,
FastpResult çıkar — I/O yok, hızlı ve deterministik test edilir."""
from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path


class FastpParseError(ValueError):
    """fastp JSON çıktısı beklenen biçimde değil."""


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
