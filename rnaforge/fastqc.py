"""FastQC çıktısını parse eder ve çalıştırır. Parser saftır: string girer,
FastQCReport çıkar — I/O yok, hızlı ve deterministik test edilir."""
from __future__ import annotations

from dataclasses import dataclass

_STATUSES = ("PASS", "WARN", "FAIL")


class FastQCParseError(ValueError):
    """FastQC çıktısı beklenen biçimde değil."""


@dataclass(frozen=True)
class FastQCReport:
    modules: dict[str, str]
    basic_stats: dict[str, str]


def parse_fastqc_report(summary_text: str, data_text: str) -> FastQCReport:
    modules: dict[str, str] = {}
    for line in summary_text.splitlines():
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) < 2:
            raise FastQCParseError(f"summary line is not tab-delimited: {line!r}")
        status, name = parts[0].strip(), parts[1].strip()
        if status not in _STATUSES:
            raise FastQCParseError(
                f"unknown FastQC status {status!r} for module {name!r} "
                f"(expected one of {_STATUSES})"
            )
        modules[name] = status

    basic_stats = _parse_basic_stats(data_text)
    return FastQCReport(modules=modules, basic_stats=basic_stats)


def _parse_basic_stats(data_text: str) -> dict[str, str]:
    stats: dict[str, str] = {}
    in_basic = False
    for line in data_text.splitlines():
        if line.startswith(">>Basic Statistics"):
            in_basic = True
            continue
        if in_basic:
            if line.startswith(">>END_MODULE"):
                return stats
            if line.startswith("#") or not line.strip():
                continue
            key, _, value = line.partition("\t")
            stats[key.strip()] = value.strip()
    raise FastQCParseError(
        "FastQC data has no '>>Basic Statistics' module — output is malformed"
    )
