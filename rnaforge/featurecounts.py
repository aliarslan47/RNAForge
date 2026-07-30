"""featureCounts çıktısını parse eder ve çalıştırır. Parserlar saftır."""
from __future__ import annotations

from dataclasses import dataclass


class FeatureCountsParseError(ValueError):
    """featureCounts çıktısı beklenen biçimde değil."""


@dataclass(frozen=True)
class FeatureCountsResult:
    gene_ids: list[str]
    counts: dict[str, list[int]]           # sütun (BAM) -> sayımlar
    assignment_rates: dict[str, float]     # sütun (BAM) -> atama oranı


def parse_counts(counts_text: str) -> tuple[list[str], dict[str, list[int]]]:
    header = None
    gene_ids: list[str] = []
    columns: list[str] = []
    counts: dict[str, list[int]] = {}
    for line in counts_text.splitlines():
        if line.startswith("#") or not line.strip():
            continue
        fields = line.split("\t")
        if header is None:
            if fields[0] != "Geneid":
                raise FeatureCountsParseError(
                    f"featureCounts counts file has no 'Geneid' header (got {fields[0]!r})"
                )
            header = fields
            columns = fields[6:]              # Geneid Chr Start End Strand Length <bam...>
            counts = {c: [] for c in columns}
            continue
        gene_ids.append(fields[0])
        for col, value in zip(columns, fields[6:]):
            counts[col].append(int(value))
    if header is None:
        raise FeatureCountsParseError("featureCounts counts file has no 'Geneid' header line")
    return gene_ids, counts


def parse_summary(summary_text: str) -> dict[str, float]:
    lines = [ln for ln in summary_text.splitlines() if ln.strip()]
    if not lines or not lines[0].startswith("Status"):
        raise FeatureCountsParseError("featureCounts summary has no 'Status' header")
    columns = lines[0].split("\t")[1:]
    assigned = {c: 0 for c in columns}
    totals = {c: 0 for c in columns}
    for line in lines[1:]:
        fields = line.split("\t")
        status = fields[0]
        for col, value in zip(columns, fields[1:]):
            v = int(value)
            totals[col] += v
            if status == "Assigned":
                assigned[col] += v
    return {c: (assigned[c] / totals[c] if totals[c] > 0 else 0.0) for c in columns}
