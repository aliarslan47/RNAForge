"""DESeq2 (R/Bioconductor) diferansiyel ekspresyon: çalıştırır ve çıktısını parse eder.
Parserlar saftır: string girer, veri çıkar (I/O yok)."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

_NUMERIC = ("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")


class DeseqParseError(ValueError):
    """DESeq2 çıktısı beklenen biçimde değil."""


@dataclass(frozen=True)
class DeseqResult:
    results: list[dict]
    metrics: dict
    results_path: Path | None = None
    normalized_path: Path | None = None


def _num(value: str):
    if value == "NA" or value == "":
        return None
    return float(value)


def parse_deseq2_results(results_text: str) -> list[dict]:
    lines = [ln for ln in results_text.splitlines() if ln.strip()]
    if not lines or lines[0].split("\t")[0] != "gene":
        raise DeseqParseError("DESeq2 results file has no 'gene' header column")
    header = lines[0].split("\t")
    rows: list[dict] = []
    for line in lines[1:]:
        fields = line.split("\t")
        row: dict = {}
        for key, val in zip(header, fields):
            row[key] = val if key == "gene" else _num(val)
        rows.append(row)
    return rows


def parse_de_metrics(metrics_text: str) -> dict:
    out: dict = {}
    for line in metrics_text.splitlines():
        if not line.strip():
            continue
        key, _, value = line.partition("\t")
        try:
            out[key] = float(value)
        except ValueError:
            out[key] = value
    return out
