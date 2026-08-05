"""featureCounts çıktısını parse eder ve çalıştırır. Parserlar saftır."""
from __future__ import annotations

import subprocess
from dataclasses import dataclass
from pathlib import Path


class FeatureCountsParseError(ValueError):
    """featureCounts çıktısı beklenen biçimde değil."""


class FeatureCountsRunError(RuntimeError):
    """featureCounts çalıştırılamadı ya da beklenen çıktıyı üretmedi."""


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


def parse_lengths(counts_text: str) -> dict[str, int]:
    """featureCounts çıktısından gen -> Length (bç). Length 6. sütun (index 5)."""
    lengths: dict[str, int] = {}
    header_seen = False
    for line in counts_text.splitlines():
        if line.startswith("#") or not line.strip():
            continue
        fields = line.split("\t")
        if not header_seen:
            if fields[0] != "Geneid":
                raise FeatureCountsParseError("featureCounts counts file has no 'Geneid' header")
            header_seen = True
            continue
        if len(fields) > 5:
            lengths[fields[0]] = int(fields[5])
    return lengths


def compute_tpm_fpkm(counts_text: str):
    """featureCounts ham çıktısından TPM ve FPKM matrisleri (gen uzunluğuyla normalize).
    Returns: (gene_ids, columns, tpm{col:[...]}, fpkm{col:[...]})."""
    gene_ids, counts = parse_counts(counts_text)
    lengths = parse_lengths(counts_text)
    kb = [max(lengths.get(g, 0), 1) / 1000.0 for g in gene_ids]     # gen uzunluğu (kb), 0-koruması
    columns = list(counts)
    tpm: dict[str, list[float]] = {}
    fpkm: dict[str, list[float]] = {}
    for col in columns:
        c = counts[col]
        total = sum(c)                                              # kütüphane büyüklüğü (atanmış okuma)
        rpk = [c[i] / kb[i] for i in range(len(gene_ids))]          # reads per kilobase
        scale = sum(rpk) / 1e6
        tpm[col] = [round(r / scale, 4) if scale > 0 else 0.0 for r in rpk]
        fpkm[col] = [round(c[i] / (kb[i] * (total / 1e6)), 4) if total > 0 else 0.0
                     for i in range(len(gene_ids))]
    return gene_ids, columns, tpm, fpkm


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


def run_featurecounts(bams: list[Path], gff: Path, out_dir: Path, feature_type: str,
                      attribute: str, paired: bool = False, threads: int = 4,
                      env: str = "rnaforge-quant-prok") -> FeatureCountsResult:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    counts_path = out_dir / "counts.txt"
    cmd = ["conda", "run", "-n", env, "featureCounts",
           "-a", str(gff), "-o", str(counts_path),
           "-t", feature_type, "-g", attribute, "-T", str(threads)]
    if paired:
        cmd += ["-p", "--countReadPairs"]
    cmd += [str(b) for b in bams]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise FeatureCountsRunError(
            f"featureCounts failed (exit {r.returncode})\ncmd: {' '.join(cmd)}\n"
            f"stderr: {r.stderr.strip()}"
        )
    summary_path = counts_path.with_name(counts_path.name + ".summary")
    if not counts_path.exists() or not summary_path.exists():
        raise FeatureCountsRunError(
            f"featureCounts reported success but output missing at {counts_path}"
        )
    gene_ids, counts = parse_counts(counts_path.read_text())
    rates = parse_summary(summary_path.read_text())
    return FeatureCountsResult(gene_ids=gene_ids, counts=counts, assignment_rates=rates)
