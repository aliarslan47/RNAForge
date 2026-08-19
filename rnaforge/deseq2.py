"""DESeq2 (R/Bioconductor) diferansiyel ekspresyon: çalıştırır ve çıktısını parse eder.
Parserlar saftır: string girer, veri çıkar (I/O yok)."""
from __future__ import annotations

import subprocess
from dataclasses import dataclass, field
from pathlib import Path

_NUMERIC = ("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")
_SCRIPT = Path(__file__).parent / "scripts" / "deseq2.R"


class DeseqParseError(ValueError):
    """DESeq2 çıktısı beklenen biçimde değil."""


class DeseqRunError(RuntimeError):
    """DESeq2 (Rscript) çalıştırılamadı ya da beklenen çıktıyı üretmedi."""


@dataclass(frozen=True)
class DeseqResult:
    results: list[dict]
    metrics: dict
    results_path: Path | None = None
    normalized_path: Path | None = None
    # Açık kontrast istendiğinde her biri için ayrı sonuç dosyası: "<test>_vs_<ref>" -> yol.
    contrast_paths: dict = field(default_factory=dict)


def format_contrasts(contrasts) -> str:
    """((test, ref), ...) -> 'test:ref;test2:ref2' (deseq2.R'nin 6. argümanı). Boş -> ''."""
    return ";".join(f"{t}:{r}" for t, r in (contrasts or ()))


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


def run_deseq2(counts_tsv: Path, coldata_tsv: Path, design: str, out_dir: Path,
               reference: str | None = None, contrasts=None,
               env: str = "rnaforge-de") -> DeseqResult:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    contrasts = tuple(contrasts or ())
    cmd = ["conda", "run", "-n", env, "Rscript", str(_SCRIPT),
           str(counts_tsv), str(coldata_tsv), design, reference or "", str(out_dir),
           format_contrasts(contrasts)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise DeseqRunError(
            f"DESeq2 failed (exit {r.returncode})\ncmd: {' '.join(cmd)}\n"
            f"stderr: {r.stderr.strip()}"
        )
    results_path = out_dir / "deseq2_results.tsv"
    metrics_path = out_dir / "de_metrics.tsv"
    normalized_path = out_dir / "normalized_counts.tsv"
    if not results_path.exists() or not metrics_path.exists():
        raise DeseqRunError(
            f"DESeq2 reported success but output missing in {out_dir}\nstderr: {r.stderr.strip()}"
        )
    contrast_paths: dict = {}
    for test, ref in contrasts:
        path = out_dir / f"deseq2_results.{test}_vs_{ref}.tsv"
        if not path.exists():
            raise DeseqRunError(
                f"DESeq2 reported success but contrast output missing: {path}\n"
                f"stderr: {r.stderr.strip()}"
            )
        contrast_paths[f"{test}_vs_{ref}"] = path
    results = parse_deseq2_results(results_path.read_text())
    metrics = parse_de_metrics(metrics_path.read_text())
    return DeseqResult(results=results, metrics=metrics,
                       results_path=results_path, normalized_path=normalized_path,
                       contrast_paths=contrast_paths)
