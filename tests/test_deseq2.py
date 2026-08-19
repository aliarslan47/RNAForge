from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from rnaforge.deseq2 import (
    DeseqParseError,
    DeseqResult,
    DeseqRunError,
    format_contrasts,
    parse_de_metrics,
    parse_deseq2_results,
    run_deseq2,
)


def test_format_contrasts_serializes_pairs():
    assert format_contrasts((("high", "control"), ("low", "control"))) == "high:control;low:control"


def test_format_contrasts_empty():
    assert format_contrasts(()) == ""
    assert format_contrasts(None) == ""


def test_run_deseq2_passes_contrasts_and_collects_paths(tmp_path, monkeypatch):
    """Wrapper contrasts'ı R'a 6. arg olarak geçmeli ve her contrast çıktı dosyasını
    toplamalı. subprocess.run mock'lanır (conda gerekmez)."""
    import subprocess as sp
    captured = {}

    def fake_run(cmd, capture_output, text):
        captured["cmd"] = cmd
        out_dir = Path(cmd[-2])          # ... out_dir, contrasts_spec
        for name in ("deseq2_results.tsv", "deseq2_results.high_vs_control.tsv",
                     "deseq2_results.low_vs_control.tsv"):
            (out_dir / name).write_text("gene\tbaseMean\tlog2FoldChange\tlfcSE\tstat\tpvalue\tpadj\n"
                                        "g0\t10\t1\t0.1\t2\t0.01\t0.02\n")
        (out_dir / "de_metrics.tsv").write_text("min_replicate_correlation\t0.9\ncontrast\thigh vs control\nn_genes\t1\n")
        (out_dir / "normalized_counts.tsv").write_text("gene\ts1\ng0\t10\n")

        class R: returncode = 0; stderr = ""
        return R()

    monkeypatch.setattr(sp, "run", fake_run)
    contrasts = (("high", "control"), ("low", "control"))
    result = run_deseq2(tmp_path / "counts.tsv", tmp_path / "coldata.tsv", "~condition",
                        tmp_path / "de", reference="control", contrasts=contrasts)
    assert captured["cmd"][-1] == "high:control;low:control"
    assert set(result.contrast_paths) == {"high_vs_control", "low_vs_control"}
    assert result.contrast_paths["high_vs_control"].exists()

_RESULTS = (
    "gene\tbaseMean\tlog2FoldChange\tlfcSE\tstat\tpvalue\tpadj\n"
    "geneA\t120.5\t2.3\t0.4\t5.7\t1e-8\t2e-7\n"
    "geneB\t80.0\t-0.1\t0.3\t-0.3\t0.7\tNA\n"
)

_METRICS = "min_replicate_correlation\t0.93\ncontrast\ttreated vs control\nn_genes\t2\n"


def test_parse_results_reads_rows_and_types():
    rows = parse_deseq2_results(_RESULTS)
    assert len(rows) == 2
    assert rows[0]["gene"] == "geneA"
    assert rows[0]["log2FoldChange"] == pytest.approx(2.3)
    assert rows[0]["padj"] == pytest.approx(2e-7)
    assert rows[1]["padj"] is None            # NA -> None


def test_parse_results_rejects_missing_header():
    with pytest.raises(DeseqParseError, match="gene"):
        parse_deseq2_results("foo\tbar\n1\t2\n")


def test_parse_metrics_typed():
    m = parse_de_metrics(_METRICS)
    assert m["min_replicate_correlation"] == pytest.approx(0.93)
    assert m["contrast"] == "treated vs control"
    assert m["n_genes"] == pytest.approx(2.0)


def _write_counts_coldata(tmp_path):
    """2 koşul × 3 replika; birkaç gende açık sinyal."""
    import random
    random.seed(4)
    samples = ["c1", "c2", "c3", "t1", "t2", "t3"]
    genes = [f"g{i}" for i in range(50)]
    counts = tmp_path / "counts.tsv"
    with counts.open("w") as f:
        f.write("gene\t" + "\t".join(samples) + "\n")
        for gi, g in enumerate(genes):
            base = random.randint(100, 300)
            row = []
            for s in samples:
                val = base + random.randint(-15, 15)
                if gi < 5 and s.startswith("t"):
                    val *= 4                          # ilk 5 gen treated'da yukari
                row.append(str(val))
            f.write(g + "\t" + "\t".join(row) + "\n")
    coldata = tmp_path / "coldata.tsv"
    with coldata.open("w") as f:
        f.write("sample\tcondition\n")
        for s in samples:
            f.write(f"{s}\t{'control' if s.startswith('c') else 'treated'}\n")
    return counts, coldata


def _write_paired_counts_coldata(tmp_path):
    """3 denek × 2 koşul (eşleştirilmiş); denek-bazlı baseline kayması + koşul sinyali."""
    import random
    random.seed(7)
    subjects = ["p1", "p2", "p3"]
    samples = [(f"{s}_{c}", c, s) for s in subjects for c in ("control", "treated")]
    genes = [f"g{i}" for i in range(40)]
    subj_offset = {"p1": 0, "p2": 120, "p3": 240}     # denekler arası baseline farkı
    counts = tmp_path / "counts.tsv"
    with counts.open("w") as f:
        f.write("gene\t" + "\t".join(sid for sid, _, _ in samples) + "\n")
        for gi, g in enumerate(genes):
            base = random.randint(100, 200)
            row = []
            for _sid, cond, subj in samples:
                val = base + subj_offset[subj] + random.randint(-10, 10)
                if gi < 5 and cond == "treated":
                    val += 400                          # ilk 5 gen treated'da yukarı
                row.append(str(val))
            f.write(g + "\t" + "\t".join(row) + "\n")
    coldata = tmp_path / "coldata.tsv"
    with coldata.open("w") as f:
        f.write("sample\tcondition\tsubject\n")
        for sid, cond, subj in samples:
            f.write(f"{sid}\t{cond}\t{subj}\n")
    return counts, coldata


@pytest.mark.skipif(shutil.which("conda") is None, reason="conda yok")
def test_run_deseq2_paired_subject_design(tmp_path):
    """Entegrasyon: '~subject + condition' coldata'da subject olduğunda çökmeden koşar
    ve koşul sinyalini bulur. Regresyon: eskiden coldata subject yazmıyordu → R'da
    'variable subject not found' ile patlıyordu. rnaforge-de yoksa skip."""
    counts, coldata = _write_paired_counts_coldata(tmp_path)
    try:
        result = run_deseq2(counts, coldata, "~subject + condition", tmp_path / "de",
                            reference="control")
    except DeseqRunError as exc:
        pytest.skip(f"DESeq2 çalıştırılamadı (env yok?): {exc}")
    by_gene = {r["gene"]: r for r in result.results}
    assert by_gene["g0"]["padj"] is not None and by_gene["g0"]["padj"] < 0.05
    assert by_gene["g0"]["log2FoldChange"] > 1


def _write_three_level(tmp_path):
    """3 seviye (control/low/high) × 3 replika; g0 high'da güçlü, low'da orta yukarı."""
    import random
    random.seed(11)
    conds = ["control", "low", "high"]
    samples = [(f"{c}{i}", c) for c in conds for i in (1, 2, 3)]
    mult = {"control": 1, "low": 2, "high": 5}
    counts = tmp_path / "counts.tsv"
    with counts.open("w") as f:
        f.write("gene\t" + "\t".join(sid for sid, _ in samples) + "\n")
        for gi in range(40):
            base = random.randint(100, 200)
            row = []
            for _sid, cond in samples:
                val = base + random.randint(-10, 10)
                if gi < 5:
                    val *= mult[cond]
                row.append(str(val))
            f.write(f"g{gi}\t" + "\t".join(row) + "\n")
    coldata = tmp_path / "coldata.tsv"
    with coldata.open("w") as f:
        f.write("sample\tcondition\n")
        for sid, cond in samples:
            f.write(f"{sid}\t{cond}\n")
    return counts, coldata


@pytest.mark.skipif(shutil.which("conda") is None, reason="conda yok")
def test_run_deseq2_multiple_contrasts(tmp_path):
    """Entegrasyon: 3-seviyeli condition'da açık kontrastlar ayrı dosyalar üretir ve
    her biri doğru karşılaştırmayı raporlar (high-vs-control > low-vs-control g0'da)."""
    counts, coldata = _write_three_level(tmp_path)
    contrasts = (("high", "control"), ("low", "control"))
    try:
        result = run_deseq2(counts, coldata, "~condition", tmp_path / "de",
                            reference="control", contrasts=contrasts)
    except DeseqRunError as exc:
        pytest.skip(f"DESeq2 çalıştırılamadı (env yok?): {exc}")
    assert set(result.contrast_paths) == {"high_vs_control", "low_vs_control"}
    high = {r["gene"]: r for r in parse_deseq2_results(
        result.contrast_paths["high_vs_control"].read_text())}
    low = {r["gene"]: r for r in parse_deseq2_results(
        result.contrast_paths["low_vs_control"].read_text())}
    assert high["g0"]["log2FoldChange"] > low["g0"]["log2FoldChange"] > 0
    assert high["g0"]["padj"] < 0.05
    # Birincil deseq2_results.tsv = ilk kontrast (high vs control)
    assert result.metrics["contrast"] == "high vs control"


@pytest.mark.skipif(shutil.which("conda") is None, reason="conda yok")
def test_run_deseq2_detects_signal(tmp_path):
    """Entegrasyon: gerçek DESeq2 sinyalli genleri anlamlı bulur. rnaforge-de yoksa skip."""
    counts, coldata = _write_counts_coldata(tmp_path)
    try:
        result = run_deseq2(counts, coldata, "~condition", tmp_path / "de", reference="control")
    except DeseqRunError as exc:
        pytest.skip(f"DESeq2 çalıştırılamadı (env yok?): {exc}")
    assert len(result.results) == 50
    by_gene = {r["gene"]: r for r in result.results}
    assert by_gene["g0"]["padj"] is not None and by_gene["g0"]["padj"] < 0.05
    assert by_gene["g0"]["log2FoldChange"] > 1
    assert result.metrics["contrast"] == "treated vs control"
    assert result.metrics["min_replicate_correlation"] > 0.5
    disp = tmp_path / "de" / "dispersions.tsv"
    assert disp.exists()
    header = disp.read_text().splitlines()[0].split("\t")
    assert header == ["gene_id", "baseMean", "dispGeneEst", "dispFit", "dispFinal"]
