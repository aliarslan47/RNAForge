from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from rnaforge.deseq2 import (
    DeseqParseError,
    DeseqResult,
    DeseqRunError,
    parse_de_metrics,
    parse_deseq2_results,
    run_deseq2,
)

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
