from __future__ import annotations

import pytest

from rnaforge.deseq2 import DeseqParseError, DeseqResult, parse_de_metrics, parse_deseq2_results

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
