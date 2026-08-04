import json
from pathlib import Path
import pytest
from rnaforge.report_html import (
    parse_de_results, load_gene_map, top_degs, embed_png, load_report_inputs,
)


def test_parse_de_results_types(tmp_path):
    p = tmp_path / "de.tsv"
    p.write_text("gene\tbaseMean\tlog2FoldChange\tlfcSE\tstat\tpvalue\tpadj\n"
                 "LT_1\t200\t3.0\t0.2\t5\t1e-9\t1e-8\n"
                 "LT_2\t50\t-2.0\t0.3\t-4\tNA\tNA\n")
    rows = parse_de_results(p)
    assert rows[0]["gene"] == "LT_1" and rows[0]["padj"] == 1e-8
    assert rows[1]["padj"] is None and rows[1]["log2FoldChange"] == -2.0


def test_load_gene_map_missing_is_empty(tmp_path):
    assert load_gene_map(tmp_path / "nope.tsv") == {}


def test_top_degs_orders_and_labels(tmp_path):
    de = [
        {"gene": "LT_1", "baseMean": 200.0, "log2FoldChange": 3.0, "padj": 1e-8},
        {"gene": "LT_2", "baseMean": 80.0, "log2FoldChange": -2.5, "padj": 1e-4},
        {"gene": "LT_3", "baseMean": 10.0, "log2FoldChange": 0.2, "padj": 1e-9},   # |lfc|<1 -> not sig
        {"gene": "LT_4", "baseMean": 10.0, "log2FoldChange": 4.0, "padj": None},   # NA -> not sig
    ]
    gm = {"LT_1": "pspA"}
    out = top_degs(de, gm, fdr=0.05, lfc=1.0, n=50)
    assert [r["gene"] for r in out] == ["pspA", "LT_2"]         # padj asc, sig only, named
    assert out[0]["direction"] == "Up" and out[1]["direction"] == "Down"


def test_top_degs_empty_when_none_significant():
    de = [{"gene": "LT_1", "baseMean": 10.0, "log2FoldChange": 0.1, "padj": 0.9}]
    assert top_degs(de, {}, fdr=0.05, lfc=1.0) == []


def test_embed_png_data_uri(tmp_path):
    p = tmp_path / "x.png"; p.write_bytes(b"\x89PNG\r\n")
    assert embed_png(p).startswith("data:image/png;base64,")


def test_load_report_inputs_loud_on_missing(tmp_path):
    with pytest.raises(FileNotFoundError):
        load_report_inputs(tmp_path / "run")   # nothing there


from rnaforge.report_html import (
    LABELS, _esc, _table, section_confidence, section_dataset, section_quality,
)


def test_labels_tr_en_differ():
    assert LABELS["tr"]["confidence"] != LABELS["en"]["confidence"]
    assert set(LABELS["tr"]) == set(LABELS["en"])   # ayni anahtar seti


def test_esc_and_table():
    assert _esc(None) == "—"
    assert _esc("<x>") == "&lt;x&gt;"
    h = _table(["A", "B"], [["1", "<i>"]])
    assert "<th>A</th>" in h and "&lt;i&gt;" in h


def test_section_confidence_banner_and_gates():
    conf = {"verdict": "SUSPECT", "counts": {"PASS": 10, "WARN": 1, "FAIL": 0},
            "profile": {"name": "prokaryote", "overrides": {}},
            "gates": [{"name": "gc_content", "status": "WARN", "measured": 0.7, "threshold": 0.6}]}
    h = section_confidence(conf, LABELS["tr"])
    assert "SUSPECT" in h and "verdict-suspect" in h and "gc_content" in h


def test_section_dataset_lists_samples():
    raw = {"organism": "E. coli", "platform": "illumina", "design": "~condition",
           "conditions": {"control": 3, "enterololin": 3},
           "samples": [{"sample_id": "c1", "condition": "control", "batch": None,
                        "paired": True, "mean_read_length": 150.0, "mean_quality": 39.0}]}
    h = section_dataset(raw, LABELS["tr"])
    assert "E. coli" in h and "c1" in h


def test_section_quality_rates():
    align = {"samples": {"c1": {"alignment_rate": 0.99}}}
    count = {"samples": {"c1": {"assignment_rate": 0.85}}, "n_genes": 4398}
    h = section_quality(align, count, {"min_length": 36, "aggressive": False}, LABELS["tr"])
    assert "c1" in h and "99" in h
