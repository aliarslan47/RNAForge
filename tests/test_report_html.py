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


from rnaforge.report_html import (
    section_de, section_figures, section_table, section_methods, section_references,
)


def test_section_de_counts():
    de = {"contrast": "t vs c", "n_genes": 4398, "n_significant": 1634,
          "fdr_threshold": 0.05, "log2fc_threshold": 1.0, "min_replicate_correlation": 0.98}
    h = section_de(de, LABELS["en"])
    assert "1634" in h and "4398" in h and "t vs c" in h


def test_section_figures_embeds(tmp_path):
    fig = tmp_path / "figures"; fig.mkdir()
    (fig / "01_pca.png").write_bytes(b"\x89PNG")
    manifest = {"figures": [{"id": "pca", "title": "PCA", "png": "01_pca.png", "svg": None}]}
    h = section_figures(manifest, fig, LABELS["en"])
    assert 'src="data:image/png;base64,' in h and "PCA" in h


def test_section_figures_loud_when_png_missing(tmp_path):
    fig = tmp_path / "figures"; fig.mkdir()
    manifest = {"figures": [{"id": "pca", "title": "PCA", "png": "missing.png", "svg": None}]}
    with pytest.raises(FileNotFoundError):
        section_figures(manifest, fig, LABELS["en"])


def test_section_table_empty_note():
    assert LABELS["tr"]["no_degs"] in section_table([], LABELS["tr"])
    h = section_table([{"gene": "pspA", "log2fc": 3.0, "padj": 1e-8,
                        "base_mean": 200.0, "direction": "Up"}], LABELS["tr"])
    assert "pspA" in h


def test_section_methods_and_references():
    from rnaforge.config import load_config
    import tempfile, os
    cfgtext = ("organism: E\norganism_type: prokaryote\nplatform: auto\n"
               "reference:\n  genome_fasta: g.fa\n  annotation_gff: g.gff\n"
               "trimming:\n  min_length: 36\n  aggressive_quality: false\n"
               "quantification:\n  feature_type: CDS\n  attribute: locus_tag\n"
               "de:\n  design: '~condition'\n  fdr_threshold: 0.05\n  log2fc_threshold: 1.0\n")
    p = tempfile.NamedTemporaryFile("w", suffix=".yaml", delete=False); p.write(cfgtext); p.close()
    cfg = load_config(p.name); os.unlink(p.name)
    hm = section_methods(cfg, LABELS["en"])
    assert "DESeq2" in hm and "bowtie2" in hm and "36" in hm
    assert "DESeq2" in section_references(LABELS["en"])


from rnaforge.report_html import render_report, N_SECTIONS


def _full_inputs(tmp_path):
    fig = tmp_path / "figures"; fig.mkdir()
    for base in ("01_pca", "02_volcano", "03_heatmap", "04_ma"):
        (fig / f"{base}.png").write_bytes(b"\x89PNG")
    return {
        "raw": {"organism": "E. coli", "platform": "illumina", "design": "~condition",
                "conditions": {"control": 2, "treated": 2},
                "samples": [{"sample_id": "c1", "condition": "control", "batch": None,
                             "paired": True, "mean_read_length": 150.0, "mean_quality": 39.0}]},
        "qc": {"samples": {}}, "trimming": {"samples": {}},
        "alignment": {"samples": {"c1": {"alignment_rate": 0.99}}},
        "count": {"samples": {"c1": {"assignment_rate": 0.85}}, "n_genes": 4398},
        "de": {"contrast": "treated vs control", "n_genes": 4398, "n_significant": 2,
               "fdr_threshold": 0.05, "log2fc_threshold": 1.0, "min_replicate_correlation": 0.98},
        "confidence": {"verdict": "TRUSTWORTHY", "counts": {"PASS": 11, "WARN": 0, "FAIL": 0},
                       "profile": {"name": "prokaryote", "overrides": {}}, "gates": []},
        "figures": {"figures": [{"id": "pca", "title": "PCA", "png": "01_pca.png", "svg": None},
                                {"id": "volcano", "title": "Volcano", "png": "02_volcano.png", "svg": None},
                                {"id": "heatmap", "title": "Heatmap", "png": "03_heatmap.png", "svg": None},
                                {"id": "ma", "title": "MA", "png": "04_ma.png", "svg": None}]},
        "de_results": [{"gene": "LT_1", "baseMean": 200.0, "log2FoldChange": 3.0, "padj": 1e-8}],
        "gene_map": {"LT_1": "pspA"},
        "figures_dir": fig,
    }


def _cfg(lang):
    from rnaforge.config import load_config
    import tempfile, os
    txt = (f"organism: E\norganism_type: prokaryote\nplatform: auto\n"
           f"reference:\n  genome_fasta: g.fa\n  annotation_gff: g.gff\n"
           f"de:\n  design: '~condition'\n  fdr_threshold: 0.05\n  log2fc_threshold: 1.0\n"
           f"report:\n  language: {lang}\n")
    p = tempfile.NamedTemporaryFile("w", suffix=".yaml", delete=False); p.write(txt); p.close()
    c = load_config(p.name); os.unlink(p.name); return c


def test_render_report_full(tmp_path):
    doc = render_report(_full_inputs(tmp_path), _cfg("tr"), version="0.1.0")
    assert doc.lstrip().startswith("<!doctype html>") or "<html" in doc
    assert doc.count('src="data:image/png;base64,') == 4        # all figures embedded
    assert "TRUSTWORTHY" in doc and "pspA" in doc and "Güvence" in doc
    assert "0.1.0" in doc


def test_render_report_language_switch(tmp_path):
    en = render_report(_full_inputs(tmp_path), _cfg("en"), version="0.1.0")
    assert "Confidence Card" in en and "Güvence" not in en


def test_render_report_includes_run_id(tmp_path):
    doc = render_report(_full_inputs(tmp_path), _cfg("tr"), version="0.1.0",
                        run_id="20260803_143036_GSE300731")
    assert "20260803_143036_GSE300731" in doc
