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
    parse_coldata, parse_normalized_counts, condition_layout, cond_mean,
)


def test_parse_coldata_order(tmp_path):
    p = tmp_path / "coldata.tsv"
    p.write_text("sample\tcondition\nc1\tcontrol\nc2\tcontrol\nt1\ttreated\n")
    assert parse_coldata(p) == [("c1", "control"), ("c2", "control"), ("t1", "treated")]


def test_parse_normalized_counts(tmp_path):
    p = tmp_path / "nc.tsv"
    p.write_text("gene\tc1\tc2\nLT_1\t10\t20\nLT_2\t5\tNA\n")
    nc = parse_normalized_counts(p)
    assert nc["LT_1"] == {"c1": 10.0, "c2": 20.0}
    assert nc["LT_2"]["c1"] == 5.0 and nc["LT_2"]["c2"] is None


def test_condition_layout_first_appearance():
    order, samples = condition_layout([("c1", "control"), ("t1", "treated"), ("c2", "control")])
    assert order == ["control", "treated"]
    assert samples == {"control": ["c1", "c2"], "treated": ["t1"]}


def test_cond_mean():
    nc = {"LT_1": {"c1": 10.0, "c2": 30.0, "t1": 100.0}}
    assert cond_mean("LT_1", ["c1", "c2"], nc) == 20.0
    assert cond_mean("LT_1", ["t1"], nc) == 100.0
    assert cond_mean("missing", ["c1"], nc) is None    # gene not in matrix -> None


def test_top_degs_keeps_gene_id():
    de = [{"gene": "LT_1", "baseMean": 100.0, "log2FoldChange": 3.0, "padj": 1e-8}]
    out = top_degs(de, {"LT_1": "pspA"}, 0.05, 1.0)
    assert out[0]["gene"] == "pspA" and out[0]["gene_id"] == "LT_1"


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
          "n_up": 807, "n_down": 827,
          "fdr_threshold": 0.05, "log2fc_threshold": 1.0, "min_replicate_correlation": 0.98}
    h = section_de(de, LABELS["en"])
    assert "1634" in h and "4398" in h and "t vs c" in h
    assert "807" in h and "827" in h    # n_up / n_down surfaced


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


from rnaforge.report_html import top_degs_by_direction


def test_top_degs_by_direction_filters():
    de = [
        {"gene": "U1", "baseMean": 100.0, "log2FoldChange": 3.0, "padj": 1e-8},
        {"gene": "U2", "baseMean": 100.0, "log2FoldChange": 2.0, "padj": 1e-4},
        {"gene": "D1", "baseMean": 100.0, "log2FoldChange": -4.0, "padj": 1e-9},
    ]
    up = top_degs_by_direction(de, {}, 0.05, 1.0, "Up", n=25)
    down = top_degs_by_direction(de, {}, 0.05, 1.0, "Down", n=25)
    assert [r["gene"] for r in up] == ["U1", "U2"]      # padj asc, Up only
    assert [r["gene"] for r in down] == ["D1"]


def test_section_table_has_up_and_down():
    de = [
        {"gene": "U1", "baseMean": 100.0, "log2FoldChange": 3.0, "padj": 1e-8},
        {"gene": "D1", "baseMean": 100.0, "log2FoldChange": -4.0, "padj": 1e-9},
    ]
    h = section_table(de, {"LT": "x"}, 0.05, 1.0, LABELS["tr"])
    assert LABELS["tr"]["up_table"] in h and LABELS["tr"]["down_table"] in h
    assert "U1" in h and "D1" in h


def test_section_table_empty_both():
    assert LABELS["tr"]["no_degs"] in section_table([], {}, 0.05, 1.0, LABELS["tr"])


def test_section_table_condition_expression_columns():
    de = [{"gene": "LT_1", "baseMean": 100.0, "log2FoldChange": 3.0, "padj": 1e-8}]
    cond_ctx = {
        "norm_counts": {"LT_1": {"c1": 10.0, "c2": 30.0, "t1": 200.0, "t2": 220.0}},
        "order": ["control", "treated"],
        "samples": {"control": ["c1", "c2"], "treated": ["t1", "t2"]},
    }
    h = section_table(de, {}, 0.05, 1.0, LABELS["tr"], cond_ctx=cond_ctx)
    assert f'control {LABELS["tr"]["mean_suffix"]}' in h     # condition header
    assert f'treated {LABELS["tr"]["mean_suffix"]}' in h
    assert "20.0" in h and "210.0" in h                       # control mean 20, treated mean 210


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
    # EN methods: proper scientific narrative, parametrized from config
    hm_en = section_methods(cfg, LABELS["en"])
    assert "DESeq2" in hm_en and "Bowtie2" in hm_en and "featureCounts" in hm_en
    assert "36" in hm_en and "Benjamini" in hm_en and "CDS" in hm_en
    # TR methods differ (bilingual)
    hm_tr = section_methods(cfg, LABELS["tr"])
    assert "medyan" in hm_tr.lower() and hm_tr != hm_en
    # References carry verified DOIs (links)
    refs = section_references(LABELS["en"])
    assert "10.1186/s13059-014-0550-8" in refs      # DESeq2
    assert "10.1038/nmeth.1923" in refs             # Bowtie2
    assert "10.1093/bioinformatics/bty560" in refs  # fastp
    assert 'href="https://doi.org/' in refs


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
        "norm_counts": {"LT_1": {"c1": 10.0, "t1": 200.0}},
        "coldata": [("c1", "control"), ("t1", "treated")],
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


from rnaforge.report_html import FIGURE_CAPTIONS, SECTION_INTRO


def test_caption_and_intro_bilingual():
    assert set(FIGURE_CAPTIONS["tr"]) == set(FIGURE_CAPTIONS["en"])
    assert FIGURE_CAPTIONS["tr"]["pca"] != FIGURE_CAPTIONS["en"]["pca"]
    assert set(SECTION_INTRO["tr"]) == set(SECTION_INTRO["en"])


def test_section_figures_shows_caption(tmp_path):
    fig = tmp_path / "figures"; fig.mkdir()
    (fig / "01_pca.png").write_bytes(b"\x89PNG")
    manifest = {"figures": [{"id": "pca", "title": "PCA", "png": "01_pca.png", "svg": None}]}
    h = section_figures(manifest, fig, LABELS["tr"], lang="tr")
    assert FIGURE_CAPTIONS["tr"]["pca"][:12] in h    # caption metni basildi


def test_render_report_has_section_intro(tmp_path):
    doc = render_report(_full_inputs(tmp_path), _cfg("tr"), version="0.1.0")
    assert 'class="intro"' in doc and SECTION_INTRO["tr"]["confidence"][:12] in doc


from rnaforge.report_html import parse_enrichment_tsv, section_enrichment, LABELS


def _enrich_rows_tsv(tmp_path):
    p = tmp_path / "enrichment_up.tsv"
    p.write_text(
        "go_id\tnamespace\tterm\tstudy_count\tstudy_n\tbg_count\tbg_n\texpected\t"
        "fold_enrichment\tp_value\tp_adj\tgenes\n"
        "GO:0009279\tCC\touter membrane\t8\t40\t20\t400\t2.0\t4.0\t1.0e-05\t2.0e-04\tpspA;ompC\n"
        "GO:0006979\tBP\tresponse to oxidative stress\t5\t40\t15\t400\t1.5\t3.3\t1.0e-03\t8.0e-01\tkatE\n"
    )
    return parse_enrichment_tsv(p)


def test_parse_enrichment_tsv_types(tmp_path):
    rows = _enrich_rows_tsv(tmp_path)
    assert rows[0]["go_id"] == "GO:0009279"
    assert rows[0]["study_count"] == 8 and rows[0]["bg_count"] == 20
    assert rows[0]["fold_enrichment"] == 4.0 and rows[0]["p_adj"] == 2.0e-04


def test_parse_enrichment_tsv_missing_empty(tmp_path):
    assert parse_enrichment_tsv(tmp_path / "nope.tsv") == []


def test_section_enrichment_present(tmp_path):
    up = _enrich_rows_tsv(tmp_path)
    inputs = {"enrichment_up": up, "enrichment_down": [], "enrichment_manifest": None,
              "enrichment_dir": tmp_path}
    html = section_enrichment(inputs, LABELS["tr"], "tr")
    assert "Fonksiyonel Zenginleştirme" in html
    assert "outer membrane" in html                      # anlamlı terim (padj<0.05) tabloda
    assert "response to oxidative stress" not in html     # padj=0.8 -> anlamlı değil, elenmiş
    assert "bulunamadı" in html                           # down boş -> "no enrichment" notu


def test_section_enrichment_not_run(tmp_path):
    inputs = {"enrichment_up": None, "enrichment_down": None}
    html = section_enrichment(inputs, LABELS["en"], "en")
    assert "was not run" in html                          # dürüst not, kırılmaz


def test_render_report_go_section_absent_note(tmp_path):
    # _full_inputs enrichment içermez -> render_report "not run" bölümü basar, kırılmaz.
    doc = render_report(_full_inputs(tmp_path), _cfg("tr"), version="0.1.0")
    assert 'id="enrichment"' in doc and "çalıştırılmadı" in doc


from rnaforge.report_html import section_methods, section_references


def test_section_enrichment_has_legend(tmp_path):
    up = _enrich_rows_tsv(tmp_path)
    inputs = {"enrichment_up": up, "enrichment_down": [], "enrichment_manifest": None,
              "enrichment_dir": tmp_path}
    html = section_enrichment(inputs, LABELS["tr"], "tr")
    assert "Biyolojik Süreç" in html and "Moleküler İşlev" in html and "Hücresel Bileşen" in html
    assert "Kat-zenginleşme" in html
    assert "<b>GO id</b>" in html          # açıklama ham HTML olarak render edilir (escape değil)


def test_methods_includes_go_only_when_enrichment_ran():
    cfg = _cfg("tr")
    without = section_methods(cfg, LABELS["tr"], enrichment_ran=False)
    withgo = section_methods(cfg, LABELS["tr"], enrichment_ran=True)
    assert "hipergeometrik" not in without
    assert "hipergeometrik" in withgo and "Benjamini" in withgo


def test_references_includes_go_only_when_enrichment_ran():
    base = section_references(LABELS["en"], enrichment_ran=False)
    withgo = section_references(LABELS["en"], enrichment_ran=True)
    assert "10.1038/75556" not in base                    # GO founding paper
    assert "10.1038/75556" in withgo and "tb02031.x" in withgo  # GO + Benjamini-Hochberg


def test_render_report_go_methods_and_refs_present(tmp_path):
    inp = _full_inputs(tmp_path)
    inp["enrichment_up"] = _enrich_rows_tsv(tmp_path)
    inp["enrichment_down"] = []
    inp["enrichment_manifest"] = None
    inp["enrichment_dir"] = tmp_path
    doc = render_report(inp, _cfg("tr"), version="0.1.0")
    assert "hipergeometrik" in doc                         # yöntemlerde GO paragrafı
    assert "10.1093/nar/gku1113" in doc                    # EBI-GOA kaynağı
    assert "Biyolojik Süreç" in doc                        # tablo açıklaması


def _kegg_rows_tsv(tmp_path):
    p = tmp_path / "kegg_up.tsv"
    p.write_text(
        "go_id\tnamespace\tterm\tstudy_count\tstudy_n\tbg_count\tbg_n\texpected\t"
        "fold_enrichment\tp_value\tp_adj\tgenes\n"
        "eco02020\tKEGG\tTwo-component system\t10\t50\t30\t500\t3.0\t3.3\t1.0e-06\t5.0e-05\tphoB;ompR\n"
    )
    return parse_enrichment_tsv(p)


def test_section_enrichment_kegg_subsection(tmp_path):
    inputs = {"enrichment_up": None, "enrichment_down": None,
              "kegg_up": _kegg_rows_tsv(tmp_path), "kegg_down": [],
              "kegg_manifest": None, "kegg_dir": tmp_path}
    html = section_enrichment(inputs, LABELS["tr"], "tr")
    assert "KEGG Yolakları" in html
    assert "Two-component system" in html          # KEGG namespace tabloda görünüyor
    assert "çalıştırılmadı" not in html            # KEGG koştu -> not-run notu yok


def test_section_enrichment_both_go_and_kegg(tmp_path):
    inputs = {"enrichment_up": _enrich_rows_tsv(tmp_path), "enrichment_down": [],
              "enrichment_manifest": None, "enrichment_dir": tmp_path,
              "kegg_up": _kegg_rows_tsv(tmp_path), "kegg_down": [],
              "kegg_manifest": None, "kegg_dir": tmp_path}
    html = section_enrichment(inputs, LABELS["en"], "en")
    assert "Gene Ontology (GO)" in html and "KEGG Pathways" in html
    assert "outer membrane" in html and "Two-component system" in html


def test_methods_and_refs_kegg_when_ran():
    cfg = _cfg("tr")
    cfg = __import__("dataclasses").replace(
        cfg, enrichment=__import__("dataclasses").replace(cfg.enrichment, kegg_organism="eco"))
    withk = section_methods(cfg, LABELS["tr"], enrichment_ran=True, kegg_ran=True)
    assert "KEGG" in withk and "eco" in withk
    refs = section_references(LABELS["en"], enrichment_ran=True, kegg_ran=True)
    assert "10.1093/nar/28.1.27" in refs           # Kanehisa & Goto 2000
    refs_no = section_references(LABELS["en"], enrichment_ran=True, kegg_ran=False)
    assert "10.1093/nar/28.1.27" not in refs_no
