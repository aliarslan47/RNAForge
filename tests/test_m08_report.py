import json
from pathlib import Path
import pytest
from rnaforge.modules import m08_report
from rnaforge.state import RunState


def _cfg(tmp_path, lang="tr"):
    from rnaforge.config import load_config
    (tmp_path / "c.yaml").write_text(
        "organism: E\norganism_type: prokaryote\nplatform: auto\n"
        "reference:\n  genome_fasta: g.fa\n  annotation_gff: g.gff\n"
        "de:\n  design: '~condition'\n  fdr_threshold: 0.05\n  log2fc_threshold: 1.0\n"
        f"report:\n  language: {lang}\n")
    return load_config(tmp_path / "c.yaml")


def _seed_run(rd):
    (rd / "statistics").mkdir(parents=True)
    (rd / "quality").mkdir(parents=True)
    (rd / "logs").mkdir(parents=True)
    de = rd / "differential_expression"; de.mkdir(parents=True)
    fig = rd / "figures"; fig.mkdir(parents=True)
    for base in ("01_pca", "02_volcano", "03_heatmap", "04_ma"):
        (fig / f"{base}.png").write_bytes(b"\x89PNG")
    (fig / "manifest.json").write_text(json.dumps({"figures": [
        {"id": "pca", "title": "PCA", "png": "01_pca.png", "svg": None},
        {"id": "volcano", "title": "Volcano", "png": "02_volcano.png", "svg": None},
        {"id": "heatmap", "title": "Heatmap", "png": "03_heatmap.png", "svg": None},
        {"id": "ma", "title": "MA", "png": "04_ma.png", "svg": None}]}))
    (fig / "gene_map.tsv").write_text("locus_tag\tgene\nLT_1\tpspA\n")
    (de / "deseq2_results.tsv").write_text(
        "gene\tbaseMean\tlog2FoldChange\tlfcSE\tstat\tpvalue\tpadj\nLT_1\t200\t3.0\t0.2\t5\t1e-9\t1e-8\n")
    (de / "normalized_counts.tsv").write_text("gene\tc1\tt1\nLT_1\t50\t400\n")
    (de / "coldata.tsv").write_text("sample\tcondition\nc1\tcontrol\nt1\ttreated\n")
    s = rd / "statistics"
    (s / "raw_statistics.json").write_text(json.dumps({"organism": "E. coli", "platform": "illumina",
        "design": "~condition", "conditions": {"control": 1, "treated": 1},
        "samples": [{"sample_id": "c1", "condition": "control", "batch": None, "paired": True,
                     "mean_read_length": 150.0, "mean_quality": 39.0}]}))
    for name in ("qc", "trimming"):
        (s / f"{name}_statistics.json").write_text(json.dumps({"samples": {}}))
    (s / "alignment_statistics.json").write_text(json.dumps({"samples": {"c1": {"alignment_rate": 0.99}}}))
    (s / "count_statistics.json").write_text(json.dumps({"samples": {"c1": {"assignment_rate": 0.85}}, "n_genes": 4398}))
    (s / "de_statistics.json").write_text(json.dumps({"contrast": "treated vs control", "n_genes": 4398,
        "n_significant": 1, "fdr_threshold": 0.05, "log2fc_threshold": 1.0, "min_replicate_correlation": 0.98}))
    (rd / "quality" / "confidence_card.json").write_text(json.dumps({"verdict": "TRUSTWORTHY",
        "counts": {"PASS": 11, "WARN": 0, "FAIL": 0}, "profile": {"name": "prokaryote", "overrides": {}}, "gates": []}))


def test_run_report_requires_m07(tmp_path):
    rd = tmp_path / "run"; _seed_run(rd)
    RunState(rd).mark_done("m06_de", [])   # m07 NOT done
    with pytest.raises(ValueError):
        m08_report.run_report(_cfg(tmp_path), tmp_path / "m.tsv", rd)


def test_run_report_writes_and_resumes(tmp_path):
    rd = tmp_path / "run"; _seed_run(rd)
    RunState(rd).mark_done("m07_figures", [])
    s = m08_report.run_report(_cfg(tmp_path), tmp_path / "m.tsv", rd)
    assert s["n_sections"] == 16
    doc = (rd / "report" / "report.html").read_text()
    assert "TRUSTWORTHY" in doc and doc.count('src="data:image/png;base64,') == 4
    assert (rd / "statistics" / "report_statistics.json").exists()
    assert RunState(rd).is_done("m08_report")
    s2 = m08_report.run_report(_cfg(tmp_path), tmp_path / "m.tsv", rd)
    assert s2.get("resumed") is True
