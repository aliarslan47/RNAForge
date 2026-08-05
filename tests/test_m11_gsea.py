"""m11 orkestrasyon testleri: ön koşul, koleksiyon/stat gürültülü hata, gate-yok, çıktı, resume."""
from __future__ import annotations

import textwrap
from pathlib import Path

import pytest

from rnaforge.config import load_config
from rnaforge.modules import m11_gsea
from rnaforge.state import RunState

GFF = (
    "NZ\tX\tCDS\t1\t9\t.\t+\t0\tOntology_term=GO:0000002;go_process=mid|0000002||IEA;gene=gA;locus_tag=LT_A\n"
    "NZ\tX\tCDS\t1\t9\t.\t+\t0\tOntology_term=GO:0000002;go_process=mid|0000002||IEA;gene=gB;locus_tag=LT_B\n"
    "NZ\tX\tCDS\t1\t9\t.\t+\t0\tOntology_term=GO:0000002;go_process=mid|0000002||IEA;gene=gC;locus_tag=LT_C\n"
)
OBO = textwrap.dedent("""\
    [Term]
    id: GO:0000001
    name: root
    namespace: biological_process

    [Term]
    id: GO:0000002
    name: mid
    namespace: biological_process
    is_a: GO:0000001 ! root
    """)
DESEQ = (
    "gene\tbaseMean\tlog2FoldChange\tlfcSE\tstat\tpvalue\tpadj\n"
    "LT_A\t100\t3.0\t0.2\t9.0\t1e-10\t1e-9\n"
    "LT_B\t100\t2.0\t0.2\t6.0\t1e-8\t1e-7\n"
    "LT_C\t100\t-2.0\t0.2\t-6.0\t1e-8\t1e-7\n"
)
DESEQ_NO_STAT = "gene\tbaseMean\tlog2FoldChange\tpadj\nLT_A\t100\t3.0\t1e-9\n"


def _setup(tmp_path, obo=True, deseq=DESEQ):
    gff = tmp_path / "g.gff"; gff.write_text(GFF)
    obo_line = ""
    if obo:
        (tmp_path / "go.obo").write_text(OBO)
        obo_line = f"  obo: {tmp_path/'go.obo'}\n"
    (tmp_path / "c.yaml").write_text(
        "organism: E\norganism_type: prokaryote\nplatform: auto\n"
        f"reference:\n  genome_fasta: g.fa\n  annotation_gff: {gff}\n"
        "de:\n  design: '~condition'\n  fdr_threshold: 0.05\n  log2fc_threshold: 1.0\n"
        f"enrichment:\n  min_term_size: 1\n  gsea_min_size: 1\n{obo_line}")
    cfg = load_config(tmp_path / "c.yaml")
    rd = tmp_path / "run"
    de = rd / "differential_expression"; de.mkdir(parents=True)
    (de / "deseq2_results.tsv").write_text(deseq)
    (rd / "statistics").mkdir(); (rd / "logs").mkdir()
    RunState(rd).mark_done("m06_de", [])
    md = tmp_path / "m.tsv"
    md.write_text("sample_id\tcondition\tfastq_1\nc1\tcontrol\t/x/a.fq\n")
    return cfg, md, rd


def _fake_r(monkeypatch):
    def fake(rnk, gmt, gene_map, out_dir, collection, min_size, max_size, title, env="rnaforge-de"):
        Path(out_dir).mkdir(parents=True, exist_ok=True)
        (Path(out_dir) / f"gsea_{collection}.tsv").write_text(
            "pathway_id\tname\tsize\tES\tNES\tpval\tpadj\tleading_edge\n"
            "GO:0000002\tmid\t3\t0.8\t1.9\t1e-3\t2e-3\tgA;gB\n")
        (Path(out_dir) / f"gsea_{collection}.png").write_bytes(b"x" * 2000)
        (Path(out_dir) / f"gsea_{collection}.svg").write_text("<svg/>")
        return "gsea.R done\n"
    monkeypatch.setattr(m11_gsea, "run_gsea_r", fake)


def test_run_gsea_requires_m06(tmp_path):
    cfg, md, rd = _setup(tmp_path)
    fresh = tmp_path / "run2"                       # m06 done DEĞİL
    (fresh / "differential_expression").mkdir(parents=True)
    (fresh / "differential_expression" / "deseq2_results.tsv").write_text(DESEQ)
    with pytest.raises(ValueError, match="m06"):
        m11_gsea.run_gsea(cfg, md, fresh)


def test_run_gsea_requires_a_collection(tmp_path):
    cfg, md, rd = _setup(tmp_path, obo=False)   # obo yok + kegg_organism yok
    with pytest.raises(ValueError, match="collection"):
        m11_gsea.run_gsea(cfg, md, rd)


def test_run_gsea_missing_stat_loud(tmp_path):
    cfg, md, rd = _setup(tmp_path, deseq=DESEQ_NO_STAT)
    with pytest.raises(ValueError, match="stat"):
        m11_gsea.run_gsea(cfg, md, rd)


def test_run_gsea_writes_outputs_and_no_gate(tmp_path, monkeypatch):
    cfg, md, rd = _setup(tmp_path)
    _fake_r(monkeypatch)
    s = m11_gsea.run_gsea(cfg, md, rd)
    assert (rd / "gsea" / "ranked.rnk").exists()
    assert (rd / "gsea" / "go.gmt").exists()
    assert (rd / "gsea" / "gsea_go.tsv").exists()
    assert (rd / "gsea" / "manifest.json").exists()
    assert (rd / "statistics" / "gsea_statistics.json").exists()
    assert s["n_ranked"] == 3
    assert s["collections"]["go"]["n_sig_pos"] == 1     # NES 1.9, padj 2e-3
    assert not (rd / "quality" / "gates.json").exists() # yeni kapı yok
    assert RunState(rd).is_done("m11_gsea")


def test_run_gsea_resume(tmp_path, monkeypatch):
    cfg, md, rd = _setup(tmp_path)
    _fake_r(monkeypatch)
    m11_gsea.run_gsea(cfg, md, rd)
    s2 = m11_gsea.run_gsea(cfg, md, rd)
    assert s2.get("resumed") is True
