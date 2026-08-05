"""m12 orkestrasyon testleri: ön koşul, obo gürültülü hata, indirgeme, gate-yok, resume."""
from __future__ import annotations

import dataclasses
import textwrap
from pathlib import Path

import pytest

from rnaforge.config import load_config
from rnaforge.modules import m12_semantic
from rnaforge.state import RunState

# GFF: kardeş GO'lar aynı gene ailesinde; obo BP zinciri.
GFF = (
    "NZ\tX\tCDS\t1\t9\t.\t+\t0\tOntology_term=GO:0000002;go_process=alpha|0000002||IEA;gene=gA;locus_tag=LT_A\n"
    "NZ\tX\tCDS\t1\t9\t.\t+\t0\tOntology_term=GO:0000003;go_process=beta|0000003||IEA;gene=gB;locus_tag=LT_B\n"
    "NZ\tX\tCDS\t1\t9\t.\t+\t0\tOntology_term=GO:0000002;go_process=alpha|0000002||IEA;gene=gC;locus_tag=LT_C\n"
)
OBO = textwrap.dedent("""\
    [Term]
    id: GO:0000001
    name: root
    namespace: biological_process

    [Term]
    id: GO:0000002
    name: alpha
    namespace: biological_process
    is_a: GO:0000001 ! root

    [Term]
    id: GO:0000003
    name: beta
    namespace: biological_process
    is_a: GO:0000001 ! root
    """)
# m09 enrichment_up.tsv (kardeş GO:0000002 ve GO:0000003, ikisi de anlamlı)
ORA_UP = (
    "go_id\tnamespace\tterm\tstudy_count\tstudy_n\tbg_count\tbg_n\texpected\tfold_enrichment\tp_value\tp_adj\tgenes\n"
    "GO:0000002\tBP\talpha\t2\t3\t2\t3\t1\t1\t1e-5\t1e-4\tgA;gC\n"
    "GO:0000003\tBP\tbeta\t1\t3\t1\t3\t1\t1\t1e-3\t2e-3\tgB\n"
)


def _setup(tmp_path, obo=True, m09=True, m11=False):
    gff = tmp_path / "g.gff"; gff.write_text(GFF)
    obo_line = ""
    if obo:
        (tmp_path / "go.obo").write_text(OBO)
        obo_line = f"  obo: {tmp_path/'go.obo'}\n"
    (tmp_path / "c.yaml").write_text(
        "organism: E\norganism_type: prokaryote\nplatform: auto\n"
        f"reference:\n  genome_fasta: g.fa\n  annotation_gff: {gff}\n"
        "de:\n  design: '~condition'\n"
        f"enrichment:\n  revigo_similarity: 0.1\n{obo_line}")
    cfg = load_config(tmp_path / "c.yaml")
    rd = tmp_path / "run"
    (rd / "statistics").mkdir(parents=True); (rd / "logs").mkdir()
    st = RunState(rd)
    if m09:
        (rd / "enrichment").mkdir()
        (rd / "enrichment" / "enrichment_up.tsv").write_text(ORA_UP)
        (rd / "enrichment" / "enrichment_down.tsv").write_text(ORA_UP.split("\n")[0] + "\n")
        st.mark_done("m09_enrichment", [])
    if m11:
        st.mark_done("m11_gsea", [])
    md = tmp_path / "m.tsv"
    md.write_text("sample_id\tcondition\tfastq_1\nc1\tcontrol\t/x/a.fq\n")
    return cfg, md, rd


def test_run_semantic_requires_go_source(tmp_path):
    cfg, md, rd = _setup(tmp_path, m09=False, m11=False)
    with pytest.raises(ValueError, match="GO enrichment source"):
        m12_semantic.run_semantic(cfg, md, rd)


def test_run_semantic_requires_obo(tmp_path):
    cfg, md, rd = _setup(tmp_path, obo=True)
    bad = dataclasses.replace(cfg.enrichment, obo=Path(tmp_path / "yok.obo"))
    cfg = dataclasses.replace(cfg, enrichment=bad)
    with pytest.raises(FileNotFoundError, match="obo"):
        m12_semantic.run_semantic(cfg, md, rd)


def test_run_semantic_reduces_ora(tmp_path):
    cfg, md, rd = _setup(tmp_path)
    s = m12_semantic.run_semantic(cfg, md, rd)
    red = rd / "semantic" / "reduced_ora_up.tsv"
    assert red.exists()
    lines = red.read_text().splitlines()
    # Uçtan uca kablo: 2 anlamlı terim okundu, temsilci(ler) yazıldı (indirgeme matematiği
    # test_semantic.py'de). Küçük fixture'da kardeşlerin tek ortak atası kök (IC=0) -> Lin=0.
    assert s["collections"]["ora_up"]["n_terms"] == 2
    assert 1 <= s["collections"]["ora_up"]["n_representatives"] <= 2
    assert lines[0].split("\t") == ["go_id", "namespace", "term", "padj", "n_collapsed", "members"]
    assert lines[1].startswith("GO:0000002\tBP\talpha")   # en iyi padj ilk temsilci


def test_run_semantic_no_gate(tmp_path):
    cfg, md, rd = _setup(tmp_path)
    m12_semantic.run_semantic(cfg, md, rd)
    assert not (rd / "quality" / "gates.json").exists()
    assert RunState(rd).is_done("m12_semantic")


def test_run_semantic_resume(tmp_path):
    cfg, md, rd = _setup(tmp_path)
    m12_semantic.run_semantic(cfg, md, rd)
    s2 = m12_semantic.run_semantic(cfg, md, rd)
    assert s2.get("resumed") is True
