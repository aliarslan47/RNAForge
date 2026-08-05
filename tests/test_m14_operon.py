"""m14 orkestrasyon testleri: ön koşul, çıktı tablosu, gate-yok, resume."""
from __future__ import annotations

import pytest

from rnaforge.config import load_config
from rnaforge.modules import m14_operon
from rnaforge.state import RunState


@pytest.fixture(autouse=True)
def _no_operon_r(monkeypatch):
    # Figür best-effort; birim testte R çağrısını atla (env-bağımsız, hızlı).
    monkeypatch.setattr(m14_operon, "run_operon_r", lambda *a, **k: "")

GFF = (
    "c1\tx\tgene\t100\t200\t.\t+\t.\tlocus_tag=A;gene=a\n"
    "c1\tx\tgene\t210\t300\t.\t+\t.\tlocus_tag=B;gene=b\n"
    "c1\tx\tgene\t500\t600\t.\t+\t.\tlocus_tag=C;gene=c\n"
)
DESEQ = (
    "gene\tbaseMean\tlog2FoldChange\tlfcSE\tstat\tpvalue\tpadj\n"
    "A\t100\t2.5\t0.2\t8\t1e-9\t1e-8\n"
    "B\t100\t2.2\t0.2\t7\t1e-8\t1e-7\n"
    "C\t100\t0.1\t0.2\t0.3\t0.7\t0.8\n"
)


def _setup(tmp_path, with_m06=True):
    gff = tmp_path / "g.gff"; gff.write_text(GFF)
    (tmp_path / "c.yaml").write_text(
        "organism: E\norganism_type: prokaryote\nplatform: auto\n"
        f"reference:\n  genome_fasta: g.fa\n  annotation_gff: {gff}\n"
        "de:\n  design: '~condition'\n  fdr_threshold: 0.05\n  log2fc_threshold: 1.0\n"
        "operon:\n  max_gap: 50\n")
    cfg = load_config(tmp_path / "c.yaml")
    rd = tmp_path / "run"
    de = rd / "differential_expression"; de.mkdir(parents=True)
    (de / "deseq2_results.tsv").write_text(DESEQ)
    (rd / "statistics").mkdir(); (rd / "logs").mkdir()
    if with_m06:
        RunState(rd).mark_done("m06_de", [])
    md = tmp_path / "m.tsv"
    md.write_text("sample_id\tcondition\tfastq_1\nc1\tcontrol\t/x/a.fq\n")
    return cfg, md, rd


def test_run_operon_requires_m06(tmp_path):
    cfg, md, rd = _setup(tmp_path, with_m06=False)
    with pytest.raises(ValueError, match="m06"):
        m14_operon.run_operon(cfg, md, rd)


def test_run_operon_writes_table_and_stats(tmp_path):
    cfg, md, rd = _setup(tmp_path)
    s = m14_operon.run_operon(cfg, md, rd)
    tsv = rd / "operon" / "operons.tsv"
    assert tsv.exists() and (rd / "statistics" / "operon_statistics.json").exists()
    assert s["n_operons"] == 2                    # [A,B] + [C]
    assert s["n_coordinated"] == 1                # A,B ikisi up
    body = tsv.read_text().splitlines()
    assert body[0].split("\t") == ["operon_id", "contig", "strand", "genes", "size",
                                    "n_tested", "n_deg", "n_up", "n_down", "mean_log2fc", "coordinated"]
    assert body[1].split("\t")[3] == "a;b" and body[1].endswith("yes")   # koordineli önce


def test_run_operon_no_gate(tmp_path):
    cfg, md, rd = _setup(tmp_path)
    m14_operon.run_operon(cfg, md, rd)
    assert not (rd / "quality" / "gates.json").exists()
    assert RunState(rd).is_done("m14_operon")


def test_run_operon_resume(tmp_path):
    cfg, md, rd = _setup(tmp_path)
    m14_operon.run_operon(cfg, md, rd)
    s2 = m14_operon.run_operon(cfg, md, rd)
    assert s2.get("resumed") is True
