"""m13 orkestrasyon testleri: ön koşul, genome gürültülü hata, çıktı tabloları, gate-yok, resume."""
from __future__ import annotations

from pathlib import Path

import pytest

from rnaforge.config import load_config
from rnaforge.modules import m13_amr
from rnaforge.state import RunState

GFF = (
    "chr1\tRefSeq\tgene\t90\t1300\t.\t+\t.\tID=g1;gene=acrB;locus_tag=LT_1\n"
    "chr1\tRefSeq\tgene\t2000\t2500\t.\t+\t.\tID=g2;gene=ompF;locus_tag=LT_2\n"
)
DESEQ = (
    "gene\tbaseMean\tlog2FoldChange\tlfcSE\tstat\tpvalue\tpadj\n"
    "LT_1\t100\t2.5\t0.2\t8\t1e-9\t1e-8\n"       # up
    "LT_2\t100\t0.1\t0.2\t0.3\t0.7\t0.8\n"        # ns
)
ABR_HEADER = ("#FILE\tSEQUENCE\tSTART\tEND\tSTRAND\tGENE\tCOVERAGE\tCOVERAGE_MAP\tGAPS\t"
              "%COVERAGE\t%IDENTITY\tDATABASE\tACCESSION\tPRODUCT\tRESISTANCE\n")


def _setup(tmp_path, with_m06=True, genome=True):
    gff = tmp_path / "g.gff"; gff.write_text(GFF)
    gpath = tmp_path / "genome.fa"
    if genome:
        gpath.write_text(">chr1\nACGT\n")
    (tmp_path / "c.yaml").write_text(
        "organism: E\norganism_type: prokaryote\nplatform: auto\n"
        f"reference:\n  genome_fasta: {gpath}\n  annotation_gff: {gff}\n"
        "de:\n  design: '~condition'\n  fdr_threshold: 0.05\n  log2fc_threshold: 1.0\n"
        "amr:\n  min_identity: 80\n  min_coverage: 80\n")
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


def _fake_abricate(monkeypatch, out_dir_holder):
    # abricate yerine: db'ye göre TSV yaz (card -> acrB hit LT_1; vfdb -> hit LT_2)
    def fake(genome_fa, db, out_tsv, env="rnaforge-amr"):
        Path(out_tsv).parent.mkdir(parents=True, exist_ok=True)
        if db == "card":
            Path(out_tsv).write_text(
                ABR_HEADER +
                "g\tchr1\t100\t1200\t+\tacrB\t1-1/1\t=\t0/0\t100.0\t99.5\tcard\tA1\tefflux\tMULTIDRUG\n")
        else:  # vfdb
            Path(out_tsv).write_text(
                ABR_HEADER +
                "g\tchr1\t2010\t2400\t+\tvfX\t1-1/1\t=\t0/0\t100.0\t95.0\tvfdb\tA2\tadhesin\t\n")
        return ""
    monkeypatch.setattr(m13_amr, "run_abricate", fake)


def test_run_amr_requires_m06(tmp_path, monkeypatch):
    cfg, md, rd = _setup(tmp_path, with_m06=False)
    _fake_abricate(monkeypatch, None)
    with pytest.raises(ValueError, match="m06"):
        m13_amr.run_amr(cfg, md, rd)


def test_run_amr_requires_genome(tmp_path, monkeypatch):
    cfg, md, rd = _setup(tmp_path, genome=False)
    _fake_abricate(monkeypatch, None)
    with pytest.raises(FileNotFoundError, match="genome_fasta"):
        m13_amr.run_amr(cfg, md, rd)


def test_run_amr_writes_tables_and_stats(tmp_path, monkeypatch):
    cfg, md, rd = _setup(tmp_path)
    _fake_abricate(monkeypatch, None)
    s = m13_amr.run_amr(cfg, md, rd)
    amr_tsv = rd / "amr" / "amr_genes.tsv"
    vir_tsv = rd / "amr" / "virulence_genes.tsv"
    assert amr_tsv.exists() and vir_tsv.exists()
    assert (rd / "statistics" / "amr_statistics.json").exists()
    assert s["n_amr_genes"] == 1 and s["n_amr_de"] == 1     # acrB -> LT_1 up
    assert s["n_vir_genes"] == 1                             # vfX -> LT_2 (ns)
    body = amr_tsv.read_text().splitlines()
    assert body[0].split("\t") == ["gene", "locus_tag", "db", "label", "pct_identity",
                                    "pct_coverage", "log2fc", "padj", "de_status"]
    assert body[1].startswith("acrB\tLT_1\tcard") and body[1].endswith("up")


def test_run_amr_no_gate(tmp_path, monkeypatch):
    cfg, md, rd = _setup(tmp_path)
    _fake_abricate(monkeypatch, None)
    m13_amr.run_amr(cfg, md, rd)
    assert not (rd / "quality" / "gates.json").exists()
    assert RunState(rd).is_done("m13_amr")


def test_run_amr_resume(tmp_path, monkeypatch):
    cfg, md, rd = _setup(tmp_path)
    _fake_abricate(monkeypatch, None)
    m13_amr.run_amr(cfg, md, rd)
    s2 = m13_amr.run_amr(cfg, md, rd)
    assert s2.get("resumed") is True
