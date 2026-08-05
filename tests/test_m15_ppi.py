"""m15 orkestrasyon testleri: ön koşul, taxid/dosya gürültülü hata, çıktı, gate-yok, resume."""
from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from rnaforge.config import load_config
from rnaforge.modules import m15_ppi
from rnaforge.state import RunState


@pytest.fixture(autouse=True)
def _no_ppi_r(monkeypatch):
    # Ağ figürü best-effort; birim testte R çağrısını atla (networkx layout saf-Python, kalır).
    monkeypatch.setattr(m15_ppi, "run_ppi_r", lambda *a, **k: "")

GFF = (
    "c1\tx\tCDS\t1\t9\t.\t+\t0\tgene=a;locus_tag=LT_1\n"
    "c1\tx\tCDS\t1\t9\t.\t+\t0\tgene=b;locus_tag=LT_2\n"
    "c1\tx\tCDS\t1\t9\t.\t+\t0\tgene=c;locus_tag=LT_3\n"
)
DESEQ = (
    "gene\tbaseMean\tlog2FoldChange\tlfcSE\tstat\tpvalue\tpadj\n"
    "LT_1\t100\t2.5\t0.2\t8\t1e-9\t1e-8\n"      # up
    "LT_2\t100\t2.2\t0.2\t7\t1e-8\t1e-7\n"      # up
    "LT_3\t100\t2.1\t0.2\t7\t1e-8\t1e-7\n"      # up
)


def _gz(path, text):
    with gzip.open(path, "wt") as f:
        f.write(text)
    return path


def _setup(tmp_path, taxid="9999", files=True):
    gff = tmp_path / "g.gff"; gff.write_text(GFF)
    sdir = tmp_path / "string"; sdir.mkdir()
    if files:
        _gz(sdir / "protein.info.txt.gz",
            "#string_protein_id\tpreferred_name\n9999.p1\ta\n9999.p2\tb\n9999.p3\tc\n")
        _gz(sdir / "protein.links.txt.gz",
            "protein1 protein2 combined_score\n9999.p1 9999.p2 900\n9999.p2 9999.p3 850\n9999.p1 9999.p3 800\n")
    taxid_line = f"  taxid: '{taxid}'\n" if taxid else ""
    (tmp_path / "c.yaml").write_text(
        "organism: E\norganism_type: prokaryote\nplatform: auto\n"
        f"reference:\n  genome_fasta: g.fa\n  annotation_gff: {gff}\n"
        "de:\n  design: '~condition'\n  fdr_threshold: 0.05\n  log2fc_threshold: 1.0\n"
        f"ppi:\n  min_score: 700\n  min_community_size: 2\n{taxid_line}  string_dir: {sdir}\n")
    cfg = load_config(tmp_path / "c.yaml")
    rd = tmp_path / "run"
    de = rd / "differential_expression"; de.mkdir(parents=True)
    (de / "deseq2_results.tsv").write_text(DESEQ)
    (rd / "statistics").mkdir(); (rd / "logs").mkdir()
    RunState(rd).mark_done("m06_de", [])
    md = tmp_path / "m.tsv"
    md.write_text("sample_id\tcondition\tfastq_1\nc1\tcontrol\t/x/a.fq\n")
    return cfg, md, rd


def test_run_ppi_requires_m06(tmp_path):
    cfg, md, rd = _setup(tmp_path)
    RunState(rd)  # m06 done; şimdi taze dizin
    fresh = tmp_path / "run2"; (fresh / "differential_expression").mkdir(parents=True)
    (fresh / "differential_expression" / "deseq2_results.tsv").write_text(DESEQ)
    with pytest.raises(ValueError, match="m06"):
        m15_ppi.run_ppi(cfg, md, fresh)


def test_run_ppi_requires_taxid(tmp_path):
    cfg, md, rd = _setup(tmp_path, taxid=None)
    with pytest.raises(ValueError, match="taxid"):
        m15_ppi.run_ppi(cfg, md, rd)


def test_run_ppi_missing_files_loud(tmp_path):
    cfg, md, rd = _setup(tmp_path, files=False)
    with pytest.raises(FileNotFoundError, match="STRING files missing"):
        m15_ppi.run_ppi(cfg, md, rd)


def test_run_ppi_writes_communities_and_stats(tmp_path):
    cfg, md, rd = _setup(tmp_path)
    s = m15_ppi.run_ppi(cfg, md, rd)
    tsv = rd / "ppi" / "communities.tsv"
    assert tsv.exists() and (rd / "statistics" / "ppi_statistics.json").exists()
    assert s["n_deg"] == 3 and s["n_deg_in_network"] == 3 and s["n_edges"] == 3
    assert s["n_communities"] == 1                 # 3 gen tek modül (üçgen)
    body = tsv.read_text().splitlines()
    assert body[0].split("\t") == ["community_id", "size", "n_up", "n_down", "dominant", "genes"]
    assert body[1].split("\t")[4] == "up" and body[1].split("\t")[5] == "a;b;c"


def test_run_ppi_no_gate(tmp_path):
    cfg, md, rd = _setup(tmp_path)
    m15_ppi.run_ppi(cfg, md, rd)
    assert not (rd / "quality" / "gates.json").exists()
    assert RunState(rd).is_done("m15_ppi")


def test_run_ppi_resume(tmp_path):
    cfg, md, rd = _setup(tmp_path)
    m15_ppi.run_ppi(cfg, md, rd)
    s2 = m15_ppi.run_ppi(cfg, md, rd)
    assert s2.get("resumed") is True
