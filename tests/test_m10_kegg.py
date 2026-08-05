"""m10 orkestrasyon testleri: ön koşul, kegg_organism/dosya gürültülü hata, gate-yok, çıktı, resume."""
from __future__ import annotations

import dataclasses
from pathlib import Path

import pytest

from rnaforge.config import load_config
from rnaforge.modules import m10_kegg
from rnaforge.state import RunState


def _cds(gene, lt):
    return f"NZ\tX\tCDS\t1\t9\t.\t+\t0\tgene={gene};locus_tag={lt}\n"


# P1 (eco00260): gA..gD (DEG'ler). P2 (eco00010): g6..g12 (ns, arka planı büyütür).
GFF = ("".join(_cds(f"g{c}", f"LT_{c}") for c in "ABCD")
       + "".join(_cds(f"g{i}", f"LT_{i}") for i in range(6, 13)))
LINKS = ("".join(f"eco:b{c}\tpath:eco00260\n" for c in "ABCD")
         + "".join(f"eco:b{i}\tpath:eco00010\n" for i in range(6, 13)))
NAMES = ("eco00260\tGlycine, serine and threonine metabolism - E. coli\n"
         "eco00010\tGlycolysis / Gluconeogenesis - E. coli\n")
GENELIST = ("".join(f"eco:b{c}\tCDS\t1..9\tg{c}; protein\n" for c in "ABCD")
            + "".join(f"eco:b{i}\tCDS\t1..9\tg{i}; protein\n" for i in range(6, 13)))
DESEQ = (
    "gene\tbaseMean\tlog2FoldChange\tlfcSE\tstat\tpvalue\tpadj\n"
    "LT_A\t100\t3.0\t0.2\t9\t1e-10\t1e-9\n"     # up
    "LT_B\t100\t2.5\t0.2\t8\t1e-9\t1e-8\n"      # up
    "LT_C\t100\t2.2\t0.2\t7\t1e-8\t1e-7\n"      # up
    "LT_D\t100\t-3.0\t0.2\t-9\t1e-10\t1e-9\n"   # down
    + "".join(f"LT_{i}\t100\t0.1\t0.2\t0.3\t0.7\t0.8\n" for i in range(6, 13))  # ns arka plan
)


def _setup(tmp_path, with_m06=True, org="eco"):
    gff = tmp_path / "g.gff"; gff.write_text(GFF)
    kdir = tmp_path / "kegg"; kdir.mkdir()
    (kdir / "pathway_links.tsv").write_text(LINKS)
    (kdir / "pathway_names.tsv").write_text(NAMES)
    (kdir / "gene_list.tsv").write_text(GENELIST)
    org_line = f"  kegg_organism: {org}\n" if org else ""
    (tmp_path / "c.yaml").write_text(
        "organism: E\norganism_type: prokaryote\nplatform: auto\n"
        f"reference:\n  genome_fasta: g.fa\n  annotation_gff: {gff}\n"
        "de:\n  design: '~condition'\n  fdr_threshold: 0.05\n  log2fc_threshold: 1.0\n"
        f"enrichment:\n  min_term_size: 3\n{org_line}  kegg_dir: {kdir}\n")
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


def _fake_r(monkeypatch):
    def fake(up_tsv, down_tsv, out_dir, top_n, title_prefix="", basename_prefix="enrichment",
             env="rnaforge-de"):
        Path(out_dir).mkdir(parents=True, exist_ok=True)
        for d in ("up", "down"):
            (Path(out_dir) / f"{basename_prefix}_{d}.png").write_bytes(b"x" * 2000)
            (Path(out_dir) / f"{basename_prefix}_{d}.svg").write_text("<svg/>")
        return "enrichment.R done\n"
    monkeypatch.setattr(m10_kegg, "run_enrichment_r", fake)


def test_run_kegg_requires_m06(tmp_path):
    cfg, md, rd = _setup(tmp_path, with_m06=False)
    with pytest.raises(ValueError, match="m06"):
        m10_kegg.run_kegg(cfg, md, rd)


def test_run_kegg_requires_organism(tmp_path):
    cfg, md, rd = _setup(tmp_path, org=None)
    with pytest.raises(ValueError, match="kegg_organism"):
        m10_kegg.run_kegg(cfg, md, rd)


def test_run_kegg_missing_files_loud(tmp_path):
    cfg, md, rd = _setup(tmp_path)
    bad = dataclasses.replace(cfg.enrichment, kegg_dir=Path(tmp_path / "yok"))
    cfg = dataclasses.replace(cfg, enrichment=bad)
    with pytest.raises(FileNotFoundError, match="KEGG files missing"):
        m10_kegg.run_kegg(cfg, md, rd)


def test_run_kegg_writes_outputs_and_stats(tmp_path, monkeypatch):
    cfg, md, rd = _setup(tmp_path)
    _fake_r(monkeypatch)
    s = m10_kegg.run_kegg(cfg, md, rd)
    assert (rd / "kegg" / "kegg_up.tsv").exists()
    assert (rd / "kegg" / "kegg_down.tsv").exists()
    assert (rd / "kegg" / "gene2pathway.tsv").exists()
    assert (rd / "kegg" / "manifest.json").exists()
    assert (rd / "statistics" / "kegg_statistics.json").exists()
    assert s["n_up_degs"] == 3 and s["n_down_degs"] == 1
    assert s["n_annotated"] == 11 and s["organism"] == "eco"
    assert s["n_sig_up"] >= 1              # up=3 gen hepsi eco00260 -> zenginleşir


def test_run_kegg_no_gate(tmp_path, monkeypatch):
    cfg, md, rd = _setup(tmp_path)
    _fake_r(monkeypatch)
    m10_kegg.run_kegg(cfg, md, rd)
    assert not (rd / "quality" / "gates.json").exists()   # yeni kapı yok
    assert RunState(rd).is_done("m10_kegg")


def test_run_kegg_resume(tmp_path, monkeypatch):
    cfg, md, rd = _setup(tmp_path)
    _fake_r(monkeypatch)
    m10_kegg.run_kegg(cfg, md, rd)
    s2 = m10_kegg.run_kegg(cfg, md, rd)
    assert s2.get("resumed") is True
