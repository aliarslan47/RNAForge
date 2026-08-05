"""m09 orkestrasyon testleri: ön koşul, gürültülü hata, gate-yokluğu, çıktı, resume."""
from __future__ import annotations

import textwrap
from pathlib import Path

import pytest

from rnaforge.config import load_config
from rnaforge.modules import m09_enrichment
from rnaforge.state import RunState

def _cds(gene, lt, go=None):
    ann = f"Ontology_term={go};go_process=t|{go[3:]}||IEA;" if go else ""
    return f"NZ\tX\tCDS\t1\t9\t.\t+\t0\t{ann}gene={gene};locus_tag={lt}\n"


# GO:0000002 ("mid"): LT_A..D (DEG'ler burada). GO:0000005 ("other"): LT_F..L (ns, arka planı büyütür).
GFF = (
    _cds("gA", "LT_A", "GO:0000002") + _cds("gB", "LT_B", "GO:0000002")
    + _cds("gC", "LT_C", "GO:0000002") + _cds("gD", "LT_D", "GO:0000002")
    + _cds("gE", "LT_E")
    + "".join(_cds(f"g{i}", f"LT_{i}", "GO:0000005") for i in range(6, 13))
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

    [Term]
    id: GO:0000005
    name: other
    namespace: biological_process
    is_a: GO:0000001 ! root
    """)
DESEQ = (
    "gene\tbaseMean\tlog2FoldChange\tlfcSE\tstat\tpvalue\tpadj\n"
    "LT_A\t100\t3.0\t0.2\t9\t1e-10\t1e-9\n"     # up
    "LT_B\t100\t2.5\t0.2\t8\t1e-9\t1e-8\n"      # up
    "LT_C\t100\t2.2\t0.2\t7\t1e-8\t1e-7\n"      # up
    "LT_D\t100\t-3.0\t0.2\t-9\t1e-10\t1e-9\n"   # down
    "LT_E\t100\t0.1\t0.2\t0.3\t0.7\t0.8\n"      # ns
    + "".join(f"LT_{i}\t100\t0.1\t0.2\t0.3\t0.7\t0.8\n" for i in range(6, 13))  # ns arka plan
)


def _setup(tmp_path, with_m06=True):
    gff = tmp_path / "g.gff"; gff.write_text(GFF)
    obo = tmp_path / "go.obo"; obo.write_text(OBO)
    (tmp_path / "c.yaml").write_text(
        "organism: E\norganism_type: prokaryote\nplatform: auto\n"
        f"reference:\n  genome_fasta: g.fa\n  annotation_gff: {gff}\n"
        "de:\n  design: '~condition'\n  fdr_threshold: 0.05\n  log2fc_threshold: 1.0\n"
        f"enrichment:\n  min_term_size: 3\n  obo: {obo}\n")
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


def _fake_r(monkeypatch, rd):
    def fake(up_tsv, down_tsv, out_dir, top_n, env="rnaforge-de"):
        Path(out_dir).mkdir(parents=True, exist_ok=True)
        for base in ("enrichment_up", "enrichment_down"):
            (Path(out_dir) / f"{base}.png").write_bytes(b"x" * 2000)
            (Path(out_dir) / f"{base}.svg").write_text("<svg/>")
        return "enrichment.R done\n"
    monkeypatch.setattr(m09_enrichment, "run_enrichment_r", fake)


def test_run_enrichment_requires_m06(tmp_path):
    cfg, md, rd = _setup(tmp_path, with_m06=False)
    with pytest.raises(ValueError, match="m06"):
        m09_enrichment.run_enrichment(cfg, md, rd)


def test_run_enrichment_missing_obo_loud(tmp_path, monkeypatch):
    cfg, md, rd = _setup(tmp_path)
    # obo yolunu var olmayan bir dosyaya çevir -> FileNotFoundError (sessiz skip değil)
    import dataclasses
    bad = dataclasses.replace(cfg.enrichment, obo=Path(tmp_path / "yok.obo"))
    cfg = dataclasses.replace(cfg, enrichment=bad)
    with pytest.raises(FileNotFoundError, match="go-basic.obo"):
        m09_enrichment.run_enrichment(cfg, md, rd)


def test_run_enrichment_writes_outputs_and_stats(tmp_path, monkeypatch):
    cfg, md, rd = _setup(tmp_path)
    _fake_r(monkeypatch, rd)
    s = m09_enrichment.run_enrichment(cfg, md, rd)
    assert (rd / "enrichment" / "enrichment_up.tsv").exists()
    assert (rd / "enrichment" / "enrichment_down.tsv").exists()
    assert (rd / "enrichment" / "gene2go.tsv").exists()
    assert (rd / "enrichment" / "manifest.json").exists()
    assert (rd / "statistics" / "enrichment_statistics.json").exists()
    assert s["n_up_degs"] == 3 and s["n_down_degs"] == 1
    assert s["n_annotated"] == 11         # LT_A..D (mid) + LT_6..12 (other); LT_E GO'suz
    # up seti (3 gen, hepsi GO:0000002) -> terim anlamlı zenginleşir
    assert s["n_sig_up"] >= 1


def test_run_enrichment_no_gate_state_clean(tmp_path, monkeypatch):
    cfg, md, rd = _setup(tmp_path)
    _fake_r(monkeypatch, rd)
    m09_enrichment.run_enrichment(cfg, md, rd)
    # m09 gates.json'a DOKUNMAZ (yeni kapı yok) — dosya oluşturmamalı.
    assert not (rd / "quality" / "gates.json").exists()
    assert RunState(rd).is_done("m09_enrichment")


def test_run_enrichment_resume(tmp_path, monkeypatch):
    cfg, md, rd = _setup(tmp_path)
    _fake_r(monkeypatch, rd)
    m09_enrichment.run_enrichment(cfg, md, rd)
    s2 = m09_enrichment.run_enrichment(cfg, md, rd)
    assert s2.get("resumed") is True
