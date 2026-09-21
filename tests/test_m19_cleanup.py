"""m19 temizlik testleri: ara dosya (trimmed + BAM) silme, sonuç koruma, keep_bam,
remove_intermediates=false no-op, resume. Gerçek dosya sistemi (araç yok)."""
from __future__ import annotations

import json
from pathlib import Path

from rnaforge.config import load_config
from rnaforge.modules import m19_cleanup
from rnaforge.state import RunState


def _cfg(tmp_path, remove=True, keep_bam=False):
    (tmp_path / "genome.fa").write_text(">chr1\nACGT\n")
    (tmp_path / "g.gff").write_text("chr1\tX\tgene\t1\t4\t.\t+\t.\tlocus_tag=a\n")
    body = (
        "organism: E\norganism_type: prokaryote\nplatform: auto\n"
        f"reference:\n  genome_fasta: {tmp_path/'genome.fa'}\n  annotation_gff: {tmp_path/'g.gff'}\n"
        "de:\n  design: '~condition'\n"
    )
    if remove is not None:
        body += "cleanup:\n"
        body += f"  remove_intermediates: {'true' if remove else 'false'}\n"
        body += f"  keep_bam: {'true' if keep_bam else 'false'}\n"
    (tmp_path / "c.yaml").write_text(body)
    return load_config(tmp_path / "c.yaml")


def _seed_run(tmp_path) -> Path:
    """trimmed/ + quantification/ (BAM + counts.tsv) olan bir koşu dizini kur."""
    rd = tmp_path / "run"
    for sid in ("c1", "t1"):
        td = rd / "trimmed" / sid
        td.mkdir(parents=True)
        (td / f"{sid}.trimmed.fastq").write_text("@r\nACGT\n+\nIIII\n" * 1000)
        qd = rd / "quantification" / sid
        qd.mkdir(parents=True)
        (qd / "aligned.sorted.bam").write_bytes(b"BAMDATA" * 1000)
        (qd / "aligned.sorted.bam.bai").write_bytes(b"BAI" * 100)
    # korunması GEREKEN sonuçlar
    (rd / "quantification" / "counts.tsv").write_text("gene\tc1\tt1\na\t5\t9\n")
    (rd / "statistics").mkdir(parents=True)
    (rd / "differential_expression").mkdir()
    (rd / "differential_expression" / "deseq2_results.tsv").write_text("gene\tlog2FC\na\t2.1\n")
    return rd


def test_cleanup_removes_trimmed_and_bam_keeps_results(tmp_path):
    cfg = _cfg(tmp_path, remove=True)
    rd = _seed_run(tmp_path)
    summary = m19_cleanup.run_cleanup(cfg, tmp_path / "m.tsv", rd)

    assert summary["removed"] is True
    assert summary["freed_bytes"] > 0
    # ara dosyalar GİTTİ
    assert not (rd / "trimmed").exists()
    assert not (rd / "quantification" / "c1" / "aligned.sorted.bam").exists()
    assert not (rd / "quantification" / "t1" / "aligned.sorted.bam.bai").exists()
    # sonuçlar DURUYOR
    assert (rd / "quantification" / "counts.tsv").exists()
    assert (rd / "differential_expression" / "deseq2_results.tsv").exists()
    # stats + state
    stats = json.loads((rd / "statistics" / "cleanup_statistics.json").read_text())
    assert stats["freed_bytes"] == summary["freed_bytes"]
    assert RunState(rd).is_done("m19_cleanup")


def test_cleanup_keep_bam_deletes_only_trimmed(tmp_path):
    cfg = _cfg(tmp_path, remove=True, keep_bam=True)
    rd = _seed_run(tmp_path)
    summary = m19_cleanup.run_cleanup(cfg, tmp_path / "m.tsv", rd)

    assert summary["removed"] is True
    assert summary["keep_bam"] is True
    assert not (rd / "trimmed").exists()                                   # trimmed silindi
    assert (rd / "quantification" / "c1" / "aligned.sorted.bam").exists()  # BAM KORUNDU


def test_cleanup_disabled_is_noop(tmp_path):
    cfg = _cfg(tmp_path, remove=False)
    rd = _seed_run(tmp_path)
    summary = m19_cleanup.run_cleanup(cfg, tmp_path / "m.tsv", rd)

    assert summary["removed"] is False
    assert summary["freed_bytes"] == 0
    # HİÇBİR ŞEY silinmedi
    assert (rd / "trimmed").exists()
    assert (rd / "quantification" / "c1" / "aligned.sorted.bam").exists()


def test_cleanup_resume_is_idempotent(tmp_path):
    cfg = _cfg(tmp_path, remove=True)
    rd = _seed_run(tmp_path)
    first = m19_cleanup.run_cleanup(cfg, tmp_path / "m.tsv", rd)
    assert "resumed" not in first
    # ikinci çağrı: state done → resume, tekrar silmeye çalışmaz (dosyalar zaten yok)
    second = m19_cleanup.run_cleanup(cfg, tmp_path / "m.tsv", rd)
    assert second.get("resumed") is True


def test_cleanup_default_enabled_when_no_config_section(tmp_path):
    """cleanup bölümü olmayan config → remove_intermediates VARSAYILAN True (otomatik)."""
    cfg = _cfg(tmp_path, remove=None)  # cleanup bölümü yazılmaz
    assert cfg.cleanup.remove_intermediates is True
    assert cfg.cleanup.keep_bam is False
