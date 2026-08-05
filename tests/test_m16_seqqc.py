"""m16 orkestrasyon testleri: ön koşul, rRNA/strand kapıları (WARN), stats, resume. Araçlar monkeypatch."""
from __future__ import annotations

import gzip
import json
from pathlib import Path

import pytest

from rnaforge.config import load_config
from rnaforge.modules import m16_seqqc
from rnaforge.modules.m03_trim import trimmed_reads
from rnaforge.metadata import load_metadata
from rnaforge.state import RunState

GENOME = ">chr1\n" + "ACGTACGTAC" * 20 + "\n"
GFF = (
    "chr1\tRefSeq\trRNA\t1\t20\t.\t+\t.\tlocus_tag=r16S;product=16S ribosomal RNA\n"
    "chr1\tRefSeq\tgene\t30\t120\t.\t+\t.\tlocus_tag=LT_1;gene=a\n"
)


def _setup(tmp_path, with_m04=True, strandedness="unstranded"):
    (tmp_path / "genome.fa").write_text(GENOME)
    (tmp_path / "g.gff").write_text(GFF)
    (tmp_path / "c.yaml").write_text(
        "organism: E\norganism_type: prokaryote\nplatform: auto\n"
        f"reference:\n  genome_fasta: {tmp_path/'genome.fa'}\n  annotation_gff: {tmp_path/'g.gff'}\n"
        f"library:\n  strandedness: {strandedness}\n"
        "de:\n  design: '~condition'\n")
    cfg = load_config(tmp_path / "c.yaml")
    # metadata'nın işaret ettiği ham FASTQ'lar (load_metadata varlık doğrular)
    for name in ("c1_1", "c1_2", "t1_1", "t1_2"):
        with gzip.open(tmp_path / f"{name}.fq.gz", "wt") as f:
            f.write("@r\nACGT\n+\nIIII\n")
    md = tmp_path / "m.tsv"
    md.write_text("sample_id\tcondition\tfastq_1\tfastq_2\n"
                  f"c1\tcontrol\t{tmp_path/'c1_1.fq.gz'}\t{tmp_path/'c1_2.fq.gz'}\n"
                  f"t1\ttreated\t{tmp_path/'t1_1.fq.gz'}\t{tmp_path/'t1_2.fq.gz'}\n")
    rd = tmp_path / "run"
    (rd / "statistics").mkdir(parents=True); (rd / "logs").mkdir()
    # trimlenmiş okumalar + BAM'ler (içerik önemsiz; araçlar monkeypatch)
    for s in load_metadata(md):
        t1, t2 = trimmed_reads(rd, s)
        t1.parent.mkdir(parents=True, exist_ok=True)
        for t in (t1, t2):
            with gzip.open(t, "wt") as f:
                f.write("@r\nACGT\n+\nIIII\n")
        bam = rd / "quantification" / s.sample_id / "aligned.sorted.bam"
        bam.parent.mkdir(parents=True, exist_ok=True); bam.write_bytes(b"BAM")
    if with_m04:
        RunState(rd).mark_done("m04_quant", [])
    return cfg, md, rd


def _fake_tools(monkeypatch, rrna_pct=8.5, infer_strand="unstranded"):
    def fake_smr(reads, ref_fasta, workdir, threads=8, env="rnaforge-seqqc"):
        wd = Path(workdir); (wd / "out").mkdir(parents=True, exist_ok=True)
        log = wd / "out" / "aligned.log"
        log.write_text(f"Total reads = 1000\n"
                       f"Total reads passing E-value threshold = {int(rrna_pct*10)} ({rrna_pct:.2f}%)\n")
        return str(log)
    monkeypatch.setattr(m16_seqqc, "run_sortmerna", fake_smr)
    fwd = 0.93 if infer_strand == "stranded" else (0.04 if infer_strand == "reverse" else 0.47)
    rev = 0.93 if infer_strand == "reverse" else (0.04 if infer_strand == "stranded" else 0.48)
    out = (f'Fraction of reads explained by "1++,1--,2+-,2-+": {fwd}\n'
           f'Fraction of reads explained by "1+-,1-+,2++,2--": {rev}\n')
    monkeypatch.setattr(m16_seqqc, "run_infer_experiment", lambda bam, bed, env="rnaforge-seqqc": out)


def test_run_seqqc_requires_m04(tmp_path):
    cfg, md, rd = _setup(tmp_path, with_m04=False)
    with pytest.raises(ValueError, match="m04"):
        m16_seqqc.run_seqqc(cfg, md, rd)


def test_run_seqqc_all_pass(tmp_path, monkeypatch):
    cfg, md, rd = _setup(tmp_path, strandedness="unstranded")
    _fake_tools(monkeypatch, rrna_pct=8.5, infer_strand="unstranded")
    s = m16_seqqc.run_seqqc(cfg, md, rd)
    assert abs(s["mean_rrna_fraction"] - 0.085) < 1e-3
    assert s["inferred_strandedness"] == "unstranded" and s["strandedness_match"] is True
    assert s["gate_counts"].get("WARN", 0) == 0        # rRNA düşük + strand uyumlu -> WARN yok
    assert (rd / "quality" / "gates.json").exists()    # kapılar güvence kartına yazıldı


def test_run_seqqc_rrna_warn(tmp_path, monkeypatch):
    cfg, md, rd = _setup(tmp_path)
    _fake_tools(monkeypatch, rrna_pct=60.0, infer_strand="unstranded")   # %60 > %20 eşik
    s = m16_seqqc.run_seqqc(cfg, md, rd)
    assert s["mean_rrna_fraction"] > 0.20
    assert s["gate_counts"].get("WARN", 0) >= 1        # rRNA yüksek -> WARN


def test_run_seqqc_strand_mismatch_warn(tmp_path, monkeypatch):
    cfg, md, rd = _setup(tmp_path, strandedness="unstranded")
    _fake_tools(monkeypatch, rrna_pct=5.0, infer_strand="stranded")      # çıkarım != beyan
    s = m16_seqqc.run_seqqc(cfg, md, rd)
    assert s["inferred_strandedness"] == "stranded" and s["strandedness_match"] is False
    assert s["gate_counts"].get("WARN", 0) >= 1        # strand uyuşmazlığı -> WARN


def test_run_seqqc_resume(tmp_path, monkeypatch):
    cfg, md, rd = _setup(tmp_path)
    _fake_tools(monkeypatch)
    m16_seqqc.run_seqqc(cfg, md, rd)
    s2 = m16_seqqc.run_seqqc(cfg, md, rd)
    assert s2.get("resumed") is True
