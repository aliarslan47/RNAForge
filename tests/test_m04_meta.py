"""m04 metatranscriptome dalı testi: gen kataloğuna Bowtie2 hizalama (rRNA'sız
okumalar). bowtie2 monkeypatch'lenir (gerçek Bowtie2 gerektirmez). Odak: FAIL kapısı
YOK (permissive metatranscriptome profili — düşük/sıfır hizalama DİAGNOSTİK, asla
GateFailure fırlatmaz), BAM sözleşme yolu, m_rrna_deplete guard, rRNA'sız okuma yoksa
m03 trimmed'e düşme + yüksek sesle log."""
from __future__ import annotations

import json
import textwrap
from pathlib import Path

import pytest

from rnaforge.bowtie2 import AlignmentResult
from rnaforge.config import load_config
from rnaforge.gates import GateFailure
from rnaforge.metadata import load_metadata
from rnaforge.modules import m04_quant
from rnaforge.modules.m01_validate import run_validation
from rnaforge.modules.m03_trim import trimmed_reads
from rnaforge.modules.m04_quant import run_quant
from rnaforge.state import RunState
from tests.conftest import write_fastq


def _setup(tmp_path):
    cat = tmp_path / "catalog.fa"
    cat.write_text(">g1\nACGT\n")
    ann = tmp_path / "catalog.gff"
    ann.write_text("g1\t.\tgene\t1\t4\t.\t+\t.\tID=g1\n")
    for n in ("c1.fastq", "c2.fastq", "t1.fastq", "t2.fastq"):
        write_fastq(tmp_path / n, 50, 100, "I")
    config_path = tmp_path / "config.yaml"
    config_path.write_text(textwrap.dedent(f"""
        organism: "gut community"
        organism_type: "metatranscriptome"
        platform: "illumina"
        reference:
          gene_catalog_fasta: "{cat}"
          catalog_annotation: "{ann}"
        rrna:
          db_fasta: "{cat}"
          env: "rnaforge-seqqc"
    """))
    metadata_path = tmp_path / "samples.tsv"
    metadata_path.write_text(
        "sample_id\tcondition\tfastq_1\n"
        "s1\tcontrol\tc1.fastq\n" "s2\tcontrol\tc2.fastq\n"
        "s3\ttreated\tt1.fastq\n" "s4\ttreated\tt2.fastq\n"
    )
    return config_path, metadata_path


def _mark_m03_done(run_dir, metadata_path):
    """Gerçek run_trim yerine (metatranscriptome profili m03'ün survival_rate kapısı
    için eşik tanımlamıyor — Task 4/6 ile aynı, bu görevin kapsamı dışında bir profil
    boşluğu) sözleşme yoluna doğrudan sahte trimlenmiş FASTQ yazıp state'i işaretliyoruz."""
    samples = load_metadata(metadata_path)
    state = RunState(run_dir)
    for sample in samples:
        out1, out2 = trimmed_reads(run_dir, sample)
        out1.parent.mkdir(parents=True, exist_ok=True)
        out1.write_text("@r\nACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIII\n")
        if out2 is not None:
            out2.write_text("@r\nACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIII\n")
        state.mark_item_done("m03_trim", sample.sample_id, {"survival_rate": 0.98})
    state.mark_done("m03_trim", [])


def _mark_m_rrna_deplete_done(run_dir, metadata_path, write_reads=True):
    """m_rrna_deplete tamamlanmışlığını işaretler; write_reads=False ise sözleşme
    dizini oluşturulmaz (rrna_depleted_reads() boş döner) — fallback testinde kullanılır."""
    samples = load_metadata(metadata_path)
    state = RunState(run_dir)
    for sample in samples:
        if write_reads:
            out_dir = run_dir / "rrna_depleted" / sample.sample_id
            out_dir.mkdir(parents=True, exist_ok=True)
            (out_dir / "other_1.fastq.gz").write_bytes(b"")
        state.mark_item_done("m_rrna_deplete", sample.sample_id, {"depletion_rate": 0.5})
    state.mark_done("m_rrna_deplete", [])


def _prep(config_path, metadata_path, run_dir, write_depleted_reads=True):
    run_validation(load_config(config_path), metadata_path, run_dir)
    _mark_m03_done(run_dir, metadata_path)
    _mark_m_rrna_deplete_done(run_dir, metadata_path, write_reads=write_depleted_reads)


def _fake_bowtie2(monkeypatch, rate=0.95):
    monkeypatch.setattr(m04_quant, "build_index",
                        lambda ref, index_dir, env="rnaforge-quant-prok": Path(index_dir) / "catalog")
    calls = []

    def fake_align(index_prefix, out_dir, fastq_1, fastq_2=None, threads=4,
                   env="rnaforge-quant-prok"):
        calls.append(Path(fastq_1))
        out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
        bam = out_dir / "aligned.sorted.bam"; bam.write_bytes(b"BAM")
        return AlignmentResult(bam=bam, alignment_rate=rate)

    monkeypatch.setattr(m04_quant, "run_bowtie2", fake_align)
    return calls


def test_quant_meta_requires_rrna_deplete(tmp_path, monkeypatch):
    _fake_bowtie2(monkeypatch)
    config_path, metadata_path = _setup(tmp_path)
    run_dir = tmp_path / "run"
    run_validation(load_config(config_path), metadata_path, run_dir)
    _mark_m03_done(run_dir, metadata_path)   # m01+m03 done, m_rrna_deplete NOT done
    with pytest.raises(ValueError, match="m_rrna_deplete"):
        run_quant(load_config(config_path), metadata_path, run_dir)


def test_quant_meta_no_fail_gate_on_low_alignment(tmp_path, monkeypatch):
    """KRİTİK: alignment_rate çok düşük (0.01 < 0.05 eşik) olsa bile GateFailure
    YÜKSELMEMELİ — metatranscriptome profili permissive/diagnostik."""
    config_path, metadata_path = _setup(tmp_path)
    run_dir = tmp_path / "run"
    _prep(config_path, metadata_path, run_dir)
    _fake_bowtie2(monkeypatch, rate=0.01)

    summary = run_quant(load_config(config_path), metadata_path, run_dir)   # must NOT raise

    assert summary["organism_type"] == "metatranscriptome"
    assert summary["read_type"] == "short"
    assert summary["n_samples"] == 4
    assert (run_dir / "quantification" / "s1" / "aligned.sorted.bam").exists()
    assert (run_dir / "quantification" / "s2" / "aligned.sorted.bam").exists()

    stats = json.loads((run_dir / "statistics" / "alignment_statistics.json").read_text())
    assert stats["organism_type"] == "metatranscriptome"
    assert stats["samples"]["s1"]["alignment_rate"] == 0.01

    gates = json.loads((run_dir / "quality" / "gates.json").read_text())["gates"]
    align_gates = [g for g in gates if g["module"] == "m04_quant" and g["name"] == "alignment_rate"]
    assert align_gates
    assert all(g["status"] != "FAIL" for g in align_gates)   # diagnostik: asla FAIL


def test_quant_meta_zero_alignment_never_raises(tmp_path, monkeypatch):
    config_path, metadata_path = _setup(tmp_path)
    run_dir = tmp_path / "run"
    _prep(config_path, metadata_path, run_dir)
    _fake_bowtie2(monkeypatch, rate=0.0)
    try:
        run_quant(load_config(config_path), metadata_path, run_dir)
    except GateFailure:
        pytest.fail("metatranscriptome alignment_rate gate must never FAIL (diagnostic only)")


def test_quant_meta_uses_gene_catalog_not_genome(tmp_path, monkeypatch):
    config_path, metadata_path = _setup(tmp_path)
    run_dir = tmp_path / "run"
    _prep(config_path, metadata_path, run_dir)
    seen_refs = []
    monkeypatch.setattr(m04_quant, "build_index",
                        lambda ref, index_dir, env="rnaforge-quant-prok":
                            (seen_refs.append(Path(ref)), Path(index_dir) / "catalog")[1])

    def fake_align(index_prefix, out_dir, fastq_1, fastq_2=None, threads=4,
                   env="rnaforge-quant-prok"):
        out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
        bam = out_dir / "aligned.sorted.bam"; bam.write_bytes(b"BAM")
        return AlignmentResult(bam=bam, alignment_rate=0.5)
    monkeypatch.setattr(m04_quant, "run_bowtie2", fake_align)

    config = load_config(config_path)
    run_quant(config, metadata_path, run_dir)
    assert seen_refs == [config.reference.gene_catalog_fasta]


def test_quant_meta_uses_rrna_depleted_reads(tmp_path, monkeypatch):
    """Girdi rRNA'sız okumalar olmalı (trimlenmiş DEĞİL) — depleted dosyalar mevcutken."""
    config_path, metadata_path = _setup(tmp_path)
    run_dir = tmp_path / "run"
    _prep(config_path, metadata_path, run_dir, write_depleted_reads=True)
    calls = _fake_bowtie2(monkeypatch, rate=0.5)
    run_quant(load_config(config_path), metadata_path, run_dir)
    assert all("rrna_depleted" in str(p) for p in calls)


def test_quant_meta_falls_back_to_trimmed_when_depleted_missing_and_logs_loudly(tmp_path, monkeypatch):
    """rRNA'sız okuma bulunamazsa (dosya yok) m03 trimmed'e düş, sessizce ham veriye DEĞİL,
    ve bunu logs/quant.log'a yüksek sesle yaz."""
    config_path, metadata_path = _setup(tmp_path)
    run_dir = tmp_path / "run"
    _prep(config_path, metadata_path, run_dir, write_depleted_reads=False)
    calls = _fake_bowtie2(monkeypatch, rate=0.5)
    run_quant(load_config(config_path), metadata_path, run_dir)
    assert all("trimmed" in str(p) for p in calls)
    log_text = (run_dir / "logs" / "quant.log").read_text()
    assert "s1" in log_text and ("trimmed" in log_text.lower() or "yedek" in log_text.lower()
                                  or "fallback" in log_text.lower() or "düş" in log_text.lower())


def test_quant_meta_resumes_cached_sample(tmp_path, monkeypatch):
    config_path, metadata_path = _setup(tmp_path)
    run_dir = tmp_path / "run"
    _prep(config_path, metadata_path, run_dir)
    _fake_bowtie2(monkeypatch, rate=0.3)
    run_quant(load_config(config_path), metadata_path, run_dir)

    calls = []
    monkeypatch.setattr(m04_quant, "run_bowtie2", lambda *a, **k: calls.append(1))
    summary = run_quant(load_config(config_path), metadata_path, run_dir)
    assert summary.get("resumed") is True
    assert calls == []
