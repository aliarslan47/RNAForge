from __future__ import annotations

import json
import textwrap
from pathlib import Path

import pytest

from rnaforge.bowtie2 import AlignmentResult
from rnaforge.config import load_config
from rnaforge.gates import FAIL, PASS, GateFailure
from rnaforge.modules import m04_quant
from rnaforge.modules.m04_quant import build_alignment_gates, run_quant
from rnaforge.quality import load_profile
from tests.conftest import write_fastq


def _res(rate: float) -> AlignmentResult:
    return AlignmentResult(bam=Path("x.bam"), alignment_rate=rate)


def test_all_above_threshold_passes():
    profile = load_profile("prokaryote")  # alignment_rate = 0.70
    gates = build_alignment_gates({"s1": _res(0.98), "s2": _res(0.85)}, profile)
    assert len(gates) == 1
    assert gates[0].name == "alignment_rate"
    assert gates[0].module == "m04_quant"
    assert gates[0].status == PASS


def test_below_threshold_fails_and_lists_offenders():
    profile = load_profile("prokaryote")
    gates = build_alignment_gates({"s1": _res(0.98), "s2": _res(0.30)}, profile)
    g = gates[0]
    assert g.status == FAIL
    assert g.samples == ("s2",)
    assert g.measured == 0.30
    assert g.threshold == 0.70


def test_override_marks_gate_overridden():
    profile = load_profile("prokaryote", {"alignment_rate": 0.20})
    gates = build_alignment_gates({"s1": _res(0.30)}, profile)
    assert gates[0].status == PASS
    assert gates[0].overridden is True
    assert gates[0].threshold == 0.20


def _setup(tmp_path, organism_type="prokaryote"):
    (tmp_path / "ref").mkdir()
    (tmp_path / "ref" / "genome.fa").write_text(">c1\n" + "ACGT" * 25 + "\n")
    (tmp_path / "ref" / "genes.gff").write_text("##gff-version 3\n")
    ref = ("genome_fasta" if organism_type == "prokaryote" else "transcriptome_fasta")
    extra = ("annotation_gff" if organism_type == "prokaryote" else "tx2gene")
    for n in ("c1.fastq", "c2.fastq", "t1.fastq", "t2.fastq"):
        write_fastq(tmp_path / n, 200, 150, "I")
    config_path = tmp_path / "config.yaml"
    config_path.write_text(textwrap.dedent(f"""
        organism: "E. coli"
        organism_type: "{organism_type}"
        reference:
          {ref}: "{tmp_path / 'ref' / 'genome.fa'}"
          {extra}: "{tmp_path / 'ref' / 'genes.gff'}"
    """))
    metadata_path = tmp_path / "samples.tsv"
    metadata_path.write_text(
        "sample_id\tcondition\tfastq_1\n"
        "s1\tcontrol\tc1.fastq\n" "s2\tcontrol\tc2.fastq\n"
        "s3\ttreated\tt1.fastq\n" "s4\ttreated\tt2.fastq\n"
    )
    return config_path, metadata_path


def _prep_m01_m03(config_path, metadata_path, run_dir, monkeypatch, survival=0.98):
    """m01 (gerçek) + m03 (fastp monkeypatch) hazırlar — m04 için trimlenmiş okuma
    ve done state gerekir."""
    from rnaforge.modules.m01_validate import run_validation
    from rnaforge.modules import m03_trim
    from rnaforge.modules.m03_trim import run_trim, trimmed_name
    from rnaforge.fastp import FastpResult
    run_validation(load_config(config_path), metadata_path, run_dir)

    def fake_fastp(fastq_1, out_dir, min_length, fastq_2=None,
                   aggressive_quality=False, env="rnaforge-qc"):
        out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
        out1 = out_dir / trimmed_name(Path(fastq_1))
        out1.write_text("@r\nACGT\n+\nIIII\n")
        (out_dir / "fastp.json").write_text("{}")
        return FastpResult(reads_before=200, reads_after=int(200 * survival),
                           survival_rate=survival, out1=out1)
    monkeypatch.setattr(m03_trim, "run_fastp", fake_fastp)
    run_trim(load_config(config_path), metadata_path, run_dir)


def _fake_bowtie2(monkeypatch, rate=0.95):
    monkeypatch.setattr(m04_quant, "build_index",
                        lambda genome, index_dir, env="rnaforge-quant-prok": Path(index_dir) / "genome")
    def fake_align(index_prefix, out_dir, fastq_1, fastq_2=None, threads=4,
                   env="rnaforge-quant-prok"):
        out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
        bam = out_dir / "aligned.sorted.bam"; bam.write_bytes(b"BAM")
        return AlignmentResult(bam=bam, alignment_rate=rate)
    monkeypatch.setattr(m04_quant, "run_bowtie2", fake_align)


def test_run_quant_requires_m03_done(tmp_path, monkeypatch):
    _fake_bowtie2(monkeypatch)
    config_path, metadata_path = _setup(tmp_path)
    from rnaforge.modules.m01_validate import run_validation
    run_dir = tmp_path / "run"
    run_validation(load_config(config_path), metadata_path, run_dir)  # yalniz m01
    with pytest.raises(ValueError, match="m03"):
        run_quant(load_config(config_path), metadata_path, run_dir)


def test_run_quant_writes_bam_and_passes(tmp_path, monkeypatch):
    config_path, metadata_path = _setup(tmp_path)
    run_dir = tmp_path / "run"
    _prep_m01_m03(config_path, metadata_path, run_dir, monkeypatch)
    _fake_bowtie2(monkeypatch, rate=0.95)
    summary = run_quant(load_config(config_path), metadata_path, run_dir)

    assert summary["n_samples"] == 4
    assert summary["read_type"] == "short"
    assert (run_dir / "quantification" / "s1" / "aligned.sorted.bam").exists()
    stats = json.loads((run_dir / "statistics" / "alignment_statistics.json").read_text())
    assert set(stats["samples"]) == {"s1", "s2", "s3", "s4"}
    gates = json.loads((run_dir / "quality" / "gates.json").read_text())["gates"]
    assert any(g["module"] == "m04_quant" and g["status"] == "PASS" for g in gates)
    assert any(g["module"] == "m03_trim" for g in gates)   # onceki kapilar korundu


def test_run_quant_low_alignment_fails(tmp_path, monkeypatch):
    config_path, metadata_path = _setup(tmp_path)
    run_dir = tmp_path / "run"
    _prep_m01_m03(config_path, metadata_path, run_dir, monkeypatch)
    _fake_bowtie2(monkeypatch, rate=0.10)
    with pytest.raises(GateFailure):
        run_quant(load_config(config_path), metadata_path, run_dir)
    gates = json.loads((run_dir / "quality" / "gates.json").read_text())["gates"]
    assert any(g["module"] == "m04_quant" and g["status"] == "FAIL" for g in gates)


def test_run_quant_resumes(tmp_path, monkeypatch):
    config_path, metadata_path = _setup(tmp_path)
    run_dir = tmp_path / "run"
    _prep_m01_m03(config_path, metadata_path, run_dir, monkeypatch)
    _fake_bowtie2(monkeypatch, rate=0.95)
    run_quant(load_config(config_path), metadata_path, run_dir)
    calls = []
    monkeypatch.setattr(m04_quant, "run_bowtie2", lambda *a, **k: calls.append(1))
    summary = run_quant(load_config(config_path), metadata_path, run_dir)
    assert summary.get("resumed") is True
    assert calls == []


def _seed_long(tmp_path, platform="ont", chemistry="cdna"):
    """Uzun-okuma run dizini doğrudan tohumlanır (m01+m03 done, raw_statistics)."""
    from rnaforge.state import RunState
    run_dir = tmp_path / "run"
    (run_dir / "statistics").mkdir(parents=True)
    (run_dir / "statistics" / "raw_statistics.json").write_text(
        json.dumps({"platform": platform, "read_type": "long", "chemistry": chemistry})
    )
    st = RunState(run_dir)
    st.mark_done("m01_validate", [])
    st.mark_done("m03_trim", [])
    return run_dir


def _long_config(tmp_path, chemistry="cdna"):
    (tmp_path / "ref").mkdir(exist_ok=True)
    (tmp_path / "ref" / "genome.fa").write_text(">c1\n" + "ACGT" * 25 + "\n")
    (tmp_path / "ref" / "genes.gff").write_text("##gff-version 3\n")
    for n in ("s1.fastq", "s2.fastq"):
        write_fastq(tmp_path / n, 50, 600, "I")   # long ONT-like reads
    config_path = tmp_path / "config.yaml"
    config_path.write_text(textwrap.dedent(f"""
        organism: "E. coli"
        organism_type: "prokaryote"
        library:
          chemistry: "{chemistry}"
        reference:
          genome_fasta: "{tmp_path / 'ref' / 'genome.fa'}"
          annotation_gff: "{tmp_path / 'ref' / 'genes.gff'}"
    """))
    metadata_path = tmp_path / "samples.tsv"
    metadata_path.write_text(
        "sample_id\tcondition\tfastq_1\n"
        "s1\tcontrol\ts1.fastq\n" "s2\ttreated\ts2.fastq\n"
    )
    return config_path, metadata_path


def _fake_minimap2(monkeypatch, rate=0.92):
    from rnaforge.minimap2 import AlignmentResult
    calls = []

    def fake_align(genome_fasta, out_dir, fastq, preset, threads=4,
                   env="rnaforge-longread"):
        calls.append(preset)
        out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
        bam = out_dir / "aligned.sorted.bam"; bam.write_bytes(b"BAM")
        return AlignmentResult(bam=bam, alignment_rate=rate)

    monkeypatch.setattr(m04_quant, "run_minimap2", fake_align)
    return calls


def test_run_quant_long_dispatches_to_minimap2(tmp_path, monkeypatch):
    config_path, metadata_path = _long_config(tmp_path, "cdna")
    run_dir = _seed_long(tmp_path, platform="ont")
    calls = _fake_minimap2(monkeypatch, rate=0.92)

    def boom(*a, **k):
        raise AssertionError("bowtie2 must not run on the long branch")
    monkeypatch.setattr(m04_quant, "run_bowtie2", boom)

    summary = run_quant(load_config(config_path), metadata_path, run_dir)
    assert summary["read_type"] == "long"
    assert summary["platform"] == "ont"
    assert calls == ["map-ont", "map-ont"]        # preset per sample, from platform
    assert summary["samples"]["s1"]["alignment_rate"] == 0.92
    assert (run_dir / "quantification" / "s1" / "aligned.sorted.bam").exists()
    # Step 6: uzun-okuma alignment FAIL kapısı yazılır (prokaryote_long, eşik 0.50);
    # 0.92 > 0.50 → PASS. Kart prokaryote_long (permissive) damgalanır.
    gates = json.loads((run_dir / "quality" / "gates.json").read_text())["gates"]
    align = [g for g in gates if g["module"] == "m04_quant" and g["name"] == "alignment_rate"]
    assert align and align[0]["status"] == "PASS"
    stats = json.loads((run_dir / "statistics" / "alignment_statistics.json").read_text())
    assert stats["read_type"] == "long"


def test_run_quant_long_catastrophic_alignment_fails(tmp_path, monkeypatch):
    """Yanlış referans → çok düşük hizalama long profilde de FAIL (sonuç geçersiz)."""
    from rnaforge.gates import GateFailure
    config_path, metadata_path = _long_config(tmp_path, "cdna")
    run_dir = _seed_long(tmp_path, platform="ont")
    _fake_minimap2(monkeypatch, rate=0.10)      # 0.10 < 0.50 long floor
    with pytest.raises(GateFailure):
        run_quant(load_config(config_path), metadata_path, run_dir)
    gates = json.loads((run_dir / "quality" / "gates.json").read_text())["gates"]
    assert any(g["module"] == "m04_quant" and g["status"] == "FAIL" for g in gates)


def test_run_quant_long_pacbio_uses_hifi_preset(tmp_path, monkeypatch):
    config_path, metadata_path = _long_config(tmp_path, "cdna")
    run_dir = _seed_long(tmp_path, platform="pacbio_hifi")
    calls = _fake_minimap2(monkeypatch)
    summary = run_quant(load_config(config_path), metadata_path, run_dir)
    assert calls == ["map-hifi", "map-hifi"]
    assert summary["platform"] == "pacbio_hifi"


def test_cli_quant_returns_zero_and_prints_verdict(tmp_path, monkeypatch, capsys):
    from rnaforge.cli import main
    from rnaforge.modules import m03_trim
    from rnaforge.modules.m03_trim import trimmed_name
    from rnaforge.fastp import FastpResult
    config_path, metadata_path = _setup(tmp_path)
    common = ["--config", str(config_path), "--metadata", str(metadata_path),
              "--runs-dir", str(tmp_path / "runs"), "--run-id", "demo"]

    def fake_fastp(fastq_1, out_dir, min_length, fastq_2=None,
                   aggressive_quality=False, env="rnaforge-qc"):
        out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
        out1 = out_dir / trimmed_name(Path(fastq_1)); out1.write_text("@r\nACGT\n+\nIIII\n")
        (out_dir / "fastp.json").write_text("{}")
        return FastpResult(200, int(200 * 0.98), 0.98, out1=out1)
    monkeypatch.setattr(m03_trim, "run_fastp", fake_fastp)
    _fake_bowtie2(monkeypatch, rate=0.95)

    assert main(["validate", *common]) == 0
    assert main(["trim", *common]) == 0
    capsys.readouterr()
    assert main(["quant", *common]) == 0
    assert "quality verdict" in capsys.readouterr().out


def test_run_quant_eukaryote_runs_salmon_and_gates(tmp_path, monkeypatch):
    import rnaforge.modules.m04_quant as m04
    from rnaforge.salmon import SalmonQuant
    from rnaforge.state import RunState
    run_dir = tmp_path / "run"; (run_dir / "statistics").mkdir(parents=True)
    (run_dir / "statistics" / "raw_statistics.json").write_text(
        '{"read_type":"short","platform":"illumina"}')
    st = RunState(run_dir); st.mark_done("m01_validate", []); st.mark_done("m03_trim", [])
    fq = tmp_path / "s1.fastq"; fq.write_text("@r\nACGT\n+\nIIII\n")
    meta = tmp_path / "m.tsv"; meta.write_text(f"sample_id\tcondition\tfastq_1\ns1\tctrl\t{fq}\n")
    monkeypatch.setattr(m04, "trimmed_reads", lambda rd, s: (fq, None))
    monkeypatch.setattr(m04, "build_salmon_index", lambda *a, **k: tmp_path / "idx")
    monkeypatch.setattr(m04, "run_salmon_quant",
        lambda *a, **k: SalmonQuant(quant_sf=tmp_path / "q.sf", mapping_rate=0.9))
    from rnaforge.config import (Config, Reference, Library, Trimming, DE, Report, Resources)
    cfg = Config(organism="human", organism_type="eukaryote", platform="illumina",
        reference=Reference(transcriptome_fasta=tmp_path / "tx.fa", tx2gene=tmp_path / "t2g.tsv"),
        library=Library(), trimming=Trimming(), de=DE(), report=Report(), resources=Resources())
    summary = m04.run_quant(cfg, meta, run_dir)
    assert summary["organism_type"] == "eukaryote"
    assert summary["samples"]["s1"]["mapping_rate"] == 0.9
    gates = json.loads((run_dir / "quality" / "gates.json").read_text())["gates"]
    assert any(g["name"] == "alignment_rate" for g in gates if g["module"] == "m04_quant")


def test_run_quant_eukaryote_long_minimap2_diagnostic(tmp_path, monkeypatch):
    import rnaforge.modules.m04_quant as m04
    from rnaforge.bowtie2 import AlignmentResult
    from rnaforge.state import RunState
    run_dir = tmp_path / "run"; (run_dir / "statistics").mkdir(parents=True)
    (run_dir / "statistics" / "raw_statistics.json").write_text(
        '{"read_type":"long","platform":"ont","chemistry":"cdna"}')
    st = RunState(run_dir); st.mark_done("m01_validate", []); st.mark_done("m03_trim", [])
    fq = tmp_path / "s1.fastq"; fq.write_text("@r\nACGT\n+\nIIII\n")
    meta = tmp_path / "m.tsv"; meta.write_text(f"sample_id\tcondition\tfastq_1\ns1\tctrl\t{fq}\n")
    monkeypatch.setattr(m04, "trimmed_reads", lambda rd, s: (fq, None))
    bam = tmp_path / "b.bam"; bam.write_text("")
    monkeypatch.setattr(m04, "run_minimap2",
        lambda *a, **k: AlignmentResult(alignment_rate=0.55, bam=bam))
    from rnaforge.config import (Config, Reference, Library, Trimming, DE, Report, Resources)
    cfg = Config(organism="human", organism_type="eukaryote", platform="ont",
        reference=Reference(transcriptome_fasta=tmp_path / "tx.fa", tx2gene=tmp_path / "t2g.tsv"),
        library=Library(chemistry="cdna"), trimming=Trimming(), de=DE(), report=Report(),
        resources=Resources())
    summary = m04.run_quant(cfg, meta, run_dir)
    assert summary["organism_type"] == "eukaryote" and summary["read_type"] == "long"
    assert summary["samples"]["s1"]["mapping_rate"] == 0.55
    gpath = run_dir / "quality" / "gates.json"
    if gpath.exists():
        gates = json.loads(gpath.read_text())["gates"]
        assert all(g.get("status") != "FAIL" for g in gates)   # diagnostik, FAIL yok


def test_run_quant_skips_already_done_samples(tmp_path, monkeypatch):
    """Faz 3 örnek-başı resume: önceden hizalanmış örnekler (işaretçi + BAM var) yeniden
    bowtie2'ye SOKULMAMALI; özet + kapı yine tüm örnekleri kapsar."""
    from rnaforge.state import RunState
    from rnaforge.metadata import load_metadata
    from rnaforge.modules.m04_quant import MODULE_NAME

    config_path, metadata_path = _setup(tmp_path)
    run_dir = tmp_path / "run"
    _prep_m01_m03(config_path, metadata_path, run_dir, monkeypatch)

    monkeypatch.setattr(m04_quant, "build_index",
                        lambda genome, index_dir, env="rnaforge-quant-prok": Path(index_dir) / "genome")
    calls = []

    def fake_align(index_prefix, out_dir, fastq_1, fastq_2=None, threads=4,
                   env="rnaforge-quant-prok"):
        out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
        calls.append(out_dir.name)
        bam = out_dir / "aligned.sorted.bam"; bam.write_bytes(b"BAM")
        return AlignmentResult(bam=bam, alignment_rate=0.95)

    monkeypatch.setattr(m04_quant, "run_bowtie2", fake_align)

    # s1, s2 hizalanmış gibi tohumla: BAM + işaretçi + payload; modül 'done' DEĞİL.
    samples = load_metadata(metadata_path)
    state = RunState(run_dir)
    for s in samples[:2]:
        bam = run_dir / "quantification" / s.sample_id / "aligned.sorted.bam"
        bam.parent.mkdir(parents=True, exist_ok=True)
        bam.write_bytes(b"BAM")
        state.mark_item_done(MODULE_NAME, s.sample_id,
                             {"alignment_rate": 0.91, "bam": str(bam)})

    summary = run_quant(load_config(config_path), metadata_path, run_dir)
    assert sorted(calls) == ["s3", "s4"]                        # yalnız kalanlar hizalandı
    assert set(summary["samples"]) == {"s1", "s2", "s3", "s4"}
    assert summary["samples"]["s1"]["alignment_rate"] == 0.91   # cached değer korundu
    gates = json.loads((run_dir / "quality" / "gates.json").read_text())["gates"]
    assert any(g["module"] == "m04_quant" and g["status"] == "PASS" for g in gates)
