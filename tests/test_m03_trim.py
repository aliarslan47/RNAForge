from __future__ import annotations

import json
import textwrap
from pathlib import Path

import pytest

from rnaforge.config import load_config
from rnaforge.fastp import FastpResult
from rnaforge.gates import FAIL, PASS, WARN, GateFailure
from rnaforge.modules import m03_trim
from rnaforge.modules.m03_trim import build_trim_gates, run_trim
from rnaforge.quality import load_profile
from tests.conftest import write_fastq


def test_all_above_threshold_passes():
    profile = load_profile("prokaryote")   # survival_rate = 0.50
    gates = build_trim_gates({"s1": 0.98, "s2": 0.90}, profile)
    assert len(gates) == 1
    g = gates[0]
    assert g.name == "survival_rate"
    assert g.module == "m03_trim"
    assert g.status == PASS


def test_below_threshold_fails_and_lists_offenders():
    profile = load_profile("prokaryote")
    gates = build_trim_gates({"s1": 0.98, "s2": 0.20}, profile)
    g = gates[0]
    assert g.status == FAIL
    assert g.samples == ("s2",)
    assert g.measured == 0.20          # en düşük survival
    assert g.threshold == 0.50


def test_below_threshold_warn_only_is_warn_not_fail():
    """Uzun-okuma: düşük survival ŞÜPHELİ (WARN), geçersiz DEĞİL (Pychopper doğal düşük)."""
    profile = load_profile("prokaryote_long")   # survival_rate = 0.20
    gates = build_trim_gates({"s1": 0.10}, profile, warn_only=True)
    assert gates[0].status == WARN
    assert gates[0].samples == ("s1",)


def test_threshold_override_marks_gate_overridden():
    profile = load_profile("prokaryote", {"survival_rate": 0.10})
    gates = build_trim_gates({"s1": 0.20}, profile)
    g = gates[0]
    assert g.status == PASS             # 0.20 >= ezilmiş 0.10
    assert g.overridden is True
    assert g.threshold == 0.10


def _setup(tmp_path):
    (tmp_path / "ref").mkdir()
    (tmp_path / "ref" / "genome.fa").write_text(">c1\nACGT\n")
    (tmp_path / "ref" / "genes.gff").write_text("##gff-version 3\n")
    for n in ("c1.fastq", "c2.fastq", "t1.fastq", "t2.fastq"):
        write_fastq(tmp_path / n, 200, 150, "I")
    config_path = tmp_path / "config.yaml"
    config_path.write_text(textwrap.dedent(f"""
        organism: "Escherichia coli"
        organism_type: "prokaryote"
        reference:
          genome_fasta: "{tmp_path / 'ref' / 'genome.fa'}"
          annotation_gff: "{tmp_path / 'ref' / 'genes.gff'}"
    """))
    metadata_path = tmp_path / "samples.tsv"
    metadata_path.write_text(
        "sample_id\tcondition\tfastq_1\n"
        "s1\tcontrol\tc1.fastq\n" "s2\tcontrol\tc2.fastq\n"
        "s3\ttreated\tt1.fastq\n" "s4\ttreated\tt2.fastq\n"
    )
    return config_path, metadata_path


def _run_m01(config_path, metadata_path, run_dir):
    from rnaforge.modules.m01_validate import run_validation
    run_validation(load_config(config_path), metadata_path, run_dir)


def _fake_fastp(monkeypatch, survival=0.98):
    def fake_run(fastq_1, out_dir, min_length, fastq_2=None,
                 aggressive_quality=False, env="rnaforge-qc"):
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        out1 = out_dir / (Path(fastq_1).stem + ".trimmed.fastq")
        out1.write_text("@r\nACGT\n+\nIIII\n")
        (out_dir / "fastp.json").write_text("{}")
        (out_dir / "fastp.html").write_text("<html></html>")
        return FastpResult(reads_before=200, reads_after=int(200 * survival),
                           survival_rate=survival, out1=out1)
    monkeypatch.setattr(m03_trim, "run_fastp", fake_run)


def test_run_trim_requires_m01_done(tmp_path, monkeypatch):
    _fake_fastp(monkeypatch)
    config_path, metadata_path = _setup(tmp_path)
    with pytest.raises(ValueError, match="m01"):
        run_trim(load_config(config_path), metadata_path, tmp_path / "run")


def test_run_trim_writes_outputs_and_passes(tmp_path, monkeypatch):
    _fake_fastp(monkeypatch, survival=0.98)
    config_path, metadata_path = _setup(tmp_path)
    run_dir = tmp_path / "run"
    _run_m01(config_path, metadata_path, run_dir)
    summary = run_trim(load_config(config_path), metadata_path, run_dir)

    assert summary["n_samples"] == 4
    assert summary["read_type"] == "short"   # dispatch: kısa-okuma fastp yolu
    assert (run_dir / "trimmed" / "s1").is_dir()
    stats = json.loads((run_dir / "statistics" / "trimming_statistics.json").read_text())
    assert set(stats["samples"]) == {"s1", "s2", "s3", "s4"}
    gates = json.loads((run_dir / "quality" / "gates.json").read_text())["gates"]
    assert any(g["module"] == "m03_trim" and g["status"] == "PASS" for g in gates)
    assert any(g["module"] == "m01" for g in gates)   # m01 kapıları korundu


def _seed_long(tmp_path, chemistry="cdna"):
    from rnaforge.state import RunState
    run_dir = tmp_path / "run"
    (run_dir / "statistics").mkdir(parents=True)
    (run_dir / "statistics" / "raw_statistics.json").write_text(
        json.dumps({"read_type": "long", "chemistry": chemistry})
    )
    RunState(run_dir).mark_done("m01_validate", [])
    return run_dir


def _long_cfg(chemistry="cdna"):
    from rnaforge.config import (
        Config, Reference, Library, Trimming, DE, Report, Resources,
    )
    return Config(
        organism="E. coli", organism_type="prokaryote", platform="auto",
        reference=Reference(), library=Library(chemistry=chemistry),
        trimming=Trimming(), de=DE(), report=Report(), resources=Resources(),
    )


def test_run_trim_long_cdna_pychopper_then_chopper(tmp_path, monkeypatch):
    import rnaforge.modules.m03_trim as m03
    from rnaforge.pychopper import PychopperStats
    run_dir = _seed_long(tmp_path, "cdna")
    fq = tmp_path / "s1.fastq"
    fq.write_text("@r\n" + "ACGT" * 50 + "\n+\n" + "I" * 200 + "\n")
    meta = tmp_path / "m.tsv"
    meta.write_text(f"sample_id\tcondition\tfastq_1\ns1\tctrl\t{fq}\n")

    calls = []

    def fake_pychopper(in_fastq, out_fastq, stats_tsv, **k):
        calls.append("pychopper")
        Path(out_fastq).parent.mkdir(parents=True, exist_ok=True)
        Path(out_fastq).write_text("@r\nACGT\n+\nIIII\n")
        return PychopperStats(pass_reads=100, primers_found=60,
                              rescue=5, unusable=35, len_fail=3)

    def fake_chopper(in_fastq, out_fastq, **k):
        calls.append("chopper")
        Path(out_fastq).parent.mkdir(parents=True, exist_ok=True)
        Path(out_fastq).write_text("@r\nACGT\n+\nIIII\n")
        return 57

    monkeypatch.setattr(m03, "run_pychopper", fake_pychopper)
    monkeypatch.setattr(m03, "run_chopper", fake_chopper)

    summary = m03.run_trim(_long_cfg("cdna"), meta, run_dir)
    assert summary["read_type"] == "long"
    assert summary["chemistry"] == "cdna"
    assert calls == ["pychopper", "chopper"]      # sıra önemli
    assert summary["samples"]["s1"]["reads_after"] == 57
    assert not (run_dir / "quality" / "gates.json").exists()  # diagnostik, FAIL kapısı yok
    from rnaforge.modules.m03_trim import trimmed_reads
    from rnaforge.metadata import load_metadata
    out1, out2 = trimmed_reads(run_dir, load_metadata(meta)[0])
    assert out1.exists() and out2 is None


def test_run_trim_long_direct_rna_chopper_only(tmp_path, monkeypatch):
    import rnaforge.modules.m03_trim as m03
    run_dir = _seed_long(tmp_path, "direct_rna")
    fq = tmp_path / "s1.fastq"
    fq.write_text("@r\n" + "ACGT" * 50 + "\n+\n" + "I" * 200 + "\n")
    meta = tmp_path / "m.tsv"
    meta.write_text(f"sample_id\tcondition\tfastq_1\ns1\tctrl\t{fq}\n")

    calls = []
    monkeypatch.setattr(m03, "run_pychopper",
                        lambda *a, **k: calls.append("pychopper"))

    def fake_chopper(in_fastq, out_fastq, **k):
        calls.append("chopper")
        Path(out_fastq).parent.mkdir(parents=True, exist_ok=True)
        Path(out_fastq).write_text("@r\nACGT\n+\nIIII\n")
        return 90

    monkeypatch.setattr(m03, "run_chopper", fake_chopper)
    summary = m03.run_trim(_long_cfg("direct_rna"), meta, run_dir)
    assert calls == ["chopper"]                   # direct-RNA'da pychopper yok
    assert summary["samples"]["s1"]["reads_after"] == 90


def test_run_trim_low_survival_fails_and_records_gate(tmp_path, monkeypatch):
    _fake_fastp(monkeypatch, survival=0.10)   # eşik 0.50'nin altında
    config_path, metadata_path = _setup(tmp_path)
    run_dir = tmp_path / "run"
    _run_m01(config_path, metadata_path, run_dir)
    with pytest.raises(GateFailure):
        run_trim(load_config(config_path), metadata_path, run_dir)
    # FAIL'de bile stats + gates.json diskte kalmalı (teşhis)
    assert (run_dir / "statistics" / "trimming_statistics.json").exists()
    gates = json.loads((run_dir / "quality" / "gates.json").read_text())["gates"]
    failed = [g for g in gates if g["module"] == "m03_trim" and g["status"] == "FAIL"]
    assert failed and failed[0]["name"] == "survival_rate"


def test_run_trim_resumes_without_rerunning(tmp_path, monkeypatch):
    _fake_fastp(monkeypatch, survival=0.98)
    config_path, metadata_path = _setup(tmp_path)
    run_dir = tmp_path / "run"
    _run_m01(config_path, metadata_path, run_dir)
    run_trim(load_config(config_path), metadata_path, run_dir)
    calls = []
    monkeypatch.setattr(m03_trim, "run_fastp", lambda *a, **k: calls.append(1))
    summary = run_trim(load_config(config_path), metadata_path, run_dir)
    assert summary.get("resumed") is True
    assert calls == []


def test_cli_trim_returns_zero_and_prints_verdict(tmp_path, monkeypatch, capsys):
    from rnaforge.cli import main
    _fake_fastp(monkeypatch, survival=0.98)
    config_path, metadata_path = _setup(tmp_path)
    common = ["--config", str(config_path), "--metadata", str(metadata_path),
              "--runs-dir", str(tmp_path / "runs"), "--run-id", "demo"]
    assert main(["validate", *common]) == 0
    capsys.readouterr()
    assert main(["trim", *common]) == 0
    assert "quality verdict" in capsys.readouterr().out


def test_trimmed_reads_single_end(tmp_path):
    from rnaforge.metadata import Sample
    from rnaforge.modules.m03_trim import trimmed_reads
    sample = Sample("s1", "control", tmp_path / "c1.fastq")
    out1, out2 = trimmed_reads(tmp_path / "run", sample)
    assert out1 == tmp_path / "run" / "trimmed" / "s1" / "c1.trimmed.fastq"
    assert out2 is None


def test_trimmed_reads_paired_end(tmp_path):
    from rnaforge.metadata import Sample
    from rnaforge.modules.m03_trim import trimmed_reads
    sample = Sample("s1", "control", tmp_path / "c1_R1.fastq", tmp_path / "c1_R2.fastq")
    out1, out2 = trimmed_reads(tmp_path / "run", sample)
    assert out1.name == "c1_R1.trimmed.fastq"
    assert out2.name == "c1_R2.trimmed.fastq"


def test_cli_trim_returns_one_on_low_survival(tmp_path, monkeypatch, capsys):
    from rnaforge.cli import main
    _fake_fastp(monkeypatch, survival=0.10)
    config_path, metadata_path = _setup(tmp_path)
    common = ["--config", str(config_path), "--metadata", str(metadata_path),
              "--runs-dir", str(tmp_path / "runs"), "--run-id", "demo"]
    assert main(["validate", *common]) == 0
    capsys.readouterr()
    assert main(["trim", *common]) == 1
    assert "quality gate" in capsys.readouterr().err.lower()
