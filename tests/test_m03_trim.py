from __future__ import annotations

import json
import textwrap
from pathlib import Path

import pytest

from rnaforge.config import load_config
from rnaforge.fastp import FastpResult
from rnaforge.gates import FAIL, PASS, GateFailure
from rnaforge.modules import m03_trim
from rnaforge.modules.m03_trim import build_trim_gates, run_trim
from rnaforge.quality import load_profile
from tests.conftest import write_fastq


def _res(survival: float) -> FastpResult:
    return FastpResult(reads_before=1000, reads_after=int(1000 * survival),
                       survival_rate=survival)


def test_all_above_threshold_passes():
    profile = load_profile("prokaryote")   # survival_rate = 0.50
    gates = build_trim_gates({"s1": _res(0.98), "s2": _res(0.90)}, profile)
    assert len(gates) == 1
    g = gates[0]
    assert g.name == "survival_rate"
    assert g.module == "m03_trim"
    assert g.status == PASS


def test_below_threshold_fails_and_lists_offenders():
    profile = load_profile("prokaryote")
    gates = build_trim_gates({"s1": _res(0.98), "s2": _res(0.20)}, profile)
    g = gates[0]
    assert g.status == FAIL
    assert g.samples == ("s2",)
    assert g.measured == 0.20          # en düşük survival
    assert g.threshold == 0.50


def test_threshold_override_marks_gate_overridden():
    profile = load_profile("prokaryote", {"survival_rate": 0.10})
    gates = build_trim_gates({"s1": _res(0.20)}, profile)
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
    assert (run_dir / "trimmed" / "s1").is_dir()
    stats = json.loads((run_dir / "statistics" / "trimming_statistics.json").read_text())
    assert set(stats["samples"]) == {"s1", "s2", "s3", "s4"}
    gates = json.loads((run_dir / "quality" / "gates.json").read_text())["gates"]
    assert any(g["module"] == "m03_trim" and g["status"] == "PASS" for g in gates)
    assert any(g["module"] == "m01" for g in gates)   # m01 kapıları korundu


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
