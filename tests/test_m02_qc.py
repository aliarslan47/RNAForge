from __future__ import annotations

import json
import textwrap
from pathlib import Path

import pytest

from rnaforge.config import load_config
from rnaforge.fastqc import FastQCReport
from rnaforge.gates import PASS, WARN, FAIL
from rnaforge.modules import m02_qc
from rnaforge.modules.m02_qc import build_qc_gates, run_qc
from tests.conftest import write_fastq


def _report(**flags) -> FastQCReport:
    base = {
        "Per base sequence quality": "PASS",
        "Adapter Content": "PASS",
        "Overrepresented sequences": "PASS",
        "Per sequence GC content": "PASS",
    }
    base.update(flags)
    return FastQCReport(modules=base, basic_stats={"Total Sequences": "500"})


def test_all_pass_gives_pass_gates():
    gates = build_qc_gates({"s1": _report()})
    assert {g.name for g in gates} == {
        "per_base_quality", "adapter_content", "overrepresented", "gc_content"
    }
    assert all(g.status == PASS for g in gates)
    assert all(g.module == "m02_qc" for g in gates)


def test_fastqc_fail_maps_to_warn_never_fail():
    gates = build_qc_gates({"s1": _report(**{"Per base sequence quality": "FAIL"})})
    pbq = next(g for g in gates if g.name == "per_base_quality")
    assert pbq.status == WARN            # FAIL degil!
    assert "s1" in pbq.samples
    assert all(g.status != FAIL for g in gates)


def test_worst_flag_across_samples_wins_and_lists_offenders():
    gates = build_qc_gates({
        "s1": _report(),
        "s2": _report(**{"Adapter Content": "WARN"}),
        "s3": _report(**{"Adapter Content": "FAIL"}),
    })
    adapter = next(g for g in gates if g.name == "adapter_content")
    assert adapter.status == WARN
    assert set(adapter.samples) == {"s2", "s3"}   # PASS olan s1 yok


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
    # Replikalı tasarım: m01 replication kapısı koşul başına ≥2 örnek ister.
    metadata_path.write_text(
        "sample_id\tcondition\tfastq_1\n"
        "s1\tcontrol\tc1.fastq\n"
        "s2\tcontrol\tc2.fastq\n"
        "s3\ttreated\tt1.fastq\n"
        "s4\ttreated\tt2.fastq\n"
    )
    return config_path, metadata_path


def _run_m01(config_path, metadata_path, run_dir):
    """m02'nin ön koşulu: m01 bu run_dir'de tamamlanmış olmalı. Gerçek
    biyoinformatik araç GEREKMEZ (m01 saf Python: platform tespiti + kapılar)."""
    from rnaforge.modules.m01_validate import run_validation
    run_validation(load_config(config_path), metadata_path, run_dir)


def _fake_fastqc(monkeypatch, flags_by_stem=None):
    """run_fastqc'yi sahte bir zip yolu döndürecek, parse_fastqc_zip'i sabit
    rapor döndürecek şekilde değiştirir; gerçek FastQC çağrılmaz."""
    def fake_run(fastq, out_dir, env="rnaforge-qc"):
        p = Path(out_dir) / (Path(fastq).stem + "_fastqc.zip")
        p.write_bytes(b"stub")
        return p
    def fake_parse(zip_path):
        return FastQCReport(
            modules={"Per base sequence quality": "PASS", "Adapter Content": "PASS",
                     "Overrepresented sequences": "PASS", "Per sequence GC content": "PASS"},
            basic_stats={"Total Sequences": "200", "%GC": "50"},
        )
    monkeypatch.setattr(m02_qc, "run_fastqc", fake_run)
    monkeypatch.setattr(m02_qc, "parse_fastqc_zip", fake_parse)


def test_run_qc_requires_m01_done(tmp_path, monkeypatch):
    """Kural 7 + tespit≠destek: m01 koşmadan ham FASTQ'ya FastQC koşma —
    platform reddi (ONT vb.) atlanmış olabilir. Net hata ver, sessiz koşma."""
    _fake_fastqc(monkeypatch)
    config_path, metadata_path = _setup(tmp_path)
    with pytest.raises(ValueError, match="m01"):
        run_qc(load_config(config_path), metadata_path, tmp_path / "run")


def test_run_qc_writes_outputs_and_gates(tmp_path, monkeypatch):
    _fake_fastqc(monkeypatch)
    config_path, metadata_path = _setup(tmp_path)
    run_dir = tmp_path / "run"
    _run_m01(config_path, metadata_path, run_dir)
    summary = run_qc(load_config(config_path), metadata_path, run_dir)

    assert summary["n_samples"] == 4
    assert (run_dir / "raw_qc" / "s1").is_dir()
    assert (run_dir / "raw_qc" / "s4").is_dir()
    stats = json.loads((run_dir / "statistics" / "qc_statistics.json").read_text())
    assert set(stats["samples"]) == {"s1", "s2", "s3", "s4"}
    gates = json.loads((run_dir / "quality" / "gates.json").read_text())["gates"]
    assert any(g["module"] == "m02_qc" for g in gates)
    assert any(g["module"] == "m01" for g in gates)  # m01 kapıları KORUNDU
    assert (run_dir / "logs" / "qc.log").exists()


def test_run_qc_resumes_without_rerunning(tmp_path, monkeypatch):
    _fake_fastqc(monkeypatch)
    config_path, metadata_path = _setup(tmp_path)
    run_dir = tmp_path / "run"
    _run_m01(config_path, metadata_path, run_dir)
    run_qc(load_config(config_path), metadata_path, run_dir)

    calls = []
    monkeypatch.setattr(m02_qc, "run_fastqc",
                        lambda *a, **k: calls.append(1))  # çağrılırsa patlar (yol yok)
    summary = run_qc(load_config(config_path), metadata_path, run_dir)
    assert summary.get("resumed") is True
    assert calls == []  # FastQC tekrar çağrılmadı


def test_run_qc_gate_counts_reported(tmp_path, monkeypatch):
    _fake_fastqc(monkeypatch)
    config_path, metadata_path = _setup(tmp_path)
    run_dir = tmp_path / "run"
    _run_m01(config_path, metadata_path, run_dir)
    summary = run_qc(load_config(config_path), metadata_path, run_dir)
    assert summary["gate_counts"]["PASS"] >= 1
    assert summary["gate_counts"].get("FAIL", 0) == 0


def test_cli_qc_returns_zero_and_prints_verdict(tmp_path, monkeypatch, capsys):
    from rnaforge.cli import main
    _fake_fastqc(monkeypatch)
    config_path, metadata_path = _setup(tmp_path)
    common = [
        "--config", str(config_path),
        "--metadata", str(metadata_path),
        "--runs-dir", str(tmp_path / "runs"),
        "--run-id", "demo",
    ]
    # Önce m01: resolve_run_dir aynı --run-id'yi aynı run dizinine çözer,
    # böylece qc, validate'in bıraktığı state'i (m01 done) görür.
    assert main(["validate", *common]) == 0
    capsys.readouterr()  # validate çıktısını temizle
    assert main(["qc", *common]) == 0
    out = capsys.readouterr().out
    assert "quality verdict" in out
