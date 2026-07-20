from __future__ import annotations

import json
import textwrap

import pytest

from rnaforge.cli import main
from rnaforge.config import load_config
from rnaforge.modules.m01_validate import run_validation
from rnaforge.platform import UnsupportedPlatformError
from tests.conftest import write_fastq


def _setup(tmp_path, fastq_maker) -> tuple:
    """Geçerli config + metadata + FASTQ üretir; (config_path, metadata_path) döner."""
    (tmp_path / "ref").mkdir()
    (tmp_path / "ref" / "genome.fa").write_text(">c1\nACGT\n")
    (tmp_path / "ref" / "genes.gff").write_text("##gff-version 3\n")

    names = ["c1.fastq", "c2.fastq", "t1.fastq", "t2.fastq"]
    for n in names:
        fastq_maker(tmp_path / n)

    config_path = tmp_path / "config.yaml"
    config_path.write_text(textwrap.dedent(f"""
        organism: "Escherichia coli"
        organism_type: "prokaryote"
        reference:
          genome_fasta: "{tmp_path / 'ref' / 'genome.fa'}"
          annotation_gff: "{tmp_path / 'ref' / 'genes.gff'}"
        de:
          design: "~condition"
    """))

    metadata_path = tmp_path / "samples.tsv"
    metadata_path.write_text(
        "sample_id\tcondition\tfastq_1\n"
        "s1\tcontrol\tc1.fastq\n"
        "s2\tcontrol\tc2.fastq\n"
        "s3\ttreated\tt1.fastq\n"
        "s4\ttreated\tt2.fastq\n"
    )
    return config_path, metadata_path


def _illumina(path):
    write_fastq(path, 50, 150, "I")


def _ont(path):
    write_fastq(path, 50, (1000, 20000), "+")


def test_validation_succeeds_on_illumina(tmp_path):
    config_path, metadata_path = _setup(tmp_path, _illumina)
    run_dir = tmp_path / "run"
    summary = run_validation(load_config(config_path), metadata_path, run_dir)

    assert summary["n_samples"] == 4
    assert summary["platform"] == "illumina"
    assert summary["organism_type"] == "prokaryote"
    assert summary["conditions"] == {"control": 2, "treated": 2}


def test_validation_writes_log_and_statistics(tmp_path):
    config_path, metadata_path = _setup(tmp_path, _illumina)
    run_dir = tmp_path / "run"
    run_validation(load_config(config_path), metadata_path, run_dir)

    assert (run_dir / "logs" / "validation.log").exists()
    stats = json.loads((run_dir / "statistics" / "raw_statistics.json").read_text())
    assert stats["n_samples"] == 4
    assert len(stats["samples"]) == 4
    assert stats["samples"][0]["mean_read_length"] == pytest.approx(150.0)


def test_validation_rejects_ont(tmp_path):
    config_path, metadata_path = _setup(tmp_path, _ont)
    with pytest.raises(UnsupportedPlatformError, match="ont"):
        run_validation(load_config(config_path), metadata_path, tmp_path / "run")


def test_validation_marks_module_done_for_resume(tmp_path):
    from rnaforge.state import RunState

    config_path, metadata_path = _setup(tmp_path, _illumina)
    run_dir = tmp_path / "run"
    run_validation(load_config(config_path), metadata_path, run_dir)
    assert RunState(run_dir).is_done("m01_validate") is True


def test_cli_validate_returns_zero_on_success(tmp_path, capsys):
    config_path, metadata_path = _setup(tmp_path, _illumina)
    code = main([
        "validate",
        "--config", str(config_path),
        "--metadata", str(metadata_path),
        "--runs-dir", str(tmp_path / "runs"),
        "--run-id", "demo",
    ])
    assert code == 0
    assert "illumina" in capsys.readouterr().out


def test_cli_validate_returns_one_on_ont(tmp_path, capsys):
    config_path, metadata_path = _setup(tmp_path, _ont)
    code = main([
        "validate",
        "--config", str(config_path),
        "--metadata", str(metadata_path),
        "--runs-dir", str(tmp_path / "runs"),
        "--run-id", "demo",
    ])
    assert code == 1
    assert "not supported" in capsys.readouterr().err


from rnaforge.gates import FAIL, PASS, GateFailure


def _illumina(path):
    return write_fastq(path, 200, 150, "I")


def test_m01_writes_gate_results(tmp_path):
    """Kapilar gorunur olmali: PASS alan kosuda da neyin kontrol edildigi yazilir."""
    config_path, metadata_path = _setup(tmp_path, _illumina)
    run_dir = tmp_path / "run"
    run_validation(load_config(config_path), metadata_path, run_dir)
    data = json.loads((run_dir / "quality" / "gates.json").read_text())
    names = {g["name"] for g in data["gates"]}
    assert {"design_rank", "replication", "paired_declared"} <= names
    assert all(g["status"] == PASS for g in data["gates"])
    assert all(g["module"] == "m01" for g in data["gates"])


def test_m01_summary_records_quality_profile(tmp_path):
    config_path, metadata_path = _setup(tmp_path, _illumina)
    run_dir = tmp_path / "run"
    summary = run_validation(load_config(config_path), metadata_path, run_dir)
    assert summary["quality_profile"] == "prokaryote"
    assert summary["permissive_profile"] is False


def test_m01_failing_gate_stops_the_run_and_is_recorded(tmp_path):
    """FAIL: kosu DURUR, ama dusen kapi gates.json'a YAZILMIS olmali —
    teshis raporunun gosterecek verisi buradan gelir (spec 3.5)."""
    config_path, metadata_path = _setup(tmp_path, _illumina)
    # replikasiz tasarim: her condition tek ornek
    metadata_path.write_text(
        "sample_id\tcondition\tfastq_1\n"
        "s1\tcontrol\tc1.fastq\n"
        "s2\ttreated\tt1.fastq\n"
    )
    run_dir = tmp_path / "run"
    with pytest.raises(GateFailure):
        run_validation(load_config(config_path), metadata_path, run_dir)

    data = json.loads((run_dir / "quality" / "gates.json").read_text())
    failed = [g for g in data["gates"] if g["status"] == FAIL]
    assert [g["name"] for g in failed] == ["replication"]
    assert failed[0]["remedy"]


def test_m01_does_not_write_statistics_when_a_gate_fails(tmp_path):
    """FAIL = sonuc GECERSIZ: biyolojik cikti URETILMEZ (spec karar 4)."""
    config_path, metadata_path = _setup(tmp_path, _illumina)
    metadata_path.write_text(
        "sample_id\tcondition\tfastq_1\n"
        "s1\tcontrol\tc1.fastq\n"
        "s2\ttreated\tt1.fastq\n"
    )
    run_dir = tmp_path / "run"
    with pytest.raises(GateFailure):
        run_validation(load_config(config_path), metadata_path, run_dir)
    assert not (run_dir / "statistics" / "raw_statistics.json").exists()
