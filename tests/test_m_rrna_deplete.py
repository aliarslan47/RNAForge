"""m_rrna_deplete testi: run_sortmerna_deplete monkeypatch'lenir (SortMeRNA gerektirmez)."""
from __future__ import annotations

import json
import textwrap
from pathlib import Path

import pytest

from rnaforge.config import load_config
from rnaforge.gates import PASS, WARN
from rnaforge.modules import m_rrna_deplete
from rnaforge.modules.m01_validate import run_validation
from rnaforge.modules.m_rrna_deplete import (
    build_rrna_gates, rrna_depleted_reads, run_rrna_deplete,
)
from rnaforge.quality import load_profile
from tests.conftest import write_fastq


def test_gate_passes_above_threshold():
    profile = load_profile("metatranscriptome")   # rrna_depletion_rate = 0.30
    gates = build_rrna_gates({"s1": 0.40, "s2": 0.90}, profile)
    assert len(gates) == 1
    g = gates[0]
    assert g.name == "rrna_depletion_rate"
    assert g.module == "m_rrna_deplete"
    assert g.status == PASS


def test_gate_below_threshold_is_warn_never_fail():
    profile = load_profile("metatranscriptome")
    gates = build_rrna_gates({"s1": 0.05, "s2": 0.90}, profile)
    g = gates[0]
    assert g.status == WARN                # asla FAIL, damgalı-permissive profil
    assert g.samples == ("s1",)
    assert g.measured == 0.05


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


def _run_m01(config_path, metadata_path, run_dir):
    run_validation(load_config(config_path), metadata_path, run_dir)


def _fake_sortmerna(monkeypatch, depletion_rate=0.40):
    def fake_run(reads, rrna_db, workdir, paired, threads=8, env="rnaforge-seqqc"):
        workdir = Path(workdir)
        workdir.mkdir(parents=True, exist_ok=True)
        other = workdir / "other.fastq"
        other.write_text("@r\nACGT\n+\nIIII\n")
        log = workdir / "aligned.log"
        log.write_text("Total reads = 100\nTotal reads passing E-value threshold = "
                       f"{int((1 - depletion_rate) * 100)}\n")
        return {"other": [other], "depletion_rate": depletion_rate, "aligned_log": log}
    monkeypatch.setattr(m_rrna_deplete, "run_sortmerna_deplete", fake_run)


def test_run_rrna_deplete_requires_m01(tmp_path, monkeypatch):
    _fake_sortmerna(monkeypatch)
    config_path, metadata_path = _setup(tmp_path)
    with pytest.raises(ValueError, match="m01"):
        run_rrna_deplete(load_config(config_path), metadata_path, tmp_path / "run")


def test_run_rrna_deplete_writes_stats_and_output(tmp_path, monkeypatch):
    _fake_sortmerna(monkeypatch, depletion_rate=0.40)
    config_path, metadata_path = _setup(tmp_path)
    run_dir = tmp_path / "run"
    _run_m01(config_path, metadata_path, run_dir)
    config = load_config(config_path)
    summary = run_rrna_deplete(config, metadata_path, run_dir)

    assert summary["n_samples"] == 4
    stats = json.loads((run_dir / "statistics" / "rrna_depletion.json").read_text())
    # Sözleşme: {sid: {"depletion_rate": float}} — TAM olarak bu şekil.
    assert set(stats) == {"s1", "s2", "s3", "s4"}
    for sid, entry in stats.items():
        assert entry == {"depletion_rate": 0.4}

    # other fastq sözleşme yolunda: rrna_depleted/<sid>/other_*.fastq.gz
    for sid in stats:
        out_dir = run_dir / "rrna_depleted" / sid
        others = sorted(out_dir.glob("other_*.fastq.gz"))
        assert len(others) == 1
        from rnaforge.metadata import load_metadata
        sample = next(s for s in load_metadata(metadata_path) if s.sample_id == sid)
        assert rrna_depleted_reads(run_dir, sample) == others

    gates = json.loads((run_dir / "quality" / "gates.json").read_text())["gates"]
    rrna_gate = next(g for g in gates if g["module"] == "m_rrna_deplete")
    # 0.40 >= eşik 0.30 -> PASS; verdict FAIL YAPMAZ (bu modülde FAIL zaten üretilemez).
    assert rrna_gate["status"] == PASS
    assert not any(g["status"] == "FAIL" for g in gates if g["module"] == "m_rrna_deplete")

    from rnaforge.state import RunState
    state = RunState(run_dir)
    assert state.is_done("m_rrna_deplete")
    for sid in stats:
        assert state.is_item_done("m_rrna_deplete", sid)


def test_run_rrna_deplete_resumes_without_rerunning(tmp_path, monkeypatch):
    _fake_sortmerna(monkeypatch, depletion_rate=0.40)
    config_path, metadata_path = _setup(tmp_path)
    run_dir = tmp_path / "run"
    _run_m01(config_path, metadata_path, run_dir)
    config = load_config(config_path)
    run_rrna_deplete(config, metadata_path, run_dir)

    calls = []
    def boom(*a, **k):
        calls.append(1)
        raise AssertionError("run_sortmerna_deplete should not be called again on resume")
    monkeypatch.setattr(m_rrna_deplete, "run_sortmerna_deplete", boom)

    summary = run_rrna_deplete(config, metadata_path, run_dir)
    assert summary["resumed"] is True
    assert summary["n_samples"] == 4
    assert not calls


def test_run_rrna_deplete_low_depletion_warns_not_fails(tmp_path, monkeypatch):
    """Düşük depletion oranı bile pipeline'ı durdurmamalı (WARN-only kapı, asla FAIL)."""
    _fake_sortmerna(monkeypatch, depletion_rate=0.01)
    config_path, metadata_path = _setup(tmp_path)
    run_dir = tmp_path / "run"
    _run_m01(config_path, metadata_path, run_dir)
    config = load_config(config_path)
    summary = run_rrna_deplete(config, metadata_path, run_dir)   # raise ETMEMELİ

    assert summary["gate_counts"].get("FAIL", 0) == 0
    gates = json.loads((run_dir / "quality" / "gates.json").read_text())["gates"]
    rrna_gate = next(g for g in gates if g["module"] == "m_rrna_deplete")
    assert rrna_gate["status"] == WARN
