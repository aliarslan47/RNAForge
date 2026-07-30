# tests/test_confidence.py
from __future__ import annotations

import json

import pytest

from rnaforge.cli import main
from rnaforge.gates import FAIL, PASS, WARN, GateResult, write_gate_results
from rnaforge.quality import load_profile
from rnaforge.report.confidence import UNKNOWN, build_confidence_card, write_confidence_card
from tests.test_m01_validate import _illumina, _setup


def _gate(name, status, module="m01"):
    return GateResult(
        name=name, module=module, status=status,
        message="check", remedy="do something",
    )


def test_card_summarises_gate_counts(tmp_path):
    write_gate_results(tmp_path, [
        _gate("design_rank", PASS),
        _gate("rrna_fraction", WARN, module="m04"),
    ])
    card = build_confidence_card(tmp_path, load_profile("prokaryote"))
    assert card["counts"] == {"PASS": 1, "WARN": 1, "FAIL": 0}
    assert card["verdict"] == "SUSPECT"


def test_verdict_is_trustworthy_when_all_pass(tmp_path):
    write_gate_results(tmp_path, [_gate("design_rank", PASS)])
    card = build_confidence_card(tmp_path, load_profile("prokaryote"))
    assert card["verdict"] == "TRUSTWORTHY"


def test_verdict_is_invalid_on_any_fail(tmp_path):
    write_gate_results(tmp_path, [_gate("design_rank", PASS), _gate("alignment_rate", FAIL)])
    card = build_confidence_card(tmp_path, load_profile("prokaryote"))
    assert card["verdict"] == "INVALID"


def test_card_records_permissive_profile_and_overrides(tmp_path):
    """Gevsetilen esik ve gevsek profil GORUNMEK zorunda (spec §3.2)."""
    write_gate_results(tmp_path, [_gate("design_rank", PASS)])
    profile = load_profile("eukaryote", overrides={"alignment_rate": 0.2})
    card = build_confidence_card(tmp_path, profile)
    assert card["profile"]["permissive"] is True
    assert card["profile"]["overrides"] == {"alignment_rate": 0.2}


def test_write_confidence_card_creates_file(tmp_path):
    write_gate_results(tmp_path, [_gate("design_rank", PASS)])
    path = write_confidence_card(tmp_path, load_profile("prokaryote"))
    assert path.name == "confidence_card.json"
    assert json.loads(path.read_text())["verdict"] == "TRUSTWORTHY"


# --- Finding 1: karti FAIL yolunda da yazilmali ------------------------------

def test_cli_validate_writes_confidence_card_on_gate_failure(tmp_path, capsys):
    """CLI'da bir kapi FAIL verdiginde de confidence_card.json diskte kalmali:
    teshis raporu kosu basarisiz oldugunda tam da bu veriye ihtiyac duyar."""
    config_path, metadata_path = _setup(tmp_path, _illumina)
    # replikasiz tasarim -> "replication" kapisi FAIL verir (bkz. test_m01_validate.py)
    metadata_path.write_text(
        "sample_id\tcondition\tfastq_1\n"
        "s1\tcontrol\tc1.fastq\n"
        "s2\ttreated\tt1.fastq\n"
    )
    code = main([
        "validate",
        "--config", str(config_path),
        "--metadata", str(metadata_path),
        "--runs-dir", str(tmp_path / "runs"),
        "--run-id", "demo",
    ])
    captured = capsys.readouterr()

    assert code == 1
    assert "quality gate(s) failed" in captured.err
    assert "replication" in captured.err

    # resolve_run_dir zaman damgasi ekler (ör. 20260720_120000_demo); tam adi bilmiyoruz.
    run_dirs = list((tmp_path / "runs").glob("*_demo"))
    assert len(run_dirs) == 1
    run_dir = run_dirs[0]
    card_path = run_dir / "quality" / "confidence_card.json"
    assert card_path.exists()
    card = json.loads(card_path.read_text())
    assert card["verdict"] == "INVALID"
    failed_names = {g["name"] for g in card["gates"] if g["status"] == FAIL}
    assert "replication" in failed_names


# --- Finding 2: hic kapi yoksa TRUSTWORTHY degil, UNKNOWN --------------------

def test_verdict_is_unknown_when_gates_file_is_missing(tmp_path):
    card = build_confidence_card(tmp_path, load_profile("prokaryote"))
    assert card["verdict"] == UNKNOWN
    assert card["counts"] == {"PASS": 0, "WARN": 0, "FAIL": 0}


def test_verdict_is_unknown_when_gates_list_is_empty(tmp_path):
    write_gate_results(tmp_path, [])
    card = build_confidence_card(tmp_path, load_profile("prokaryote"))
    assert card["verdict"] == UNKNOWN


# --- Finding 4: WARN + FAIL birlikteyken FAIL kazanir (INVALID) -------------

def test_verdict_is_invalid_when_warn_and_fail_both_present(tmp_path):
    write_gate_results(tmp_path, [
        _gate("rrna_fraction", WARN, module="m04"),
        _gate("alignment_rate", FAIL, module="m04"),
    ])
    card = build_confidence_card(tmp_path, load_profile("prokaryote"))
    assert card["verdict"] == "INVALID"
    assert card["counts"] == {"PASS": 0, "WARN": 1, "FAIL": 1}


# --- Finding 3: atomik yazim + bozuk gates.json'a dayaniklilik ---------------

def test_write_confidence_card_leaves_no_tmp_file_behind(tmp_path):
    write_gate_results(tmp_path, [_gate("design_rank", PASS)])
    write_confidence_card(tmp_path, load_profile("prokaryote"))
    leftovers = list((tmp_path / "quality").glob("confidence_card.json.tmp*"))
    assert leftovers == []


def test_corrupt_gates_file_does_not_crash_card_build(tmp_path):
    """Bozuk gates.json build_confidence_card'i JSONDecodeError ile patlatmamali;
    gates.py'deki desenle tutarli: kenara al, yuksek sesle uyar, calismaya devam et."""
    gates_dir = tmp_path / "quality"
    gates_dir.mkdir()
    (gates_dir / "gates.json").write_text('{"gates": [{"name": "design_rank"')

    with pytest.warns(UserWarning, match="corrupt"):
        card = build_confidence_card(tmp_path, load_profile("prokaryote"))

    assert card["verdict"] == UNKNOWN
    corrupt_candidates = list(gates_dir.glob("gates.json.corrupt*"))
    assert len(corrupt_candidates) == 1
