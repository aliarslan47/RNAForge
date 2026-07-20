# tests/test_confidence.py
from __future__ import annotations

import json

from rnaforge.gates import FAIL, PASS, WARN, GateResult, write_gate_results
from rnaforge.quality import load_profile
from rnaforge.report.confidence import build_confidence_card, write_confidence_card


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
