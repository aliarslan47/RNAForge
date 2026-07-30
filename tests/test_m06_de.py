from __future__ import annotations

from rnaforge.gates import FAIL, PASS, WARN
from rnaforge.modules.m06_de import build_de_gates
from rnaforge.quality import load_profile


def test_high_correlation_passes():
    profile = load_profile("prokaryote")  # replicate_correlation = 0.85
    gates = build_de_gates(0.95, profile)
    assert len(gates) == 1
    assert gates[0].name == "replicate_correlation"
    assert gates[0].module == "m06_de"
    assert gates[0].status == PASS


def test_low_correlation_warns_never_fails():
    profile = load_profile("prokaryote")
    gates = build_de_gates(0.40, profile)
    assert gates[0].status == WARN          # FAIL degil!
    assert gates[0].measured == 0.40
    assert gates[0].threshold == 0.85
    assert all(g.status != FAIL for g in gates)


def test_override_marks_overridden():
    profile = load_profile("prokaryote", {"replicate_correlation": 0.20})
    gates = build_de_gates(0.40, profile)
    assert gates[0].status == PASS
    assert gates[0].overridden is True
