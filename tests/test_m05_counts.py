from __future__ import annotations

from rnaforge.gates import FAIL, PASS
from rnaforge.modules.m05_counts import build_count_gates
from rnaforge.quality import load_profile


def test_all_above_threshold_passes():
    profile = load_profile("prokaryote")  # assignment_rate = 0.50
    gates = build_count_gates({"s1": 0.95, "s2": 0.80}, profile)
    assert len(gates) == 1
    assert gates[0].name == "assignment_rate"
    assert gates[0].module == "m05_counts"
    assert gates[0].status == PASS


def test_below_threshold_fails():
    profile = load_profile("prokaryote")
    gates = build_count_gates({"s1": 0.95, "s2": 0.10}, profile)
    g = gates[0]
    assert g.status == FAIL
    assert g.samples == ("s2",)
    assert g.measured == 0.10
    assert g.threshold == 0.50


def test_override_marks_overridden():
    profile = load_profile("prokaryote", {"assignment_rate": 0.05})
    gates = build_count_gates({"s1": 0.10}, profile)
    assert gates[0].status == PASS
    assert gates[0].overridden is True
