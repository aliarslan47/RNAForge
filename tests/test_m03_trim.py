from __future__ import annotations

from rnaforge.fastp import FastpResult
from rnaforge.gates import FAIL, PASS
from rnaforge.modules.m03_trim import build_trim_gates
from rnaforge.quality import load_profile


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
