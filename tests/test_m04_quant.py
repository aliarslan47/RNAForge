from __future__ import annotations

from pathlib import Path

from rnaforge.bowtie2 import AlignmentResult
from rnaforge.gates import FAIL, PASS
from rnaforge.modules.m04_quant import build_alignment_gates
from rnaforge.quality import load_profile


def _res(rate: float) -> AlignmentResult:
    return AlignmentResult(bam=Path("x.bam"), alignment_rate=rate)


def test_all_above_threshold_passes():
    profile = load_profile("prokaryote")  # alignment_rate = 0.70
    gates = build_alignment_gates({"s1": _res(0.98), "s2": _res(0.85)}, profile)
    assert len(gates) == 1
    assert gates[0].name == "alignment_rate"
    assert gates[0].module == "m04_quant"
    assert gates[0].status == PASS


def test_below_threshold_fails_and_lists_offenders():
    profile = load_profile("prokaryote")
    gates = build_alignment_gates({"s1": _res(0.98), "s2": _res(0.30)}, profile)
    g = gates[0]
    assert g.status == FAIL
    assert g.samples == ("s2",)
    assert g.measured == 0.30
    assert g.threshold == 0.70


def test_override_marks_gate_overridden():
    profile = load_profile("prokaryote", {"alignment_rate": 0.20})
    gates = build_alignment_gates({"s1": _res(0.30)}, profile)
    assert gates[0].status == PASS
    assert gates[0].overridden is True
    assert gates[0].threshold == 0.20
