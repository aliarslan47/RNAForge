from __future__ import annotations

from rnaforge.fastqc import FastQCReport
from rnaforge.gates import PASS, WARN, FAIL
from rnaforge.modules.m02_qc import build_qc_gates


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
