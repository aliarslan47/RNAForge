"""`rnaforge run` orkestratörü: tam pipeline'ı tek komutta, sırayla, stop-on-FAIL."""
from __future__ import annotations

import pytest

from rnaforge import cli
from rnaforge.cli import build_run_sequence, main
from rnaforge.gates import FAIL, GateFailure, GateResult


def test_build_run_sequence_default_is_core_m01_to_report():
    assert build_run_sequence() == [
        "validate", "qc", "trim", "quant", "counts", "de", "figures", "report"
    ]


def test_build_run_sequence_from_to_slices_core():
    assert build_run_sequence(start="trim", end="counts") == ["trim", "quant", "counts"]


def test_build_run_sequence_include_runs_optionals_before_report():
    seq = build_run_sequence(include=["kegg", "enrich"])
    assert seq[-1] == "report"                    # report her zaman en son
    # opsiyoneller report'tan ÖNCE ve kanonik sırada (enrich < kegg)
    assert seq[seq.index("figures") + 1:seq.index("report")] == ["enrich", "kegg"]


def test_build_run_sequence_unknown_include_raises():
    with pytest.raises(ValueError, match="unknown"):
        build_run_sequence(include=["nonsense"])


def test_build_run_sequence_from_after_to_raises():
    with pytest.raises(ValueError, match="after"):
        build_run_sequence(start="report", end="validate")


def _fake_dispatch(monkeypatch, fail_at=None):
    calls = []

    def make(name):
        def fn(args):
            calls.append(name)
            if name == fail_at:
                raise GateFailure([GateResult(
                    name="demo_gate", module=name, status=FAIL,
                    message=f"{name} failed", remedy="fix it")])
            return 0
        return fn

    monkeypatch.setattr(cli, "_STAGE_DISPATCH",
                        {name: make(name) for name in cli._STAGE_DISPATCH})
    return calls


def _run_args(tmp_path, *extra):
    return ["run", "--config", str(tmp_path / "c.yaml"),
            "--metadata", str(tmp_path / "m.tsv"),
            "--runs-dir", str(tmp_path / "runs"), "--run-id", "demo", *extra]


def test_run_executes_all_core_stages_in_order(tmp_path, monkeypatch):
    calls = _fake_dispatch(monkeypatch)
    assert main(_run_args(tmp_path)) == 0
    assert calls == ["validate", "qc", "trim", "quant", "counts", "de", "figures", "report"]


def test_run_stops_on_first_gate_failure(tmp_path, monkeypatch):
    calls = _fake_dispatch(monkeypatch, fail_at="counts")
    assert main(_run_args(tmp_path)) == 1                 # GateFailure → exit 1
    assert calls == ["validate", "qc", "trim", "quant", "counts"]   # de/figures/report KOŞMAZ
