"""Kurulum ön-uçuş kontrolü: gerekli conda env'leri var mı, actionable rapor."""
from __future__ import annotations

from rnaforge import cli
from rnaforge.cli import main
from rnaforge.preflight import ENVIRONMENTS, environment_report


def test_environment_report_marks_present_and_missing():
    available = {"rnaforge-core", "rnaforge-de"}
    report = environment_report(available)
    assert {name for name, _desc, present in report if present} == available
    # her bilinen env raporda bir kez
    assert {name for name, _desc, _present in report} == set(ENVIRONMENTS)


def test_environment_report_all_present():
    report = environment_report(set(ENVIRONMENTS))
    assert all(present for _n, _d, present in report)


def test_doctor_exits_zero_when_all_present(monkeypatch, capsys):
    monkeypatch.setattr(cli, "list_conda_envs", lambda: set(ENVIRONMENTS))
    assert main(["doctor"]) == 0
    assert "OK" in capsys.readouterr().out


def test_doctor_exits_one_and_names_missing(monkeypatch, capsys):
    monkeypatch.setattr(cli, "list_conda_envs", lambda: {"rnaforge-core"})
    assert main(["doctor"]) == 1
    out = capsys.readouterr().out
    assert "rnaforge-de" in out and "missing" in out.lower()
