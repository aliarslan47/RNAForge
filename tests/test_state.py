from __future__ import annotations

import json
from datetime import datetime

from rnaforge.state import RunState, new_run_dir


def test_new_run_dir_has_timestamp_and_id(tmp_path):
    now = datetime(2026, 7, 16, 14, 30, 22)
    run_dir = new_run_dir(tmp_path, "demo", now=now)
    assert run_dir.name == "20260716_143022_demo"
    assert run_dir.exists()


def test_module_not_done_initially(tmp_path):
    assert RunState(tmp_path).is_done("m01_validate") is False


def test_mark_done_persists_across_instances(tmp_path):
    RunState(tmp_path).mark_done("m01_validate", ["a.json"])
    # yeni instance = süreç yeniden başladı
    assert RunState(tmp_path).is_done("m01_validate") is True


def test_completed_modules_listed_in_order(tmp_path):
    state = RunState(tmp_path)
    state.mark_done("m01_validate", [])
    state.mark_done("m02_qc", [])
    assert state.completed_modules() == ["m01_validate", "m02_qc"]


def test_state_file_is_valid_json(tmp_path):
    RunState(tmp_path).mark_done("m01_validate", ["out/a.json"])
    data = json.loads((tmp_path / "state.json").read_text())
    assert data["modules"]["m01_validate"]["outputs"] == ["out/a.json"]
    assert "completed_at" in data["modules"]["m01_validate"]


def test_heartbeat_writes_valid_timestamp(tmp_path):
    state = RunState(tmp_path)
    state.heartbeat()
    path = tmp_path / "heartbeat.txt"
    assert path.exists()
    # İçerik geçerli bir ISO timestamp olmalı; bozuksa fromisoformat raise eder.
    datetime.fromisoformat(path.read_text().strip())


def test_corrupt_state_file_does_not_crash(tmp_path):
    (tmp_path / "state.json").write_text("{ bozuk json")
    # bozuk durum = hiç ilerleme yok kabul edilir, çökme YOK
    assert RunState(tmp_path).is_done("m01_validate") is False
