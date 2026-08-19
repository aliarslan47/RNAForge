from __future__ import annotations

import json
from datetime import datetime

from rnaforge.state import RunState, resolve_run_dir


def test_resolve_run_dir_has_timestamp_and_id(tmp_path):
    now = datetime(2026, 7, 16, 14, 30, 22)
    run_dir = resolve_run_dir(tmp_path, "demo", now=now)
    assert run_dir.name == "20260716_143022_demo"
    assert run_dir.exists()


def test_resolve_run_dir_reuses_existing_run_id(tmp_path):
    """Resume'un dayanağı: aynı run_id ikinci kez yeni dizin AÇMAMALI.
    Açarsa state.json erişilemez kalır ve 'kaldığı yerden devam' sessizce çöker."""
    first = resolve_run_dir(tmp_path, "demo", now=datetime(2026, 7, 16, 14, 30, 22))
    second = resolve_run_dir(tmp_path, "demo", now=datetime(2026, 7, 16, 15, 0, 0))
    assert second == first
    assert len(list(tmp_path.iterdir())) == 1


def test_resolve_run_dir_does_not_confuse_similar_ids(tmp_path):
    resolve_run_dir(tmp_path, "a_run", now=datetime(2026, 7, 16, 14, 30, 22))
    other = resolve_run_dir(tmp_path, "run", now=datetime(2026, 7, 16, 15, 0, 0))
    assert other.name == "20260716_150000_run"


def test_heartbeat_is_throttled_to_interval(tmp_path):
    """HEARTBEAT_INTERVAL_SECONDS gerçekten uygulanmalı; aksi halde sabit
    tanımlı ama karşılıksız kalır (PLAN §15 iddiası boşa çıkar)."""
    state = RunState(tmp_path)
    state.heartbeat()
    first = (tmp_path / "heartbeat.txt").read_text()
    (tmp_path / "heartbeat.txt").write_text("MARKER\n")
    state.heartbeat()  # interval dolmadı -> yazmamalı
    assert (tmp_path / "heartbeat.txt").read_text() == "MARKER\n"
    state.heartbeat(force=True)  # force -> yazmalı
    assert (tmp_path / "heartbeat.txt").read_text() != "MARKER\n"
    assert first  # ilk çağrı her zaman yazar


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


# --- Faz 3: örnek-başı checkpoint ---------------------------------------------

def test_item_not_done_initially(tmp_path):
    assert RunState(tmp_path).is_item_done("m03_trim", "s1") is False
    assert RunState(tmp_path).item_payload("m03_trim", "s1") is None


def test_mark_item_done_persists_across_instances(tmp_path):
    """Örnek-başı işaretçi + payload süreç yeniden başlasa da kalıcı olmalı (resume)."""
    RunState(tmp_path).mark_item_done("m03_trim", "s1", {"survival_rate": 0.9})
    reloaded = RunState(tmp_path)
    assert reloaded.is_item_done("m03_trim", "s1") is True
    assert reloaded.item_payload("m03_trim", "s1") == {"survival_rate": 0.9}


def test_mark_item_done_does_not_mark_module_done(tmp_path):
    """KRİTİK: örnek-başı işaretçi modülü 'done' YAPMAMALI — yoksa downstream bağımlılık
    guard'ı (m04 is_done('m03_trim')) yarıda kalmış aşamayı tamamlanmış sanar."""
    state = RunState(tmp_path)
    state.mark_item_done("m03_trim", "s1", {})
    assert state.is_done("m03_trim") is False
    assert RunState(tmp_path).is_done("m03_trim") is False


def test_mark_done_and_item_markers_coexist(tmp_path):
    """Aşama-düzeyi mark_done ile örnek-düzeyi işaretçiler ayrı yaşar (birbirini silmez)."""
    state = RunState(tmp_path)
    state.mark_item_done("m03_trim", "s1", {"survival_rate": 0.8})
    state.mark_done("m03_trim", ["out.json"])
    assert state.is_done("m03_trim") is True
    assert state.is_item_done("m03_trim", "s1") is True


def test_item_state_file_is_valid_json(tmp_path):
    RunState(tmp_path).mark_item_done("m04_quant", "s1", {"alignment_rate": 0.95})
    data = json.loads((tmp_path / "state.json").read_text())
    assert data["items"]["m04_quant"]["s1"] == {"alignment_rate": 0.95}
