"""Run durumu: resume + heartbeat (PLAN §15 — kapatma dayanıklılığı).

Durum atomik yazılır (geçici dosya + os.replace): yarıda kapanma
bozuk state.json bırakmaz.
"""
from __future__ import annotations

import json
import os
from datetime import datetime
from pathlib import Path

STATE_FILE = "state.json"
HEARTBEAT_FILE = "heartbeat.txt"
HEARTBEAT_INTERVAL_SECONDS = 10


def new_run_dir(base: Path | str, run_id: str, now: datetime | None = None) -> Path:
    stamp = (now or datetime.now()).strftime("%Y%m%d_%H%M%S")
    run_dir = Path(base) / f"{stamp}_{run_id}"
    run_dir.mkdir(parents=True, exist_ok=True)
    return run_dir


class RunState:
    def __init__(self, run_dir: Path | str):
        self.run_dir = Path(run_dir)
        self.run_dir.mkdir(parents=True, exist_ok=True)

    @property
    def _path(self) -> Path:
        return self.run_dir / STATE_FILE

    def _read(self) -> dict:
        if not self._path.exists():
            return {"modules": {}}
        try:
            data = json.loads(self._path.read_text())
        except (json.JSONDecodeError, OSError):
            # Bozuk durum dosyası = ilerleme yok say. Çökmek yerine baştan koş.
            return {"modules": {}}
        data.setdefault("modules", {})
        return data

    def _write(self, data: dict) -> None:
        tmp = self._path.with_suffix(".json.tmp")
        tmp.write_text(json.dumps(data, indent=2))
        os.replace(tmp, self._path)  # atomik

    def mark_done(self, module: str, outputs: list[str]) -> None:
        data = self._read()
        data["modules"][module] = {
            "completed_at": datetime.now().isoformat(timespec="seconds"),
            "outputs": [str(o) for o in outputs],
        }
        self._write(data)

    def is_done(self, module: str) -> bool:
        return module in self._read()["modules"]

    def completed_modules(self) -> list[str]:
        return list(self._read()["modules"].keys())

    def heartbeat(self) -> None:
        (self.run_dir / HEARTBEAT_FILE).write_text(
            datetime.now().isoformat(timespec="seconds") + "\n"
        )
