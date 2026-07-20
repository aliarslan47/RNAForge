"""Run durumu: resume + heartbeat (PLAN §15 — kapatma dayanıklılığı).

Durum atomik yazılır (geçici dosya + os.replace): yarıda kapanma
bozuk state.json bırakmaz.
"""
from __future__ import annotations

import json
import os
import time
from datetime import datetime
from pathlib import Path

STATE_FILE = "state.json"
HEARTBEAT_FILE = "heartbeat.txt"
HEARTBEAT_INTERVAL_SECONDS = 10


def resolve_run_dir(base: Path | str, run_id: str, now: datetime | None = None) -> Path:
    """Aynı run_id için VAR OLAN dizini yeniden kullan; yoksa yenisini aç.

    Resume'un dayanağı budur: yeni zaman damgalı dizin açmak state.json'ı
    erişilemez kılar ve "kaldığı yerden devam" iddiasını sessizce boşa çıkarır
    (PLAN §15). Aynı run_id ile birden çok eski dizin varsa en yenisi seçilir.
    """
    base = Path(base)
    # ????????_?????? = %Y%m%d_%H%M%S; gevşek `*_id` deseni "a_run"u "run" sanardı.
    existing = sorted(p for p in base.glob(f"????????_??????_{run_id}") if p.is_dir())
    if existing:
        return existing[-1]
    stamp = (now or datetime.now()).strftime("%Y%m%d_%H%M%S")
    run_dir = base / f"{stamp}_{run_id}"
    run_dir.mkdir(parents=True, exist_ok=True)
    return run_dir


class RunState:
    def __init__(self, run_dir: Path | str):
        self.run_dir = Path(run_dir)
        self.run_dir.mkdir(parents=True, exist_ok=True)
        self._last_heartbeat: float | None = None

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
        # Sabit .tmp adı yerine pid'li: resume ile aynı dizine iki süreç yazarsa çakışmasın.
        tmp = self._path.with_suffix(f".json.tmp.{os.getpid()}")
        with tmp.open("w") as handle:
            json.dump(data, handle, indent=2)
            handle.flush()
            os.fsync(handle.fileno())  # replace'ten önce diske insin, yoksa atomiklik iddiası eksik
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

    def heartbeat(self, force: bool = False) -> None:
        """En fazla HEARTBEAT_INTERVAL_SECONDS'ta bir yaz (PLAN §15: 10 sn kadans).

        Çağıranlar her birimde çağırır; kadansı uygulamak burasının işi. Aksi halde
        sabit tanımlı ama karşılıksız kalır ve binlerce örnekte gereksiz I/O olur.
        """
        now = time.monotonic()
        if not force and self._last_heartbeat is not None:
            if now - self._last_heartbeat < HEARTBEAT_INTERVAL_SECONDS:
                return
        self._last_heartbeat = now
        (self.run_dir / HEARTBEAT_FILE).write_text(
            datetime.now().isoformat(timespec="seconds") + "\n"
        )
