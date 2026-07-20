"""Kalite kapıları: bir koşunun sonucuna güvenilip güvenilemeyeceğini ölçer.

Sözleşme (spec 2026-07-20):
  FAIL = sonuç GEÇERSİZ -> pipeline durur, biyolojik çıktı üretilmez
  WARN = sonuç ŞÜPHELİ  -> üretilir ama damgalanır
  PASS = kapı geçildi
"""
from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from datetime import datetime
from pathlib import Path

PASS = "PASS"
WARN = "WARN"
FAIL = "FAIL"
STATUSES = (PASS, WARN, FAIL)

GATES_FILE = "gates.json"
QUALITY_DIR = "quality"


@dataclass(frozen=True)
class GateResult:
    name: str
    module: str
    status: str
    message: str
    remedy: str
    measured: float | None = None
    threshold: float | None = None
    overridden: bool = False
    samples: list[str] = field(default_factory=list)

    def __post_init__(self) -> None:
        if self.status not in STATUSES:
            raise ValueError(f"status must be one of {STATUSES}, got {self.status!r}")
        if not self.message.strip():
            raise ValueError(f"gate {self.name}: message must not be empty")
        if not self.remedy.strip():
            # Ne yapılacağını söylemeyen kapı, müşteriyi çıkmaza sokar.
            raise ValueError(f"gate {self.name}: remedy must not be empty")


class GateFailure(Exception):
    """Bir veya daha fazla kapı FAIL verdi; pipeline burada durur."""

    def __init__(self, failures: list[GateResult]):
        self.failures = failures
        names = ", ".join(f.name for f in failures)
        detail = "\n".join(f"  - {f.name}: {f.message} -> {f.remedy}" for f in failures)
        super().__init__(f"quality gate(s) failed: {names}\n{detail}")


def raise_if_failed(results: list[GateResult]) -> None:
    failures = [r for r in results if r.status == FAIL]
    if failures:
        raise GateFailure(failures)


def write_gate_results(run_dir: Path | str, results: list[GateResult]) -> Path:
    """Sonuçları quality/gates.json'a EKLE. Aynı modülün eski sonuçları değiştirilir
    (--force ile yeniden koşma), diğer modüllerinkine dokunulmaz (resume uyumu)."""
    quality_dir = Path(run_dir) / QUALITY_DIR
    quality_dir.mkdir(parents=True, exist_ok=True)
    path = quality_dir / GATES_FILE

    existing: list[dict] = []
    if path.exists():
        try:
            existing = json.loads(path.read_text()).get("gates", [])
        except (json.JSONDecodeError, OSError):
            existing = []

    modules = {r.module for r in results}
    kept = [g for g in existing if g.get("module") not in modules]
    payload = {
        "updated_at": datetime.now().isoformat(timespec="seconds"),
        "gates": kept + [asdict(r) for r in results],
    }
    path.write_text(json.dumps(payload, indent=2))
    return path
