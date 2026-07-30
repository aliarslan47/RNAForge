"""Kalite kapıları: bir koşunun sonucuna güvenilip güvenilemeyeceğini ölçer.

Sözleşme (spec 2026-07-20):
  FAIL = sonuç GEÇERSİZ -> pipeline durur, biyolojik çıktı üretilmez
  WARN = sonuç ŞÜPHELİ  -> üretilir ama damgalanır
  PASS = kapı geçildi
"""
from __future__ import annotations

import json
import os
import warnings
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
    samples: tuple[str, ...] = field(default_factory=tuple)

    def __post_init__(self) -> None:
        if self.status not in STATUSES:
            raise ValueError(f"status must be one of {STATUSES}, got {self.status!r}")
        if not self.message.strip():
            raise ValueError(f"gate {self.name}: message must not be empty")
        if not self.remedy.strip():
            # Ne yapılacağını söylemeyen kapı, müşteriyi çıkmaza sokar.
            raise ValueError(f"gate {self.name}: remedy must not be empty")
        # frozen dataclass yerinde mutasyonu engeller ama alan liste ise
        # result.samples.append(...) hala calisir; tuple'a cevirerek kapatiyoruz.
        # __setattr__ frozen'da bloklu, object.__setattr__ ile bypass ediyoruz.
        if not isinstance(self.samples, tuple):
            object.__setattr__(self, "samples", tuple(self.samples))


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
    (--force ile yeniden koşma), diğer modüllerinkine dokunulmaz (resume uyumu).

    Yazım atomiktir (pid'li geçici dosya + fsync + os.replace): yarıda kesilen
    bir süreç asla yarım/bozuk gates.json bırakmaz. gates.json çökme sonrası
    tanı raporu için okunduğundan bu, sessiz veri kaybına karşı bağlayıcıdır.
    """
    quality_dir = Path(run_dir) / QUALITY_DIR
    quality_dir.mkdir(parents=True, exist_ok=True)
    path = quality_dir / GATES_FILE

    existing: list[dict] = []
    if path.exists():
        try:
            existing = json.loads(path.read_text()).get("gates", [])
        except (json.JSONDecodeError, OSError) as exc:
            # Bozuk dosyayı SESSİZCE atmak, önceki modüllerin sonuçlarını
            # geri getirilemez şekilde siler (bkz. spec: crash-survival).
            # Onun yerine kenara alıp yüksek sesle bildiriyoruz; adli iz kalsın.
            corrupt_path = path.with_name(f"{path.name}.corrupt.{os.getpid()}")
            os.replace(path, corrupt_path)
            warnings.warn(
                f"gates.json was corrupt and could not be parsed ({exc}); "
                f"the damaged file was preserved at {corrupt_path} and a fresh "
                "gates.json is being started. Any gate results recorded only in "
                "the damaged file are lost — inspect it manually if needed.",
                stacklevel=2,
            )
            existing = []

    modules = {r.module for r in results}
    kept = [g for g in existing if g.get("module") not in modules]
    payload = {
        "updated_at": datetime.now().isoformat(timespec="seconds"),
        "gates": kept + [asdict(r) for r in results],
    }

    # state.py'deki _write() ile aynı desen: pid'li tmp dosya + fsync + os.replace.
    # Sabit adlı .tmp yerine pid'li: aynı dizine iki süreç yazarsa çakışmasın.
    tmp = path.with_name(f"{path.name}.tmp.{os.getpid()}")
    with tmp.open("w") as handle:
        json.dump(payload, handle, indent=2)
        handle.flush()
        os.fsync(handle.fileno())  # replace'ten önce diske insin
    os.replace(tmp, path)  # atomik: yarıda kesilme olursa eski dosya olduğu gibi kalır
    return path
