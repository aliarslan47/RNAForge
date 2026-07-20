"""Güvence kartı: koşunun sonucuna ne kadar güvenilebileceğinin tek sayfalık özeti.

PASS alan koşuda da üretilir — müşteri NEYİN kontrol edildiğini görmelidir.
Görünmeyen güvence, güvence değildir (spec 2026-07-20 §3.4).
"""
from __future__ import annotations

import json
import os
import warnings
from datetime import datetime
from pathlib import Path

from rnaforge.gates import FAIL, GATES_FILE, PASS, QUALITY_DIR, WARN
from rnaforge.quality import Profile

CARD_FILE = "confidence_card.json"

TRUSTWORTHY = "TRUSTWORTHY"
SUSPECT = "SUSPECT"
INVALID = "INVALID"
UNKNOWN = "UNKNOWN"  # hic kapi kontrol edilmemis -> bu bir garanti DEGIL (finding 2)


def _load_gates(gates_path: Path) -> list[dict]:
    if not gates_path.exists():
        return []
    try:
        return json.loads(gates_path.read_text()).get("gates", [])
    except (json.JSONDecodeError, OSError) as exc:
        # gates.py'deki write_gate_results ile AYNI desen: bozuk dosyayi sessizce
        # yok saymak yerine kenara al ve yuksek sesle bildir (adli iz kalsin).
        corrupt_path = gates_path.with_name(f"{gates_path.name}.corrupt.{os.getpid()}")
        os.replace(gates_path, corrupt_path)
        warnings.warn(
            f"gates.json was corrupt and could not be parsed ({exc}); "
            f"the damaged file was preserved at {corrupt_path} and the confidence "
            "card is being built as if no gates were recorded (verdict=UNKNOWN). "
            "Inspect the damaged file manually if needed.",
            stacklevel=2,
        )
        return []


def build_confidence_card(run_dir: Path | str, profile: Profile) -> dict:
    gates_path = Path(run_dir) / QUALITY_DIR / GATES_FILE
    gates = _load_gates(gates_path)

    counts = {
        PASS: sum(1 for g in gates if g["status"] == PASS),
        WARN: sum(1 for g in gates if g["status"] == WARN),
        FAIL: sum(1 for g in gates if g["status"] == FAIL),
    }
    if not gates:
        # Sifir kapi kontrol edildi -> bu bir GARANTI degil, sadece bilgi eksikligi.
        # TRUSTWORTHY ile ayni gorunmemeli (finding 2).
        verdict = UNKNOWN
    elif counts[FAIL]:
        verdict = INVALID
    elif counts[WARN]:
        verdict = SUSPECT
    else:
        verdict = TRUSTWORTHY

    return {
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "verdict": verdict,
        "counts": counts,
        "profile": {
            "name": profile.name,
            "permissive": profile.permissive,
            "description": profile.description,
            "overrides": profile.overrides(),
        },
        "gates": gates,
    }


def write_confidence_card(run_dir: Path | str, profile: Profile) -> Path:
    path = Path(run_dir) / QUALITY_DIR / CARD_FILE
    path.parent.mkdir(parents=True, exist_ok=True)
    card = build_confidence_card(run_dir, profile)

    # gates.py'deki write_gate_results ile AYNI desen: pid'li tmp dosya + fsync +
    # os.replace. Bu dosya cokme sonrasi tani icin okunur; yariya yazilmis bir
    # kart, hicbir kart olmamasindan daha kotu bir yaniltici olurdu.
    tmp = path.with_name(f"{path.name}.tmp.{os.getpid()}")
    with tmp.open("w") as handle:
        json.dump(card, handle, indent=2)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(tmp, path)
    return path
