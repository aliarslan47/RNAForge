"""Güvence kartı: koşunun sonucuna ne kadar güvenilebileceğinin tek sayfalık özeti.

PASS alan koşuda da üretilir — müşteri NEYİN kontrol edildiğini görmelidir.
Görünmeyen güvence, güvence değildir (spec 2026-07-20 §3.4).
"""
from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path

from rnaforge.gates import FAIL, GATES_FILE, PASS, QUALITY_DIR, WARN
from rnaforge.quality import Profile

CARD_FILE = "confidence_card.json"

TRUSTWORTHY = "TRUSTWORTHY"
SUSPECT = "SUSPECT"
INVALID = "INVALID"


def build_confidence_card(run_dir: Path | str, profile: Profile) -> dict:
    gates_path = Path(run_dir) / QUALITY_DIR / GATES_FILE
    gates = json.loads(gates_path.read_text())["gates"] if gates_path.exists() else []

    counts = {
        PASS: sum(1 for g in gates if g["status"] == PASS),
        WARN: sum(1 for g in gates if g["status"] == WARN),
        FAIL: sum(1 for g in gates if g["status"] == FAIL),
    }
    if counts[FAIL]:
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
    path.write_text(json.dumps(build_confidence_card(run_dir, profile), indent=2))
    return path
