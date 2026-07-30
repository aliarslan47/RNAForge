"""m03 — Nazik read trimming (fastp).

fastp adapter'ı temizler ve çok kısa okumaları eler AMA nazikçe (PLAN §4.2,
Williams 2016). m03'ün veri kapısı `survival_rate`'tir: trimming sonrası okumaların
çok azı hayatta kaldıysa kütüphane bozuktur → FAIL, koşu durur (m02'nin aksine)."""
from __future__ import annotations

from rnaforge.fastp import FastpResult
from rnaforge.gates import FAIL, PASS, GateResult
from rnaforge.quality import Profile

MODULE_NAME = "m03_trim"
_GATE = "survival_rate"


def build_trim_gates(results: dict[str, FastpResult], profile: Profile) -> list[GateResult]:
    threshold = profile.threshold(_GATE)
    offenders = sorted(sid for sid, r in results.items() if r.survival_rate < threshold)
    lowest = min((r.survival_rate for r in results.values()), default=1.0)
    overridden = _GATE in profile.overrides()
    if offenders:
        status = FAIL
        message = (
            f"trimming sonrası survival eşiğin altında ({len(offenders)} örnek: "
            f"{', '.join(offenders)}); en düşük {lowest:.2f} < {threshold:.2f}. "
            "Okumaların çoğu min_length filtresini geçemedi."
        )
    else:
        status = PASS
        message = f"tüm örnekler survival ≥ {threshold:.2f} (en düşük {lowest:.2f})."
    return [GateResult(
        name=_GATE, module=MODULE_NAME, status=status, message=message,
        remedy=("Bu örneklerde okumaların çoğu çok kısa: okuma uzunluğunu ve platformu "
                "doğrulayın, yanlış/bozuk veri ya da beklenenden kısa okumalar olabilir."),
        measured=lowest, threshold=threshold, overridden=overridden,
        samples=tuple(offenders),
    )]
