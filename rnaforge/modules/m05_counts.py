"""m05 — Count Matrix (prokaryot: featureCounts).

m04 BAM'lerini anotasyona göre sayıp gen×örnek count matrisi (ortak sözleşme,
PLAN §5) üretir. Veri kapısı `assignment_rate`: featureCounts'un gene atadığı
okuma oranı profil eşiğinin altındaysa FAIL — çok düşük atama yanlış anotasyon/
tür demektir, sayımlar güvenilmez."""
from __future__ import annotations

from rnaforge.gates import FAIL, PASS, GateResult
from rnaforge.quality import Profile

MODULE_NAME = "m05_counts"
_GATE = "assignment_rate"


def build_count_gates(assignment_rates: dict[str, float],
                      profile: Profile) -> list[GateResult]:
    threshold = profile.threshold(_GATE)
    offenders = sorted(sid for sid, r in assignment_rates.items() if r < threshold)
    lowest = min(assignment_rates.values(), default=1.0)
    overridden = _GATE in profile.overrides()
    if offenders:
        status = FAIL
        message = (
            f"gene atama oranı eşiğin altında ({len(offenders)} örnek: "
            f"{', '.join(offenders)}); en düşük {lowest:.2f} < {threshold:.2f}. "
            "Düşük atama yanlış anotasyon/tür → güvenilmez sayımlar."
        )
    else:
        status = PASS
        message = f"tüm örnekler assignment ≥ {threshold:.2f} (en düşük {lowest:.2f})."
    return [GateResult(
        name=_GATE, module=MODULE_NAME, status=status, message=message,
        remedy=("Anotasyon (GFF/GTF) ile referans genomun eşleştiğini ve feature_type/"
                "attribute config'inin anotasyon formatına uyduğunu doğrulayın."),
        measured=lowest, threshold=threshold, overridden=overridden,
        samples=tuple(offenders),
    )]
