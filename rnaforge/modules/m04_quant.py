"""m04 — Quantification ROUTER (prokaryot: bowtie2 genom hizalama).

organism_type'a göre dallanır (PLAN §5). Şimdilik yalnız prokaryot yolu bağlı;
eukaryote (Salmon) net NotImplementedError verir. Veri kapısı `alignment_rate`:
bowtie2 overall alignment rate profil eşiğinin altındaysa FAIL — düşük hizalama
güvenilmez count matrisi üretir (PLAN §3)."""
from __future__ import annotations

from rnaforge.bowtie2 import AlignmentResult
from rnaforge.gates import FAIL, PASS, GateResult
from rnaforge.quality import Profile

MODULE_NAME = "m04_quant"
_GATE = "alignment_rate"


def build_alignment_gates(results: dict[str, AlignmentResult],
                          profile: Profile) -> list[GateResult]:
    threshold = profile.threshold(_GATE)
    offenders = sorted(sid for sid, r in results.items() if r.alignment_rate < threshold)
    lowest = min((r.alignment_rate for r in results.values()), default=1.0)
    overridden = _GATE in profile.overrides()
    if offenders:
        status = FAIL
        message = (
            f"hizalama oranı eşiğin altında ({len(offenders)} örnek: "
            f"{', '.join(offenders)}); en düşük {lowest:.2f} < {threshold:.2f}. "
            "Düşük hizalama güvenilmez count matrisi üretir."
        )
    else:
        status = PASS
        message = f"tüm örnekler alignment ≥ {threshold:.2f} (en düşük {lowest:.2f})."
    return [GateResult(
        name=_GATE, module=MODULE_NAME, status=status, message=message,
        remedy=("Referans genomu, kütüphane kimyasını ve tür kimliğini doğrulayın; "
                "yanlış referans veya kontaminasyon hizalamayı düşürür."),
        measured=lowest, threshold=threshold, overridden=overridden,
        samples=tuple(offenders),
    )]
