"""m06 — Differential Expression (DESeq2).

m05 count matrisini R/Bioconductor DESeq2 ile diferansiyel ekspresyona çevirir —
pipeline'ın biyolojik çıktısı. İlk ORTAK (organizma-agnostik) analiz adımı.
Veri kapısı `replicate_correlation`: koşul-içi replikalar zayıf korele ise WARN
(sonuç ŞÜPHELİ damgalanır ama ÜRETİLİR — düşük korelasyon DE'yi geçersiz kılmaz,
gücü düşürür). m06 asla FAIL üretmez."""
from __future__ import annotations

from rnaforge.gates import PASS, WARN, GateResult
from rnaforge.quality import Profile

MODULE_NAME = "m06_de"
_GATE = "replicate_correlation"


def build_de_gates(min_correlation: float, profile: Profile) -> list[GateResult]:
    threshold = profile.threshold(_GATE)
    overridden = _GATE in profile.overrides()
    if min_correlation < threshold:
        status = WARN
        message = (
            f"koşul-içi replika korelasyonu düşük (min {min_correlation:.2f} < "
            f"{threshold:.2f}). DE üretildi ama ŞÜPHELİ: replikalar zayıf kümeleniyor "
            "(olası aykırı örnek / batch etkisi)."
        )
    else:
        status = PASS
        message = f"replika korelasyonu yeterli (min {min_correlation:.2f} ≥ {threshold:.2f})."
    return [GateResult(
        name=_GATE, module=MODULE_NAME, status=status, message=message,
        remedy=("PCA/heatmap ile aykırı örnek arayın; batch/covariate varsa design formülüne "
                "ekleyin (`~batch + condition`). Düşük korelasyon DE gücünü düşürür."),
        measured=min_correlation, threshold=threshold, overridden=overridden,
    )]
