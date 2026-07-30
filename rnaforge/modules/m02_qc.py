"""m02 — Ham okuma kalite kontrolü (FastQC).

Ham QC diagnostiktir: kötü ham kalite BEKLENEN ve m03 trimming'in düzelttiği
şeydir. Bu yüzden m02 asla koşuyu durdurmaz — FastQC FAIL bile bizde WARN olur
(sonuç ŞÜPHELİ damgalanır, GEÇERSİZ değil). Bkz. spec 2026-07-30."""
from __future__ import annotations

from rnaforge.fastqc import FastQCReport
from rnaforge.gates import FAIL, PASS, WARN, GateResult

MODULE_NAME = "m02_qc"

# curated set: FastQC modül adı -> (kapı adı, WARN remedy)
_GATE_MAP = {
    "Per base sequence quality": (
        "per_base_quality",
        "m03 trimming (fastp) kalite profilini iyileştirir; sonra tekrar bak.",
    ),
    "Adapter Content": (
        "adapter_content",
        "m03 fastp adapter'ı otomatik temizler; adapter kontaminasyonu bekleniyor olabilir.",
    ),
    "Overrepresented sequences": (
        "overrepresented",
        "rRNA/adapter kontaminasyonu olabilir; kütüphane kimyasını (selection) doğrula.",
    ),
    "Per sequence GC content": (
        "gc_content",
        "beklenmedik GC → tür karışımı/kontaminasyon olabilir; girdiyi doğrula.",
    ),
}

# FastQC bayragi -> bizim durum. FAIL BILEREK WARN'a eslenir (m02 durdurmaz).
_FLAG_TO_STATUS = {"PASS": PASS, "WARN": WARN, "FAIL": WARN}
_SEVERITY = {PASS: 0, WARN: 1, FAIL: 2}


def build_qc_gates(reports: dict[str, FastQCReport]) -> list[GateResult]:
    gates: list[GateResult] = []
    for fastqc_module, (gate_name, remedy) in _GATE_MAP.items():
        offenders: list[str] = []
        status = PASS
        for sample_id, report in reports.items():
            flag = report.modules.get(fastqc_module, "PASS")
            mapped = _FLAG_TO_STATUS.get(flag, WARN)
            if mapped != PASS:
                offenders.append(sample_id)
            if _SEVERITY[mapped] > _SEVERITY[status]:
                status = mapped
        if status == PASS:
            message = f"FastQC '{fastqc_module}': tüm örnekler PASS."
        else:
            message = (
                f"FastQC '{fastqc_module}': {len(offenders)} örnek işaretlendi "
                f"({', '.join(sorted(offenders))}). Ham okuma; sonuç ŞÜPHELİ damgalanır."
            )
        gates.append(GateResult(
            name=gate_name, module=MODULE_NAME, status=status,
            message=message, remedy=remedy, samples=tuple(sorted(offenders)),
        ))
    assert all(g.status != FAIL for g in gates)  # sozlesme: m02 asla FAIL uretmez
    return gates
