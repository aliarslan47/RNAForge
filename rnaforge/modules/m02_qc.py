"""m02 — Ham okuma kalite kontrolü (FastQC).

Ham QC diagnostiktir: kötü ham kalite BEKLENEN ve m03 trimming'in düzelttiği
şeydir. Bu yüzden m02 asla koşuyu durdurmaz — FastQC FAIL bile bizde WARN olur
(sonuç ŞÜPHELİ damgalanır, GEÇERSİZ değil). Bkz. spec 2026-07-30."""
from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

from rnaforge.config import Config
from rnaforge.fastqc import FastQCReport, parse_fastqc_zip, run_fastqc
from rnaforge.gates import FAIL, PASS, WARN, GateResult, write_gate_results
from rnaforge.metadata import load_metadata
from rnaforge.state import RunState

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


def run_qc(config: Config, metadata_path: Path, run_dir: Path,
           force: bool = False) -> dict:
    run_dir = Path(run_dir)
    raw_qc_dir = run_dir / "raw_qc"
    stats_dir = run_dir / "statistics"
    logs_dir = run_dir / "logs"
    for d in (raw_qc_dir, stats_dir, logs_dir):
        d.mkdir(parents=True, exist_ok=True)
    state = RunState(run_dir)
    stats_path = stats_dir / "qc_statistics.json"

    if not force and state.is_done(MODULE_NAME) and stats_path.exists():
        summary = json.loads(stats_path.read_text())
        summary["resumed"] = True
        return summary

    # Ön koşul: m01 bu run_dir'de tamamlanmış olmalı. Aksi halde platform reddi
    # (ONT/PacBio) atlanmış olabilir ve FastQC'yi desteklenmeyen girdiye koşardık.
    if not state.is_done("m01_validate"):
        raise ValueError(
            "m02 (qc) requires m01 (validate) to have completed in this run "
            f"directory first: {run_dir}. Run `rnaforge validate` with the same "
            "--run-id, then re-run qc."
        )

    log_path = logs_dir / "qc.log"
    with log_path.open("w") as log_file:
        def log(msg: str) -> None:
            log_file.write(msg + "\n")
            log_file.flush()

        samples = load_metadata(metadata_path)
        log(f"m02 FastQC: {len(samples)} sample(s)")
        reports = {}
        module_flags = {}
        for sample in samples:
            state.heartbeat()
            sample_out = raw_qc_dir / sample.sample_id
            sample_out.mkdir(parents=True, exist_ok=True)
            zip_path = run_fastqc(sample.fastq_1, sample_out)
            report = parse_fastqc_zip(zip_path)
            reports[sample.sample_id] = report
            module_flags[sample.sample_id] = report.modules
            log(f"{sample.sample_id}: FastQC OK ({zip_path.name})")

        gates = build_qc_gates(reports)
        write_gate_results(run_dir, gates)
        for g in gates:
            log(f"gate {g.name}: {g.status} — {g.message}")

        gate_counts = dict(Counter(g.status for g in gates))
        summary = {
            "n_samples": len(samples),
            "samples": {sid: r.basic_stats for sid, r in reports.items()},
            "module_flags": module_flags,
            "gate_counts": gate_counts,
        }
        stats_path.write_text(json.dumps(summary, indent=2))
        log(f"qc statistics written: {stats_path}")

    state.mark_done(MODULE_NAME, [str(stats_path), str(log_path)])
    return summary
