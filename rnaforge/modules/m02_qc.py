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
from rnaforge.qcplots import QCPlotError, render_qc_figure
from rnaforge.quality import Profile, load_profile
from rnaforge.routing import require_short_read
from rnaforge.state import RunState

MODULE_NAME = "m02_qc"
_DEDUP_GATE = "dedup_fraction"

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


def build_dedup_gate(dedup_by_sample: dict[str, float | None],
                     profile: Profile) -> GateResult:
    """Duplikasyon WARN kapısı. `deduplication` = benzersiz okuma yüzdesi (FastQC);
    fraksiyon eşiğin altındaysa (aşırı duplikasyon) WARN — m02 asla FAIL üretmez.
    RNA-seq'te bir miktar duplikasyon beklenir; eşik bilinçle gevşektir."""
    thr = profile.threshold(_DEDUP_GATE)
    fracs = {sid: v / 100.0 for sid, v in dedup_by_sample.items() if v is not None}
    mean = sum(fracs.values()) / len(fracs) if fracs else 1.0
    offenders = sorted(sid for sid, f in fracs.items() if f < thr)
    overridden = profile.is_overridden(_DEDUP_GATE)
    if offenders:
        status = WARN
        message = (f"benzersiz-okuma fraksiyonu eşiğin altında ({len(offenders)} örnek: "
                   f"{', '.join(offenders)}); ortalama {mean:.1%} < {thr:.0%}. Yüksek "
                   "duplikasyon; sonuç ŞÜPHELİ damgalanabilir.")
    else:
        status = PASS
        message = (f"benzersiz-okuma fraksiyonu tüm örneklerde ≥ {thr:.0%} "
                   f"(ortalama {mean:.1%}).")
    return GateResult(
        name=_DEDUP_GATE, module=MODULE_NAME, status=status, message=message,
        remedy=("Yüksek duplikasyon PCR fazlası veya düşük kütüphane karmaşıklığı olabilir; "
                "girdi miktarı/PCR döngülerini gözden geçirin. RNA-seq'te yüksek ekspresyon "
                "nedeniyle bir miktar duplikasyon normaldir."),
        measured=round(mean, 4), threshold=thr, overridden=overridden,
        samples=tuple(offenders))


def mean_per_base_composition(
        reports: dict[str, FastQCReport]) -> tuple[list[str], dict[str, list[float]]]:
    """Örnekler arası pozisyon-başına ortalama A/T/G/C. En kısa profil uzunluğuna
    hizalanır (aynı okuma uzunluğunda FastQC binning'i aynıdır). Veri yoksa ([],{})."""
    profiles = [r.per_base_content for r in reports.values() if r.per_base_content]
    if not profiles:
        return [], {}
    n = min(len(p) for p in profiles)
    labels = [profiles[0][i][0] for i in range(n)]
    bases = ["A", "T", "G", "C"]
    means = {b: [] for b in bases}
    for i in range(n):
        for b in bases:
            vals = [p[i][1].get(b, 0.0) for p in profiles]
            means[b].append(round(sum(vals) / len(vals), 3))
    return labels, means


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
    require_short_read(run_dir, "qc")  # long-read QC (NanoPlot) not built yet

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

        profile = load_profile(config.organism_type, config.quality)
        gates = build_qc_gates(reports)
        dedup_by_sample = {sid: r.deduplication for sid, r in reports.items()}
        gates.append(build_dedup_gate(dedup_by_sample, profile))
        write_gate_results(run_dir, gates)
        for g in gates:
            log(f"gate {g.name}: {g.status} — {g.message}")

        # Per-base baz kompozisyonu figürü (F1) — best-effort; başarısızlık SESSİZCE
        # yutulmaz: log'a ve summary'ye yazılır (m02 diagnostik, bunun için durmaz).
        labels, comp = mean_per_base_composition(reports)
        composition_figure = None
        composition_error = None
        if labels:
            figures_dir = run_dir / "figures"
            fig_path = figures_dir / "qc_per_base_composition.png"
            spec = {
                "type": "lines", "title": "Per-base baz kompozisyonu (örnek ortalaması)",
                "xlabel": "Okuma pozisyonu (bç)", "ylabel": "%",
                "x": labels, "series": comp,
            }
            try:
                render_qc_figure(spec, fig_path)
                composition_figure = fig_path.name
                log(f"per-base composition figure written: {fig_path}")
            except QCPlotError as exc:
                composition_error = str(exc)
                log(f"WARNING: per-base composition figure FAILED (diagnostik atlandı): {exc}")

        gate_counts = dict(Counter(g.status for g in gates))
        summary = {
            "n_samples": len(samples),
            "samples": {sid: r.basic_stats for sid, r in reports.items()},
            "module_flags": module_flags,
            "gate_counts": gate_counts,
            "deduplication": {sid: v for sid, v in dedup_by_sample.items()},
            "mean_dedup_fraction": round(
                sum(v for v in dedup_by_sample.values() if v is not None)
                / max(1, sum(1 for v in dedup_by_sample.values() if v is not None)) / 100.0, 4),
            "per_base_composition_figure": composition_figure,
            "per_base_composition_error": composition_error,
        }
        stats_path.write_text(json.dumps(summary, indent=2))
        log(f"qc statistics written: {stats_path}")

    state.mark_done(MODULE_NAME, [str(stats_path), str(log_path)])
    return summary
