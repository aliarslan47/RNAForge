"""m00 — Basecalling (ONT ham sinyal FAST5/POD5 → FASTQ).

Pipeline'ın girdisi FASTQ'dur. Ham sinyal geldiğinde (metadata fastq_1 bir POD5/FAST5
dosyası ya da dizini) m00 bunu dorado (GPU) ile basecall eder, fastq_1'i üretilen FASTQ'ya
yeniden yönlendiren çözülmüş metadata yazar; m01 varsa bunu tercih eder. FASTQ örnekleri
aynen geçer (passthrough). Diagnostik — FAIL kapısı yok (basecalling ya read üretir ya da
yüksek sesle hata verir). m00 OPSİYONEL: yalnız ham sinyal varsa `rnaforge basecall` koşulur;
saf FASTQ akışında hiç çağrılmaz (o zaman m01 kullanıcının metadata'sını kullanır)."""
from __future__ import annotations

import json
from pathlib import Path

from rnaforge.basecall import (
    basecalled_metadata_path,
    convert_fast5_to_pod5,
    is_signal_input,
    run_dorado,
)
from rnaforge.config import Config
from rnaforge.metadata import Sample, load_metadata
from rnaforge.state import RunState

MODULE_NAME = "m00_basecall"

_COLUMNS = ("sample_id", "condition", "fastq_1", "fastq_2", "batch", "subject")


def _write_resolved_metadata(path: Path, rows: list[Sample]) -> None:
    """Çözülmüş metadata MUTLAK yollarla yazılır: load_metadata yolları metadata
    dosyasının dizinine göre çözer; göreli yol yazmak (run_dir göreli olduğunda)
    yolu ikilerdi (canlı e2e'de yakalandı).

    Keyfi kovaryat sütunları (sex, lane, genotype...) korunur (Faz 3) — düşürmek ONT
    yolunda kovaryat design'larını sessizce bozardı. Kovaryat sütun birleşimi tüm
    örneklerden toplanır; bir örnekte yoksa boş yazılır."""
    covariate_cols: list[str] = []
    for s in rows:
        for key in s.covariates:
            if key not in covariate_cols:
                covariate_cols.append(key)
    columns = list(_COLUMNS) + covariate_cols
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fh:
        fh.write("\t".join(columns) + "\n")
        for s in rows:
            fh.write("\t".join([
                s.sample_id, s.condition, str(Path(s.fastq_1).resolve()),
                str(Path(s.fastq_2).resolve()) if s.fastq_2 else "",
                s.batch or "", s.subject or "",
                *[s.covariates.get(c, "") for c in covariate_cols],
            ]) + "\n")


def run_basecall(config: Config, metadata_path: Path, run_dir: Path,
                 force: bool = False) -> dict:
    run_dir = Path(run_dir)
    bc_dir = run_dir / "basecalled"
    stats_dir = run_dir / "statistics"
    logs_dir = run_dir / "logs"
    for d in (bc_dir, stats_dir, logs_dir):
        d.mkdir(parents=True, exist_ok=True)
    state = RunState(run_dir)
    stats_path = stats_dir / "basecall_statistics.json"
    resolved_path = basecalled_metadata_path(run_dir)

    if not force and state.is_done(MODULE_NAME) and stats_path.exists() and resolved_path.exists():
        summary = json.loads(stats_path.read_text())
        summary["resumed"] = True
        return summary

    bc = config.basecall
    models_dir = Path(bc.models_dir) if bc.models_dir else bc_dir / "_models"
    log_path = logs_dir / "basecall.log"
    with log_path.open("w") as log_file:
        def log(msg: str) -> None:
            log_file.write(msg + "\n")
            log_file.flush()

        samples = load_metadata(metadata_path)
        log(f"m00 basecall: {len(samples)} sample(s), model={bc.model}, device={bc.device}")
        per_sample: dict[str, dict] = {}
        resolved: list[Sample] = []
        for sample in samples:
            state.heartbeat()
            kind = is_signal_input(sample.fastq_1)
            if kind is None:
                # FASTQ passthrough — basecalling gerekmez.
                per_sample[sample.sample_id] = {"input_kind": "fastq", "reads": None}
                resolved.append(sample)
                log(f"{sample.sample_id}: FASTQ passthrough ({sample.fastq_1})")
                continue

            sample_dir = bc_dir / sample.sample_id
            pod5 = sample.fastq_1
            if kind == "fast5":
                pod5 = convert_fast5_to_pod5(
                    sample.fastq_1, sample_dir / "converted.pod5", env=bc.env)
                log(f"{sample.sample_id}: FAST5 → POD5 dönüştürüldü")
            out_fastq = sample_dir / f"{sample.sample_id}.fastq"
            reads = run_dorado(pod5, out_fastq, dorado_bin=bc.dorado_bin,
                               model=bc.model, device=bc.device, models_dir=models_dir)
            per_sample[sample.sample_id] = {"input_kind": kind, "reads": reads}
            log(f"{sample.sample_id}: {kind} → basecalled {reads} reads ({out_fastq})")
            # ONT tek-uçlu; fastq_1 basecall çıktısına yönlendirilir, fastq_2 yok.
            resolved.append(Sample(sample.sample_id, sample.condition, out_fastq,
                                   None, sample.batch, sample.subject,
                                   sample.covariates))

        _write_resolved_metadata(resolved_path, resolved)
        summary = {
            "n_samples": len(samples),
            "model": bc.model,
            "device": bc.device,
            "samples": per_sample,
            "resolved_metadata": str(resolved_path),
        }
        stats_path.write_text(json.dumps(summary, indent=2))
        log(f"resolved metadata written: {resolved_path}")

    state.mark_done(MODULE_NAME, [str(stats_path), str(resolved_path), str(log_path)])
    return summary
