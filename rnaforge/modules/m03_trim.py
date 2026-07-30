"""m03 — Nazik read trimming (fastp).

fastp adapter'ı temizler ve çok kısa okumaları eler AMA nazikçe (PLAN §4.2,
Williams 2016). m03'ün veri kapısı `survival_rate`'tir: trimming sonrası okumaların
çok azı hayatta kaldıysa kütüphane bozuktur → FAIL, koşu durur (m02'nin aksine)."""
from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

from rnaforge.config import Config
from rnaforge.fastp import FastpResult, run_fastp, trimmed_name
from rnaforge.gates import FAIL, PASS, GateResult, raise_if_failed, write_gate_results
from rnaforge.metadata import Sample, load_metadata
from rnaforge.quality import Profile, load_profile
from rnaforge.state import RunState

MODULE_NAME = "m03_trim"
_GATE = "survival_rate"


def trimmed_reads(run_dir, sample: Sample) -> tuple[Path, Path | None]:
    """m03'ün bir örnek için ürettiği trimlenmiş FASTQ yol(lar)ı. Adlandırma
    kuralının TEK kaynağı: m03 buraya yazar, m04 buradan okur (drift önlenir)."""
    d = Path(run_dir) / "trimmed" / sample.sample_id
    out1 = d / trimmed_name(sample.fastq_1)
    out2 = d / trimmed_name(sample.fastq_2) if sample.fastq_2 else None
    return out1, out2


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


def run_trim(config: Config, metadata_path: Path, run_dir: Path,
             force: bool = False) -> dict:
    run_dir = Path(run_dir)
    trimmed_dir = run_dir / "trimmed"
    stats_dir = run_dir / "statistics"
    logs_dir = run_dir / "logs"
    for d in (trimmed_dir, stats_dir, logs_dir):
        d.mkdir(parents=True, exist_ok=True)
    state = RunState(run_dir)
    stats_path = stats_dir / "trimming_statistics.json"

    if not force and state.is_done(MODULE_NAME) and stats_path.exists():
        summary = json.loads(stats_path.read_text())
        summary["resumed"] = True
        return summary

    if not state.is_done("m01_validate"):
        raise ValueError(
            "m03 (trim) requires m01 (validate) to have completed in this run "
            f"directory first: {run_dir}. Run `rnaforge validate` with the same "
            "--run-id, then re-run trim."
        )

    profile = load_profile(config.organism_type, config.quality)
    log_path = logs_dir / "trim.log"
    with log_path.open("w") as log_file:
        def log(msg: str) -> None:
            log_file.write(msg + "\n")
            log_file.flush()

        samples = load_metadata(metadata_path)
        log(f"m03 fastp: {len(samples)} sample(s), min_length={config.trimming.min_length}, "
            f"aggressive_quality={config.trimming.aggressive_quality}")
        results = {}
        per_sample = {}
        for sample in samples:
            state.heartbeat()
            sample_out = trimmed_dir / sample.sample_id
            result = run_fastp(
                sample.fastq_1, sample_out, min_length=config.trimming.min_length,
                fastq_2=sample.fastq_2,
                aggressive_quality=config.trimming.aggressive_quality,
            )
            results[sample.sample_id] = result
            per_sample[sample.sample_id] = {
                "reads_before": result.reads_before,
                "reads_after": result.reads_after,
                "survival_rate": result.survival_rate,
            }
            log(f"{sample.sample_id}: survival={result.survival_rate:.3f} "
                f"({result.reads_after}/{result.reads_before})")

        gates = build_trim_gates(results, profile)
        # Sıra (m01 deseni): stats yaz → gates yaz → EN SON raise.
        summary = {
            "n_samples": len(samples),
            "samples": per_sample,
            "gate_counts": dict(Counter(g.status for g in gates)),
        }
        stats_path.write_text(json.dumps(summary, indent=2))
        write_gate_results(run_dir, gates)
        for g in gates:
            log(f"gate {g.name}: {g.status} — {g.message}")
        raise_if_failed(gates)   # FAIL varsa BURADA durur (stats+gates zaten diskte)
        log(f"trimming statistics written: {stats_path}")

    state.mark_done(MODULE_NAME, [str(stats_path), str(log_path)])
    return summary
