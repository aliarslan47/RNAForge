"""m03 — Nazik read trimming (fastp).

fastp adapter'ı temizler ve çok kısa okumaları eler AMA nazikçe (PLAN §4.2,
Williams 2016). m03'ün veri kapısı `survival_rate`'tir: trimming sonrası okumaların
çok azı hayatta kaldıysa kütüphane bozuktur → FAIL, koşu durur (m02'nin aksine)."""
from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

from rnaforge.chopper import run_chopper
from rnaforge.config import CHEMISTRY, Config
from rnaforge.fastp import FastpResult, run_fastp, trimmed_name
from rnaforge.gates import FAIL, PASS, WARN, GateResult, raise_if_failed, write_gate_results
from rnaforge.metadata import Sample, load_metadata
from rnaforge.pychopper import run_pychopper
from rnaforge.quality import Profile, load_profile, profile_name_for
from rnaforge.routing import resolve_read_type
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


def build_trim_gates(survival_rates: dict[str, float], profile: Profile,
                     warn_only: bool = False) -> list[GateResult]:
    """survival_rate kapısı. warn_only=True → eşiğin altı WARN (uzun-okuma: Pychopper
    tam-boy olmayanı atar, düşük survival şüpheli ama geçersiz değil)."""
    threshold = profile.threshold(_GATE)
    offenders = sorted(sid for sid, r in survival_rates.items() if r < threshold)
    lowest = min(survival_rates.values(), default=1.0)
    overridden = _GATE in profile.overrides()
    if offenders:
        status = WARN if warn_only else FAIL
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
    # read_type yönlendirmesi (m04 router deseni): kısa → fastp, uzun → Pychopper/chopper.
    read_type = resolve_read_type(run_dir)
    if read_type == "long":
        summary = _trim_long(config, metadata_path, run_dir,
                             trimmed_dir, stats_dir, logs_dir, state)
    else:
        summary = _trim_short(config, metadata_path, run_dir,
                             trimmed_dir, stats_dir, logs_dir, state)

    state.mark_done(MODULE_NAME, [str(stats_path), str(logs_dir / "trim.log")])
    return summary


def _trim_short(config: Config, metadata_path: Path, run_dir: Path,
                trimmed_dir: Path, stats_dir: Path, logs_dir: Path,
                state: RunState) -> dict:
    """Kısa-okuma trimming (fastp). survival_rate FAIL kapısı korunur."""
    stats_path = stats_dir / "trimming_statistics.json"
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

        gates = build_trim_gates(
            {sid: r.survival_rate for sid, r in results.items()}, profile)
        # Sıra (m01 deseni): stats yaz → gates yaz → EN SON raise.
        summary = {
            "read_type": "short",
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
    return summary


def _count_fastx(path) -> int:
    import gzip as _gz
    opener = _gz.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        return sum(1 for _ in fh) // 4


def _trim_long(config: Config, metadata_path: Path, run_dir: Path,
               trimmed_dir: Path, stats_dir: Path, logs_dir: Path,
               state: RunState) -> dict:
    """Uzun-okuma ön-işleme. cdna: Pychopper (yönlendir/kes) + chopper (filtre);
    direct_rna: yalnız chopper. Diagnostik — FAIL kapısı yok (long profil Step 6)."""
    chemistry = config.library.chemistry
    if chemistry not in CHEMISTRY:
        raise ValueError(
            f"long-read trim requires library.chemistry in {CHEMISTRY}, "
            f"got {chemistry!r} (m01 should have enforced this for ONT)."
        )
    stats_path = stats_dir / "trimming_statistics.json"
    log_path = logs_dir / "trim.log"
    min_len = config.trimming.min_length
    with log_path.open("w") as log_file:
        def log(msg: str) -> None:
            log_file.write(msg + "\n")
            log_file.flush()

        samples = load_metadata(metadata_path)
        log(f"m03 long-read ({chemistry}): {len(samples)} sample(s)")
        per_sample: dict[str, dict] = {}
        for sample in samples:
            state.heartbeat()
            out1, _ = trimmed_reads(run_dir, sample)
            out1.parent.mkdir(parents=True, exist_ok=True)
            reads_before = _count_fastx(sample.fastq_1)

            if chemistry == "cdna":
                work = out1.parent / "pychopper_full_length.fastq"
                stats_tsv = out1.parent / "pychopper_stats.tsv"
                ps = run_pychopper(sample.fastq_1, work, stats_tsv)
                reads_after = run_chopper(work, out1, min_len=min_len)
                log(f"{sample.sample_id}: pychopper primers={ps.primers_found} "
                    f"rescue={ps.rescue} unusable={ps.unusable}; chopper kept {reads_after}")
            else:  # direct_rna
                reads_after = run_chopper(sample.fastq_1, out1, min_len=min_len)
                log(f"{sample.sample_id}: chopper kept {reads_after}")

            survival = reads_after / reads_before if reads_before else 0.0
            per_sample[sample.sample_id] = {
                "reads_before": reads_before,
                "reads_after": reads_after,
                "survival_rate": round(survival, 4),
            }

        # Step 6: uzun-okuma survival WARN kapısı (prokaryote_long; asla FAIL —
        # Pychopper tam-boy olmayanı atar, düşük survival şüpheli ama geçersiz değil).
        profile = load_profile(profile_name_for(config.organism_type, "long"),
                               config.quality)
        survivals = {sid: v["survival_rate"] for sid, v in per_sample.items()}
        gates = build_trim_gates(survivals, profile, warn_only=True)
        summary = {
            "read_type": "long",
            "chemistry": chemistry,
            "n_samples": len(samples),
            "samples": per_sample,
            "gate_counts": dict(Counter(g.status for g in gates)),
        }
        stats_path.write_text(json.dumps(summary, indent=2))
        write_gate_results(run_dir, gates)
        for g in gates:
            log(f"gate {g.name}: {g.status} — {g.message}")
        log(f"trimming statistics written: {stats_path}")
    return summary
