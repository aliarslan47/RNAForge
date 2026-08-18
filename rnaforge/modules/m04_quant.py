"""m04 — Quantification ROUTER.

organism_type'a göre dallanır (PLAN §5): prokaryot → bowtie2 (kısa) / minimap2 (uzun);
eukaryote → Salmon (decoy-aware transkriptom niceleme). Veri kapısı `alignment_rate`:
prokaryot-kısa bowtie2 overall rate / ökaryot salmon mapping_rate profil eşiğinin
altındaysa FAIL — düşük hizalama güvenilmez count matrisi üretir (PLAN §3)."""
from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

from rnaforge.bowtie2 import AlignmentResult, build_index, run_bowtie2
from rnaforge.config import Config
from rnaforge.gates import FAIL, PASS, GateResult, raise_if_failed, write_gate_results
from rnaforge.metadata import load_metadata
from rnaforge.minimap2 import minimap2_preset, run_minimap2
from rnaforge.modules.m03_trim import trimmed_reads
from rnaforge.quality import Profile, load_profile, profile_name_for
from rnaforge.routing import resolve_platform, resolve_read_type
from rnaforge.salmon import build_salmon_index, run_salmon_quant
from rnaforge.state import RunState
from dataclasses import dataclass as _dataclass

MODULE_NAME = "m04_quant"
_GATE = "alignment_rate"


@_dataclass
class _MappingAdapter:
    """build_alignment_gates .alignment_rate bekler; salmon mapping_rate'i uyarlar."""
    alignment_rate: float


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


def run_quant(config: Config, metadata_path: Path, run_dir: Path,
              force: bool = False) -> dict:
    run_dir = Path(run_dir)
    quant_dir = run_dir / "quantification"
    stats_dir = run_dir / "statistics"
    logs_dir = run_dir / "logs"
    for d in (quant_dir, stats_dir, logs_dir):
        d.mkdir(parents=True, exist_ok=True)
    state = RunState(run_dir)
    stats_path = stats_dir / "alignment_statistics.json"

    if not force and state.is_done(MODULE_NAME) and stats_path.exists():
        summary = json.loads(stats_path.read_text())
        summary["resumed"] = True
        return summary

    if not state.is_done("m03_trim"):
        raise ValueError(
            "m04 (quant) requires m03 (trim) to have completed in this run directory "
            f"first: {run_dir}. Run `rnaforge trim` with the same --run-id, then re-run quant."
        )
    # ROUTER: önce organism_type (ökaryot → Salmon), sonra read_type (prokaryot: kısa
    # bowtie2 / uzun minimap2). Ayrım YALNIZ m04/m05'te; m06+ agnostik. Step-1'in
    # `require_short_read` muhafızı read_type dispatch'le değişti; muhafız m05'te kalır.
    if config.organism_type == "eukaryote":
        summary = _quant_euk(config, metadata_path, run_dir,
                             quant_dir, stats_dir, logs_dir, state)
    else:
        read_type = resolve_read_type(run_dir)
        if read_type == "long":
            summary = _quant_long(config, metadata_path, run_dir,
                                  quant_dir, stats_dir, logs_dir, state)
        else:
            summary = _quant_short(config, metadata_path, run_dir,
                                   quant_dir, stats_dir, logs_dir, state)

    state.mark_done(MODULE_NAME, [str(stats_path), str(logs_dir / "quant.log")])
    return summary


def _quant_short(config: Config, metadata_path: Path, run_dir: Path,
                 quant_dir: Path, stats_dir: Path, logs_dir: Path,
                 state: RunState) -> dict:
    """Kısa-okuma hizalama (bowtie2). alignment_rate FAIL kapısı korunur."""
    stats_path = stats_dir / "alignment_statistics.json"
    profile = load_profile(config.organism_type, config.quality)
    log_path = logs_dir / "quant.log"
    with log_path.open("w") as log_file:
        def log(msg: str) -> None:
            log_file.write(msg + "\n")
            log_file.flush()

        samples = load_metadata(metadata_path)
        index_prefix = build_index(config.reference.genome_fasta, quant_dir / "_index")
        log(f"m04 bowtie2: index built, {len(samples)} sample(s)")
        results = {}
        per_sample = {}
        for sample in samples:
            state.heartbeat()
            t1, t2 = trimmed_reads(run_dir, sample)
            result = run_bowtie2(index_prefix, quant_dir / sample.sample_id, t1,
                                 fastq_2=t2, threads=config.resources.threads)
            results[sample.sample_id] = result
            per_sample[sample.sample_id] = {
                "alignment_rate": result.alignment_rate, "bam": str(result.bam),
            }
            log(f"{sample.sample_id}: alignment_rate={result.alignment_rate:.3f}")

        gates = build_alignment_gates(results, profile)
        summary = {
            "read_type": "short",
            "n_samples": len(samples), "samples": per_sample,
            "gate_counts": dict(Counter(g.status for g in gates)),
        }
        stats_path.write_text(json.dumps(summary, indent=2))
        write_gate_results(run_dir, gates)
        for g in gates:
            log(f"gate {g.name}: {g.status} — {g.message}")
        raise_if_failed(gates)
        log(f"alignment statistics written: {stats_path}")
    return summary


def _quant_euk(config: Config, metadata_path: Path, run_dir: Path,
               quant_dir: Path, stats_dir: Path, logs_dir: Path,
               state: RunState) -> dict:
    """Ökaryot niceleme (Salmon, decoy-aware). mapping_rate → alignment_rate FAIL kapısı
    (eukaryote.yml permissive). Trimlenmiş okuma → salmon quant -l A."""
    stats_path = stats_dir / "alignment_statistics.json"
    profile = load_profile(config.organism_type, config.quality)
    log_path = logs_dir / "quant.log"
    with log_path.open("w") as log_file:
        def log(msg: str) -> None:
            log_file.write(msg + "\n")
            log_file.flush()

        samples = load_metadata(metadata_path)
        index_dir = build_salmon_index(
            config.reference.transcriptome_fasta, quant_dir / "_index",
            genome_fasta=config.reference.genome_fasta,
            threads=config.resources.threads, log=log)
        log(f"m04 salmon: index ready, {len(samples)} sample(s)")
        results = {}
        per_sample = {}
        for sample in samples:
            state.heartbeat()
            t1, t2 = trimmed_reads(run_dir, sample)
            q = run_salmon_quant(index_dir, quant_dir / sample.sample_id, t1,
                                 fastq_2=t2, threads=config.resources.threads)
            results[sample.sample_id] = _MappingAdapter(q.mapping_rate)
            per_sample[sample.sample_id] = {
                "mapping_rate": q.mapping_rate, "quant_sf": str(q.quant_sf)}
            log(f"{sample.sample_id}: mapping_rate={q.mapping_rate:.3f}")

        gates = build_alignment_gates(results, profile)
        summary = {
            "read_type": "short", "organism_type": "eukaryote",
            "n_samples": len(samples), "samples": per_sample,
            "gate_counts": dict(Counter(g.status for g in gates)),
        }
        stats_path.write_text(json.dumps(summary, indent=2))
        write_gate_results(run_dir, gates)
        for g in gates:
            log(f"gate {g.name}: {g.status} — {g.message}")
        raise_if_failed(gates)
        log(f"alignment statistics written: {stats_path}")
    return summary


def _quant_long(config: Config, metadata_path: Path, run_dir: Path,
                quant_dir: Path, stats_dir: Path, logs_dir: Path,
                state: RunState) -> dict:
    """Uzun-okuma hizalama (minimap2). Preset platformdan (ont→map-ont,
    pacbio_hifi→map-hifi). Diagnostik — FAIL kapısı yok (long profil Step 6);
    alignment_rate yalnız istatistik olarak kaydedilir."""
    platform = resolve_platform(run_dir)
    preset = minimap2_preset(platform)
    stats_path = stats_dir / "alignment_statistics.json"
    log_path = logs_dir / "quant.log"
    with log_path.open("w") as log_file:
        def log(msg: str) -> None:
            log_file.write(msg + "\n")
            log_file.flush()

        samples = load_metadata(metadata_path)
        log(f"m04 minimap2 ({platform} → -ax {preset}): {len(samples)} sample(s)")
        results = {}
        per_sample = {}
        for sample in samples:
            state.heartbeat()
            t1, _ = trimmed_reads(run_dir, sample)   # ONT tek-uçlu
            result = run_minimap2(config.reference.genome_fasta,
                                  quant_dir / sample.sample_id, t1,
                                  preset=preset, threads=config.resources.threads)
            results[sample.sample_id] = result
            per_sample[sample.sample_id] = {
                "alignment_rate": result.alignment_rate, "bam": str(result.bam),
            }
            log(f"{sample.sample_id}: alignment_rate={result.alignment_rate:.3f}")

        # Step 6: uzun-okuma alignment FAIL kapısı (prokaryote_long, eşik 0.50).
        # Yanlış referans/kontaminasyon her platformda GEÇERSİZ → FAIL (survival/assignment
        # WARN'ken bu FAIL: katastrofik hizalama gerçek bir geçersizlik sinyali).
        profile = load_profile(profile_name_for(config.organism_type, "long"),
                               config.quality)
        gates = build_alignment_gates(results, profile)
        summary = {
            "read_type": "long",
            "platform": platform,
            "n_samples": len(samples),
            "samples": per_sample,
            "gate_counts": dict(Counter(g.status for g in gates)),
        }
        stats_path.write_text(json.dumps(summary, indent=2))
        write_gate_results(run_dir, gates)
        for g in gates:
            log(f"gate {g.name}: {g.status} — {g.message}")
        raise_if_failed(gates)   # FAIL varsa BURADA durur (stats+gates zaten diskte)
        log(f"alignment statistics written: {stats_path}")
    return summary
