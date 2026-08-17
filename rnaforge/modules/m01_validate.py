"""m01 — Girdi doğrulama + platform tespiti.

Bu modül pipeline'ın kapısıdır: config, metadata ve FASTQ'lar burada
doğrulanır. Hata varsa BURADA durulur, sessiz devam yoktur (PLAN §13, Kural 7).
"""
from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

from rnaforge.config import REQUIRED_REFERENCE, Config
from rnaforge.gates import raise_if_failed, write_gate_results
from rnaforge.basecall import basecalled_metadata_path
from rnaforge.metadata import Sample, load_metadata, validate_design
from rnaforge.platform import PlatformInfo, detect_platform, read_type_for, require_supported
from rnaforge.quality import load_profile
from rnaforge.state import RunState

MODULE_NAME = "m01_validate"


def _check_reference(config: Config) -> None:
    # Yönlendirme sözleşmesi TEK kaynaktan gelir (config.REQUIRED_REFERENCE);
    # buradaki ikinci bir kopya sessizce sürüklenirdi.
    fields = REQUIRED_REFERENCE[config.organism_type]
    for field in fields:
        path = getattr(config.reference, field)
        if not path.exists():
            raise FileNotFoundError(
                f"reference.{field} does not exist: {path} "
                f"(required for organism_type={config.organism_type})"
            )


def _sample_stats(sample: Sample, info: PlatformInfo) -> dict:
    return {
        "sample_id": sample.sample_id,
        "condition": sample.condition,
        "batch": sample.batch,
        "paired": sample.fastq_2 is not None,
        "platform": info.platform,
        "mean_read_length": info.mean_read_length,
        "mean_quality": info.mean_quality,
        "n_reads_sampled": info.n_reads_sampled,
    }


def run_validation(
    config: Config, metadata_path: Path, run_dir: Path, force: bool = False
) -> dict:
    run_dir = Path(run_dir)
    # m00 (basecall) ham sinyali FASTQ'ya çevirdiyse çözülmüş metadata'yı tercih et
    # (ham-sinyal → FASTQ handoff sözleşmesi). Yoksa kullanıcının metadata'sı (değişmez).
    _resolved = basecalled_metadata_path(run_dir)
    if _resolved.exists():
        metadata_path = _resolved
    logs_dir = run_dir / "logs"
    stats_dir = run_dir / "statistics"
    logs_dir.mkdir(parents=True, exist_ok=True)
    stats_dir.mkdir(parents=True, exist_ok=True)
    state = RunState(run_dir)
    stats_path = stats_dir / "raw_statistics.json"
    log_path = logs_dir / "validation.log"

    if not force and state.is_done(MODULE_NAME) and stats_path.exists():
        # Resume: bu modül bu run dizininde zaten bitmiş. İşi tekrarlamak yerine
        # önceki özeti geri ver (PLAN §15). --force ile bilerek ezilebilir.
        summary = json.loads(stats_path.read_text())
        summary["resumed"] = True
        return summary

    # Log satır satır AÇIK dosyaya yazılır ve flush edilir: koşu hata verip düşerse
    # nedenini gösteren satırlar diskte kalmalı. Sadece başarıda yazmak, en çok
    # ihtiyaç duyulan anda log'u yok ediyordu.
    with log_path.open("w") as log_file:

        def log(message: str) -> None:
            log_file.write(message + "\n")
            log_file.flush()

        log(f"organism={config.organism} organism_type={config.organism_type}")
        _check_reference(config)
        log("reference files: OK")

        samples = load_metadata(metadata_path)
        log(f"metadata: {len(samples)} sample(s) loaded from {metadata_path}")

        profile = load_profile(config.organism_type, config.quality)
        log(f"quality profile: {profile.name} (permissive={profile.permissive})")

        design_gates = validate_design(samples, config.de.design, paired=config.paired)
        write_gate_results(run_dir, design_gates)
        for gate in design_gates:
            log(f"gate {gate.name}: {gate.status} — {gate.message}")
        # Kapılar ÖNCE yazılır, SONRA zorlanır: FAIL'de de gates.json diskte kalmalı,
        # teşhis raporu onu okuyacak (spec §3.5).
        raise_if_failed(design_gates)
        log(f"design formula {config.de.design!r}: OK")

        per_sample: list[dict] = []
        platforms: set[str] = set()
        for sample in samples:
            state.heartbeat()
            info = detect_platform(sample.fastq_1)
            require_supported(info, sample.fastq_1)  # desteklenmiyorsa BURADA durur
            platforms.add(info.platform)
            per_sample.append(_sample_stats(sample, info))
            log(f"{sample.sample_id}: platform={info.platform} "
                f"mean_read_length={info.mean_read_length}")

        detected = sorted(platforms)
        if config.platform != "auto":
            # Kullanıcı platformu AÇIKÇA belirtti → ona güvenilir (tespit yalnız tanısal).
            # Uzunluk-tabanlı tespit kısa ONT cDNA'yı (ör. Nano3P-seq 3'-uç yakalama)
            # illumina sanabilir; açık config bu yanlış-sınıflamayı ezer (config = otorite,
            # tespit = kolaylık). require_supported yine tanımlanamayanı örnek-başı reddeder.
            platform = config.platform
            if detected != [platform]:
                log(f"platform: config={platform!r} (tespit {detected} — config'e güvenildi, "
                    "kısa cDNA/amplikon okuma uzunluk-tabanlı tespiti yanıltabilir)")
        else:
            if len(platforms) > 1:
                raise ValueError(
                    f"samples come from mixed platforms: {', '.join(detected)}. "
                    "A single run must use one platform, or set platform explicitly in config."
                )
            platform = platforms.pop()

        read_type = read_type_for(platform)
        chemistry = config.library.chemistry
        # cDNA vs direct-RNA is undetectable from FASTQ but changes m03 (Pychopper).
        # ONT long reads must declare it; HiFi is inferred, short reads don't care.
        if read_type == "long" and platform == "ont" and chemistry is None:
            raise ValueError(
                "ONT long-read input requires library.chemistry to be set "
                "('cdna' or 'direct_rna'): it cannot be detected from the FASTQ "
                "and it changes preprocessing (cDNA needs Pychopper). "
                "Set library.chemistry in the config and re-run validate."
            )

        conditions = dict(Counter(s.condition for s in samples))
        summary = {
            "organism": config.organism,
            "organism_type": config.organism_type,
            "platform": platform,
            "read_type": read_type,
            "chemistry": chemistry,
            "n_samples": len(samples),
            "conditions": conditions,
            "design": config.de.design,
            "samples": per_sample,
            "quality_profile": profile.name,
            "permissive_profile": profile.permissive,
        }

        stats_path.write_text(json.dumps(summary, indent=2))
        log(f"raw statistics written: {stats_path}")

    state.mark_done(MODULE_NAME, [str(stats_path), str(log_path)])
    return summary
