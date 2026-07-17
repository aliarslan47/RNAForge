"""m01 — Girdi doğrulama + platform tespiti.

Bu modül pipeline'ın kapısıdır: config, metadata ve FASTQ'lar burada
doğrulanır. Hata varsa BURADA durulur, sessiz devam yoktur (PLAN §13, Kural 7).
"""
from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

from rnaforge.config import Config
from rnaforge.metadata import Sample, load_metadata, validate_design
from rnaforge.platform import PlatformInfo, detect_platform, require_supported
from rnaforge.state import RunState

MODULE_NAME = "m01_validate"


def _check_reference(config: Config) -> None:
    fields = {
        "prokaryote": ("genome_fasta", "annotation_gff"),
        "eukaryote": ("transcriptome_fasta", "tx2gene"),
    }[config.organism_type]
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


def run_validation(config: Config, metadata_path: Path, run_dir: Path) -> dict:
    run_dir = Path(run_dir)
    logs_dir = run_dir / "logs"
    stats_dir = run_dir / "statistics"
    logs_dir.mkdir(parents=True, exist_ok=True)
    stats_dir.mkdir(parents=True, exist_ok=True)
    state = RunState(run_dir)
    lines: list[str] = []

    def log(message: str) -> None:
        lines.append(message)

    log(f"organism={config.organism} organism_type={config.organism_type}")
    _check_reference(config)
    log("reference files: OK")

    samples = load_metadata(metadata_path)
    log(f"metadata: {len(samples)} sample(s) loaded from {metadata_path}")

    validate_design(samples, config.de.design)
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

    if len(platforms) > 1:
        raise ValueError(
            f"samples come from mixed platforms: {', '.join(sorted(platforms))}. "
            "A single run must use one platform."
        )
    platform = platforms.pop()

    if config.platform != "auto" and config.platform != platform:
        raise ValueError(
            f"config says platform={config.platform!r} but the FASTQ files look like "
            f"{platform!r}. Fix the config, or set platform: auto."
        )

    conditions = dict(Counter(s.condition for s in samples))
    summary = {
        "organism": config.organism,
        "organism_type": config.organism_type,
        "platform": platform,
        "n_samples": len(samples),
        "conditions": conditions,
        "design": config.de.design,
        "samples": per_sample,
    }

    stats_path = stats_dir / "raw_statistics.json"
    stats_path.write_text(json.dumps(summary, indent=2))
    log(f"raw statistics written: {stats_path}")
    (logs_dir / "validation.log").write_text("\n".join(lines) + "\n")

    state.mark_done(MODULE_NAME, [str(stats_path)])
    return summary
