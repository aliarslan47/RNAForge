"""read_type dispatch seam. m01 detects read_type and writes it into the run's
raw_statistics.json (single source of truth); the per-stage routers read it back
here. Long-read stages are not built yet — they must fail loudly, never run the
short-read tool on long reads (Rule 7 / feedback_gurultulu_hata)."""
from __future__ import annotations

import json
from pathlib import Path


def resolve_read_type(run_dir: Path | str) -> str:
    """Return the read_type m01 recorded for this run."""
    stats = Path(run_dir) / "statistics" / "raw_statistics.json"
    if not stats.exists():
        raise ValueError(
            f"cannot resolve read_type: {stats} not found. "
            "Run `rnaforge validate` (m01) with the same --run-id first."
        )
    data = json.loads(stats.read_text())
    read_type = data.get("read_type")
    if read_type is None:
        raise ValueError(
            f"cannot resolve read_type: no 'read_type' key in {stats}. "
            "Re-run `rnaforge validate` (m01) to regenerate it."
        )
    return read_type


def resolve_platform(run_dir: Path | str) -> str:
    """Return the platform m01 recorded for this run (ont|pacbio_hifi|illumina).

    Used by the long-read aligner (m04) to pick the minimap2 preset. Loud on a
    missing file / missing key — never guess a platform (Rule 7)."""
    stats = Path(run_dir) / "statistics" / "raw_statistics.json"
    if not stats.exists():
        raise ValueError(
            f"cannot resolve platform: {stats} not found. "
            "Run `rnaforge validate` (m01) with the same --run-id first."
        )
    data = json.loads(stats.read_text())
    platform = data.get("platform")
    if platform is None:
        raise ValueError(
            f"cannot resolve platform: no 'platform' key in {stats}. "
            "Re-run `rnaforge validate` (m01) to regenerate it."
        )
    return platform


def require_short_read(run_dir: Path | str, stage: str) -> None:
    """Guard at a short-read-only stage. No-op for short; loud stop for long."""
    read_type = resolve_read_type(run_dir)
    if read_type == "long":
        raise NotImplementedError(
            f"long-read {stage} is not implemented yet (long-read arm Step 1 only "
            "routes; NanoPlot/Pychopper/minimap2/featureCounts -L come in later "
            f"steps). This run is read_type={read_type!r}."
        )
