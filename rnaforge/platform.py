"""FASTQ'dan platform tespiti (saf stdlib — harici araç gerekmez).

Eşik mantığı ali-wgs-pipeline/ali_wgs/detect.py'den uyarlandı.
DİKKAT: Kütüphane kimyası (stranded / rRNA-polyA) FASTQ'da YOKTUR;
tespit edilemez, config'ten gelir (PLAN §4.1).
"""
from __future__ import annotations

import gzip
from dataclasses import dataclass
from pathlib import Path

SUPPORTED_PLATFORMS = ("illumina",)


class UnsupportedPlatformError(RuntimeError):
    """Girdi tespit edildi ama MVP'de desteklenmiyor (PLAN Kural 7)."""


@dataclass(frozen=True)
class PlatformInfo:
    platform: str  # illumina | ont | pacbio_hifi | unknown
    mean_read_length: float
    n50: int
    mean_quality: float
    n_reads_sampled: int


def _open(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt")


def _sample_fastq(path: Path, max_reads: int) -> tuple[list[int], list[float]]:
    lengths: list[int] = []
    quals: list[float] = []
    with _open(path) as fh:
        for i, line in enumerate(fh):
            position = i % 4
            if position == 1:
                lengths.append(len(line.strip()))
            elif position == 3:
                q = line.strip()
                if q:
                    quals.append(sum(ord(c) - 33 for c in q) / len(q))
                if len(lengths) >= max_reads:
                    break
    return lengths, quals


def _n50(lengths: list[int]) -> int:
    if not lengths:
        return 0
    ordered = sorted(lengths, reverse=True)
    half, acc = sum(ordered) / 2, 0
    for length in ordered:
        acc += length
        if acc >= half:
            return length
    return ordered[-1]


def detect_platform(
    fastq: Path,
    short_read_max_len: int = 350,
    hifi_min_qual: float = 25.0,
    max_reads: int = 5000,
) -> PlatformInfo:
    lengths, quals = _sample_fastq(Path(fastq), max_reads)
    if not lengths:
        return PlatformInfo("unknown", 0.0, 0, 0.0, 0)

    mean_len = sum(lengths) / len(lengths)
    mean_q = sum(quals) / len(quals) if quals else 0.0
    n50 = _n50(lengths)

    if mean_len <= short_read_max_len:
        platform = "illumina"
    elif mean_q >= hifi_min_qual and n50 >= 5000:
        platform = "pacbio_hifi"
    else:
        platform = "ont"

    return PlatformInfo(
        platform=platform,
        mean_read_length=round(mean_len, 1),
        n50=n50,
        mean_quality=round(mean_q, 1),
        n_reads_sampled=len(lengths),
    )


def require_supported(info: PlatformInfo, fastq: Path) -> None:
    """Desteklenmeyen platformu net mesajla reddet. Sessiz devam YOK."""
    if info.platform in SUPPORTED_PLATFORMS:
        return
    raise UnsupportedPlatformError(
        f"detected platform {info.platform!r} is not supported in the MVP "
        f"(supported: {', '.join(SUPPORTED_PLATFORMS)}).\n"
        f"  file: {fastq}\n"
        f"  mean read length: {info.mean_read_length}, N50: {info.n50}, "
        f"mean quality: {info.mean_quality}, reads sampled: {info.n_reads_sampled}\n"
        f"Long-read support (ONT/PacBio) needs a different route (minimap2) "
        f"and is planned for a later phase. Running the Illumina route on this "
        f"input would produce wrong results, so it is refused."
    )
