"""FASTQ'dan platform tespiti (saf stdlib — harici araç gerekmez).

Eşik mantığı ali-wgs-pipeline/ali_wgs/detect.py'den uyarlandı.
DİKKAT: Kütüphane kimyası (stranded / rRNA-polyA) FASTQ'da YOKTUR;
tespit edilemez, config'ten gelir (PLAN §4.1).
"""
from __future__ import annotations

import gzip
from dataclasses import dataclass
from pathlib import Path

SUPPORTED_PLATFORMS = ("illumina", "ont", "pacbio_hifi")

READ_TYPES = ("short", "long")

_PLATFORM_READ_TYPE = {
    "illumina": "short",
    "ont": "long",
    "pacbio_hifi": "long",
}


def read_type_for(platform: str) -> str:
    """Detected platform -> read_type. 'unknown' has no route (Rule 7)."""
    try:
        return _PLATFORM_READ_TYPE[platform]
    except KeyError:
        raise ValueError(
            f"cannot derive read_type for platform {platform!r}; "
            f"known platforms: {', '.join(_PLATFORM_READ_TYPE)}"
        ) from None


class UnsupportedPlatformError(RuntimeError):
    """Girdi tespit edildi ama MVP'de desteklenmiyor (PLAN Kural 7)."""


@dataclass(frozen=True)
class PlatformInfo:
    platform: str  # illumina | ont | pacbio_hifi | unknown
    mean_read_length: float
    n50: int
    mean_quality: float
    n_reads_sampled: int


def _is_gzip(path: Path) -> bool:
    """gzip'i UZANTIDAN değil İÇERİKTEN (magic byte 1f 8b) tespit et: yanlış adlandırılmış
    dosya (gz içerik/.fastq ad, ya da düz metin/.gz ad) sessiz çöp ya da BadGzipFile
    çökmesi üretmesin. Boş dosya → gzip değil (metin olarak okunur → 'unknown' → reddedilir)."""
    with open(path, "rb") as fh:
        return fh.read(2) == b"\x1f\x8b"


def _open(path: Path):
    return gzip.open(path, "rt") if _is_gzip(path) else open(path, "rt")


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
    """Refuse only input we cannot identify. Long reads (ONT/PacBio) are routed
    by read_type downstream; unidentifiable input has no safe route (Rule 7)."""
    if info.platform in SUPPORTED_PLATFORMS:
        return
    raise UnsupportedPlatformError(
        f"could not identify the sequencing platform for this input "
        f"(detected {info.platform!r}; supported: {', '.join(SUPPORTED_PLATFORMS)}).\n"
        f"  file: {fastq}\n"
        f"  mean read length: {info.mean_read_length}, N50: {info.n50}, "
        f"mean quality: {info.mean_quality}, reads sampled: {info.n_reads_sampled}\n"
        f"Running any route on unidentifiable reads would produce wrong results, "
        f"so it is refused."
    )
