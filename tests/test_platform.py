from __future__ import annotations

import pytest

from rnaforge.platform import (
    UnsupportedPlatformError,
    detect_platform,
    require_supported,
)
from tests.conftest import write_fastq


def test_short_reads_detected_as_illumina(illumina_fastq):
    info = detect_platform(illumina_fastq)
    assert info.platform == "illumina"
    assert info.mean_read_length == pytest.approx(150.0)
    assert info.n_reads_sampled == 200


def test_long_noisy_reads_detected_as_ont(ont_fastq):
    assert detect_platform(ont_fastq).platform == "ont"


def test_long_high_quality_reads_detected_as_pacbio(pacbio_fastq):
    assert detect_platform(pacbio_fastq).platform == "pacbio_hifi"


def test_gzipped_fastq_supported(tmp_path):
    path = write_fastq(tmp_path / "reads.fastq.gz", 100, 150, "I", gzipped=True)
    assert detect_platform(path).platform == "illumina"


def test_illumina_is_supported(illumina_fastq):
    require_supported(detect_platform(illumina_fastq), illumina_fastq)  # raise etmemeli


def test_ont_rejected_with_actionable_message(ont_fastq):
    """PLAN Kural 7: tespit etmek != desteklemek. Sessizce yanlis araçla koşulmaz."""
    info = detect_platform(ont_fastq)
    with pytest.raises(UnsupportedPlatformError) as exc:
        require_supported(info, ont_fastq)
    message = str(exc.value)
    assert "ont" in message.lower()
    assert "illumina" in message.lower()
    assert str(ont_fastq) in message


def test_pacbio_rejected(pacbio_fastq):
    with pytest.raises(UnsupportedPlatformError):
        require_supported(detect_platform(pacbio_fastq), pacbio_fastq)


def test_empty_fastq_is_unknown_and_rejected(tmp_path):
    path = tmp_path / "empty.fastq"
    path.write_text("")
    info = detect_platform(path)
    assert info.platform == "unknown"
    with pytest.raises(UnsupportedPlatformError):
        require_supported(info, path)


def test_read_type_for_maps_platforms():
    from rnaforge.platform import read_type_for
    assert read_type_for("illumina") == "short"
    assert read_type_for("ont") == "long"
    assert read_type_for("pacbio_hifi") == "long"


def test_read_type_for_rejects_unknown():
    from rnaforge.platform import read_type_for
    with pytest.raises(ValueError):
        read_type_for("unknown")
