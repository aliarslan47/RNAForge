"""rRNA depletion testi: SortMeRNA --other ile rRNA'sız FASTQ üretimi."""
from __future__ import annotations

import gzip
import shutil
from pathlib import Path

import pytest

from rnaforge.rrna_deplete import parse_depletion_rate, run_sortmerna_deplete


def test_parse_depletion_rate(tmp_path):
    """rRNA fraksiyonunu aligned.log'tan çıkar."""
    log = tmp_path / "aligned.log"
    log.write_text(
        "Total reads = 1000\n"
        "Total reads passing E-value threshold = 120 (12.00)\n")
    # depletion_rate = rRNA fraksiyonu = çıkarılan pay
    assert abs(parse_depletion_rate(log) - 0.12) < 1e-6


def test_parse_depletion_rate_missing(tmp_path):
    """Eksik log dosyası → 0.0 döndürür."""
    log = tmp_path / "x.log"
    log.write_text("garbage")
    assert parse_depletion_rate(log) == 0.0


# Integration test with skip guard
_HAS_CONDA = shutil.which("conda") is not None


@pytest.mark.skipif(not _HAS_CONDA, reason="conda/sortmerna not available")
def test_run_sortmerna_deplete_tiny_fastq(tmp_path):
    """Sentetik rRNA ve non-rRNA FASTQ ile depletion'u doğrula."""
    # Sentetik rrna_db FASTA (min 19bp SortMeRNA seed uzunluğu)
    rrna_db = tmp_path / "rrna.fa"
    rrna_db.write_text(">rrna_16s\nACGTACGTACGTACGTACGTACGTACGT\n")

    # Sentetik okuma: %50 rRNA, %50 non-rRNA (min 19bp)
    reads_fq = tmp_path / "reads.fq"
    with open(reads_fq, "w") as f:
        # rRNA okuma (rRNA'ya benzer)
        f.write("@rrna_read\nACGTACGTACGTACGTACGTACGTACGT\n+\nI" * 28 + "\n")
        # non-rRNA okuma (rRNA'dan farklı)
        f.write("@other_read\nTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n+\nI" * 28 + "\n")

    workdir = tmp_path / "sortmerna_work"
    result = run_sortmerna_deplete([reads_fq], rrna_db, workdir, paired=False, threads=2)

    # Kontrol: other FASTQ var ve depletion_rate [0,1] arasında
    assert "other" in result
    assert len(result["other"]) > 0
    assert all(p.exists() for p in result["other"])
    assert 0.0 <= result["depletion_rate"] <= 1.0
    assert "aligned_log" in result
    assert result["aligned_log"].exists()
