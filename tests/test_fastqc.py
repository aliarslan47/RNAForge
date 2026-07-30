from __future__ import annotations

import shutil
import zipfile
from pathlib import Path

import pytest

from rnaforge.fastqc import (
    FastQCParseError,
    FastQCReport,
    FastQCRunError,
    parse_fastqc_report,
    parse_fastqc_zip,
    run_fastqc,
)
from tests.conftest import write_fastq

SUMMARY = (
    "PASS\tBasic Statistics\tt.fastq\n"
    "FAIL\tPer base sequence quality\tt.fastq\n"
    "WARN\tPer sequence GC content\tt.fastq\n"
    "PASS\tAdapter Content\tt.fastq\n"
)

DATA = (
    ">>Basic Statistics\tpass\n"
    "#Measure\tValue\n"
    "Filename\tt.fastq\n"
    "Encoding\tSanger / Illumina 1.9\n"
    "Total Sequences\t500\n"
    "Sequences flagged as poor quality\t0\n"
    "Sequence length\t150\n"
    "%GC\t50\n"
    ">>END_MODULE\n"
    ">>Per base sequence quality\tfail\n"
    "#Base\tMean\n"
    "1\t30.0\n"
    ">>END_MODULE\n"
)


def test_parse_reads_module_flags():
    report = parse_fastqc_report(SUMMARY, DATA)
    assert report.modules["Per base sequence quality"] == "FAIL"
    assert report.modules["Per sequence GC content"] == "WARN"
    assert report.modules["Basic Statistics"] == "PASS"


def test_parse_reads_basic_stats():
    report = parse_fastqc_report(SUMMARY, DATA)
    assert report.basic_stats["Total Sequences"] == "500"
    assert report.basic_stats["%GC"] == "50"
    assert report.basic_stats["Sequence length"] == "150"
    assert report.basic_stats["Encoding"] == "Sanger / Illumina 1.9"


def test_parse_rejects_summary_with_unknown_status():
    bad = "GOOD\tBasic Statistics\tt.fastq\n"
    with pytest.raises(FastQCParseError, match="status"):
        parse_fastqc_report(bad, DATA)


def test_parse_rejects_missing_basic_statistics_module():
    data_without_basic = ">>Per base sequence quality\tfail\n>>END_MODULE\n"
    with pytest.raises(FastQCParseError, match="Basic Statistics"):
        parse_fastqc_report(SUMMARY, data_without_basic)


def _make_zip(tmp_path: Path, name: str = "t_fastqc") -> Path:
    inner = tmp_path / name
    inner.mkdir()
    (inner / "summary.txt").write_text(SUMMARY)
    (inner / "fastqc_data.txt").write_text(DATA)
    zip_path = tmp_path / f"{name}.zip"
    with zipfile.ZipFile(zip_path, "w") as zf:
        zf.write(inner / "summary.txt", f"{name}/summary.txt")
        zf.write(inner / "fastqc_data.txt", f"{name}/fastqc_data.txt")
    return zip_path


def test_parse_fastqc_zip_reads_inner_files(tmp_path):
    report = parse_fastqc_zip(_make_zip(tmp_path))
    assert report.modules["Per base sequence quality"] == "FAIL"
    assert report.basic_stats["Total Sequences"] == "500"


@pytest.mark.skipif(shutil.which("conda") is None, reason="conda yok")
def test_run_fastqc_produces_parsable_zip(tmp_path):
    """Entegrasyon: gerçek FastQC'nin çıktısı parser'ımızla okunabilmeli.
    rnaforge-qc env / fastqc yoksa skip (yanlış alarm vermesin)."""
    fastq = write_fastq(tmp_path / "s.fastq", 500, 150, "I")
    out_dir = tmp_path / "raw_qc"
    out_dir.mkdir()
    try:
        zip_path = run_fastqc(fastq, out_dir)
    except FastQCRunError as exc:
        pytest.skip(f"FastQC çalıştırılamadı (env yok?): {exc}")
    assert zip_path.exists()
    report = parse_fastqc_zip(zip_path)
    assert "Basic Statistics" in report.modules
    assert int(report.basic_stats["Total Sequences"]) == 500
