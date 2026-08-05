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


# --- per-base kompozisyon + duplikasyon parser'ları (F1) ---
from rnaforge.fastqc import parse_per_base_content, parse_deduplication  # noqa: E402

_PBC = (
    ">>Per base sequence content\tfail\n"
    "#Base\tG\tA\tT\tC\n"
    "1\t22.5\t28.0\t19.7\t29.7\n"
    "2-3\t27.0\t25.1\t24.1\t23.7\n"
    ">>END_MODULE\n"
    ">>Sequence Duplication Levels\twarn\n"
    "#Total Deduplicated Percentage\t56.966\n"
    "#Duplication Level\tPercentage of deduplicated\tPercentage of total\n"
    ">>END_MODULE\n"
)


def test_parse_per_base_content_reads_header_order():
    rows = parse_per_base_content(_PBC)
    assert rows[0][0] == "1"
    assert rows[0][1] == {"G": 22.5, "A": 28.0, "T": 19.7, "C": 29.7}
    assert rows[1][0] == "2-3"
    assert rows[1][1]["A"] == 25.1


def test_parse_per_base_content_absent_returns_empty():
    assert parse_per_base_content(">>Basic Statistics\tpass\n>>END_MODULE\n") == []


def test_parse_deduplication_reads_value():
    assert parse_deduplication(_PBC) == 56.966


def test_parse_deduplication_absent_returns_none():
    assert parse_deduplication(">>Basic Statistics\tpass\n>>END_MODULE\n") is None


def test_parse_fastqc_report_fills_new_fields():
    from rnaforge.fastqc import parse_fastqc_report
    full_data = DATA + _PBC
    report = parse_fastqc_report(SUMMARY, full_data)
    assert report.deduplication == 56.966
    assert len(report.per_base_content) == 2
