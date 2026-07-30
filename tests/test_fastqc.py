from __future__ import annotations

import pytest

from rnaforge.fastqc import FastQCReport, FastQCParseError, parse_fastqc_report

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
