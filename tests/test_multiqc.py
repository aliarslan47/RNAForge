from __future__ import annotations

from rnaforge.multiqc import count_modules, parse_general_stats

_GS = (
    "Sample\tFastQC_percent_gc\tfeatureCounts_percent_assigned\n"
    "ctrl_rep1\t50\t72.1\n"
    "ctrl_rep2\t51\t70.3\n"
)


def test_parse_general_stats_rows():
    rows = parse_general_stats(_GS)
    assert len(rows) == 2
    assert rows[0]["Sample"] == "ctrl_rep1"
    assert rows[1]["featureCounts_percent_assigned"] == "70.3"


def test_parse_general_stats_empty():
    assert parse_general_stats("") == []
    assert parse_general_stats("   \n") == []


def test_count_modules_ignores_meta(tmp_path):
    d = tmp_path / "multiqc_report_data"
    d.mkdir()
    for n in ("multiqc_general_stats.txt", "multiqc_sources.txt", "multiqc_citations.txt",
              "multiqc_fastqc.txt", "multiqc_featureCounts.txt", "multiqc_samtools_stats.txt"):
        (d / n).write_text("x")
    assert count_modules(d) == 3   # fastqc, featureCounts, samtools_stats


def test_count_modules_missing_dir(tmp_path):
    assert count_modules(tmp_path / "nope") == 0
