from __future__ import annotations

import pytest

from rnaforge.featurecounts import (
    FeatureCountsParseError, FeatureCountsResult, parse_counts, parse_summary,
)

_COUNTS = """# Program:featureCounts
Geneid\tChr\tStart\tEnd\tStrand\tLength\ts1.bam\ts2.bam
geneA\tchr1\t101\t1100\t+\t1000\t150\t80
geneB\tchr1\t2101\t3100\t+\t1000\t150\t60
"""

_SUMMARY = """Status\ts1.bam\ts2.bam
Assigned\t300\t140
Unassigned_Unmapped\t0\t0
Unassigned_NoFeatures\t0\t60
"""


def test_parse_counts_reads_genes_and_columns():
    genes, counts = parse_counts(_COUNTS)
    assert genes == ["geneA", "geneB"]
    assert counts["s1.bam"] == [150, 150]
    assert counts["s2.bam"] == [80, 60]
    assert list(counts.keys()) == ["s1.bam", "s2.bam"]   # insertion order = BAM sirasi


def test_parse_counts_rejects_missing_header():
    with pytest.raises(FeatureCountsParseError, match="Geneid"):
        parse_counts("# only a comment\n")


def test_parse_summary_computes_assignment_rate():
    rates = parse_summary(_SUMMARY)
    assert rates["s1.bam"] == pytest.approx(1.0)          # 300/300
    assert rates["s2.bam"] == pytest.approx(140 / 200)    # 140/(140+0+60)


def test_parse_summary_zero_total_is_zero():
    assert parse_summary("Status\tx.bam\nAssigned\t0\n")["x.bam"] == 0.0
