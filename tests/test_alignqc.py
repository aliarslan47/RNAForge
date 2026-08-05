from __future__ import annotations

import pytest

from rnaforge.alignqc import (
    AlignQCParseError,
    aggregate_histograms,
    parse_read_distribution,
    parse_samtools_coverage,
    parse_samtools_stats,
)

_STATS = (
    "SN\traw total sequences:\t11134902\t# excluding ...\n"
    "SN\treads paired:\t11134902\t# paired-end technology bit set\n"
    "SN\tinsert size average:\t271.2\n"
    "SN\tinsert size standard deviation:\t98.2\n"
    "IS\t0\t40531\t31880\t6739\t1912\n"
    "IS\t1\t0\t0\t0\t0\n"
    "IS\t270\t5000\t4000\t900\t100\n"
)

_COVERAGE = (
    "#rname\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq\n"
    "NZ_CP009273.1\t1\t4631469\t11052654\t4381813\t94.6096\t352.861\t38.7\t32.5\n"
    "contig2\t1\t1000\t50\t900\t90.0\t5.5\t38\t30\n"
)

_READ_DIST = (
    "Total Reads                   11052654\n"
    "Total Tags                    19949160\n"
    "Total Assigned Tags           18287406\n"
    "=====================================================================\n"
    "Group               Total_bases         Tag_count           Tags/Kb\n"
    "CDS_Exons           3989128             18286062            4583.97\n"
    "5'UTR_Exons         0                   0                   0.00\n"
    "3'UTR_Exons         0                   0                   0.00\n"
    "Introns             0                   0                   0.00\n"
    "TSS_up_1kb          518460              764                 1.47\n"
    "=====================================================================\n"
)


def test_parse_samtools_stats_insert_size_and_hist():
    d = parse_samtools_stats(_STATS)
    assert d["reads_paired"] == 11134902
    assert d["insert_size_average"] == 271.2
    assert d["insert_size_sd"] == 98.2
    assert (270, 5000) in d["histogram"]


def test_parse_samtools_stats_single_end_zero_insert():
    se = "SN\treads paired:\t0\t#\nSN\tinsert size average:\t0.0\n"
    d = parse_samtools_stats(se)
    assert d["reads_paired"] == 0
    assert d["insert_size_average"] == 0.0


def test_parse_samtools_coverage_per_contig():
    rows = parse_samtools_coverage(_COVERAGE)
    assert rows[0]["contig"] == "NZ_CP009273.1"
    assert rows[0]["meandepth"] == 352.861
    assert rows[1]["contig"] == "contig2"
    assert rows[1]["numreads"] == 50


def test_parse_read_distribution_percentages():
    d = parse_read_distribution(_READ_DIST)
    assert d["total_tags"] == 19949160
    assert d["assigned_tags"] == 18287406
    # CDS büyük çoğunluk
    assert d["percentages"]["CDS_Exons"] > 90.0
    # intergenik = atanmamış / total
    assert d["percentages"]["Intergenic"] == pytest.approx(8.33, abs=0.1)
    assert d["percentages"]["Introns"] == 0.0


def test_parse_read_distribution_no_total_raises():
    with pytest.raises(AlignQCParseError):
        parse_read_distribution("Group Total_bases Tag_count Tags/Kb\n")


def test_aggregate_histograms_sums_and_caps():
    h1 = [(100, 5), (200, 10), (5000, 99)]
    h2 = [(100, 3), (200, 2)]
    labels, totals = aggregate_histograms([h1, h2], max_insert=1000)
    assert labels == ["100", "200"]      # 5000 kesildi
    assert totals == [8, 12]
