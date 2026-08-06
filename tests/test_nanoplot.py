from __future__ import annotations

import pytest

from rnaforge.nanoplot import NanoStats, NanoStatsParseError, parse_nanostats

_REAL = """Metrics\tdataset
number_of_reads\t89032
number_of_bases\t74816374.0
median_read_length\t724.0
mean_read_length\t840.3
read_length_stdev\t786.1
n50\t1371.0
mean_qual\t6.7
median_qual\t9.3
longest_read_(with_Q):1\t30092 (5.5)
highest_Q_read_(with_length):1\t18.3 (195)
Reads >Q10:\t38965 (43.8%) 50.4Mb
Reads >Q15:\t309 (0.3%) 0.3Mb
"""


def test_parse_nanostats_core_fields():
    s = parse_nanostats(_REAL)
    assert isinstance(s, NanoStats)
    assert s.number_of_reads == 89032
    assert s.number_of_bases == 74816374
    assert s.mean_read_length == pytest.approx(840.3)
    assert s.median_read_length == pytest.approx(724.0)
    assert s.n50 == pytest.approx(1371.0)
    assert s.mean_qual == pytest.approx(6.7)
    assert s.median_qual == pytest.approx(9.3)
    assert s.reads_above_q10_pct == pytest.approx(43.8)


def test_parse_nanostats_missing_core_raises():
    with pytest.raises(NanoStatsParseError):
        parse_nanostats("Metrics\tdataset\nnumber_of_reads\t10\n")
