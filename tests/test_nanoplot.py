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


import shutil  # noqa: E402

from rnaforge.nanoplot import run_nanoplot  # noqa: E402

_HAS_ENV = shutil.which("conda") is not None


@pytest.mark.skipif(not _HAS_ENV, reason="conda/rnaforge-longread not available")
def test_run_nanoplot_on_tiny_fastq(tmp_path):
    fq = tmp_path / "r.fastq"
    fq.write_text(
        "@r1\n" + "ACGT" * 60 + "\n+\n" + "I" * 240 + "\n"
        "@r2\n" + "ACGT" * 80 + "\n+\n" + "I" * 320 + "\n"
        "@r3\n" + "ACGT" * 50 + "\n+\n" + "I" * 200 + "\n"
    )
    stats_path = run_nanoplot(fq, tmp_path / "out")
    assert stats_path.name == "NanoStats.txt"
    assert stats_path.exists()
    s = parse_nanostats(stats_path.read_text())
    assert s.number_of_reads == 3
