from __future__ import annotations

import shutil

import pytest

from rnaforge.pychopper import (
    PychopperParseError, PychopperStats, parse_pychopper_stats, run_pychopper,
)

_HAS = shutil.which("conda") is not None

_STATS = (
    "Category\tName\tValue\n"
    "ReadStats\tPassReads\t2000.0\n"
    "ReadStats\tLenFail\t63.0\n"
    "ReadStats\tQcFail\t0.0\n"
    "Classification\tPrimers_found\t1166.0\n"
    "Classification\tRescue\t28.0\n"
    "Classification\tUnusable\t820.0\n"
    "Strand\t+\t419.0\n"
)


def test_parse_pychopper_stats():
    s = parse_pychopper_stats(_STATS)
    assert isinstance(s, PychopperStats)
    assert s.pass_reads == 2000
    assert s.primers_found == 1166
    assert s.rescue == 28
    assert s.unusable == 820
    assert s.len_fail == 63


def test_parse_pychopper_stats_missing_raises():
    with pytest.raises(PychopperParseError):
        parse_pychopper_stats("Category\tName\tValue\nReadStats\tPassReads\t10.0\n")


@pytest.mark.skipif(not _HAS, reason="conda/rnaforge-longread not available")
def test_run_pychopper_tolerates_plot_crash(tmp_path):
    fq = tmp_path / "in.fastq"
    body = "".join(
        f"@r{i}\n" + "ACGT" * 60 + "\n+\n" + "I" * 240 + "\n" for i in range(200)
    )
    fq.write_text(body)
    stats = run_pychopper(fq, tmp_path / "fl.fastq", tmp_path / "stats.tsv")
    assert isinstance(stats, PychopperStats)
    assert stats.pass_reads == 200
    assert (tmp_path / "fl.fastq").exists()
