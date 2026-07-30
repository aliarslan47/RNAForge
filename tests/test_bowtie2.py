from __future__ import annotations

import pytest

from rnaforge.bowtie2 import AlignmentResult, Bowtie2ParseError, parse_bowtie2_summary

_STDERR = """300 reads; of these:
  300 (100.00%) were unpaired; of these:
    12 (4.00%) aligned 0 times
    288 (96.00%) aligned exactly 1 time
96.00% overall alignment rate
"""


def test_parse_reads_overall_rate():
    assert parse_bowtie2_summary(_STDERR) == pytest.approx(0.96)


def test_parse_zero_rate():
    assert parse_bowtie2_summary("0.00% overall alignment rate\n") == 0.0


def test_parse_rejects_missing_summary():
    with pytest.raises(Bowtie2ParseError, match="overall alignment rate"):
        parse_bowtie2_summary("some unrelated bowtie2 chatter\n")
