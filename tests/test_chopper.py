from __future__ import annotations

import shutil

import pytest

from rnaforge.chopper import ChopperRunError, parse_kept, run_chopper

_HAS = shutil.which("conda") is not None


def test_parse_kept_reads_number():
    assert parse_kept("Kept 1097 reads out of 1109 reads\n") == 1097


def test_parse_kept_absent_is_none():
    assert parse_kept("no summary here") is None


@pytest.mark.skipif(not _HAS, reason="conda/rnaforge-longread not available")
def test_run_chopper_filters(tmp_path):
    fq = tmp_path / "in.fastq"
    fq.write_text(
        "@a\n" + "ACGT" * 40 + "\n+\n" + "I" * 160 + "\n"
        "@b\n" + "ACGT" * 40 + "\n+\n" + "I" * 160 + "\n"
        "@c\nACGT\n+\nIIII\n"
    )
    out = tmp_path / "out.fastq"
    kept = run_chopper(fq, out, min_qual=7, min_len=50)
    assert kept == 2
    assert out.exists()
