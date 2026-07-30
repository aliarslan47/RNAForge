from __future__ import annotations

import json
import shutil

import pytest

from rnaforge.fastp import (
    FastpParseError,
    FastpResult,
    FastpRunError,
    parse_fastp_json,
    run_fastp,
)
from tests.conftest import write_fastq

_GOOD = json.dumps({
    "summary": {
        "before_filtering": {"total_reads": 1000},
        "after_filtering": {"total_reads": 900},
    }
})


def test_parse_reads_before_after_and_survival():
    r = parse_fastp_json(_GOOD)
    assert r.reads_before == 1000
    assert r.reads_after == 900
    assert r.survival_rate == pytest.approx(0.9)


def test_parse_survival_zero_when_no_reads_before():
    text = json.dumps({"summary": {
        "before_filtering": {"total_reads": 0},
        "after_filtering": {"total_reads": 0}}})
    assert parse_fastp_json(text).survival_rate == 0.0


def test_parse_rejects_missing_summary():
    with pytest.raises(FastpParseError, match="summary"):
        parse_fastp_json(json.dumps({"filtering_result": {}}))


@pytest.mark.skipif(shutil.which("conda") is None, reason="conda yok")
def test_run_fastp_trims_and_reports(tmp_path):
    """Entegrasyon: gerçek fastp nazik trimming; kısa-olmayan okumalar korunur
    (survival ~1.0). rnaforge-qc / fastp yoksa skip."""
    fastq = write_fastq(tmp_path / "s.fastq", 500, 150, "I")
    out_dir = tmp_path / "trimmed"
    try:
        result = run_fastp(fastq, out_dir, min_length=36)
    except FastpRunError as exc:
        pytest.skip(f"fastp çalıştırılamadı (env yok?): {exc}")
    assert result.out1.exists()
    assert (out_dir / "fastp.json").exists()
    assert (out_dir / "fastp.html").exists()
    assert result.reads_before == 500
    assert result.survival_rate > 0.95   # nazik: 150 bp okumalar 36 filtresini geçer
