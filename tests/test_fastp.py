from __future__ import annotations

import json

import pytest

from rnaforge.fastp import FastpParseError, FastpResult, parse_fastp_json

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
