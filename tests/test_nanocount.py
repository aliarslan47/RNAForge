from __future__ import annotations

import rnaforge.nanocount as nanocount
from rnaforge.nanocount import parse_nanocount, run_nanocount


# Gerçek NanoCount v1.x çıktı başlığı: transcript_name, raw, est_count, tpm (TAB).
_NC = (
    "transcript_name\traw\test_count\ttpm\n"
    "ENST001\t0.90\t120.4\t5000.0\n"
    "ENST002\t0.10\t13.6\t800.0\n"
)


def test_parse_nanocount_transcript_to_est_count():
    assert parse_nanocount(_NC) == {"ENST001": 120.4, "ENST002": 13.6}


def test_parse_nanocount_empty_or_header_only():
    assert parse_nanocount("transcript_name\traw\test_count\ttpm\n") == {}
    assert parse_nanocount("") == {}


def test_run_nanocount_invokes_tool_and_returns_parsed(tmp_path, monkeypatch):
    """run_nanocount NanoCount'u çağırır, çıktı TSV'sini parse edip döndürür."""
    bam = tmp_path / "aligned.sorted.bam"; bam.write_bytes(b"BAM")
    out = tmp_path / "nanocount.tsv"
    captured = {}

    def fake_run(cmd, capture_output, text):
        captured["cmd"] = cmd
        out.write_text(_NC)
        class R:
            returncode = 0; stderr = ""; stdout = ""
        return R()

    monkeypatch.setattr(nanocount.subprocess, "run", fake_run)
    got = run_nanocount(bam, out)
    assert got == {"ENST001": 120.4, "ENST002": 13.6}
    assert "NanoCount" in captured["cmd"]
    assert str(bam) in captured["cmd"] and str(out) in captured["cmd"]


def test_run_nanocount_nonzero_exit_raises(tmp_path, monkeypatch):
    bam = tmp_path / "a.bam"; bam.write_bytes(b"BAM")

    def fake_run(cmd, capture_output, text):
        class R:
            returncode = 1; stderr = "boom"; stdout = ""
        return R()

    monkeypatch.setattr(nanocount.subprocess, "run", fake_run)
    import pytest
    with pytest.raises(nanocount.NanoCountRunError):
        run_nanocount(bam, tmp_path / "o.tsv")
