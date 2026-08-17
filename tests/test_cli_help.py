from __future__ import annotations

from rnaforge.cli import build_parser


def test_build_parser_top_level_help_renders():
    """`rnaforge --help` argparse'ta çökmemeli (literal % escape edilmeli)."""
    parser = build_parser()
    text = parser.format_help()          # bare % olsaydı ValueError fırlatırdı
    assert "seqqc" in text


def test_each_subparser_help_renders():
    parser = build_parser()
    for action in parser._actions:
        for sub in getattr(action, "choices", {} or {}).values() if getattr(action, "choices", None) else []:
            sub.format_help()            # her alt-komut yardımı da render olmalı


def test_seqqc_help_shows_literal_percent():
    """%% doğru şekilde tek % olarak render olmalı (çift kalmamalı)."""
    parser = build_parser()
    text = parser.format_help()
    assert "rRNA%" in text
    assert "rRNA%%" not in text
