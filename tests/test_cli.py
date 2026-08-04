from __future__ import annotations

import pytest

from rnaforge.cli import main


def test_version_flag_prints_version(capsys):
    with pytest.raises(SystemExit) as exc:
        main(["--version"])
    assert exc.value.code == 0
    assert "0.1.0" in capsys.readouterr().out


def test_no_command_returns_usage_error():
    assert main([]) == 2


def test_figures_subcommand_parses(capsys):
    from rnaforge.cli import build_parser
    args = build_parser().parse_args(["figures","--config","c.yaml","--metadata","m.tsv","--run-id","r"])
    assert args.command == "figures"
    assert args.run_id == "r"
