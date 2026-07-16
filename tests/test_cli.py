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
