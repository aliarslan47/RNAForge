from __future__ import annotations

import json

import pytest

from rnaforge.routing import require_short_read, resolve_read_type


def _write_stats(run_dir, read_type):
    stats = run_dir / "statistics"
    stats.mkdir(parents=True)
    (stats / "raw_statistics.json").write_text(json.dumps({"read_type": read_type}))


def test_resolve_read_type_reads_short(tmp_path):
    _write_stats(tmp_path, "short")
    assert resolve_read_type(tmp_path) == "short"


def test_resolve_read_type_missing_file_raises(tmp_path):
    with pytest.raises(ValueError):
        resolve_read_type(tmp_path)


def test_require_short_read_passes_for_short(tmp_path):
    _write_stats(tmp_path, "short")
    require_short_read(tmp_path, "qc")  # must NOT raise


def test_require_short_read_blocks_long(tmp_path):
    _write_stats(tmp_path, "long")
    with pytest.raises(NotImplementedError) as exc:
        require_short_read(tmp_path, "qc")
    assert "qc" in str(exc.value)
    assert "long" in str(exc.value).lower()
