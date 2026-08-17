from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from rnaforge.basecall import is_signal_input


def test_is_signal_input_pod5_file(tmp_path):
    p = tmp_path / "reads.pod5"; p.write_bytes(b"\x00")
    assert is_signal_input(p) == "pod5"


def test_is_signal_input_fast5_file(tmp_path):
    p = tmp_path / "reads.fast5"; p.write_bytes(b"\x00")
    assert is_signal_input(p) == "fast5"


def test_is_signal_input_directory_of_pod5(tmp_path):
    d = tmp_path / "signal"; d.mkdir()
    (d / "a.pod5").write_bytes(b"\x00")
    (d / "b.pod5").write_bytes(b"\x00")
    assert is_signal_input(d) == "pod5"


def test_is_signal_input_directory_of_fast5(tmp_path):
    d = tmp_path / "signal"; d.mkdir()
    (d / "a.fast5").write_bytes(b"\x00")
    assert is_signal_input(d) == "fast5"


def test_is_signal_input_fastq_is_none(tmp_path):
    p = tmp_path / "reads.fastq"; p.write_text("@r\nACGT\n+\nIIII\n")
    assert is_signal_input(p) is None


def test_is_signal_input_fastq_gz_is_none(tmp_path):
    p = tmp_path / "reads.fastq.gz"; p.write_bytes(b"\x1f\x8b")
    assert is_signal_input(p) is None


def test_is_signal_input_empty_dir_is_none(tmp_path):
    d = tmp_path / "empty"; d.mkdir()
    assert is_signal_input(d) is None


# --- Integration: real dorado GPU basecall (env/tool-gated) ---
_DORADO = Path("/home/ali/tools/dorado-2.1.1-linux-x64/bin/dorado")
_POD5 = Path("/home/ali/tools/pod5_test/r10.pod5")


@pytest.mark.skipif(not (_DORADO.exists() and _POD5.exists() and shutil.which("nvidia-smi")),
                    reason="dorado/test-pod5/GPU not available")
def test_run_dorado_basecalls_pod5_to_fastq(tmp_path):
    from rnaforge.basecall import run_dorado
    out = tmp_path / "out.fastq"
    n = run_dorado(_POD5, out, dorado_bin=_DORADO,
                   models_dir=Path("/home/ali/tools/models"))
    assert out.exists()
    assert n >= 1
    assert out.read_text().startswith("@")
