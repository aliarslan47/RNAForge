"""Sentetik FASTQ fixture'ları. Gerçek/müşteri verisi ASLA kullanılmaz (PLAN Kural 8)."""
from __future__ import annotations

import gzip
import random
from pathlib import Path

import pytest


def _record(name: str, seq_len: int, qual_char: str) -> str:
    seq = "".join(random.choice("ACGT") for _ in range(seq_len))
    return f"@{name}\n{seq}\n+\n{qual_char * seq_len}\n"


def write_fastq(path: Path, n_reads: int, seq_len, qual_char: str, gzipped: bool = False) -> Path:
    """seq_len: int (sabit) veya (min, max) tuple (değişken uzunluk)."""
    def length() -> int:
        return seq_len if isinstance(seq_len, int) else random.randint(*seq_len)

    body = "".join(_record(f"read{i}", length(), qual_char) for i in range(n_reads))
    if gzipped:
        with gzip.open(path, "wt") as fh:
            fh.write(body)
    else:
        path.write_text(body)
    return path


@pytest.fixture(autouse=True)
def _seed():
    random.seed(1337)


@pytest.fixture
def illumina_fastq(tmp_path) -> Path:
    # 150 bp sabit, Q40 ('I')
    return write_fastq(tmp_path / "illumina.fastq", 200, 150, "I")


@pytest.fixture
def ont_fastq(tmp_path) -> Path:
    # uzun ve gürültülü: 1-20 kb, Q10 ('+')
    return write_fastq(tmp_path / "ont.fastq", 200, (1000, 20000), "+")


@pytest.fixture
def pacbio_fastq(tmp_path) -> Path:
    # uzun ve yüksek kaliteli: 8-20 kb, Q40 ('I')
    return write_fastq(tmp_path / "pacbio.fastq", 200, (8000, 20000), "I")
