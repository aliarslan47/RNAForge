from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from rnaforge.bowtie2 import (
    AlignmentResult,
    Bowtie2ParseError,
    Bowtie2RunError,
    build_index,
    parse_bowtie2_summary,
    run_bowtie2,
)

_STDERR = """300 reads; of these:
  300 (100.00%) were unpaired; of these:
    12 (4.00%) aligned 0 times
    288 (96.00%) aligned exactly 1 time
96.00% overall alignment rate
"""


def test_parse_reads_overall_rate():
    assert parse_bowtie2_summary(_STDERR) == pytest.approx(0.96)


def test_parse_zero_rate():
    assert parse_bowtie2_summary("0.00% overall alignment rate\n") == 0.0


def test_parse_rejects_missing_summary():
    with pytest.raises(Bowtie2ParseError, match="overall alignment rate"):
        parse_bowtie2_summary("some unrelated bowtie2 chatter\n")


def _synthetic_genome_and_reads(tmp_path, aligned=True):
    import random
    random.seed(11)
    genome = "".join(random.choice("ACGT") for _ in range(5000))
    (tmp_path / "genome.fa").write_text(">chr1\n" + genome + "\n")
    reads = tmp_path / "reads.fastq"
    with reads.open("w") as f:
        for i in range(200):
            if aligned:
                p = random.randint(0, len(genome) - 100)
                seq = genome[p:p + 100]
            else:
                seq = "".join(random.choice("ACGT") for _ in range(100))
            f.write(f"@r{i}\n{seq}\n+\n{'I' * 100}\n")
    return tmp_path / "genome.fa", reads


@pytest.mark.skipif(shutil.which("conda") is None, reason="conda yok")
def test_build_index_and_align_reports_high_rate(tmp_path):
    """Entegrasyon: genomdan türetilmiş okumalar yüksek hizalanır. rnaforge-quant-prok
    yoksa skip."""
    genome, reads = _synthetic_genome_and_reads(tmp_path, aligned=True)
    try:
        prefix = build_index(genome, tmp_path / "idx")
        result = run_bowtie2(prefix, tmp_path / "aln", reads)
    except Bowtie2RunError as exc:
        pytest.skip(f"bowtie2 çalıştırılamadı (env yok?): {exc}")
    assert result.bam.exists()
    assert Path(str(result.bam) + ".bai").exists()   # samtools index -> aligned.sorted.bam.bai
    assert result.alignment_rate > 0.95


@pytest.mark.skipif(shutil.which("conda") is None, reason="conda yok")
def test_random_reads_align_poorly(tmp_path):
    genome, reads = _synthetic_genome_and_reads(tmp_path, aligned=False)
    try:
        prefix = build_index(genome, tmp_path / "idx")
        result = run_bowtie2(prefix, tmp_path / "aln", reads)
    except Bowtie2RunError as exc:
        pytest.skip(f"bowtie2 çalıştırılamadı: {exc}")
    assert result.alignment_rate < 0.70   # genom-dışı okumalar eşiğin altında
