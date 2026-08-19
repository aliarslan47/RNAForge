from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from rnaforge.minimap2 import (
    AlignmentResult,
    Minimap2ParseError,
    Minimap2RunError,
    minimap2_preset,
    parse_flagstat_mapped,
    run_minimap2,
)

# Real `samtools flagstat` (1.24) output shape: primary=2000, primary mapped=1900.
_FLAGSTAT = """2000 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
1900 + 0 mapped (95.00% : N/A)
2000 + 0 primary
1900 + 0 primary mapped (95.00% : N/A)
0 + 0 paired in sequencing
0 + 0 read1
0 + 0 read2
0 + 0 properly paired (N/A : N/A)
0 + 0 with itself and mate mapped
0 + 0 singletons (N/A : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
"""


def test_parse_flagstat_mapped_fraction():
    assert parse_flagstat_mapped(_FLAGSTAT) == pytest.approx(1900 / 2000)


def test_parse_flagstat_all_mapped():
    text = "500 + 0 primary\n500 + 0 primary mapped (100.00% : N/A)\n"
    assert parse_flagstat_mapped(text) == 1.0


def test_parse_flagstat_missing_raises():
    with pytest.raises(Minimap2ParseError):
        parse_flagstat_mapped("some unrelated samtools chatter\n")


def test_parse_flagstat_zero_primary_raises():
    # primary=0 is not a valid alignment rate (empty input) — fail loudly.
    with pytest.raises(Minimap2ParseError):
        parse_flagstat_mapped("0 + 0 primary\n0 + 0 primary mapped (N/A : N/A)\n")


def test_preset_ont():
    assert minimap2_preset("ont") == "map-ont"


def test_preset_pacbio_hifi():
    assert minimap2_preset("pacbio_hifi") == "map-hifi"


def test_preset_unknown_platform_raises():
    with pytest.raises(Minimap2RunError):
        minimap2_preset("illumina")


def _synthetic_genome_and_reads(tmp_path, aligned=True):
    import random
    random.seed(11)
    genome = "".join(random.choice("ACGT") for _ in range(20000))
    (tmp_path / "genome.fa").write_text(">chr1\n" + genome + "\n")
    reads = tmp_path / "reads.fastq"
    with reads.open("w") as f:
        for i in range(200):
            if aligned:
                p = random.randint(0, len(genome) - 600)
                seq = genome[p:p + 600]        # long ONT-like reads
            else:
                seq = "".join(random.choice("ACGT") for _ in range(600))
            f.write(f"@r{i}\n{seq}\n+\n{'I' * len(seq)}\n")
    return tmp_path / "genome.fa", reads


@pytest.mark.skipif(shutil.which("conda") is None, reason="conda yok")
def test_run_minimap2_aligns_derived_reads(tmp_path):
    """Entegrasyon: genomdan türetilmiş uzun okumalar yüksek hizalanır.
    rnaforge-longread yoksa skip."""
    genome, reads = _synthetic_genome_and_reads(tmp_path, aligned=True)
    try:
        result = run_minimap2(genome, tmp_path / "aln", reads, preset="map-ont")
    except Minimap2RunError as exc:
        pytest.skip(f"minimap2 çalıştırılamadı (env yok?): {exc}")
    assert isinstance(result, AlignmentResult)
    assert result.bam.exists()
    assert Path(str(result.bam) + ".bai").exists()
    assert result.alignment_rate > 0.90


def _has_lr():
    import subprocess
    return subprocess.run(["conda", "run", "-n", "rnaforge-longread", "samtools", "--version"],
                          capture_output=True, text=True).returncode == 0


@__import__("pytest").mark.skipif(not _has_lr(), reason="rnaforge-longread yok")
def test_count_primary_alignments(tmp_path):
    import subprocess
    from rnaforge.minimap2 import count_primary_alignments
    sam = tmp_path / "a.sam"
    sam.write_text(
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:tx1\tLN:100\n@SQ\tSN:tx2\tLN:100\n"
        "r1\t0\ttx1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n"
        "r2\t0\ttx1\t5\t60\t4M\t*\t0\t0\tACGT\tIIII\n"
        "r3\t0\ttx2\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n"
        "r3\t256\ttx1\t1\t0\t4M\t*\t0\t0\tACGT\tIIII\n"
        "r4\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n")
    bam = tmp_path / "a.bam"
    subprocess.run(["conda", "run", "-n", "rnaforge-longread", "bash", "-c",
                    f"samtools sort -o {bam} {sam} && samtools index {bam}"], check=True)
    assert count_primary_alignments(bam) == {"tx1": 2, "tx2": 1}


def test_run_minimap2_secondary_n_adds_flag(tmp_path, monkeypatch):
    """İzoform EM (NanoCount) için minimap2 ikincil hizalamaları saklamalı: secondary_n
    verilince komuta '-N 10' eklenir; verilmeyince eklenmez (gen-yolu regresyon-güvenli)."""
    import rnaforge.minimap2 as mm
    genome = tmp_path / "tx.fa"; genome.write_text(">t1\nACGT\n")
    reads = tmp_path / "r.fq"; reads.write_text("@r\nACGT\n+\nIIII\n")
    captured = {}

    def fake_run(cmd):
        if "minimap2" in cmd:
            captured["mm"] = cmd
            Path(cmd[cmd.index("-o") + 1]).write_text("")          # sam
        elif "sort" in cmd:
            Path(cmd[cmd.index("-o") + 1]).write_bytes(b"BAM")     # bam
        class R:
            returncode = 0; stderr = ""; stdout = _FLAGSTAT
        return R()

    monkeypatch.setattr(mm, "_run", fake_run)

    run_minimap2(genome, tmp_path / "a", reads, preset="map-ont", secondary_n=10)
    assert "-N" in captured["mm"] and captured["mm"][captured["mm"].index("-N") + 1] == "10"

    run_minimap2(genome, tmp_path / "b", reads, preset="map-ont")
    assert "-N" not in captured["mm"]                              # default: ikincil yok
