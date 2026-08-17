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
