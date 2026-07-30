from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from rnaforge.bowtie2 import Bowtie2RunError
from rnaforge.featurecounts import (
    FeatureCountsParseError,
    FeatureCountsResult,
    FeatureCountsRunError,
    parse_counts,
    parse_summary,
    run_featurecounts,
)

_COUNTS = """# Program:featureCounts
Geneid\tChr\tStart\tEnd\tStrand\tLength\ts1.bam\ts2.bam
geneA\tchr1\t101\t1100\t+\t1000\t150\t80
geneB\tchr1\t2101\t3100\t+\t1000\t150\t60
"""

_SUMMARY = """Status\ts1.bam\ts2.bam
Assigned\t300\t140
Unassigned_Unmapped\t0\t0
Unassigned_NoFeatures\t0\t60
"""


def test_parse_counts_reads_genes_and_columns():
    genes, counts = parse_counts(_COUNTS)
    assert genes == ["geneA", "geneB"]
    assert counts["s1.bam"] == [150, 150]
    assert counts["s2.bam"] == [80, 60]
    assert list(counts.keys()) == ["s1.bam", "s2.bam"]   # insertion order = BAM sirasi


def test_parse_counts_rejects_missing_header():
    with pytest.raises(FeatureCountsParseError, match="Geneid"):
        parse_counts("# only a comment\n")


def test_parse_summary_computes_assignment_rate():
    rates = parse_summary(_SUMMARY)
    assert rates["s1.bam"] == pytest.approx(1.0)          # 300/300
    assert rates["s2.bam"] == pytest.approx(140 / 200)    # 140/(140+0+60)


def test_parse_summary_zero_total_is_zero():
    assert parse_summary("Status\tx.bam\nAssigned\t0\n")["x.bam"] == 0.0


def _genome_gtf_bam(tmp_path):
    """Sentetik genom + 2-gen GTF + o bölgelerden okuma → BAM (gerçek bowtie2+samtools)."""
    import random
    from rnaforge.bowtie2 import build_index, run_bowtie2
    random.seed(8)
    genome = "".join(random.choice("ACGT") for _ in range(6000))
    (tmp_path / "genome.fa").write_text(">chr1\n" + genome + "\n")
    gtf = tmp_path / "genes.gtf"
    gtf.write_text(
        'chr1\tsrc\texon\t101\t1100\t.\t+\t.\tgene_id "geneA";\n'
        'chr1\tsrc\texon\t2101\t3100\t.\t+\t.\tgene_id "geneB";\n'
    )
    reads = tmp_path / "reads.fastq"
    with reads.open("w") as f:
        for i in range(300):
            lo, hi = (101, 1100) if i % 2 == 0 else (2101, 3100)
            p = random.randint(lo - 1, hi - 100)
            f.write(f"@r{i}\n{genome[p:p+100]}\n+\n{'I'*100}\n")
    prefix = build_index(tmp_path / "genome.fa", tmp_path / "idx")
    result = run_bowtie2(prefix, tmp_path / "aln", reads)
    return gtf, result.bam


@pytest.mark.skipif(shutil.which("conda") is None, reason="conda yok")
def test_run_featurecounts_counts_genes(tmp_path):
    """Entegrasyon: gerçek featureCounts 2 geni sayar, yüksek atama. Env yoksa skip."""
    try:
        gtf, bam = _genome_gtf_bam(tmp_path)
        result = run_featurecounts([bam], gtf, tmp_path / "fc",
                                   feature_type="exon", attribute="gene_id")
    except (FeatureCountsRunError, Bowtie2RunError) as exc:
        pytest.skip(f"featureCounts/bowtie2 çalıştırılamadı: {exc}")
    assert set(result.gene_ids) == {"geneA", "geneB"}
    col = list(result.counts.keys())[0]
    assert sum(result.counts[col]) > 0
    assert result.assignment_rates[col] > 0.9
