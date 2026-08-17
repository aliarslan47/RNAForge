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


def _fake_fc_subprocess(monkeypatch, captured):
    """subprocess.run'ı yakalar; parse'ın geçmesi için minimal counts.txt+.summary yazar."""
    import rnaforge.featurecounts as fc

    def fake_run(cmd, capture_output=True, text=True):
        captured["cmd"] = cmd
        out = Path(cmd[cmd.index("-o") + 1])
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(
            "# prog\nGeneid\tChr\tStart\tEnd\tStrand\tLength\ta.bam\n"
            "g1\tc1\t1\t80\t+\t80\t5\n"
        )
        out.with_name(out.name + ".summary").write_text(
            "Status\ta.bam\nAssigned\t5\nUnassigned_NoFeatures\t1\n"
        )

        class _R:
            returncode = 0
            stderr = ""
            stdout = ""
        return _R()

    monkeypatch.setattr(fc.subprocess, "run", fake_run)


def test_run_featurecounts_long_read_adds_L(tmp_path, monkeypatch):
    from rnaforge.featurecounts import run_featurecounts
    captured: dict = {}
    _fake_fc_subprocess(monkeypatch, captured)
    run_featurecounts([Path("a.bam")], Path("g.gtf"), tmp_path / "fc",
                      feature_type="exon", attribute="gene_id",
                      paired=True, long_read=True)   # long → -L, paired ignored
    assert "-L" in captured["cmd"]
    assert "-p" not in captured["cmd"]


def test_run_featurecounts_short_has_no_L_but_paired(tmp_path, monkeypatch):
    from rnaforge.featurecounts import run_featurecounts
    captured: dict = {}
    _fake_fc_subprocess(monkeypatch, captured)
    run_featurecounts([Path("a.bam")], Path("g.gtf"), tmp_path / "fc",
                      feature_type="exon", attribute="gene_id",
                      paired=True, long_read=False)
    assert "-L" not in captured["cmd"]
    assert "-p" in captured["cmd"]


from rnaforge.featurecounts import compute_tpm_fpkm, parse_lengths

_FC = (
    "# Program:featureCounts\n"
    "Geneid\tChr\tStart\tEnd\tStrand\tLength\ts1.bam\ts2.bam\n"
    "gA\tc1\t1\t1000\t+\t1000\t100\t0\n"
    "gB\tc1\t2000\t2999\t+\t1000\t300\t50\n"
    "gC\tc1\t4000\t4999\t+\t2000\t200\t50\n"
)


def test_parse_lengths():
    assert parse_lengths(_FC) == {"gA": 1000, "gB": 1000, "gC": 2000}


def test_compute_tpm_fpkm_shapes_and_values():
    genes, cols, tpm, fpkm = compute_tpm_fpkm(_FC)
    assert genes == ["gA", "gB", "gC"] and cols == ["s1.bam", "s2.bam"]
    # s1: RPK = 100/1, 300/1, 200/2 = 100,300,100 -> toplam 500; TPM = /500*1e6
    assert tpm["s1.bam"][0] == 200000.0 and tpm["s1.bam"][1] == 600000.0
    # TPM sütun toplamı ~1e6
    assert abs(sum(tpm["s1.bam"]) - 1_000_000) < 1
    # FPKM: gA s1 = 100/(1 * (600/1e6)) çünkü s1 toplam=600
    assert abs(fpkm["s1.bam"][0] - 100 / (1.0 * (600/1e6))) < 1


def test_compute_tpm_empty_column():
    fc = _FC.replace("100\t0", "0\t0").replace("300\t50", "0\t0").replace("200\t50", "0\t0")
    _, _, tpm, _ = compute_tpm_fpkm(fc)
    assert all(v == 0.0 for v in tpm["s1.bam"])   # sıfır sayım -> 0, çökme yok
