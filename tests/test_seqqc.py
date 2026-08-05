"""m16 seqqc helper testleri: rRNA fasta çıkarımı, subsample, sortmerna/infer parser, BED."""
from __future__ import annotations

import gzip

from rnaforge.seqqc import (
    gff_to_bed, parse_infer_experiment, parse_sortmerna_log, rrna_fasta_from_reference,
    subsample_fastq,
)

GENOME = ">chr1\n" + "ACGTACGTAC" * 5 + "\n"     # 50 bp
GFF = (
    "chr1\tRefSeq\trRNA\t1\t10\t.\t+\t.\tlocus_tag=r16S;product=16S ribosomal RNA\n"
    "chr1\tRefSeq\trRNA\t21\t30\t.\t-\t.\tlocus_tag=r23S;product=23S ribosomal RNA\n"
    "chr1\tRefSeq\tgene\t1\t10\t.\t+\t.\tlocus_tag=LT_1;gene=a\n"
)


def test_rrna_fasta_from_reference(tmp_path):
    (tmp_path / "g.fa").write_text(GENOME)
    (tmp_path / "g.gff").write_text(GFF)
    n = rrna_fasta_from_reference(tmp_path / "g.fa", tmp_path / "g.gff", tmp_path / "rrna.fa")
    assert n == 2                                    # iki rRNA
    txt = (tmp_path / "rrna.fa").read_text()
    assert ">r16S" in txt and ">r23S" in txt
    assert "ACGTACGTAC" in txt                        # + strand ilk 10 bp
    # − strand: 21-30 = "ACGTACGTAC" -> ters-tümler
    assert "GTACGTACGT" in txt


def test_subsample_fastq(tmp_path):
    src = tmp_path / "r.fastq.gz"
    with gzip.open(src, "wt") as f:
        for i in range(10):
            f.write(f"@r{i}\nACGT\n+\nIIII\n")
    n = subsample_fastq(src, 3, tmp_path / "sub.fastq")
    assert n == 3
    assert (tmp_path / "sub.fastq").read_text().count("@r") == 3


def test_parse_sortmerna_log(tmp_path):
    log = tmp_path / "aligned.log"
    log.write_text("Total reads = 100000\n"
                   "Total reads passing E-value threshold = 8500 (8.50%)\n")
    assert abs(parse_sortmerna_log(log) - 0.085) < 1e-6


def test_parse_sortmerna_log_missing(tmp_path):
    assert parse_sortmerna_log(tmp_path / "nope.log") == 0.0


def test_gff_to_bed(tmp_path):
    (tmp_path / "g.gff").write_text(GFF)
    n = gff_to_bed(tmp_path / "g.gff", tmp_path / "g.bed")
    assert n == 1                                    # yalnız gene feature
    cols = (tmp_path / "g.bed").read_text().splitlines()[0].split("\t")
    assert cols[0] == "chr1" and cols[1] == "0" and cols[2] == "10"   # 0-tabanlı start
    assert cols[3] == "LT_1" and cols[5] == "+" and cols[9] == "1"    # tek blok


INFER_FWD = ('This is PairEnd Data\n'
             'Fraction of reads failed to determine: 0.02\n'
             'Fraction of reads explained by "1++,1--,2+-,2-+": 0.93\n'
             'Fraction of reads explained by "1+-,1-+,2++,2--": 0.05\n')
INFER_UNSTR = ('Fraction of reads failed to determine: 0.05\n'
               'Fraction of reads explained by "1++,1--,2+-,2-+": 0.47\n'
               'Fraction of reads explained by "1+-,1-+,2++,2--": 0.48\n')
INFER_REV = ('Fraction of reads explained by "1++,1--,2+-,2-+": 0.04\n'
             'Fraction of reads explained by "1+-,1-+,2++,2--": 0.92\n')


def test_parse_infer_experiment_forward():
    s, fwd, rev = parse_infer_experiment(INFER_FWD)
    assert s == "stranded" and fwd == 0.93


def test_parse_infer_experiment_unstranded():
    s, _, _ = parse_infer_experiment(INFER_UNSTR)
    assert s == "unstranded"


def test_parse_infer_experiment_reverse():
    s, _, rev = parse_infer_experiment(INFER_REV)
    assert s == "reverse" and rev == 0.92
