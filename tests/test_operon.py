"""m14 operon testleri: tahmin (mesafe/strand/contig) + DE koordinasyon agregasyonu."""
from __future__ import annotations

from rnaforge.operon import aggregate_operon_de, predict_operons

GFF = (
    "c1\tx\tgene\t100\t200\t.\t+\t.\tlocus_tag=A;gene=a\n"
    "c1\tx\tgene\t210\t300\t.\t+\t.\tlocus_tag=B;gene=b\n"    # gap 9 -> A ile aynı operon
    "c1\tx\tgene\t500\t600\t.\t+\t.\tlocus_tag=C;gene=c\n"    # gap 199 -> yeni operon
    "c1\tx\tgene\t610\t700\t.\t-\t.\tlocus_tag=D;gene=d\n"    # ters strand -> yeni operon
    "c2\tx\tgene\t120\t210\t.\t+\t.\tlocus_tag=E;gene=e\n"    # farklı contig -> yeni operon
)


def _gff(tmp_path):
    p = tmp_path / "g.gff"; p.write_text(GFF); return p


def test_predict_operons_groups_and_breaks(tmp_path):
    ops = predict_operons(_gff(tmp_path), max_gap=50)
    assert [o["locus_tags"] for o in ops] == [["A", "B"], ["C"], ["D"], ["E"]]
    assert ops[0]["strand"] == "+" and ops[0]["size"] == 2
    assert ops[0]["symbols"] == ["a", "b"]


def test_predict_operons_gap_threshold(tmp_path):
    # max_gap 5 -> A(100-200) ve B(210-300) arası gap 9 > 5 -> ayrı operon
    ops = predict_operons(_gff(tmp_path), max_gap=5)
    assert ops[0]["locus_tags"] == ["A"] and ops[1]["locus_tags"] == ["B"]


DESEQ = (
    "gene\tbaseMean\tlog2FoldChange\tlfcSE\tstat\tpvalue\tpadj\n"
    "A\t100\t2.5\t0.2\t8\t1e-9\t1e-8\n"       # up
    "B\t100\t2.2\t0.2\t7\t1e-8\t1e-7\n"       # up
    "C\t100\t-3.0\t0.2\t-9\t1e-9\t1e-8\n"     # down
    "D\t100\t0.1\t0.2\t0.3\t0.7\t0.8\n"       # ns
)


def test_aggregate_coordinated_operon(tmp_path):
    tsv = tmp_path / "de.tsv"; tsv.write_text(DESEQ)
    ops = predict_operons(_gff(tmp_path), max_gap=50)
    agg = aggregate_operon_de(ops, tsv, fdr=0.05, lfc=1.0)
    ab = next(o for o in agg if o["locus_tags"] == ["A", "B"])
    assert ab["coordinated"] is True and ab["n_up"] == 2 and ab["n_deg"] == 2
    assert agg[0]["coordinated"] is True                 # koordineli önce sıralı


def test_aggregate_mixed_not_coordinated(tmp_path):
    # A up, B down aynı operonda olsaydı coordinated=False; burada kurgu:
    deseq = DESEQ.replace("B\t100\t2.2\t0.2\t7\t1e-8\t1e-7\n", "B\t100\t-2.2\t0.2\t-7\t1e-8\t1e-7\n")
    tsv = tmp_path / "de.tsv"; tsv.write_text(deseq)
    ops = predict_operons(_gff(tmp_path), max_gap=50)
    agg = aggregate_operon_de(ops, tsv, fdr=0.05, lfc=1.0)
    ab = next(o for o in agg if o["locus_tags"] == ["A", "B"])
    assert ab["coordinated"] is False and ab["n_up"] == 1 and ab["n_down"] == 1


def test_aggregate_single_gene_not_coordinated(tmp_path):
    tsv = tmp_path / "de.tsv"; tsv.write_text(DESEQ)
    ops = predict_operons(_gff(tmp_path), max_gap=50)
    agg = aggregate_operon_de(ops, tsv, fdr=0.05, lfc=1.0)
    c = next(o for o in agg if o["locus_tags"] == ["C"])
    assert c["coordinated"] is False and c["n_down"] == 1   # tek gen -> koordineli değil
