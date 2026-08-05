"""m13 abricate testleri: parser, koordinat eşleme, DE overlay."""
from __future__ import annotations

from rnaforge.abricate import (
    gene_coords, map_hits_to_genes, overlay_de, parse_abricate,
)

HEADER = ("#FILE\tSEQUENCE\tSTART\tEND\tSTRAND\tGENE\tCOVERAGE\tCOVERAGE_MAP\tGAPS\t"
          "%COVERAGE\t%IDENTITY\tDATABASE\tACCESSION\tPRODUCT\tRESISTANCE\n")


def test_parse_abricate_filters_and_types(tmp_path):
    tsv = tmp_path / "a.tsv"
    tsv.write_text(
        HEADER +
        "g\tchr1\t100\t1200\t+\tacrB\t1-1100/1100\t===\t0/0\t100.0\t99.5\tcard\tACC1\tefflux pump\tMULTIDRUG\n"
        "g\tchr1\t5000\t5300\t+\tlowid\t1-300/300\t===\t0/0\t100.0\t60.0\tcard\tACC2\tx\t\n")
    hits = parse_abricate(tsv, 80, 80)
    assert len(hits) == 1                       # düşük %id (60) elendi
    h = hits[0]
    assert h["gene"] == "acrB" and h["contig"] == "chr1" and h["start"] == 100
    assert h["pct_id"] == 99.5 and h["db"] == "card" and h["resistance"] == "MULTIDRUG"


GFF = (
    "chr1\tRefSeq\tgene\t90\t1300\t.\t+\t.\tID=gene-LT1;gene=acrB;locus_tag=LT_1\n"
    "chr1\tRefSeq\tgene\t2000\t2500\t.\t+\t.\tID=gene-LT2;gene=ompF;locus_tag=LT_2\n"
    "chr1\tProtein\tCDS\t90\t1300\t.\t+\t0\tlocus_tag=LT_1\n"   # CDS -> yok sayılır
)


def test_gene_coords_from_gff(tmp_path):
    gff = tmp_path / "g.gff"; gff.write_text(GFF)
    genes = gene_coords(gff)
    assert len(genes) == 2                       # yalnız gene feature
    assert genes[0]["locus_tag"] == "LT_1" and genes[0]["symbol"] == "acrB"
    assert genes[0]["start"] == 90 and genes[0]["end"] == 1300


def test_map_hits_best_overlap():
    genes = [{"contig": "chr1", "start": 90, "end": 1300, "locus_tag": "LT_1", "symbol": "acrB"},
             {"contig": "chr1", "start": 2000, "end": 2500, "locus_tag": "LT_2", "symbol": "ompF"}]
    hits = [{"contig": "chr1", "start": 100, "end": 1200, "gene": "acrB", "pct_id": 99.0, "pct_cov": 100,
             "db": "card", "product": "", "resistance": ""},
            {"contig": "chr9", "start": 1, "end": 50, "gene": "orphan", "pct_id": 95.0, "pct_cov": 100,
             "db": "card", "product": "", "resistance": ""}]        # farklı contig -> eşleşmez
    mapped, n_unmapped = map_hits_to_genes(hits, genes)
    assert n_unmapped == 1
    assert len(mapped) == 1 and mapped[0]["locus_tag"] == "LT_1"


def test_map_dedup_keeps_best_identity():
    genes = [{"contig": "chr1", "start": 90, "end": 1300, "locus_tag": "LT_1", "symbol": "acrB"}]
    hits = [{"contig": "chr1", "start": 100, "end": 1200, "gene": "a", "pct_id": 90.0, "pct_cov": 100,
             "db": "card", "product": "", "resistance": ""},
            {"contig": "chr1", "start": 95, "end": 1250, "gene": "a", "pct_id": 99.0, "pct_cov": 100,
             "db": "card", "product": "", "resistance": ""}]
    mapped, _ = map_hits_to_genes(hits, genes)
    assert len(mapped) == 1 and mapped[0]["pct_id"] == 99.0   # en yüksek %id


DESEQ = (
    "gene\tbaseMean\tlog2FoldChange\tlfcSE\tstat\tpvalue\tpadj\n"
    "LT_1\t100\t2.5\t0.2\t8\t1e-9\t1e-8\n"       # up
    "LT_2\t100\t-3.0\t0.2\t-9\t1e-9\t1e-8\n"      # down
    "LT_3\t100\t0.1\t0.2\t0.3\t0.7\t0.8\n"        # ns
    "LT_4\t0\tNA\tNA\tNA\tNA\tNA\n"               # untested
)


def test_overlay_de_status(tmp_path):
    tsv = tmp_path / "de.tsv"; tsv.write_text(DESEQ)
    mapped = [{"locus_tag": lt, "symbol": lt, "pct_id": 99, "pct_cov": 100, "db": "card",
               "gene": lt, "product": "", "resistance": "", "contig": "c", "start": 1, "end": 2}
              for lt in ("LT_1", "LT_2", "LT_3", "LT_4", "LT_9")]
    out = overlay_de(mapped, tsv, fdr=0.05, lfc=1.0)
    st = {r["locus_tag"]: r["de_status"] for r in out}
    assert st["LT_1"] == "up" and st["LT_2"] == "down"
    assert st["LT_3"] == "ns" and st["LT_4"] == "untested"
    assert st["LT_9"] == "untested"              # matriste yok
    assert out[0]["de_status"] in ("up", "down") # up/down önce sıralı
