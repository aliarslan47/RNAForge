"""m10 KEGG annotation testleri: REST parse, sembol join, belirsizlik atma, global-map filtresi."""
from __future__ import annotations

from rnaforge.kegg_annotation import build_gene2pathway, parse_kegg

LINKS = (
    "eco:b0002\tpath:eco00260\n"
    "eco:b0003\tpath:eco00260\n"
    "eco:b0002\tpath:eco01100\n"   # global map -> filtrelenir
    "eco:b0114\tpath:eco00010\n"
)
NAMES = (
    "eco00260\tGlycine, serine and threonine metabolism - Escherichia coli K-12 MG1655\n"
    "eco00010\tGlycolysis / Gluconeogenesis - Escherichia coli K-12 MG1655\n"
    "eco01100\tMetabolic pathways - Escherichia coli K-12 MG1655\n"
)
GENELIST = (
    "eco:b0002\tCDS\t337..2799\tthrA; fused aspartate kinase/homoserine dehydrogenase 1\n"
    "eco:b0003\tCDS\t2801..3733\tthrB; homoserine kinase\n"
    "eco:b0114\tCDS\t100..200\taceE; pyruvate dehydrogenase\n"
)


def _write3(tmp_path):
    (tmp_path / "links.tsv").write_text(LINKS)
    (tmp_path / "names.tsv").write_text(NAMES)
    (tmp_path / "genes.tsv").write_text(GENELIST)
    return tmp_path / "links.tsv", tmp_path / "names.tsv", tmp_path / "genes.tsv"


def test_parse_kegg_symbol_pathways_and_global_filter(tmp_path):
    links, names, genes = _write3(tmp_path)
    s2p, meta = parse_kegg(links, names, genes)
    assert s2p["thrA"] == {"eco00260"}          # global map eco01100 hariç
    assert s2p["aceE"] == {"eco00010"}
    assert meta["eco00260"] == ("KEGG", "Glycine, serine and threonine metabolism")  # org eki kırpıldı
    assert "eco01100" not in meta               # global map meta'ya girmez


def _gff(tmp_path, extra=""):
    body = (
        "NZ\tX\tCDS\t1\t9\t.\t+\t0\tgene=thrA;locus_tag=LT_A\n"
        "NZ\tX\tCDS\t1\t9\t.\t+\t0\tgene=thrB;locus_tag=LT_B\n"
        "NZ\tX\tCDS\t1\t9\t.\t+\t0\tgene=aceE;locus_tag=LT_C\n"
    ) + extra
    p = tmp_path / "g.gff"; p.write_text(body)
    return p


def test_build_gene2pathway_symbol_join(tmp_path):
    links, names, genes = _write3(tmp_path)
    gff = _gff(tmp_path)
    g2p, meta, sym, stats = build_gene2pathway(gff, links, names, genes)
    assert g2p["LT_A"] == {"eco00260"}          # thrA -> LT_A
    assert g2p["LT_C"] == {"eco00010"}          # aceE -> LT_C
    assert stats["n_annotated"] == 3
    assert stats["n_pathways"] == 2


def test_build_gene2pathway_drops_ambiguous_symbol(tmp_path):
    # aceE iki locus'a (LT_C, LT_C2) -> belirsiz sembol join'de ATILIR (yalancı yok)
    links, names, genes = _write3(tmp_path)
    gff = _gff(tmp_path, extra="NZ\tX\tCDS\t1\t9\t.\t+\t0\tgene=aceE;locus_tag=LT_C2\n")
    g2p, meta, sym, stats = build_gene2pathway(gff, links, names, genes)
    assert "LT_C" not in g2p and "LT_C2" not in g2p   # aceE belirsiz -> atıldı
    assert g2p["LT_A"] == {"eco00260"}                # thrA hâlâ benzersiz
