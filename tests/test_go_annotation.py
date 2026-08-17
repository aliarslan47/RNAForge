"""m09 annotation birleştirme testleri: GFF parse, obo/propagation, GAF doldurma."""
from __future__ import annotations

import textwrap

from rnaforge.go_annotation import (
    build_gene2go, fill_from_gaf, parse_gff_go, parse_obo, propagate,
)

# thrA GO'lu (Ontology_term + go_process/function); thrL GO'suz (map'e girmemeli).
GFF_2CDS = (
    "NZ\tRefSeq\tCDS\t190\t255\t.\t+\t0\t"
    "ID=cds1;gene=thrL;locus_tag=BW_RS00005;product=leader\n"
    "NZ\tProtein Homology\tCDS\t337\t2799\t.\t+\t0\t"
    "ID=cds2;Ontology_term=GO:0008652,GO:0004072;gene=thrA;"
    "go_function=aspartate kinase activity|0004072||IEA;"
    "go_process=amino acid biosynthetic process|0008652||IEA;locus_tag=BW_RS00010\n"
)


def _write(tmp_path, name, body):
    p = tmp_path / name
    p.write_text(body)
    return p


def test_parse_gff_go_extracts_ids_namespace_name(tmp_path):
    gff = _write(tmp_path, "g.gff", GFF_2CDS)
    g2go, meta, sym = parse_gff_go(gff)
    assert g2go["BW_RS00010"] == {"GO:0008652", "GO:0004072"}
    assert meta["GO:0008652"] == ("BP", "amino acid biosynthetic process")
    assert meta["GO:0004072"] == ("MF", "aspartate kinase activity")
    assert sym["BW_RS00010"] == "thrA"


def test_parse_gff_go_skips_genes_without_go(tmp_path):
    gff = _write(tmp_path, "g.gff", GFF_2CDS)
    g2go, _, sym = parse_gff_go(gff)
    assert "BW_RS00005" not in g2go   # GO'suz gen map'te yok
    assert sym["BW_RS00005"] == "thrL"  # ama sembolü bilinir


OBO = textwrap.dedent("""\
    format-version: 1.2

    [Term]
    id: GO:0000001
    name: root
    namespace: biological_process

    [Term]
    id: GO:0000002
    name: mid
    namespace: biological_process
    is_a: GO:0000001 ! root

    [Term]
    id: GO:0000003
    name: leaf
    namespace: biological_process
    is_a: GO:0000002 ! mid
    relationship: part_of GO:0000001

    [Term]
    id: GO:0000009
    name: dead
    namespace: biological_process
    is_obsolete: true

    [Typedef]
    id: part_of
    name: part of
    """)


def test_parse_obo_parents_and_obsolete(tmp_path):
    obo = parse_obo(_write(tmp_path, "go.obo", OBO))
    assert obo["GO:0000003"]["parents"] == {"GO:0000002", "GO:0000001"}
    assert obo["GO:0000003"]["name"] == "leaf"
    assert obo["GO:0000003"]["namespace"] == "BP"
    assert obo["GO:0000009"]["obsolete"] is True
    assert obo["GO:0000002"]["parents"] == {"GO:0000001"}


def test_propagate_adds_ancestors(tmp_path):
    obo = parse_obo(_write(tmp_path, "go.obo", OBO))
    out = propagate({"geneX": {"GO:0000003"}}, obo)
    assert out["geneX"] == {"GO:0000003", "GO:0000002", "GO:0000001"}


def test_propagate_drops_obsolete(tmp_path):
    obo = parse_obo(_write(tmp_path, "go.obo", OBO))
    out = propagate({"geneY": {"GO:0000009", "GO:0000002"}}, obo)
    assert "GO:0000009" not in out["geneY"]           # obsolete atılır
    assert out["geneY"] == {"GO:0000002", "GO:0000001"}


def test_propagate_cycle_safe(tmp_path):
    # Yapay döngü: A->B, B->A. propagate takılmamalı.
    obo = {
        "GO:A": {"name": "a", "namespace": "BP", "parents": {"GO:B"}, "obsolete": False},
        "GO:B": {"name": "b", "namespace": "BP", "parents": {"GO:A"}, "obsolete": False},
    }
    out = propagate({"g": {"GO:A"}}, obo)
    assert out["g"] == {"GO:A", "GO:B"}


# --- GAF doldurma ---
# GFF: geneA GO'lu (thrA), geneB GO'suz (mysB), geneC GO'suz ama sembol çift (dupS).
GFF_FOR_GAF = (
    "NZ\tX\tCDS\t1\t9\t.\t+\t0\tOntology_term=GO:0000002;gene=thrA;locus_tag=LT_A\n"
    "NZ\tX\tCDS\t1\t9\t.\t+\t0\tgene=mysB;locus_tag=LT_B\n"
    "NZ\tX\tCDS\t1\t9\t.\t+\t0\tgene=dupS;locus_tag=LT_C1\n"
    "NZ\tX\tCDS\t1\t9\t.\t+\t0\tgene=dupS;locus_tag=LT_C2\n"
)
# GAF: thrA (GFF'te GO var -> eklenmez), mysB (benzersiz -> eklenir), dupS (belirsiz -> atılır).
GAF = (
    "!gaf-version: 2.2\n"
    "DB\tU1\tthrA\t\tGO:0009999\tref\tIEA\t\tP\t\t\tprotein\ttaxon\t2020\tDB\t\t\n"
    "DB\tU2\tmysB\t\tGO:0000003\tref\tIDA\t\tF\t\t\tprotein\ttaxon\t2020\tDB\t\t\n"
    "DB\tU3\tdupS\t\tGO:0000001\tref\tIEA\t\tC\t\t\tprotein\ttaxon\t2020\tDB\t\t\n"
)


def test_fill_from_gaf_only_ungapped_unique_symbol(tmp_path):
    gff = _write(tmp_path, "g.gff", GFF_FOR_GAF)
    gaf = _write(tmp_path, "e.gaf", GAF)
    gene2go, meta, gene_symbol = parse_gff_go(gff)
    additions, sources = fill_from_gaf(gene2go, gene_symbol, gaf, meta)
    assert additions == {"LT_B": {"GO:0000003"}}      # yalnız mysB
    assert sources[("LT_B", "GO:0000003")] == "GOA"
    assert "LT_A" not in additions                    # GFF otorite
    assert "LT_C1" not in additions and "LT_C2" not in additions  # belirsiz sembol atıldı
    assert meta["GO:0000003"][0] == "MF"              # aspect F -> MF


def test_build_gene2go_integrates_and_stamps(tmp_path):
    gff = _write(tmp_path, "g.gff", GFF_FOR_GAF)
    gaf = _write(tmp_path, "e.gaf", GAF)
    obo = parse_obo(_write(tmp_path, "go.obo", OBO))
    g2go, meta, direct, sources, stats, gene_symbol = build_gene2go(gff, obo, gaf_path=gaf)
    # build_gene2go artık gene_symbol'ü de döndürür (m09 GFF'i iki kez parse etmesin)
    assert isinstance(gene_symbol, dict)
    assert gene_symbol == parse_gff_go(gff)[2]
    # LT_A: GO:0000002 (GFF) -> propagate -> +GO:0000001
    assert g2go["LT_A"] == {"GO:0000002", "GO:0000001"}
    assert sources[("LT_A", "GO:0000002")] == "GFF"
    # LT_B: GO:0000003 (GAF) -> propagate -> +GO:0000002, GO:0000001
    assert g2go["LT_B"] == {"GO:0000003", "GO:0000002", "GO:0000001"}
    assert sources[("LT_B", "GO:0000003")] == "GOA"
    # direct propagation öncesi
    assert direct["LT_B"] == {"GO:0000003"}
    assert stats["n_gff"] == 1 and stats["n_goa"] == 1


def test_build_gene2go_without_gaf_logs(tmp_path):
    gff = _write(tmp_path, "g.gff", GFF_FOR_GAF)
    obo = parse_obo(_write(tmp_path, "go.obo", OBO))
    msgs = []
    g2go, meta, direct, sources, stats, gene_symbol = build_gene2go(gff, obo, gaf_path=None, log=msgs.append)
    assert stats["n_goa"] == 0
    assert any("GAF" in m for m in msgs)              # sessiz değil
