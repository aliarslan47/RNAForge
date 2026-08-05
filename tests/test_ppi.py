"""m15 PPI testleri: STRING parser, sembol-join, alt-ağ, community, özet."""
from __future__ import annotations

import gzip

from rnaforge.ppi import (
    build_deg_network, detect_communities, parse_string_info, parse_string_links,
    string_to_locus, summarize_communities,
)


def _gz(path, text):
    with gzip.open(path, "wt") as f:
        f.write(text)
    return path


def test_parse_string_info(tmp_path):
    info = _gz(tmp_path / "info.gz",
              "#string_protein_id\tpreferred_name\tsize\tannotation\n"
              "511145.b0001\tthrL\t21\tleader\n511145.b0002\tthrA\t820\tkinase\n")
    d = parse_string_info(info)
    assert d == {"511145.b0001": "thrL", "511145.b0002": "thrA"}


def test_parse_string_links_filter(tmp_path):
    links = _gz(tmp_path / "links.gz",
               "protein1 protein2 combined_score\n"
               "511145.b0001 511145.b0002 850\n511145.b0001 511145.b0003 500\n")
    edges = parse_string_links(links, min_score=700)
    assert edges == [("511145.b0001", "511145.b0002", 850)]   # 500 elendi


def test_string_to_locus_symbol_join(tmp_path):
    info = {"511145.b0001": "thrL", "511145.b0002": "thrA", "511145.b9": "dupS"}
    gene_symbol = {"LT_1": "thrL", "LT_2": "thrA", "LT_3": "dupS", "LT_4": "dupS"}  # dupS belirsiz
    s2l = string_to_locus(info, gene_symbol)
    assert s2l["511145.b0001"] == "LT_1" and s2l["511145.b0002"] == "LT_2"
    assert "511145.b9" not in s2l                             # belirsiz sembol atıldı


def test_build_network_only_deg_deg():
    edges = [("s1", "s2", 900), ("s1", "s3", 800), ("s2", "sX", 900)]
    string2lt = {"s1": "LT_1", "s2": "LT_2", "s3": "LT_3", "sX": "LT_X"}
    deg = {"LT_1", "LT_2", "LT_3"}                            # LT_X DEG değil
    g = build_deg_network(deg, edges, string2lt)
    assert set(g.nodes) == {"LT_1", "LT_2", "LT_3"}
    assert g.has_edge("LT_1", "LT_2") and not g.has_edge("LT_2", "LT_X")
    assert g["LT_1"]["LT_2"]["weight"] == 0.9


def test_detect_communities_two_triangles():
    import networkx as nx
    g = nx.Graph()
    for a, b in [("A", "B"), ("B", "C"), ("A", "C")]:        # üçgen 1
        g.add_edge(a, b, weight=1.0)
    for a, b in [("D", "E"), ("E", "F"), ("D", "F")]:        # üçgen 2 (ayrık)
        g.add_edge(a, b, weight=1.0)
    comms = detect_communities(g, seed=42)
    assert len(comms) == 2
    assert {"A", "B", "C"} in [set(c) for c in comms]


def test_detect_communities_empty():
    import networkx as nx
    assert detect_communities(nx.Graph()) == []


def test_summarize_direction_and_filter():
    comms = [["LT_1", "LT_2", "LT_3"], ["LT_9"]]              # ikinci < min_size
    gene_symbol = {"LT_1": "a", "LT_2": "b", "LT_3": "c", "LT_9": "z"}
    de = {"LT_1": (2.0, 1e-5), "LT_2": (1.5, 1e-4), "LT_3": (-0.5, 1e-3)}
    out = summarize_communities(comms, gene_symbol, de, min_size=3)
    assert len(out) == 1                                     # tek gen modül elendi
    assert out[0]["size"] == 3 and out[0]["n_up"] == 2 and out[0]["dominant"] == "up"
    assert out[0]["genes"] == ["a", "b", "c"]
