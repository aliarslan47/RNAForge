"""m15 — STRING PPI alt-ağı + Louvain community. STRING parser + sembol-join + networkx modül tespiti.

STRING protein id (`<taxid>.<b-number>`) → preferred_name (sembol) → bizim locus_tag (tam+benzersiz).
DEG-DEG kenarları alt-ağı → louvain_communities (deterministik seed). Elle Louvain yok (networkx güvenilir).
"""
from __future__ import annotations

import gzip
from pathlib import Path

import networkx as nx

from rnaforge.go_annotation import _symbol_to_locus


def parse_string_info(info_gz: Path) -> dict[str, str]:
    """STRING info (tab, `#string_protein_id preferred_name …`) -> {string_id: symbol}."""
    out: dict[str, str] = {}
    with gzip.open(Path(info_gz), "rt") as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) >= 2 and cols[1]:
                out[cols[0]] = cols[1]
    return out


def parse_string_links(links_gz: Path, min_score: int) -> list[tuple[str, str, int]]:
    """STRING links (BOŞLUKLA ayrılmış: `p1 p2 combined_score`) -> eşik üstü kenarlar."""
    edges = []
    with gzip.open(Path(links_gz), "rt") as f:
        header = True
        for line in f:
            if header:
                header = False
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            try:
                score = int(parts[2])
            except ValueError:
                continue
            if score >= min_score:
                edges.append((parts[0], parts[1], score))
    return edges


def string_to_locus(info: dict[str, str], gene_symbol: dict[str, str]) -> dict[str, str]:
    """string_id → symbol (info) + symbol → locus_tag (GFF, tam+benzersiz) → string_id → locus_tag."""
    sym2lt = _symbol_to_locus(gene_symbol)      # belirsiz sembol atılır
    out: dict[str, str] = {}
    for sid, sym in info.items():
        lt = sym2lt.get(sym)
        if lt is not None:
            out[sid] = lt
    return out


def build_deg_network(deg_ids: set[str], edges: list[tuple[str, str, int]],
                      string2lt: dict[str, str]) -> nx.Graph:
    """STRING kenarlarından DEG-DEG alt-ağı. İki ucu da DEG olan kenarlar; weight=score/1000."""
    g = nx.Graph()
    for a, b, score in edges:
        la, lb = string2lt.get(a), string2lt.get(b)
        if la is None or lb is None or la == lb:
            continue
        if la in deg_ids and lb in deg_ids:
            g.add_edge(la, lb, weight=score / 1000.0)
    return g


def detect_communities(g: nx.Graph, seed: int = 42) -> list[list[str]]:
    """Louvain community (ağırlıklı, deterministik seed). Boş graf -> []."""
    if g.number_of_edges() == 0:
        return []
    comms = nx.community.louvain_communities(g, weight="weight", seed=seed)
    return [sorted(c) for c in comms]


def summarize_communities(communities: list[list[str]], gene_symbol: dict[str, str],
                          de: dict[str, tuple[float | None, float | None]],
                          min_size: int = 3) -> list[dict]:
    """Modül başına üye sembol, boyut, n_up/n_down, dominant yön. `< min_size` elenir; boyuta göre sıralı."""
    out = []
    for i, members in enumerate(communities, 1):
        if len(members) < min_size:
            continue
        n_up = n_down = 0
        for lt in members:
            l2fc, _ = de.get(lt, (None, None))
            if l2fc is None:
                continue
            if l2fc > 0:
                n_up += 1
            elif l2fc < 0:
                n_down += 1
        dominant = "up" if n_up > n_down else ("down" if n_down > n_up else "mixed")
        symbols = sorted(gene_symbol.get(lt, lt) for lt in members)
        out.append({"community_id": f"module_{i}", "size": len(members),
                    "n_up": n_up, "n_down": n_down, "dominant": dominant, "genes": symbols})
    out.sort(key=lambda r: -r["size"])
    return out
