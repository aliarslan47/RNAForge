"""m15 — STRING PPI alt-ağı + Louvain community. STRING parser + sembol-join + networkx modül tespiti.

STRING protein id (`<taxid>.<b-number>`) → preferred_name (sembol) → bizim locus_tag (tam+benzersiz).
DEG-DEG kenarları alt-ağı → louvain_communities (deterministik seed). Elle Louvain yok (networkx güvenilir).
"""
from __future__ import annotations

import gzip
import subprocess
from pathlib import Path

import networkx as nx

from rnaforge.go_annotation import _symbol_to_locus

_SCRIPT = Path(__file__).parent / "scripts" / "ppi.R"


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


def network_layout(g: nx.Graph, communities: list[list[str]], gene_symbol: dict[str, str],
                   de: dict[str, tuple[float | None, float | None]], top_modules: int = 8,
                   min_size: int = 3, seed: int = 42):
    """En büyük modüllerin alt-ağını konumla (spring layout) → düğüm + kenar kayıtları (ggplot için).
    Hairball değil: yalnız en büyük `top_modules` modül. Boşsa ([],[])."""
    sized = sorted((c for c in communities if len(c) >= min_size), key=len, reverse=True)[:top_modules]
    node2mod = {lt: i for i, comm in enumerate(sized, 1) for lt in comm}
    if not node2mod:
        return [], []
    h = g.subgraph(node2mod.keys())
    pos = nx.spring_layout(h, seed=seed, weight="weight")
    nodes = []
    for n in h.nodes:
        x, y = pos[n]
        l2fc, _ = de.get(n, (None, None))
        direction = "up" if (l2fc or 0) > 0 else ("down" if (l2fc or 0) < 0 else "ns")
        nodes.append({"locus_tag": n, "symbol": gene_symbol.get(n, n),
                      "x": x, "y": y, "module": f"M{node2mod[n]}", "direction": direction,
                      "degree": h.degree(n)})
    edges = []
    for a, b in h.edges:
        edges.append({"x1": pos[a][0], "y1": pos[a][1], "x2": pos[b][0], "y2": pos[b][1]})
    return nodes, edges


def write_network_tsv(nodes: list[dict], edges: list[dict], nodes_path: Path, edges_path: Path) -> None:
    with Path(nodes_path).open("w") as f:
        f.write("locus_tag\tsymbol\tx\ty\tmodule\tdirection\tdegree\n")
        for n in nodes:
            f.write(f'{n["locus_tag"]}\t{n["symbol"]}\t{n["x"]:.5f}\t{n["y"]:.5f}\t'
                    f'{n["module"]}\t{n["direction"]}\t{n["degree"]}\n')
    with Path(edges_path).open("w") as f:
        f.write("x1\ty1\tx2\ty2\n")
        for e in edges:
            f.write(f'{e["x1"]:.5f}\t{e["y1"]:.5f}\t{e["x2"]:.5f}\t{e["y2"]:.5f}\n')


def run_ppi_r(nodes_tsv: Path, edges_tsv: Path, out_dir: Path, env: str = "rnaforge-de") -> str:
    """ppi.R (modül-renkli ağ figürü). stdout/stderr döndür, hatada gürültülü yüksel."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = ["conda", "run", "-n", env, "Rscript", str(_SCRIPT),
           str(nodes_tsv), str(edges_tsv), str(out_dir)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"ppi.R failed (exit {r.returncode}):\n{r.stderr}")
    return (r.stdout or "") + (r.stderr or "")
