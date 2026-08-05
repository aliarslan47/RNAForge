"""m10 — KEGG pathway annotation (pure Python). Gene(locus_tag) -> KEGG pathway set.

KEGG REST 3 dosyasından (link/pathway, list/pathway, list/<org>) gen→pathway kurar; KEGG
b-number → gen sembolü → bizim locus_tag (TAM+BENZERSİZ join, m09 deseni — belirsiz ATILIR).
Global/overview haritaları (çok geniş) ORA'dan hariç tutulur. Organizma-agnostik (org config'ten).
"""
from __future__ import annotations

import re
from pathlib import Path

from rnaforge.go_annotation import _symbol_to_locus, parse_gff_go

# KEGG "Global and overview maps" — ORA'yı bozacak kadar geniş, dışlanır (numerik id, agnostik).
GLOBAL_MAPS = frozenset({
    "01100", "01110", "01120", "01200", "01210", "01212",
    "01220", "01230", "01232", "01240", "01250",
})

_NUM = re.compile(r"(\d+)$")


def _pathway_num(pid: str) -> str:
    m = _NUM.search(pid)
    return m.group(1) if m else ""


def _symbol_of(namefield: str) -> str:
    """KEGG list <org> 4. sütunundan gen sembolü: 'thrA; fused ...' -> 'thrA' (ilk `;`/`,` öncesi)."""
    return namefield.split(";", 1)[0].split(",", 1)[0].strip()


def parse_kegg(links_path: Path, names_path: Path, genelist_path: Path):
    """KEGG REST dosyaları -> (symbol2pathways, pathway_meta).

    Returns:
        symbol2pathways: dict[symbol, set[pathway_id]]   (global map'ler hariç)
        pathway_meta: dict[pathway_id, ("KEGG", name)]   (org eki kırpılmış ad)
    """
    # gene_id -> symbol
    id2sym: dict[str, str] = {}
    for line in Path(genelist_path).read_text().splitlines():
        if not line:
            continue
        cols = line.split("\t")
        if len(cols) < 4:
            continue
        sym = _symbol_of(cols[3])
        if sym:
            id2sym[cols[0]] = sym

    # pathway_id -> name (org eki kırp: "Glycolysis - Escherichia coli ..." -> "Glycolysis")
    names: dict[str, str] = {}
    for line in Path(names_path).read_text().splitlines():
        if not line:
            continue
        cols = line.split("\t")
        if len(cols) < 2:
            continue
        pid = cols[0].replace("path:", "")
        names[pid] = cols[1].split(" - ", 1)[0].strip()

    # gene_id -> pathways (global hariç); sonra symbol'e taşı
    symbol2pathways: dict[str, set[str]] = {}
    pathway_meta: dict[str, tuple[str, str]] = {}
    for line in Path(links_path).read_text().splitlines():
        if not line:
            continue
        cols = line.split("\t")
        if len(cols) < 2:
            continue
        gene_id = cols[0]
        pid = cols[1].replace("path:", "")
        if _pathway_num(pid) in GLOBAL_MAPS:
            continue
        sym = id2sym.get(gene_id)
        if not sym:
            continue
        symbol2pathways.setdefault(sym, set()).add(pid)
        pathway_meta.setdefault(pid, ("KEGG", names.get(pid, pid)))
    return symbol2pathways, pathway_meta


def build_gene2pathway(gff_path: Path, links_path: Path, names_path: Path,
                       genelist_path: Path):
    """GFF sembolleriyle KEGG'i birleştir: gene(locus_tag) -> KEGG pathway set.

    Returns:
        gene2pathway: dict[locus_tag, set[pathway_id]]
        pathway_meta: dict[pathway_id, ("KEGG", name)]
        gene_symbol: dict[locus_tag, symbol]   (run_ora 'genes' sütunu için)
        stats: dict
    """
    symbol2pathways, pathway_meta = parse_kegg(links_path, names_path, genelist_path)
    _, _, gene_symbol = parse_gff_go(gff_path)
    sym2lt = _symbol_to_locus(gene_symbol)      # belirsiz sembol ATILIR (yalancı yok)

    gene2pathway: dict[str, set[str]] = {}
    for sym, pathways in symbol2pathways.items():
        lt = sym2lt.get(sym)
        if lt is not None and pathways:
            gene2pathway.setdefault(lt, set()).update(pathways)

    stats = {
        "n_annotated": len(gene2pathway),
        "n_pathways": len({p for ps in gene2pathway.values() for p in ps}),
    }
    return gene2pathway, pathway_meta, gene_symbol, stats
