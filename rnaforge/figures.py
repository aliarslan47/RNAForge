"""m07 — Visualization helpers + R runner. Pure Python side; plotting lives in figures.R."""
from __future__ import annotations
from pathlib import Path


def gene_name_map(gff_path: Path) -> dict[str, str]:
    """CDS satırlarından locus_tag -> gene adı. gene= yoksa map'e girmez."""
    out: dict[str, str] = {}
    for line in Path(gff_path).read_text().splitlines():
        if not line or line.startswith("#") or "\tCDS\t" not in line:
            continue
        attrs = line.rstrip().split("\t")[8]
        d = dict(kv.split("=", 1) for kv in attrs.split(";") if "=" in kv)
        lt, gene = d.get("locus_tag"), d.get("gene")
        if lt and gene and lt not in out:
            out[lt] = gene
    return out
