"""m07 — Visualization helpers + R runner. Pure Python side; plotting lives in figures.R."""
from __future__ import annotations
import json
from pathlib import Path

FIGURE_SPECS: list[tuple[str, str, str]] = [
    ("pca", "01_pca", "PCA"),
    ("volcano", "02_volcano", "Volcano"),
    ("heatmap", "03_heatmap", "Heatmap"),
    ("ma", "04_ma", "MA plot"),
]


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


def write_gene_map(gff_path: Path, out_path: Path) -> None:
    m = gene_name_map(gff_path)
    with Path(out_path).open("w") as f:
        f.write("locus_tag\tgene\n")
        for lt, g in m.items():
            f.write(f"{lt}\t{g}\n")


def build_manifest(fig_dir: Path) -> dict:
    fig_dir = Path(fig_dir)
    figures = []
    for _id, base, title in FIGURE_SPECS:
        png, svg = fig_dir / f"{base}.png", fig_dir / f"{base}.svg"
        if not png.exists():
            raise FileNotFoundError(f"m07 figure not rendered: {png} (figures.R failed?)")
        figures.append({"id": _id, "title": title, "png": png.name,
                        "svg": svg.name if svg.exists() else None})
    return {"figures": figures}


def write_manifest(fig_dir: Path) -> Path:
    p = Path(fig_dir) / "manifest.json"
    p.write_text(json.dumps(build_manifest(fig_dir), indent=2))
    return p
