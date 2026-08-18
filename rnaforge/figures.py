"""m07 — Visualization helpers + R runner. Pure Python side; plotting lives in figures.R."""
from __future__ import annotations
import json
import subprocess
from pathlib import Path

FIGURE_SPECS: list[tuple[str, str, str]] = [
    ("pca", "01_pca", "PCA"),
    ("sample_correlation", "02_sample_correlation", "Örnek korelasyonu"),
    ("expression_dist", "03_expression_dist", "Ekspresyon dağılımı"),
    ("dispersion", "04_dispersion", "Dispersiyon"),
    ("pval_histogram", "05_pval_histogram", "p-değeri dağılımı"),
    ("volcano", "06_volcano", "Volcano"),
    ("ma", "07_ma", "MA plot"),
    ("heatmap", "08_heatmap", "Heatmap"),
]


def gene_name_map(gff_path: Path | None) -> dict[str, str]:
    """CDS satırlarından locus_tag -> gene adı. gene= yoksa map'e girmez.
    gff_path None ise (ökaryot: GFF yok, transcriptome_fasta+tx2gene) boş map
    döner → figürler gen kimliğiyle (ör. ENSG) etiketler."""
    out: dict[str, str] = {}
    if gff_path is None:
        return out
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


_SCRIPT = Path(__file__).parent / "scripts" / "figures.R"


def run_figures_r(de_dir: Path, gene_map: Path, fdr: float, lfc: float,
                  out_dir: Path, env: str = "rnaforge-de") -> str:
    """Run figures.R; return its combined stdout/stderr for logging. Raise loudly on failure."""
    de_dir, out_dir = Path(de_dir), Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = ["conda", "run", "-n", env, "Rscript", str(_SCRIPT),
           str(de_dir / "normalized_counts.tsv"), str(de_dir / "deseq2_results.tsv"),
           str(de_dir / "coldata.tsv"), str(gene_map), str(de_dir / "dispersions.tsv"),
           str(fdr), str(lfc), str(out_dir)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"figures.R failed (exit {r.returncode}):\n{r.stderr}")
    return (r.stdout or "") + (r.stderr or "")
