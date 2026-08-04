"""m07 — Visualization. m06 DE ciktisindan 4 statik figur (PCA/Volcano/Heatmap/MA).
Organizma-agnostik; YENI veri-kapisi YOK (gorsel biyolojiyi gecersiz kilmaz)."""
from __future__ import annotations
import json
from pathlib import Path

from rnaforge.config import Config
from rnaforge.figures import run_figures_r, write_gene_map, write_manifest, build_manifest
from rnaforge.state import RunState

MODULE_NAME = "m07_figures"


def run_figures(config: Config, metadata_path: Path, run_dir: Path,
                force: bool = False) -> dict:
    run_dir = Path(run_dir)
    de_dir = run_dir / "differential_expression"
    fig_dir = run_dir / "figures"
    stats_dir = run_dir / "statistics"; logs_dir = run_dir / "logs"
    for d in (fig_dir, stats_dir, logs_dir):
        d.mkdir(parents=True, exist_ok=True)
    state = RunState(run_dir)
    stats_path = stats_dir / "figures_statistics.json"

    if not force and state.is_done(MODULE_NAME) and stats_path.exists():
        summary = json.loads(stats_path.read_text()); summary["resumed"] = True
        return summary
    if not state.is_done("m06_de"):
        raise ValueError(
            "m07 (figures) requires m06 (de) to have completed in this run directory "
            f"first: {run_dir}. Run `rnaforge de` with the same --run-id, then re-run figures.")

    log_path = logs_dir / "figures.log"
    with log_path.open("w") as log_file:
        gene_map = fig_dir / "gene_map.tsv"
        write_gene_map(config.reference.annotation_gff, gene_map)
        state.heartbeat()
        run_figures_r(de_dir, gene_map, config.de.fdr_threshold,
                      config.de.log2fc_threshold, fig_dir)
        manifest = build_manifest(fig_dir)   # eksik PNG -> yuksek sesle
        write_manifest(fig_dir)
        summary = {
            "n_figures": len(manifest["figures"]),
            "figures": {f["id"]: f["png"] for f in manifest["figures"]},
            "formats": ["png", "svg"],
        }
        stats_path.write_text(json.dumps(summary, indent=2))
        log_file.write(f"m07 figures done: {summary['n_figures']} figures\n")

    state.mark_done(MODULE_NAME, [str(stats_path), str(log_path)])
    return summary
