"""m10 — KEGG pathway zenginleştirme (ORA). m06 DEG'lerinden KEGG pathway ORA + dot-plot.
m09 motorunu (enrichment.py) DEĞİŞTİRMEDEN kullanır. Gate YOK; verdict m06/m07'den taşınır."""
from __future__ import annotations

import json
from pathlib import Path

from rnaforge.config import Config
from rnaforge.enrichment import (
    all_tested_genes, build_enrichment_manifest, deg_sets, run_enrichment_r,
    run_ora, write_enrichment_manifest, write_ora_tsv,
)
from rnaforge.kegg_annotation import build_gene2pathway
from rnaforge.state import RunState

MODULE_NAME = "m10_kegg"
_KEGG_FILES = ("pathway_links.tsv", "pathway_names.tsv", "gene_list.tsv")


def _kegg_dir(config: Config) -> Path:
    org = config.enrichment.kegg_organism
    if not org:
        raise ValueError(
            "m10 (kegg) requires config.enrichment.kegg_organism (KEGG organism code, e.g. 'eco'). "
            "Without it KEGG cannot be mapped; set enrichment.kegg_organism.")
    d = config.enrichment.kegg_dir or Path("references/kegg") / org
    return Path(d)


def _require_kegg_files(kegg_dir: Path, org: str) -> tuple[Path, Path, Path]:
    paths = tuple(kegg_dir / f for f in _KEGG_FILES)
    missing = [p for p in paths if not p.exists()]
    if missing:
        raise FileNotFoundError(
            f"m10 (kegg): KEGG files missing in {kegg_dir}: {', '.join(p.name for p in missing)}. "
            f"Download once from KEGG REST (academic use):\n"
            f"  curl -s https://rest.kegg.jp/link/pathway/{org} > {kegg_dir/'pathway_links.tsv'}\n"
            f"  curl -s https://rest.kegg.jp/list/pathway/{org} > {kegg_dir/'pathway_names.tsv'}\n"
            f"  curl -s https://rest.kegg.jp/list/{org} > {kegg_dir/'gene_list.tsv'}")
    return paths


def run_kegg(config: Config, metadata_path: Path, run_dir: Path,
             force: bool = False) -> dict:
    run_dir = Path(run_dir)
    de_dir = run_dir / "differential_expression"
    out_dir = run_dir / "kegg"
    stats_dir = run_dir / "statistics"
    logs_dir = run_dir / "logs"
    for d in (out_dir, stats_dir, logs_dir):
        d.mkdir(parents=True, exist_ok=True)
    state = RunState(run_dir)
    stats_path = stats_dir / "kegg_statistics.json"

    if not force and state.is_done(MODULE_NAME) and stats_path.exists():
        summary = json.loads(stats_path.read_text()); summary["resumed"] = True
        return summary
    if not state.is_done("m06_de"):
        raise ValueError(
            "m10 (kegg) requires m06 (de) to have completed in this run directory "
            f"first: {run_dir}. Run `rnaforge de` with the same --run-id, then re-run kegg.")

    kegg_dir = _kegg_dir(config)                        # kegg_organism yoksa gürültülü hata
    links, names, genelist = _require_kegg_files(kegg_dir, config.enrichment.kegg_organism)
    gff = config.reference.annotation_gff
    deseq_tsv = de_dir / "deseq2_results.tsv"

    log_path = logs_dir / "kegg.log"
    with log_path.open("w") as log_file:
        gene2pathway, pathway_meta, gene_symbol, ann_stats = build_gene2pathway(
            gff, links, names, genelist)
        _write_gene2pathway_tsv(out_dir / "gene2pathway.tsv", gene2pathway, pathway_meta)
        state.heartbeat()

        up, down = deg_sets(deseq_tsv, config.de.fdr_threshold, config.de.log2fc_threshold)
        background = all_tested_genes(deseq_tsv)
        mts = config.enrichment.min_term_size
        up_rows = run_ora(up, background, gene2pathway, pathway_meta, gene_symbol, mts)
        down_rows = run_ora(down, background, gene2pathway, pathway_meta, gene_symbol, mts)
        write_ora_tsv(up_rows, out_dir / "kegg_up.tsv")
        write_ora_tsv(down_rows, out_dir / "kegg_down.tsv")
        log_file.write(f"m10: up_degs={len(up)} down_degs={len(down)} "
                       f"pathways_up={len(up_rows)} pathways_down={len(down_rows)}\n")
        state.heartbeat()

        r_out = run_enrichment_r(out_dir / "kegg_up.tsv", out_dir / "kegg_down.tsv",
                                 out_dir, config.enrichment.top_n,
                                 title_prefix="KEGG Pathway", basename_prefix="kegg")
        if r_out:
            log_file.write(r_out if r_out.endswith("\n") else r_out + "\n")
        manifest = build_enrichment_manifest(out_dir, basename_prefix="kegg")
        write_enrichment_manifest(out_dir, basename_prefix="kegg")

        n_sig = lambda rows: sum(1 for r in rows if r["p_adj"] < 0.05)
        summary = {
            "n_up_degs": len(up), "n_down_degs": len(down),
            "n_pathways_up": len(up_rows), "n_pathways_down": len(down_rows),
            "n_sig_up": n_sig(up_rows), "n_sig_down": n_sig(down_rows),
            "background_size": len(background),
            "n_annotated": ann_stats["n_annotated"], "n_pathways": ann_stats["n_pathways"],
            "n_figures": len(manifest["figures"]),
            "organism": config.enrichment.kegg_organism,
        }
        stats_path.write_text(json.dumps(summary, indent=2))
        log_file.write(f"m10 kegg done: {summary}\n")

    # GATE YOK — verdict m06/m07'den değişmeden taşınır.
    state.mark_done(MODULE_NAME, [str(stats_path), str(log_path)])
    return summary


def _write_gene2pathway_tsv(path, gene2pathway, pathway_meta):
    """Denetim izi: gene, pathway_id, name."""
    with Path(path).open("w") as f:
        f.write("gene\tpathway_id\tname\n")
        for gene in sorted(gene2pathway):
            for pid in sorted(gene2pathway[gene]):
                _, name = pathway_meta.get(pid, ("KEGG", ""))
                f.write(f"{gene}\t{pid}\t{name}\n")
