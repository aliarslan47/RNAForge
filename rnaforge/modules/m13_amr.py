"""m13 — AMR + virülans overlay. abricate (CARD/VFDB) genome taraması → koordinatla locus_tag →
DE durumu overlay. Prokaryot-odaklı; gate YOK (verdict m06/m07'den taşınır)."""
from __future__ import annotations

import json
from pathlib import Path

from rnaforge.abricate import (
    gene_coords, map_hits_to_genes, overlay_de, parse_abricate, run_abricate,
)
from rnaforge.config import Config
from rnaforge.state import RunState

MODULE_NAME = "m13_amr"
_TSV_HEADER = ["gene", "locus_tag", "db", "label", "pct_identity", "pct_coverage",
               "log2fc", "padj", "de_status"]


def _write_amr_tsv(rows: list[dict], path: Path) -> None:
    with Path(path).open("w") as f:
        f.write("\t".join(_TSV_HEADER) + "\n")
        for r in rows:
            label = r.get("resistance") or r.get("product") or ""
            f.write("\t".join([
                r.get("symbol") or r.get("gene", ""), r.get("locus_tag", ""), r.get("db", ""),
                label, f'{r["pct_id"]:.1f}', f'{r["pct_cov"]:.1f}',
                f'{r["log2fc"]:.3f}' if r.get("log2fc") is not None else "",
                f'{r["padj"]:.3e}' if r.get("padj") is not None else "",
                r.get("de_status", ""),
            ]) + "\n")


def run_amr(config: Config, metadata_path: Path, run_dir: Path, force: bool = False) -> dict:
    run_dir = Path(run_dir)
    de_dir = run_dir / "differential_expression"
    out_dir = run_dir / "amr"
    stats_dir = run_dir / "statistics"
    logs_dir = run_dir / "logs"
    for d in (out_dir, stats_dir, logs_dir):
        d.mkdir(parents=True, exist_ok=True)
    state = RunState(run_dir)
    stats_path = stats_dir / "amr_statistics.json"

    if not force and state.is_done(MODULE_NAME) and stats_path.exists():
        summary = json.loads(stats_path.read_text()); summary["resumed"] = True
        return summary
    if not state.is_done("m06_de"):
        raise ValueError(
            "m13 (amr) requires m06 (de) to have completed in this run directory "
            f"first: {run_dir}. Run `rnaforge de` with the same --run-id, then re-run amr.")
    genome = config.reference.genome_fasta
    if genome is None or not Path(genome).exists():
        raise FileNotFoundError(
            f"m13 (amr) requires config.reference.genome_fasta (abricate scans it); not found: {genome}.")

    gff = config.reference.annotation_gff
    deseq_tsv = de_dir / "deseq2_results.tsv"
    a = config.amr
    genes = gene_coords(gff)
    log_path = logs_dir / "amr.log"

    def _one(db: str, out_name: str, log_file) -> dict:
        raw = out_dir / f"raw_{db}.tsv"
        stderr = run_abricate(genome, db, raw, env=a.env)
        if stderr:
            log_file.write(stderr if stderr.endswith("\n") else stderr + "\n")
        hits = parse_abricate(raw, a.min_identity, a.min_coverage)
        mapped, n_unmapped = map_hits_to_genes(hits, genes)
        rows = overlay_de(mapped, deseq_tsv, config.de.fdr_threshold, config.de.log2fc_threshold)
        _write_amr_tsv(rows, out_dir / out_name)
        n_de = sum(1 for r in rows if r["de_status"] in ("up", "down"))
        log_file.write(f"m13: {db} hits={len(hits)} mapped={len(rows)} "
                       f"unmapped={n_unmapped} DE={n_de}\n")
        return {"n_genes": len(rows), "n_de": n_de, "n_unmapped": n_unmapped}

    with log_path.open("w") as log_file:
        amr = _one(a.amr_db, "amr_genes.tsv", log_file)
        state.heartbeat()
        vir = _one(a.virulence_db, "virulence_genes.tsv", log_file)
        state.heartbeat()
        summary = {
            "amr_db": a.amr_db, "virulence_db": a.virulence_db,
            "n_amr_genes": amr["n_genes"], "n_amr_de": amr["n_de"],
            "n_vir_genes": vir["n_genes"], "n_vir_de": vir["n_de"],
            "n_unmapped_amr": amr["n_unmapped"], "n_unmapped_vir": vir["n_unmapped"],
        }
        stats_path.write_text(json.dumps(summary, indent=2))
        log_file.write(f"m13 amr done: {summary}\n")

    # GATE YOK — verdict m06/m07'den değişmeden taşınır.
    state.mark_done(MODULE_NAME, [str(stats_path), str(log_path)])
    return summary
