"""m13 — AMR + virülans overlay. abricate (CARD/VFDB) + AMRFinderPlus genome taraması → koordinatla
locus_tag → DE durumu overlay. AMR tablosu iki aracı (CARD ↔ AMRFinderPlus) YAN YANA gösterir (konkordans).
Prokaryot-odaklı; gate YOK (verdict m06/m07'den taşınır)."""
from __future__ import annotations

import json
from pathlib import Path

from rnaforge.abricate import (
    gene_coords, map_hits_to_genes, overlay_de, parse_abricate, parse_amrfinder,
    run_abricate, run_amrfinder,
)
from rnaforge.config import Config
from rnaforge.state import RunState

MODULE_NAME = "m13_amr"
_AMR_HEADER = ["gene", "locus_tag", "card", "amrfinder", "pct_identity",
               "log2fc", "padj", "de_status"]
_VIR_HEADER = ["gene", "locus_tag", "db", "label", "pct_identity", "pct_coverage",
               "log2fc", "padj", "de_status"]


def _write_amr_tsv(rows: list[dict], path: Path) -> None:
    """AMR: CARD ve AMRFinderPlus sınıfları YAN YANA (konkordans)."""
    with Path(path).open("w") as f:
        f.write("\t".join(_AMR_HEADER) + "\n")
        for r in rows:
            f.write("\t".join([
                r.get("gene", ""), r.get("locus_tag", ""),
                r.get("card") or "—", r.get("amrfinder") or "—",
                f'{r["pct_id"]:.1f}' if r.get("pct_id") is not None else "",
                f'{r["log2fc"]:.3f}' if r.get("log2fc") is not None else "",
                f'{r["padj"]:.3e}' if r.get("padj") is not None else "",
                r.get("de_status", ""),
            ]) + "\n")


def _write_vir_tsv(rows: list[dict], path: Path) -> None:
    with Path(path).open("w") as f:
        f.write("\t".join(_VIR_HEADER) + "\n")
        for r in rows:
            label = r.get("resistance") or r.get("product") or ""
            f.write("\t".join([
                r.get("symbol") or r.get("gene", ""), r.get("locus_tag", ""), r.get("db", ""),
                label, f'{r["pct_id"]:.1f}', f'{r["pct_cov"]:.1f}',
                f'{r["log2fc"]:.3f}' if r.get("log2fc") is not None else "",
                f'{r["padj"]:.3e}' if r.get("padj") is not None else "",
                r.get("de_status", ""),
            ]) + "\n")


def _label(hit: dict | None) -> str:
    if not hit:
        return ""
    return hit.get("resistance") or hit.get("product") or ""


def _merge_amr(card: list[dict], afp: list[dict]) -> list[dict]:
    """CARD ∪ AMRFinderPlus'ı locus_tag'te birleştir (yan yana). Her gen için card/amrfinder sınıfı."""
    card_by = {h["locus_tag"]: h for h in card}
    afp_by = {h["locus_tag"]: h for h in afp}
    merged = []
    for lt in sorted(set(card_by) | set(afp_by)):
        c, a = card_by.get(lt), afp_by.get(lt)
        pcts = [h["pct_id"] for h in (c, a) if h]
        merged.append({
            "locus_tag": lt, "gene": (c or a).get("symbol") or (c or a).get("gene", ""),
            "card": _label(c), "amrfinder": _label(a),
            "pct_id": max(pcts) if pcts else None,
        })
    return merged


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
    fdr, lfc = config.de.fdr_threshold, config.de.log2fc_threshold
    log_path = logs_dir / "amr.log"

    with log_path.open("w") as log_file:
        # --- AMR: CARD (abricate) + AMRFinderPlus, yan yana ---
        card_raw = out_dir / f"raw_{a.amr_db}.tsv"
        stderr = run_abricate(genome, a.amr_db, card_raw, env=a.env)
        if stderr:
            log_file.write(stderr if stderr.endswith("\n") else stderr + "\n")
        card_hits = parse_abricate(card_raw, a.min_identity, a.min_coverage)
        card_mapped, _ = map_hits_to_genes(card_hits, genes)

        afp_mapped = []
        if a.amrfinder_organism:
            afp_raw = out_dir / "raw_amrfinderplus.tsv"
            run_amrfinder(genome, afp_raw, a.amrfinder_organism, env=a.amrfinder_env)
            afp_hits = parse_amrfinder(afp_raw, a.min_identity, a.min_coverage)
            afp_mapped, _ = map_hits_to_genes(afp_hits, genes)
            log_file.write(f"m13: AMRFinderPlus ({a.amrfinder_organism}) mapped={len(afp_mapped)}\n")
        else:
            log_file.write("m13: amrfinder_organism yok -> yalnız CARD\n")

        merged = _merge_amr(card_mapped, afp_mapped)
        merged = overlay_de(merged, deseq_tsv, fdr, lfc)
        _write_amr_tsv(merged, out_dir / "amr_genes.tsv")
        n_both = sum(1 for r in merged if r["card"] and r["amrfinder"])
        log_file.write(f"m13: AMR card={len(card_mapped)} amrfinder={len(afp_mapped)} "
                       f"union={len(merged)} both={n_both} "
                       f"DE={sum(1 for r in merged if r['de_status'] in ('up','down'))}\n")
        state.heartbeat()

        # --- Virülans: VFDB (tek araç) ---
        vir_raw = out_dir / f"raw_{a.virulence_db}.tsv"
        stderr = run_abricate(genome, a.virulence_db, vir_raw, env=a.env)
        if stderr:
            log_file.write(stderr if stderr.endswith("\n") else stderr + "\n")
        vir_hits = parse_abricate(vir_raw, a.min_identity, a.min_coverage)
        vir_mapped, _ = map_hits_to_genes(vir_hits, genes)
        vir_rows = overlay_de(vir_mapped, deseq_tsv, fdr, lfc)
        _write_vir_tsv(vir_rows, out_dir / "virulence_genes.tsv")
        state.heartbeat()

        summary = {
            "amr_db": a.amr_db, "virulence_db": a.virulence_db,
            "amrfinder_organism": a.amrfinder_organism,
            "n_amr_genes": len(merged),
            "n_amr_card": len(card_mapped), "n_amr_amrfinder": len(afp_mapped), "n_amr_both": n_both,
            "n_amr_de": sum(1 for r in merged if r["de_status"] in ("up", "down")),
            "n_vir_genes": len(vir_rows),
            "n_vir_de": sum(1 for r in vir_rows if r["de_status"] in ("up", "down")),
        }
        stats_path.write_text(json.dumps(summary, indent=2))
        log_file.write(f"m13 amr done: {summary}\n")

    # GATE YOK — verdict m06/m07'den değişmeden taşınır.
    state.mark_done(MODULE_NAME, [str(stats_path), str(log_path)])
    return summary
