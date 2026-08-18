"""m09 — GO fonksiyonel zenginleştirme (ORA). m06 DE çıktısından artan/azalan DEG'ler için
GO over-representation + dot-plot. Organizma-agnostik; YENİ veri-kapısı YOK (verdict m06'dan taşınır)."""
from __future__ import annotations

import json
from pathlib import Path

from rnaforge.config import Config
from rnaforge.enrichment import (
    all_tested_genes, build_enrichment_manifest, deg_sets, run_enrichment_r,
    run_ora, write_enrichment_manifest, write_ora_tsv,
)
from rnaforge.go_annotation import build_gene2go, parse_obo
from rnaforge.state import RunState

MODULE_NAME = "m09_enrichment"


def _require_obo(config: Config) -> Path:
    obo = config.enrichment.obo
    if obo is None:
        raise ValueError(
            "m09 (enrich) requires config.enrichment.obo (go-basic.obo path). "
            "Download it once: `curl -L -o references/go/go-basic.obo "
            "http://purl.obolibrary.org/obo/go/go-basic.obo` and set enrichment.obo.")
    if not Path(obo).exists():
        raise FileNotFoundError(
            f"m09 (enrich): go-basic.obo not found at {obo}. "
            "Download: `curl -L -o " + str(obo) + " "
            "http://purl.obolibrary.org/obo/go/go-basic.obo`.")
    return Path(obo)


def run_enrichment(config: Config, metadata_path: Path, run_dir: Path,
                   force: bool = False) -> dict:
    run_dir = Path(run_dir)
    de_dir = run_dir / "differential_expression"
    enrich_dir = run_dir / "enrichment"
    stats_dir = run_dir / "statistics"
    logs_dir = run_dir / "logs"
    for d in (enrich_dir, stats_dir, logs_dir):
        d.mkdir(parents=True, exist_ok=True)
    state = RunState(run_dir)
    stats_path = stats_dir / "enrichment_statistics.json"

    if not force and state.is_done(MODULE_NAME) and stats_path.exists():
        summary = json.loads(stats_path.read_text()); summary["resumed"] = True
        return summary
    if not state.is_done("m06_de"):
        raise ValueError(
            "m09 (enrich) requires m06 (de) to have completed in this run directory "
            f"first: {run_dir}. Run `rnaforge de` with the same --run-id, then re-run enrich.")

    gff = config.reference.annotation_gff              # ökaryotta None
    transcriptome = config.reference.transcriptome_fasta  # ökaryot sembol kaynağı
    obo_path = _require_obo(config)          # eksikse yüksek sesle hata (sessiz skip yok)
    deseq_tsv = de_dir / "deseq2_results.tsv"

    log_path = logs_dir / "enrichment.log"
    with log_path.open("w") as log_file:
        def log(msg):
            log_file.write(msg + "\n")

        obo = parse_obo(obo_path)
        state.heartbeat()
        gene2go, go_meta, direct, sources, ann_stats, gene_symbol = build_gene2go(
            gff, obo, gaf_path=config.enrichment.gaf,
            transcriptome_fasta=transcriptome, log=log)
        _write_gene2go_tsv(enrich_dir / "gene2go.tsv", gene2go, direct, sources, go_meta)
        state.heartbeat()

        up, down = deg_sets(deseq_tsv, config.de.fdr_threshold, config.de.log2fc_threshold)
        background = all_tested_genes(deseq_tsv)
        mts = config.enrichment.min_term_size
        up_rows = run_ora(up, background, gene2go, go_meta, gene_symbol, mts)
        down_rows = run_ora(down, background, gene2go, go_meta, gene_symbol, mts)
        write_ora_tsv(up_rows, enrich_dir / "enrichment_up.tsv")
        write_ora_tsv(down_rows, enrich_dir / "enrichment_down.tsv")
        log(f"m09: up_degs={len(up)} down_degs={len(down)} "
            f"terms_up={len(up_rows)} terms_down={len(down_rows)}")
        state.heartbeat()

        r_out = run_enrichment_r(enrich_dir / "enrichment_up.tsv",
                                 enrich_dir / "enrichment_down.tsv",
                                 enrich_dir, config.enrichment.top_n)
        if r_out:
            log_file.write(r_out if r_out.endswith("\n") else r_out + "\n")
        manifest = build_enrichment_manifest(enrich_dir)
        write_enrichment_manifest(enrich_dir)

        n_sig = lambda rows: sum(1 for r in rows if r["p_adj"] < 0.05)
        summary = {
            "n_up_degs": len(up), "n_down_degs": len(down),
            "n_terms_up": len(up_rows), "n_terms_down": len(down_rows),
            "n_sig_up": n_sig(up_rows), "n_sig_down": n_sig(down_rows),
            "background_size": len(background),
            "n_annotated": ann_stats["n_annotated"],
            "n_gff": ann_stats["n_gff"], "n_goa": ann_stats["n_goa"],
            "n_figures": len(manifest["figures"]),
        }
        stats_path.write_text(json.dumps(summary, indent=2))
        log(f"m09 enrichment done: {summary}")

    # GATE YOK — verdict m06/m07'den değişmeden taşınır (m09 dokunmaz).
    state.mark_done(MODULE_NAME, [str(stats_path), str(log_path)])
    return summary


def _write_gene2go_tsv(path, gene2go, direct, sources, go_meta):
    """Denetim izi: gene, go_id, namespace, name, source(GFF|GOA|—), direct|propagated."""
    with Path(path).open("w") as f:
        f.write("gene\tgo_id\tnamespace\tname\tsource\tkind\n")
        for gene in sorted(gene2go):
            for go_id in sorted(gene2go[gene]):
                ns, name = go_meta.get(go_id, ("?", ""))
                is_direct = go_id in direct.get(gene, set())
                src = sources.get((gene, go_id), "—") if is_direct else "—"
                kind = "direct" if is_direct else "propagated"
                f.write(f"{gene}\t{go_id}\t{ns}\t{name}\t{src}\t{kind}\n")
