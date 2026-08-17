"""m11 — GSEA. m06 ranked listesinden (DESeq2 stat) fgsea ile gen-seti zenginleştirme.
Gen setleri m09 (GO) / m10 (KEGG) kurucularından; motor fgsea. Gate YOK; verdict m06/m07'den taşınır."""
from __future__ import annotations

import json
from pathlib import Path

from rnaforge.config import Config
from rnaforge.figures import write_gene_map
from rnaforge.go_annotation import build_gene2go, parse_obo
from rnaforge.gsea import invert_to_gmt, run_gsea_r, write_rnk
from rnaforge.kegg_annotation import build_gene2pathway
from rnaforge.modules.m10_kegg import _KEGG_FILES
from rnaforge.state import RunState

MODULE_NAME = "m11_gsea"


def _collection_stats(tsv: Path) -> dict:
    """gsea_<coll>.tsv -> {n_sets, n_sig_pos, n_sig_neg}."""
    lines = Path(tsv).read_text().splitlines() if Path(tsv).exists() else []
    if len(lines) < 2:
        return {"n_sets": 0, "n_sig_pos": 0, "n_sig_neg": 0}
    header = lines[0].split("\t")
    ni, pi = header.index("NES"), header.index("padj")
    pos = neg = 0
    for line in lines[1:]:
        c = line.split("\t")
        try:
            nes, padj = float(c[ni]), float(c[pi])
        except (ValueError, IndexError):
            continue
        if padj < 0.05:
            pos += nes > 0
            neg += nes < 0
    return {"n_sets": len(lines) - 1, "n_sig_pos": pos, "n_sig_neg": neg}


def run_gsea(config: Config, metadata_path: Path, run_dir: Path,
             force: bool = False) -> dict:
    run_dir = Path(run_dir)
    de_dir = run_dir / "differential_expression"
    out_dir = run_dir / "gsea"
    stats_dir = run_dir / "statistics"
    logs_dir = run_dir / "logs"
    for d in (out_dir, stats_dir, logs_dir):
        d.mkdir(parents=True, exist_ok=True)
    state = RunState(run_dir)
    stats_path = stats_dir / "gsea_statistics.json"

    if not force and state.is_done(MODULE_NAME) and stats_path.exists():
        summary = json.loads(stats_path.read_text()); summary["resumed"] = True
        return summary
    if not state.is_done("m06_de"):
        raise ValueError(
            "m11 (gsea) requires m06 (de) to have completed in this run directory "
            f"first: {run_dir}. Run `rnaforge de` with the same --run-id, then re-run gsea.")

    gff = config.reference.annotation_gff
    deseq_tsv = de_dir / "deseq2_results.tsv"
    log_path = logs_dir / "gsea.log"
    figures = []
    collections: dict[str, dict] = {}

    with log_path.open("w") as log_file:
        def log(m):
            log_file.write(m + "\n")

        n_ranked = write_rnk(deseq_tsv, out_dir / "ranked.rnk")   # stat yoksa gürültülü hata
        gene_map = out_dir / "gene_map.tsv"
        write_gene_map(gff, gene_map)
        log(f"m11: ranked genes={n_ranked}")
        state.heartbeat()

        planned = _resolve_collections(config, log)               # obo/kegg hazır olanlar
        if not planned:
            raise ValueError(
                "m11 (gsea): no gene-set collection available. Configure enrichment.obo (GO) "
                "and/or enrichment.kegg_organism (KEGG) with their reference files.")

        for coll, gene2set, meta, title in planned:
            gmt = out_dir / f"{coll}.gmt"
            n_sets = invert_to_gmt(gene2set, meta, gmt)
            log(f"m11: {coll} gene sets in GMT={n_sets}")
            r_out = run_gsea_r(out_dir / "ranked.rnk", gmt, gene_map, out_dir, coll,
                               config.enrichment.gsea_min_size, config.enrichment.gsea_max_size, title)
            if r_out:
                log_file.write(r_out if r_out.endswith("\n") else r_out + "\n")
            collections[coll] = _collection_stats(out_dir / f"gsea_{coll}.tsv")
            png = out_dir / f"gsea_{coll}.png"
            if png.exists():
                svg = out_dir / f"gsea_{coll}.svg"
                figures.append({"id": f"gsea_{coll}", "title": coll.upper(),
                                "png": png.name, "svg": svg.name if svg.exists() else None})
            state.heartbeat()

        (out_dir / "manifest.json").write_text(json.dumps({"figures": figures}, indent=2))
        summary = {
            "n_ranked": n_ranked,
            "collections": collections,
            "n_figures": len(figures),
        }
        stats_path.write_text(json.dumps(summary, indent=2))
        log(f"m11 gsea done: {summary}")

    # GATE YOK — verdict m06/m07'den değişmeden taşınır.
    state.mark_done(MODULE_NAME, [str(stats_path), str(log_path)])
    return summary


def _resolve_collections(config: Config, log):
    """Hazır gen-seti koleksiyonlarını çöz. obo/kegg konfigüre ama dosyası eksikse gürültülü hata;
    hiç konfigüre değilse atla (log)."""
    gff = config.reference.annotation_gff
    e = config.enrichment
    planned = []

    if e.obo is not None:
        if not Path(e.obo).exists():
            raise FileNotFoundError(
                f"m11 (gsea): GO requested but go-basic.obo not found at {e.obo}. "
                "Download it (see m09) or unset enrichment.obo.")
        obo = parse_obo(e.obo)
        gene2go, go_meta, _, _, _, _ = build_gene2go(gff, obo, gaf_path=e.gaf, log=log)
        planned.append(("go", gene2go, go_meta, "GSEA — Gene Ontology"))
    else:
        log("m11: enrichment.obo yok -> GO koleksiyonu atlandı")

    if e.kegg_organism:
        kegg_dir = e.kegg_dir or Path("references/kegg") / e.kegg_organism
        missing = [f for f in _KEGG_FILES if not (kegg_dir / f).exists()]
        if missing:
            raise FileNotFoundError(
                f"m11 (gsea): KEGG requested (organism={e.kegg_organism}) but files missing in "
                f"{kegg_dir}: {', '.join(missing)} (see m10 for download).")
        g2p, p_meta, _, _ = build_gene2pathway(
            gff, kegg_dir / "pathway_links.tsv", kegg_dir / "pathway_names.tsv",
            kegg_dir / "gene_list.tsv")
        planned.append(("kegg", g2p, p_meta, "GSEA — KEGG"))
    else:
        log("m11: enrichment.kegg_organism yok -> KEGG koleksiyonu atlandı")

    return planned
