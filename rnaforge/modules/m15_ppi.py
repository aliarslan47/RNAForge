"""m15 — PPI + community detection. DEG'ler arasında STRING alt-ağı → Louvain modülleri.
Gate YOK (verdict m06/m07'den taşınır); STRING etkileşimleri kanıt-skorlu (rapora dürüst not)."""
from __future__ import annotations

import json
from pathlib import Path

from rnaforge.config import Config
from rnaforge.enrichment import deg_sets
from rnaforge.go_annotation import parse_gff_go
from rnaforge.ppi import (
    build_deg_network, detect_communities, parse_string_info, parse_string_links,
    string_to_locus, summarize_communities,
)
from rnaforge.state import RunState

MODULE_NAME = "m15_ppi"
_TSV_HEADER = ["community_id", "size", "n_up", "n_down", "dominant", "genes"]


def _write_communities_tsv(rows: list[dict], path: Path) -> None:
    with Path(path).open("w") as f:
        f.write("\t".join(_TSV_HEADER) + "\n")
        for r in rows:
            f.write("\t".join([
                r["community_id"], str(r["size"]), str(r["n_up"]), str(r["n_down"]),
                r["dominant"], ";".join(r["genes"]),
            ]) + "\n")


def _string_files(config: Config) -> tuple[Path, Path]:
    taxid = config.ppi.taxid
    if not taxid:
        raise ValueError(
            "m15 (ppi) requires config.ppi.taxid (STRING organism, e.g. '511145' for E. coli K-12). "
            "Without it the STRING network cannot be selected.")
    d = config.ppi.string_dir or Path("references/string") / taxid
    info = Path(d) / "protein.info.txt.gz"
    links = Path(d) / "protein.links.txt.gz"
    missing = [p for p in (info, links) if not p.exists()]
    if missing:
        raise FileNotFoundError(
            f"m15 (ppi): STRING files missing in {d}: {', '.join(p.name for p in missing)}. "
            f"Download once (academic use):\n"
            f"  curl -s https://stringdb-downloads.org/download/protein.info.v12.0/{taxid}.protein.info.v12.0.txt.gz -o {info}\n"
            f"  curl -s https://stringdb-downloads.org/download/protein.links.v12.0/{taxid}.protein.links.v12.0.txt.gz -o {links}")
    return info, links


def run_ppi(config: Config, metadata_path: Path, run_dir: Path, force: bool = False) -> dict:
    run_dir = Path(run_dir)
    de_dir = run_dir / "differential_expression"
    out_dir = run_dir / "ppi"
    stats_dir = run_dir / "statistics"
    logs_dir = run_dir / "logs"
    for d in (out_dir, stats_dir, logs_dir):
        d.mkdir(parents=True, exist_ok=True)
    state = RunState(run_dir)
    stats_path = stats_dir / "ppi_statistics.json"

    if not force and state.is_done(MODULE_NAME) and stats_path.exists():
        summary = json.loads(stats_path.read_text()); summary["resumed"] = True
        return summary
    if not state.is_done("m06_de"):
        raise ValueError(
            "m15 (ppi) requires m06 (de) to have completed in this run directory "
            f"first: {run_dir}. Run `rnaforge de` with the same --run-id, then re-run ppi.")

    info_gz, links_gz = _string_files(config)      # taxid/dosya yoksa gürültülü hata
    deseq_tsv = de_dir / "deseq2_results.tsv"
    p = config.ppi
    log_path = logs_dir / "ppi.log"

    with log_path.open("w") as log_file:
        _, _, gene_symbol = parse_gff_go(config.reference.annotation_gff)
        info = parse_string_info(info_gz)
        string2lt = string_to_locus(info, gene_symbol)
        up, down = deg_sets(deseq_tsv, config.de.fdr_threshold, config.de.log2fc_threshold)
        deg_ids = set(up) | set(down)
        de = {lt: (1.0, None) for lt in up}
        de.update({lt: (-1.0, None) for lt in down})
        log_file.write(f"m15: DEGs={len(deg_ids)} string_mapped={len(string2lt)}\n")
        state.heartbeat()

        edges = parse_string_links(links_gz, p.min_score)
        g = build_deg_network(deg_ids, edges, string2lt)
        comms = detect_communities(g, seed=42)
        rows = summarize_communities(comms, gene_symbol, de, p.min_community_size)
        _write_communities_tsv(rows, out_dir / "communities.tsv")
        state.heartbeat()

        summary = {
            "n_deg": len(deg_ids), "n_deg_in_network": g.number_of_nodes(),
            "n_edges": g.number_of_edges(), "n_communities": len(rows),
            "min_score": p.min_score, "taxid": p.taxid,
        }
        stats_path.write_text(json.dumps(summary, indent=2))
        log_file.write(f"m15 ppi done: {summary}\n")

    # GATE YOK — verdict m06/m07'den değişmeden taşınır.
    state.mark_done(MODULE_NAME, [str(stats_path), str(log_path)])
    return summary
