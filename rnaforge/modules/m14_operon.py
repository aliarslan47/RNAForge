"""m14 — Operon analizi. GFF'ten intergenik-mesafeyle operon tahmini + operon-düzeyi DE koordinasyonu.
Saf Python; gate YOK (verdict m06/m07'den taşınır); operonlar TAHMİN (deneysel değil)."""
from __future__ import annotations

import json
from pathlib import Path

from rnaforge.config import Config
from rnaforge.operon import aggregate_operon_de, predict_operons, run_operon_r
from rnaforge.state import RunState

MODULE_NAME = "m14_operon"
_TSV_HEADER = ["operon_id", "contig", "strand", "genes", "size", "n_tested",
               "n_deg", "n_up", "n_down", "mean_log2fc", "coordinated"]


def _write_operons_tsv(rows: list[dict], path: Path) -> None:
    with Path(path).open("w") as f:
        f.write("\t".join(_TSV_HEADER) + "\n")
        for r in rows:
            f.write("\t".join([
                r["operon_id"], r["contig"], r["strand"], ";".join(r["symbols"]),
                str(r["size"]), str(r["n_tested"]), str(r["n_deg"]),
                str(r["n_up"]), str(r["n_down"]),
                f'{r["mean_log2fc"]:.3f}' if r.get("mean_log2fc") is not None else "",
                "yes" if r["coordinated"] else "no",
            ]) + "\n")


def run_operon(config: Config, metadata_path: Path, run_dir: Path, force: bool = False) -> dict:
    run_dir = Path(run_dir)
    de_dir = run_dir / "differential_expression"
    out_dir = run_dir / "operon"
    stats_dir = run_dir / "statistics"
    logs_dir = run_dir / "logs"
    for d in (out_dir, stats_dir, logs_dir):
        d.mkdir(parents=True, exist_ok=True)
    state = RunState(run_dir)
    stats_path = stats_dir / "operon_statistics.json"

    if not force and state.is_done(MODULE_NAME) and stats_path.exists():
        summary = json.loads(stats_path.read_text()); summary["resumed"] = True
        return summary
    if not state.is_done("m06_de"):
        raise ValueError(
            "m14 (operon) requires m06 (de) to have completed in this run directory "
            f"first: {run_dir}. Run `rnaforge de` with the same --run-id, then re-run operon.")

    gff = config.reference.annotation_gff
    deseq_tsv = de_dir / "deseq2_results.tsv"
    max_gap = config.operon.max_gap
    log_path = logs_dir / "operon.log"

    with log_path.open("w") as log_file:
        operons = predict_operons(gff, max_gap)
        state.heartbeat()
        rows = aggregate_operon_de(operons, deseq_tsv, config.de.fdr_threshold,
                                   config.de.log2fc_threshold)
        _write_operons_tsv(rows, out_dir / "operons.tsv")
        n_multi = sum(1 for r in rows if r["size"] >= 2)
        n_coord = sum(1 for r in rows if r["coordinated"])
        # Figür best-effort: çekirdek tabloyu bozmasın; hata loglanır (sessiz değil).
        n_figures = 0
        try:
            r_out = run_operon_r(out_dir / "operons.tsv", out_dir)
            if r_out:
                log_file.write(r_out if r_out.endswith("\n") else r_out + "\n")
            if (out_dir / "operon_coord.png").exists():
                (out_dir / "manifest.json").write_text(json.dumps(
                    {"figures": [{"id": "operon_coord", "title": "Koordineli operonlar",
                                  "png": "operon_coord.png",
                                  "svg": "operon_coord.svg" if (out_dir / "operon_coord.svg").exists() else None}]},
                    indent=2))
                n_figures = 1
        except Exception as exc:  # noqa: BLE001 — figür opsiyonel; yüksek sesle logla, koşuyu bozma
            log_file.write(f"m14: operon figürü üretilemedi (opsiyonel): {exc}\n")
        summary = {
            "n_operons": len(rows), "n_multi_gene": n_multi,
            "n_coordinated": n_coord, "max_gap": max_gap, "n_figures": n_figures,
        }
        stats_path.write_text(json.dumps(summary, indent=2))
        log_file.write(f"m14 operon done: {summary}\n")

    # GATE YOK — verdict m06/m07'den değişmeden taşınır.
    state.mark_done(MODULE_NAME, [str(stats_path), str(log_path)])
    return summary
