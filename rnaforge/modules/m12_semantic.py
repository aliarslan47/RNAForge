"""m12 — Semantic reduction (REVIGO-benzeri). m09 ORA GO + m11 GSEA GO terimlerini Lin
benzerliğiyle temsilcilere indirger. Saf Python; figür yok; gate YOK (verdict m06'dan taşınır)."""
from __future__ import annotations

import csv
import json
from pathlib import Path

from rnaforge.config import Config
from rnaforge.go_annotation import build_gene2go, parse_obo
from rnaforge.semantic import (
    compute_ic, lin_distance_matrix, reduce_terms, run_semantic_r, write_distance_matrix,
)
from rnaforge.state import RunState

_FIG_TITLES = {"ora_up": "Anlamsal harita — Artan (GO)", "ora_down": "Anlamsal harita — Azalan (GO)",
               "gsea_go": "Anlamsal harita — GSEA (GO)"}

MODULE_NAME = "m12_semantic"
_REDUCED_HEADER = ["go_id", "namespace", "term", "padj", "n_collapsed", "members"]


def _read_sig_terms(tsv: Path, id_col: str, term_col: str, padj_col: str,
                    ns_col: str | None, go_meta: dict) -> list[dict]:
    """Bir zenginleştirme TSV'sinden anlamlı (padj<0.05) GO terimlerini çıkar."""
    if not tsv.exists():
        return []
    out = []
    with tsv.open() as f:
        for r in csv.DictReader(f, delimiter="\t"):
            gid = r.get(id_col, "")
            if not gid.startswith("GO:"):
                continue
            try:
                padj = float(r[padj_col])
            except (ValueError, KeyError, TypeError):
                continue
            if padj >= 0.05:
                continue
            ns = r.get(ns_col) if ns_col else go_meta.get(gid, ("?", ""))[0]
            out.append({"go_id": gid, "namespace": ns or "?",
                        "term": r.get(term_col, ""), "padj": padj})
    return out


def _write_reduced(reps: list[dict], path: Path) -> None:
    with Path(path).open("w") as f:
        f.write("\t".join(_REDUCED_HEADER) + "\n")
        for r in reps:
            f.write("\t".join([
                r["go_id"], r["namespace"], r["term"],
                f'{r["padj"]:.6e}' if r.get("padj") is not None else "",
                str(r["n_collapsed"]), ";".join(r["members"]),
            ]) + "\n")


def run_semantic(config: Config, metadata_path: Path, run_dir: Path,
                 force: bool = False) -> dict:
    run_dir = Path(run_dir)
    out_dir = run_dir / "semantic"
    stats_dir = run_dir / "statistics"
    logs_dir = run_dir / "logs"
    for d in (out_dir, stats_dir, logs_dir):
        d.mkdir(parents=True, exist_ok=True)
    state = RunState(run_dir)
    stats_path = stats_dir / "semantic_statistics.json"

    if not force and state.is_done(MODULE_NAME) and stats_path.exists():
        summary = json.loads(stats_path.read_text()); summary["resumed"] = True
        return summary
    if not (state.is_done("m09_enrichment") or state.is_done("m11_gsea")):
        raise ValueError(
            "m12 (semantic) requires a GO enrichment source: run `rnaforge enrich` (m09) "
            f"and/or `rnaforge gsea` (m11) with the same --run-id first: {run_dir}.")
    obo_path = config.enrichment.obo
    if obo_path is None or not Path(obo_path).exists():
        raise FileNotFoundError(
            f"m12 (semantic) requires config.enrichment.obo (go-basic.obo); not found: {obo_path}.")

    log_path = logs_dir / "semantic.log"
    with log_path.open("w") as log_file:
        def log(m):
            log_file.write(m + "\n")

        obo = parse_obo(obo_path)
        gene2go, go_meta, _, _, _ = build_gene2go(
            config.reference.annotation_gff, obo, gaf_path=config.enrichment.gaf, log=log)
        ic = compute_ic(gene2go)
        thr = config.enrichment.revigo_similarity
        state.heartbeat()

        enrich_dir = run_dir / "enrichment"
        gsea_dir = run_dir / "gsea"
        sources = [
            ("ora_up", _read_sig_terms(enrich_dir / "enrichment_up.tsv",
                                       "go_id", "term", "p_adj", "namespace", go_meta)),
            ("ora_down", _read_sig_terms(enrich_dir / "enrichment_down.tsv",
                                         "go_id", "term", "p_adj", "namespace", go_meta)),
            ("gsea_go", _read_sig_terms(gsea_dir / "gsea_go.tsv",
                                        "pathway_id", "name", "padj", None, go_meta)),
        ]
        collections: dict[str, dict] = {}
        figures = []
        for name, terms in sources:
            if not terms:
                continue
            reps = reduce_terms(terms, obo, ic, thr)
            _write_reduced(reps, out_dir / f"reduced_{name}.tsv")
            collections[name] = {"n_terms": len(terms), "n_representatives": len(reps)}
            log(f"m12: {name} {len(terms)} terms -> {len(reps)} representatives")
            # REVIGO-benzeri MDS scatter — best-effort (çekirdek tabloyu bozmasın; ≥3 temsilci gerekir)
            if len(reps) >= 3:
                rep_ids = [r["go_id"] for r in reps]
                dist = lin_distance_matrix(rep_ids, obo, ic)
                dist_path = out_dir / f"mds_{name}_dist.tsv"
                write_distance_matrix(rep_ids, dist, dist_path)
                try:
                    r_out = run_semantic_r(dist_path, out_dir / f"reduced_{name}.tsv", out_dir,
                                           f"mds_{name}", _FIG_TITLES.get(name, name))
                    if r_out:
                        log_file.write(r_out if r_out.endswith("\n") else r_out + "\n")
                    if (out_dir / f"mds_{name}.png").exists():
                        figures.append({"id": f"mds_{name}", "title": _FIG_TITLES.get(name, name),
                                        "png": f"mds_{name}.png",
                                        "svg": f"mds_{name}.svg" if (out_dir / f"mds_{name}.svg").exists() else None})
                except Exception as exc:  # noqa: BLE001 — figür opsiyonel; logla, koşuyu bozma
                    log(f"m12: {name} MDS figürü üretilemedi (opsiyonel): {exc}")
            state.heartbeat()

        if not collections:
            raise ValueError(
                "m12 (semantic): no significant GO terms found in m09/m11 outputs to reduce.")

        if figures:
            (out_dir / "manifest.json").write_text(json.dumps({"figures": figures}, indent=2))
        summary = {"collections": collections, "similarity_threshold": thr,
                   "n_figures": len(figures)}
        stats_path.write_text(json.dumps(summary, indent=2))
        log(f"m12 semantic done: {summary}")

    # GATE YOK — verdict m06/m07'den değişmeden taşınır.
    state.mark_done(MODULE_NAME, [str(stats_path), str(log_path)])
    return summary
