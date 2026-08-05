"""m17 — Hizalama-sonrası QC (diagnostik): insert-size dağılımı (samtools stats),
kontig derinliği/coverage (samtools coverage) ve okuma dağılımı (RSeQC
read_distribution). m04 sonrası koşar. KAPI ÜRETMEZ — saf tanısal figür/tablo;
verdict'i etkilemez. Figürler best-effort'tur; başarısızlık SESSİZCE yutulmaz,
log'a ve stats'e yazılır."""
from __future__ import annotations

import json
from pathlib import Path

from rnaforge.alignqc import (
    aggregate_histograms, parse_read_distribution, parse_samtools_coverage,
    parse_samtools_stats, run_read_distribution, run_samtools_coverage, run_samtools_stats,
)
from rnaforge.config import Config
from rnaforge.metadata import load_metadata
from rnaforge.qcplots import QCPlotError, render_qc_figure
from rnaforge.seqqc import gff_to_bed
from rnaforge.state import RunState

MODULE_NAME = "m17_alignqc"


def _bam(run_dir: Path, sample_id: str) -> Path:
    return run_dir / "quantification" / sample_id / "aligned.sorted.bam"


def run_alignqc(config: Config, metadata_path: Path, run_dir: Path,
                force: bool = False) -> dict:
    run_dir = Path(run_dir)
    out_dir = run_dir / "alignqc"
    stats_dir = run_dir / "statistics"
    logs_dir = run_dir / "logs"
    figures_dir = run_dir / "figures"
    for d in (out_dir, stats_dir, logs_dir, figures_dir):
        d.mkdir(parents=True, exist_ok=True)
    state = RunState(run_dir)
    stats_path = stats_dir / "alignqc_statistics.json"

    if not force and state.is_done(MODULE_NAME) and stats_path.exists():
        summary = json.loads(stats_path.read_text()); summary["resumed"] = True
        return summary
    if not state.is_done("m04_quant"):
        raise ValueError(
            "m17 (alignqc) requires m04 (quant) to have completed in this run directory "
            f"first: {run_dir}. Run `rnaforge quant` with the same --run-id, then re-run alignqc.")

    samples = load_metadata(metadata_path)
    log_path = logs_dir / "alignqc.log"
    with log_path.open("w") as log_file:
        def log(msg: str) -> None:
            log_file.write(msg + "\n"); log_file.flush()

        # --- F2 insert-size + F4 coverage: samtools ---
        insert_per_sample: dict[str, float] = {}
        sd_per_sample: dict[str, float] = {}
        histograms: list[list[tuple[int, int]]] = []
        paired = False
        contig_depths: dict[str, list[float]] = {}
        for s in samples:
            bam = _bam(run_dir, s.sample_id)
            if not bam.exists():
                log(f"{s.sample_id}: BAM yok, atlandı ({bam})")
                continue
            st = parse_samtools_stats(run_samtools_stats(bam))
            insert_per_sample[s.sample_id] = round(st["insert_size_average"], 2)
            sd_per_sample[s.sample_id] = round(st["insert_size_sd"], 2)
            if st["reads_paired"] > 0 and st["insert_size_average"] > 0:
                paired = True
                histograms.append(st["histogram"])
            for row in parse_samtools_coverage(run_samtools_coverage(bam)):
                contig_depths.setdefault(row["contig"], []).append(row["meandepth"])
            state.heartbeat()
            log(f"{s.sample_id}: insert≈{insert_per_sample[s.sample_id]} coverage OK")

        coverage_per_contig = {c: round(sum(v) / len(v), 2) for c, v in contig_depths.items()}
        genome_mean = (round(sum(coverage_per_contig.values()) / len(coverage_per_contig), 2)
                       if coverage_per_contig else 0.0)

        # --- F3 read-distribution: RSeQC (BED prokaryot CDS-only) ---
        bed = out_dir / "genes.bed"
        gff_to_bed(config.reference.annotation_gff, bed)
        rd_accum: dict[str, list[float]] = {}
        for s in samples:
            bam = _bam(run_dir, s.sample_id)
            if not bam.exists():
                continue
            rd = parse_read_distribution(run_read_distribution(bam, bed))
            for grp, pct in rd["percentages"].items():
                rd_accum.setdefault(grp, []).append(pct)
            state.heartbeat()
        read_distribution = {g: round(sum(v) / len(v), 2) for g, v in rd_accum.items()}
        log(f"read distribution (örnek ort.): {read_distribution}")

        # --- figürler (best-effort; hata sessizce yutulmaz) ---
        insert_size_figure = None
        coverage_figure = None
        figure_errors: dict[str, str] = {}
        if paired and histograms:
            labels, totals = aggregate_histograms(histograms)
            fig = figures_dir / "alignqc_insert_size.png"
            try:
                render_qc_figure({"type": "bars", "title": "Insert-size (fragment uzunluğu) dağılımı",
                                  "xlabel": "Insert-size (bç)", "ylabel": "Çift sayısı",
                                  "x": labels, "y": totals}, fig)
                insert_size_figure = fig.name
                log(f"insert-size figürü: {fig}")
            except QCPlotError as exc:
                figure_errors["insert_size"] = str(exc)
                log(f"WARNING: insert-size figürü BAŞARISIZ: {exc}")
        else:
            log("insert-size: paired-end değil, figür atlandı")

        if coverage_per_contig:
            fig = figures_dir / "alignqc_coverage.png"
            contigs = sorted(coverage_per_contig, key=coverage_per_contig.get, reverse=True)[:40]
            try:
                render_qc_figure({"type": "bars", "title": "Kontig başına ortalama derinlik (coverage)",
                                  "xlabel": "Kontig", "ylabel": "Ort. derinlik (×)",
                                  "x": contigs, "y": [coverage_per_contig[c] for c in contigs]}, fig)
                coverage_figure = fig.name
                log(f"coverage figürü: {fig}")
            except QCPlotError as exc:
                figure_errors["coverage"] = str(exc)
                log(f"WARNING: coverage figürü BAŞARISIZ: {exc}")

        summary = {
            "paired": paired,
            "insert_size_mean_per_sample": insert_per_sample,
            "insert_size_sd_per_sample": sd_per_sample,
            "insert_size_mean": (round(sum(insert_per_sample.values()) / len(insert_per_sample), 2)
                                 if insert_per_sample else 0.0),
            "coverage_per_contig": coverage_per_contig,
            "genome_mean_depth": genome_mean,
            "read_distribution": read_distribution,
            "insert_size_figure": insert_size_figure,
            "coverage_figure": coverage_figure,
            "figure_errors": figure_errors,
        }
        stats_path.write_text(json.dumps(summary, indent=2))
        log(f"m17 alignqc done: paired={paired} genome_mean_depth={genome_mean}")

    state.mark_done(MODULE_NAME, [str(stats_path), str(log_path)])
    return summary
