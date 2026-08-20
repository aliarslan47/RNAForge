"""m_taxonomy — Metatranskriptom: Kraken2/Bracken topluluk kompozisyonu (DIAGNOSTIC).

Per-sample Kraken2 taksonomik profilleme + Bracken bolluk zenginleştirmesi, rRNA'sı
çıkarılmış okumalar üzerinde (`m_rrna_deplete.rrna_depleted_reads`). Her örneğin
Bracken fraksiyonları (`kraken2.parse_bracken`) `taxonomy/<sid>.bracken` altına yazılır,
sonra tüm örnekler tek `taxonomy/abundance_matrix.tsv` bolluk matrisine birleştirilir
(satır=taxon, sütun=örnek; bir örnekte eksik taxon = 0.0).

KAPI YOK — saf tanısal topluluk kompozisyonu. Bu modül biyoloji üzerine ASLA raise
etmez; yalnız eksik önkoşul (m_rrna_deplete tamamlanmamış) veya girdi eksikliği
(rRNA'sız okuma bulunamadı) gibi durumlarda yüksek sesle durur (Kural 7)."""
from __future__ import annotations

import csv
import json
from pathlib import Path

from rnaforge.config import Config
from rnaforge.kraken2 import parse_bracken, run_bracken, run_kraken2
from rnaforge.metadata import Sample, load_metadata
from rnaforge.modules.m_rrna_deplete import rrna_depleted_reads
from rnaforge.state import RunState

MODULE_NAME = "m_taxonomy"


def bracken_path(run_dir: Path, sample: Sample) -> Path:
    """m_taxonomy'nin bir örnek için ürettiği Bracken çıktı yolu (sözleşme)."""
    return Path(run_dir) / "taxonomy" / f"{sample.sample_id}.bracken"


def build_abundance_matrix(
    per_sample: dict[str, dict[str, float]],
) -> tuple[list[str], dict[str, dict[str, float]]]:
    """Örnek başına {taxon: fraction} sözlüklerini birleştirir.

    Returns: (sorted taxon adları, {taxon: {sample_id: fraction}}) — her örnekte HER
    taxon için bir değer vardır; o örnekte görülmeyen taxon 0.0 olur (union, eksik=0)."""
    taxa = sorted({taxon for fractions in per_sample.values() for taxon in fractions})
    matrix = {
        taxon: {sid: fractions.get(taxon, 0.0) for sid, fractions in per_sample.items()}
        for taxon in taxa
    }
    return taxa, matrix


def write_abundance_matrix(
    path: Path, sample_ids: list[str], taxa: list[str],
    matrix: dict[str, dict[str, float]],
) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["taxon", *sample_ids])
        for taxon in taxa:
            writer.writerow([taxon, *(matrix[taxon][sid] for sid in sample_ids)])
    return path


def run_taxonomy(config: Config, metadata_path: Path, run_dir: Path,
                 force: bool = False) -> dict:
    run_dir = Path(run_dir)
    out_dir = run_dir / "taxonomy"
    stats_dir = run_dir / "statistics"
    logs_dir = run_dir / "logs"
    for d in (out_dir, stats_dir, logs_dir):
        d.mkdir(parents=True, exist_ok=True)
    state = RunState(run_dir)
    stats_path = stats_dir / "taxonomy_abundance.json"
    matrix_path = out_dir / "abundance_matrix.tsv"

    if not force and state.is_done(MODULE_NAME) and stats_path.exists():
        summary = json.loads(stats_path.read_text())
        summary["resumed"] = True
        return summary

    if not state.is_done("m_rrna_deplete"):
        raise ValueError(
            "m_taxonomy requires m_rrna_deplete to have completed in this run directory "
            f"first: {run_dir}. Run `rnaforge rrna-deplete` with the same --run-id, then "
            "re-run taxonomy."
        )

    samples = load_metadata(metadata_path)
    log_path = logs_dir / "taxonomy.log"
    with log_path.open("w") as log_file:
        def log(msg: str) -> None:
            log_file.write(msg + "\n")
            log_file.flush()

        log(f"m_taxonomy: {len(samples)} sample(s), kraken2_db={config.taxonomy.kraken2_db}, "
            f"env={config.taxonomy.env}")

        per_sample: dict[str, dict[str, float]] = {}
        for sample in samples:
            state.heartbeat()
            sid = sample.sample_id
            out_path = bracken_path(run_dir, sample)
            if not force and state.is_item_done(MODULE_NAME, sid) and out_path.exists():
                per_sample[sid] = parse_bracken(out_path)
                log(f"{sid}: resumed (cached bracken, {len(per_sample[sid])} taxa)")
                continue

            reads = rrna_depleted_reads(run_dir, sample)
            if not reads:
                raise ValueError(
                    f"m_taxonomy: no rRNA-depleted reads found for sample {sid} at "
                    f"{run_dir / 'rrna_depleted' / sid} (m_rrna_deplete produced none). "
                    "Run `rnaforge rrna-deplete` first."
                )
            paired = len(reads) > 1
            kraken_prefix = run_dir / "taxonomy_work" / sid / "kraken2"
            report = run_kraken2(
                reads, config.taxonomy.kraken2_db, kraken_prefix, paired=paired,
                threads=config.resources.threads, env=config.taxonomy.env,
            )
            out_path.parent.mkdir(parents=True, exist_ok=True)
            run_bracken(
                report, config.taxonomy.kraken2_db, out_path,
                read_len=config.taxonomy.bracken_read_len,
                level=config.taxonomy.bracken_level, env=config.taxonomy.env,
            )
            fractions = parse_bracken(out_path)
            per_sample[sid] = fractions
            state.mark_item_done(MODULE_NAME, sid, {"n_taxa": len(fractions)})
            log(f"{sid}: {len(fractions)} taxa (bracken={out_path})")

        sample_ids = [s.sample_id for s in samples]
        taxa, matrix = build_abundance_matrix(per_sample)
        write_abundance_matrix(matrix_path, sample_ids, taxa, matrix)
        log(f"abundance matrix written: {matrix_path} ({len(taxa)} taxa x {len(sample_ids)} samples)")

        summary = {
            "n_samples": len(samples),
            "n_taxa": len(taxa),
            "samples": {sid: len(fractions) for sid, fractions in per_sample.items()},
            "abundance_matrix": str(matrix_path),
        }
        stats_path.write_text(json.dumps(summary, indent=2))
        log(f"m_taxonomy done: {len(samples)} sample(s), {len(taxa)} taxa")

    state.mark_done(MODULE_NAME, [str(stats_path), str(matrix_path), str(log_path)])
    return summary
