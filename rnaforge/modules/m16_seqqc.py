"""m16 — Sekans QC: rRNA% (SortMeRNA) + strandedness (RSeQC) + iki WARN kapısı. m04 sonrası.
Kötü girdiyi yakalar; verdict'e WARN olarak akar (asla FAIL üretmez)."""
from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

from rnaforge.config import Config
from rnaforge.gates import PASS, WARN, GateResult, write_gate_results
from rnaforge.metadata import load_metadata
from rnaforge.modules.m03_trim import trimmed_reads
from rnaforge.quality import load_profile
from rnaforge.seqqc import (
    gff_to_bed, parse_infer_experiment, parse_sortmerna_log, rrna_fasta_from_reference,
    run_infer_experiment, run_sortmerna, subsample_fastq,
)
from rnaforge.state import RunState

MODULE_NAME = "m16_seqqc"
RRNA_SUBSAMPLE = 200000     # rRNA% tahmini için örnek başına okuma (hız; QC tahmini yeterli)


def _bam(run_dir: Path, sample_id: str) -> Path:
    return run_dir / "quantification" / sample_id / "aligned.sorted.bam"


def _build_gates(mean_rrna: float, thr: float, inferred: str, declared: str) -> list[GateResult]:
    rrna = GateResult(
        name="rrna_fraction", module=MODULE_NAME,
        status=WARN if mean_rrna > thr else PASS,
        message=(f"ortalama rRNA fraksiyonu {mean_rrna:.1%} (eşik {thr:.0%}) — depletion "
                 + ("zayıf, kütüphane rRNA açısından zengin" if mean_rrna > thr else "yeterli")),
        remedy="rRNA yüksekse etkin okuma derinliği düşer; rRNA-depletion/rRNA-removal gözden geçirilmeli.",
        measured=round(mean_rrna, 4), threshold=thr)
    strand = GateResult(
        name="strandedness_match", module=MODULE_NAME,
        status=PASS if inferred == declared else WARN,
        message=(f"veriden çıkarılan strandedness '{inferred}', config beyanı '{declared}'"
                 + (" — uyumlu" if inferred == declared else " — UYUŞMUYOR, sayımlar yanlış olabilir")),
        remedy="Uyuşmazlıkta config.library.strandedness çıkarımla eşleşecek şekilde düzeltilip yeniden koşulmalı.")
    return [rrna, strand]


def run_seqqc(config: Config, metadata_path: Path, run_dir: Path, force: bool = False) -> dict:
    run_dir = Path(run_dir)
    out_dir = run_dir / "seqqc"
    stats_dir = run_dir / "statistics"
    logs_dir = run_dir / "logs"
    for d in (out_dir, stats_dir, logs_dir):
        d.mkdir(parents=True, exist_ok=True)
    state = RunState(run_dir)
    stats_path = stats_dir / "seqqc_statistics.json"

    if not force and state.is_done(MODULE_NAME) and stats_path.exists():
        summary = json.loads(stats_path.read_text()); summary["resumed"] = True
        return summary
    if not state.is_done("m04_quant"):
        raise ValueError(
            "m16 (seqqc) requires m04 (quant) to have completed in this run directory "
            f"first: {run_dir}. Run `rnaforge quant` with the same --run-id, then re-run seqqc.")

    samples = load_metadata(metadata_path)
    profile = load_profile(config.organism_type, config.quality)
    thr = profile.threshold("rrna_fraction")
    genome, gff = config.reference.genome_fasta, config.reference.annotation_gff
    log_path = logs_dir / "seqqc.log"

    with log_path.open("w") as log_file:
        # --- rRNA% (SortMeRNA), referans genomdan çıkarılan rRNA'ya karşı ---
        ref_fa = out_dir / "rrna_ref.fasta"
        n_rrna = rrna_fasta_from_reference(genome, gff, ref_fa)
        log_file.write(f"m16: rRNA referansı {n_rrna} dizi\n")
        per_sample: dict[str, float] = {}
        for s in samples:
            t1, t2 = trimmed_reads(run_dir, s)
            reads = []
            sub1 = out_dir / f"{s.sample_id}_sub_1.fastq"
            subsample_fastq(t1, RRNA_SUBSAMPLE, sub1); reads.append(sub1)
            if t2 is not None:
                sub2 = out_dir / f"{s.sample_id}_sub_2.fastq"
                subsample_fastq(t2, RRNA_SUBSAMPLE, sub2); reads.append(sub2)
            log = run_sortmerna(reads, ref_fa, out_dir / f"smr_{s.sample_id}",
                                threads=config.resources.threads)
            frac = parse_sortmerna_log(log)
            per_sample[s.sample_id] = round(frac, 4)
            for r in reads:
                r.unlink(missing_ok=True)
            state.heartbeat()
        mean_rrna = sum(per_sample.values()) / len(per_sample) if per_sample else 0.0
        log_file.write(f"m16: rRNA per-örnek {per_sample} mean={mean_rrna:.4f}\n")

        # --- strandedness (RSeQC infer_experiment), örnek başına, çoğunluk ---
        bed = out_dir / "genes.bed"
        gff_to_bed(gff, bed)
        inferred_per: dict[str, str] = {}
        for s in samples:
            bam = _bam(run_dir, s.sample_id)
            if not bam.exists():
                continue
            strand, _, _ = parse_infer_experiment(run_infer_experiment(bam, bed))
            inferred_per[s.sample_id] = strand
            state.heartbeat()
        inferred = Counter(inferred_per.values()).most_common(1)[0][0] if inferred_per else "unknown"
        declared = config.library.strandedness
        log_file.write(f"m16: strandedness per-örnek {inferred_per} -> {inferred} (beyan {declared})\n")

        gates = _build_gates(mean_rrna, thr, inferred, declared)
        write_gate_results(run_dir, gates)          # -> güvence kartına akar (WARN)

        summary = {
            "mean_rrna_fraction": round(mean_rrna, 4), "rrna_per_sample": per_sample,
            "inferred_strandedness": inferred, "declared_strandedness": declared,
            "strandedness_match": inferred == declared,
            "rrna_threshold": thr, "gate_counts": dict(Counter(g.status for g in gates)),
        }
        stats_path.write_text(json.dumps(summary, indent=2))
        log_file.write(f"m16 seqqc done: {summary}\n")

    state.mark_done(MODULE_NAME, [str(stats_path), str(log_path)])
    return summary
