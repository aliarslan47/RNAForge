"""m05 — Count Matrix (prokaryot: featureCounts).

m04 BAM'lerini anotasyona göre sayıp gen×örnek count matrisi (ortak sözleşme,
PLAN §5) üretir. Veri kapısı `assignment_rate`: featureCounts'un gene atadığı
okuma oranı profil eşiğinin altındaysa FAIL — çok düşük atama yanlış anotasyon/
tür demektir, sayımlar güvenilmez."""
from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

from rnaforge.config import Config
from rnaforge.featurecounts import compute_tpm_fpkm, run_featurecounts
from rnaforge.gates import FAIL, PASS, GateResult, raise_if_failed, write_gate_results
from rnaforge.metadata import load_metadata
from rnaforge.quality import Profile, load_profile
from rnaforge.state import RunState

MODULE_NAME = "m05_counts"
_GATE = "assignment_rate"


def build_count_gates(assignment_rates: dict[str, float],
                      profile: Profile) -> list[GateResult]:
    threshold = profile.threshold(_GATE)
    offenders = sorted(sid for sid, r in assignment_rates.items() if r < threshold)
    lowest = min(assignment_rates.values(), default=1.0)
    overridden = _GATE in profile.overrides()
    if offenders:
        status = FAIL
        message = (
            f"gene atama oranı eşiğin altında ({len(offenders)} örnek: "
            f"{', '.join(offenders)}); en düşük {lowest:.2f} < {threshold:.2f}. "
            "Düşük atama yanlış anotasyon/tür → güvenilmez sayımlar."
        )
    else:
        status = PASS
        message = f"tüm örnekler assignment ≥ {threshold:.2f} (en düşük {lowest:.2f})."
    return [GateResult(
        name=_GATE, module=MODULE_NAME, status=status, message=message,
        remedy=("Anotasyon (GFF/GTF) ile referans genomun eşleştiğini ve feature_type/"
                "attribute config'inin anotasyon formatına uyduğunu doğrulayın."),
        measured=lowest, threshold=threshold, overridden=overridden,
        samples=tuple(offenders),
    )]


def run_counts(config: Config, metadata_path: Path, run_dir: Path,
               force: bool = False) -> dict:
    run_dir = Path(run_dir)
    quant_dir = run_dir / "quantification"
    stats_dir = run_dir / "statistics"
    logs_dir = run_dir / "logs"
    for d in (quant_dir, stats_dir, logs_dir):
        d.mkdir(parents=True, exist_ok=True)
    state = RunState(run_dir)
    stats_path = stats_dir / "count_statistics.json"

    if not force and state.is_done(MODULE_NAME) and stats_path.exists():
        summary = json.loads(stats_path.read_text())
        summary["resumed"] = True
        return summary

    if not state.is_done("m04_quant"):
        raise ValueError(
            "m05 (counts) requires m04 (quant) to have completed in this run directory "
            f"first: {run_dir}. Run `rnaforge quant` with the same --run-id, then re-run counts."
        )

    profile = load_profile(config.organism_type, config.quality)
    log_path = logs_dir / "counts.log"
    with log_path.open("w") as log_file:
        def log(msg: str) -> None:
            log_file.write(msg + "\n")
            log_file.flush()

        samples = load_metadata(metadata_path)
        bams = [quant_dir / s.sample_id / "aligned.sorted.bam" for s in samples]
        paired = any(s.fastq_2 is not None for s in samples)
        log(f"m05 featureCounts: {len(samples)} sample(s), "
            f"feature_type={config.quantification.feature_type}, "
            f"attribute={config.quantification.attribute}, paired={paired}")
        result = run_featurecounts(
            bams, config.reference.annotation_gff, quant_dir / "_featurecounts",
            feature_type=config.quantification.feature_type,
            attribute=config.quantification.attribute,
            paired=paired, threads=config.resources.threads,
        )
        if not result.gene_ids:
            raise ValueError(
                "featureCounts assigned reads to no genes (empty matrix). Likely the "
                f"feature_type ({config.quantification.feature_type!r}) or attribute "
                f"({config.quantification.attribute!r}) does not match the annotation."
            )
        state.heartbeat()

        # Sütun→sample_id KONUMLA (bams örnek sırasında verildi).
        columns = list(result.counts.keys())
        sample_ids = [s.sample_id for s in samples]
        assignment_by_sample = {
            sid: result.assignment_rates[col] for sid, col in zip(sample_ids, columns)
        }

        # counts.tsv sözleşmesi: gene\t<sample_id...>
        matrix_path = quant_dir / "counts.tsv"
        with matrix_path.open("w") as fh:
            fh.write("gene\t" + "\t".join(sample_ids) + "\n")
            for i, gene in enumerate(result.gene_ids):
                row = [str(result.counts[col][i]) for col in columns]
                fh.write(gene + "\t" + "\t".join(row) + "\n")
        log(f"count matrix written: {matrix_path} ({len(result.gene_ids)} genes)")

        # TPM / FPKM ekspresyon değerleri (gen uzunluğuyla normalize) — deliverable + rapor "Ekspresyon Düzeyi"
        _genes, _cols, tpm, fpkm = compute_tpm_fpkm((quant_dir / "_featurecounts" / "counts.txt").read_text())
        col2sid = dict(zip(_cols, sample_ids))
        for name, mat in (("tpm.tsv", tpm), ("fpkm.tsv", fpkm)):
            with (quant_dir / name).open("w") as fh:
                fh.write("gene\t" + "\t".join(sample_ids) + "\n")
                for i, gene in enumerate(_genes):
                    fh.write(gene + "\t" + "\t".join(f'{mat[c][i]:g}' for c in _cols) + "\n")
        log(f"expression matrices written: tpm.tsv, fpkm.tsv")

        gates = build_count_gates(assignment_by_sample, profile)
        summary = {
            "n_samples": len(samples), "n_genes": len(result.gene_ids),
            "samples": {sid: {"assignment_rate": assignment_by_sample[sid]} for sid in sample_ids},
            "gate_counts": dict(Counter(g.status for g in gates)),
            "expression_values": ["tpm.tsv", "fpkm.tsv"],
        }
        stats_path.write_text(json.dumps(summary, indent=2))
        write_gate_results(run_dir, gates)
        for g in gates:
            log(f"gate {g.name}: {g.status} — {g.message}")
        raise_if_failed(gates)

    state.mark_done(MODULE_NAME, [str(stats_path), str(log_path)])
    return summary
