"""m06 — Differential Expression (DESeq2).

m05 count matrisini R/Bioconductor DESeq2 ile diferansiyel ekspresyona çevirir —
pipeline'ın biyolojik çıktısı. İlk ORTAK (organizma-agnostik) analiz adımı.
Veri kapısı `replicate_correlation`: koşul-içi replikalar zayıf korele ise WARN
(sonuç ŞÜPHELİ damgalanır ama ÜRETİLİR — düşük korelasyon DE'yi geçersiz kılmaz,
gücü düşürür). m06 asla FAIL üretmez."""
from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

from rnaforge.config import Config
from rnaforge.deseq2 import run_deseq2
from rnaforge.gates import PASS, WARN, GateResult, write_gate_results
from rnaforge.metadata import load_metadata
from rnaforge.quality import Profile, load_profile
from rnaforge.state import RunState

MODULE_NAME = "m06_de"
_GATE = "replicate_correlation"


def _write_coldata(samples, path: Path) -> None:
    has_batch = any(s.batch for s in samples)
    with path.open("w") as f:
        header = "sample\tcondition" + ("\tbatch" if has_batch else "")
        f.write(header + "\n")
        for s in samples:
            row = f"{s.sample_id}\t{s.condition}"
            if has_batch:
                row += f"\t{s.batch or 'NA'}"
            f.write(row + "\n")


def build_de_gates(min_correlation: float, profile: Profile) -> list[GateResult]:
    threshold = profile.threshold(_GATE)
    overridden = _GATE in profile.overrides()
    if min_correlation < threshold:
        status = WARN
        message = (
            f"koşul-içi replika korelasyonu düşük (min {min_correlation:.2f} < "
            f"{threshold:.2f}). DE üretildi ama ŞÜPHELİ: replikalar zayıf kümeleniyor "
            "(olası aykırı örnek / batch etkisi)."
        )
    else:
        status = PASS
        message = f"replika korelasyonu yeterli (min {min_correlation:.2f} ≥ {threshold:.2f})."
    return [GateResult(
        name=_GATE, module=MODULE_NAME, status=status, message=message,
        remedy=("PCA/heatmap ile aykırı örnek arayın; batch/covariate varsa design formülüne "
                "ekleyin (`~batch + condition`). Düşük korelasyon DE gücünü düşürür."),
        measured=min_correlation, threshold=threshold, overridden=overridden,
    )]


def run_de(config: Config, metadata_path: Path, run_dir: Path,
           force: bool = False) -> dict:
    run_dir = Path(run_dir)
    de_dir = run_dir / "differential_expression"
    stats_dir = run_dir / "statistics"
    logs_dir = run_dir / "logs"
    for d in (de_dir, stats_dir, logs_dir):
        d.mkdir(parents=True, exist_ok=True)
    state = RunState(run_dir)
    stats_path = stats_dir / "de_statistics.json"

    if not force and state.is_done(MODULE_NAME) and stats_path.exists():
        summary = json.loads(stats_path.read_text())
        summary["resumed"] = True
        return summary

    if not state.is_done("m05_counts"):
        raise ValueError(
            "m06 (de) requires m05 (counts) to have completed in this run directory "
            f"first: {run_dir}. Run `rnaforge counts` with the same --run-id, then re-run de."
        )

    profile = load_profile(config.organism_type, config.quality)
    log_path = logs_dir / "de.log"
    with log_path.open("w") as log_file:
        def log(msg: str) -> None:
            log_file.write(msg + "\n")
            log_file.flush()

        samples = load_metadata(metadata_path)
        coldata_path = de_dir / "coldata.tsv"
        _write_coldata(samples, coldata_path)
        counts_tsv = run_dir / "quantification" / "counts.tsv"
        log(f"m06 DESeq2: design={config.de.design!r} reference={config.de.reference!r}")
        state.heartbeat()
        result = run_deseq2(counts_tsv, coldata_path, config.de.design, de_dir,
                            reference=config.de.reference)

        fdr = config.de.fdr_threshold
        lfc = config.de.log2fc_threshold
        n_sig = sum(
            1 for r in result.results
            if r.get("padj") is not None and r["padj"] < fdr
            and r.get("log2FoldChange") is not None and abs(r["log2FoldChange"]) >= lfc
        )
        min_corr = float(result.metrics.get("min_replicate_correlation", 1.0))
        gates = build_de_gates(min_corr, profile)

        summary = {
            "n_genes": len(result.results),
            "n_significant": n_sig,
            "contrast": result.metrics.get("contrast", ""),
            "min_replicate_correlation": min_corr,
            "fdr_threshold": fdr, "log2fc_threshold": lfc,
            "gate_counts": dict(Counter(g.status for g in gates)),
        }
        stats_path.write_text(json.dumps(summary, indent=2))
        write_gate_results(run_dir, gates)
        for g in gates:
            log(f"gate {g.name}: {g.status} — {g.message}")
        log(f"DE done: {n_sig} significant / {len(result.results)} genes, "
            f"contrast={summary['contrast']}")

    state.mark_done(MODULE_NAME, [str(stats_path), str(log_path)])
    return summary
