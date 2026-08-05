"""m18 — MultiQC toplu görünüm (diagnostik/kapstone). Run dizinini tarar ve tüm
araç çıktılarını (FastQC, fastp, featureCounts, samtools stats, RSeQC) tek toplu
HTML'de birleştirir. KAPI ÜRETMEZ. En son çalışmalı ki tüm çıktıları toplasın."""
from __future__ import annotations

import json
from pathlib import Path

from rnaforge.config import Config
from rnaforge.multiqc import count_modules, parse_general_stats
from rnaforge.multiqc import run_multiqc as run_multiqc_tool
from rnaforge.state import RunState

MODULE_NAME = "m18_multiqc"


def run_multiqc(config: Config, metadata_path: Path, run_dir: Path,
                force: bool = False) -> dict:
    run_dir = Path(run_dir)
    out_dir = run_dir / "multiqc"
    stats_dir = run_dir / "statistics"
    logs_dir = run_dir / "logs"
    for d in (out_dir, stats_dir, logs_dir):
        d.mkdir(parents=True, exist_ok=True)
    state = RunState(run_dir)
    stats_path = stats_dir / "multiqc_statistics.json"

    if not force and state.is_done(MODULE_NAME) and stats_path.exists():
        summary = json.loads(stats_path.read_text()); summary["resumed"] = True
        return summary
    if not state.is_done("m04_quant"):
        raise ValueError(
            "m18 (multiqc) requires m04 (quant) to have completed in this run directory "
            f"first: {run_dir}. Run the pipeline with the same --run-id, then re-run multiqc.")

    log_path = logs_dir / "multiqc.log"
    with log_path.open("w") as log_file:
        # run_dir taranır, out_dir'e yazılır; MultiQC -o hedefini otomatik atlar.
        report = run_multiqc_tool(run_dir, out_dir)
        data_dir = out_dir / "multiqc_report_data"
        gs_path = data_dir / "multiqc_general_stats.txt"
        n_samples = len(parse_general_stats(gs_path.read_text())) if gs_path.exists() else 0
        n_modules = count_modules(data_dir)
        # Rapor report/ altında; multiqc run_dir/multiqc altında → göreli link.
        summary = {
            "report_relpath": f"../multiqc/{report.name}", "report_name": report.name,
            "n_modules": n_modules, "n_samples": n_samples,
        }
        stats_path.write_text(json.dumps(summary, indent=2))
        log_file.write(f"m18 multiqc done: {summary}\n")

    state.mark_done(MODULE_NAME, [str(stats_path), str(log_path)])
    return summary
