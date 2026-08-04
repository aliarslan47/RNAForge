"""m08 — HTML report. m06/m07 ciktilarindan tek self-contained report.html.
Organizma-agnostik; YENI veri-kapisi YOK (rapor biyolojiyi gecersiz kilmaz)."""
from __future__ import annotations
import json
from pathlib import Path

from rnaforge import __version__
from rnaforge.config import Config
from rnaforge.report_html import load_report_inputs, render_report, N_SECTIONS
from rnaforge.state import RunState

MODULE_NAME = "m08_report"


def run_report(config: Config, metadata_path: Path, run_dir: Path,
               force: bool = False) -> dict:
    run_dir = Path(run_dir)
    report_dir = run_dir / "report"
    stats_dir = run_dir / "statistics"
    logs_dir = run_dir / "logs"
    for d in (report_dir, stats_dir, logs_dir):
        d.mkdir(parents=True, exist_ok=True)
    state = RunState(run_dir)
    stats_path = stats_dir / "report_statistics.json"

    if not force and state.is_done(MODULE_NAME) and stats_path.exists():
        summary = json.loads(stats_path.read_text()); summary["resumed"] = True
        return summary
    if not state.is_done("m07_figures"):
        raise ValueError(
            "m08 (report) requires m07 (figures) to have completed in this run directory "
            f"first: {run_dir}. Run `rnaforge figures` with the same --run-id, then re-run report.")

    log_path = logs_dir / "report.log"
    with log_path.open("w") as log_file:
        inputs = load_report_inputs(run_dir)
        state.heartbeat()
        doc = render_report(inputs, config, version=__version__)
        report_path = report_dir / "report.html"
        report_path.write_text(doc)
        summary = {
            "report": "report/report.html",
            "language": config.report.language,
            "n_sections": N_SECTIONS,
        }
        stats_path.write_text(json.dumps(summary, indent=2))
        log_file.write(f"m08 report done: {report_path} ({len(doc)} bytes)\n")

    state.mark_done(MODULE_NAME, [str(report_dir / "report.html"), str(stats_path), str(log_path)])
    return summary
