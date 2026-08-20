"""m08 — HTML report. m06/m07 ciktilarindan tek self-contained report.html.
Organizma-agnostik; YENI veri-kapisi YOK (rapor biyolojiyi gecersiz kilmaz)."""
from __future__ import annotations
import json
from pathlib import Path

from rnaforge import __version__
from rnaforge.config import Config
from rnaforge.report_html import load_report_inputs, render_report
from rnaforge.state import RunState
from rnaforge.versions import capture_tool_versions

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
        # F2.2: yazılım tablosu için gerçek kurulu sürümleri yakala (best-effort; conda
        # yoksa boş → rapor curated fallback kullanır, çökmez).
        try:
            inputs["software_versions"] = capture_tool_versions()
        except Exception as exc:  # asla rapor üretimini bloklama
            log_file.write(f"software version capture failed (curated fallback): {exc}\n")
            inputs["software_versions"] = {}
        state.heartbeat()
        doc = render_report(inputs, config, version=__version__, run_id=run_dir.name)
        report_path = report_dir / "report.html"
        report_path.write_text(doc)
        summary = {
            "report": "report/report.html",
            "language": config.report.language,
            # Render edilen bölüm sayısını doc'tan CANLI say: opsiyonel bölümler daima <section>
            # olarak gelir, taksonomi yalnız metatranskriptomda → non-meta 17, meta 18 (doğru).
            "n_sections": doc.count("<section"),
        }
        stats_path.write_text(json.dumps(summary, indent=2))
        log_file.write(f"m08 report done: {report_path} ({len(doc)} bytes)\n")

    state.mark_done(MODULE_NAME, [str(report_dir / "report.html"), str(stats_path), str(log_path)])
    return summary
