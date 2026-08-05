"""MultiQC toplu görünüm: bir run dizinini tarar (FastQC, fastp, featureCounts,
samtools stats, RSeQC…) ve tek toplu HTML üretir. Parser saftır."""
from __future__ import annotations

import subprocess
from pathlib import Path


class MultiQCRunError(RuntimeError):
    """MultiQC çalıştırılamadı ya da beklenen çıktıyı üretmedi."""


def parse_general_stats(text: str) -> list[dict]:
    """multiqc_general_stats.txt (TSV) → satır başına {sütun: değer}. İlk sütun
    örnek adı. Boş/başlıksızsa boş liste."""
    lines = [ln for ln in text.splitlines() if ln.strip()]
    if not lines:
        return []
    header = lines[0].split("\t")
    rows: list[dict] = []
    for line in lines[1:]:
        fields = line.split("\t")
        rows.append(dict(zip(header, fields)))
    return rows


def count_modules(data_dir: Path) -> int:
    """MultiQC data dizinindeki modül sayısı (multiqc_<modul>.txt dosyaları,
    general_stats ve sources hariç)."""
    data_dir = Path(data_dir)
    if not data_dir.exists():
        return 0
    n = 0
    for p in data_dir.glob("multiqc_*.txt"):
        if p.stem in ("multiqc_general_stats", "multiqc_sources", "multiqc_citations"):
            continue
        n += 1
    return n


def run_multiqc(scan_dir: Path, out_dir: Path, report_name: str = "multiqc_report",
                env: str = "rnaforge-seqqc") -> Path:
    """MultiQC'yi scan_dir üzerinde çalıştırır; out_dir/<report_name>.html döner."""
    scan_dir = Path(scan_dir)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = ["conda", "run", "-n", env, "multiqc", str(scan_dir),
           "-o", str(out_dir), "-n", report_name, "-f", "-q"]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise MultiQCRunError(
            f"multiqc failed (exit {r.returncode})\ncmd: {' '.join(cmd)}\nstderr: {r.stderr.strip()}")
    report = out_dir / f"{report_name}.html"
    if not report.exists():
        raise MultiQCRunError(
            f"multiqc reported success but no report at {report} "
            f"(out_dir: {[p.name for p in out_dir.iterdir()]})")
    return report
