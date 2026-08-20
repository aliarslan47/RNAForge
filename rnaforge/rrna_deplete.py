"""Metatranskriptom rRNA depletion: SortMeRNA ile rRNA okumalarını ÇIKAR (--other).

seqqc.run_sortmerna yalnız ÖLÇER; bu modül rRNA'sız FASTQ üretir (downstream girdisi).
"""
from __future__ import annotations

import re
import shutil
import subprocess
from pathlib import Path


def parse_depletion_rate(log_path: Path) -> float:
    """aligned.log'tan rRNA fraksiyonu (0-1) = çıkarılan pay. seqqc.parse_sortmerna_log deseni."""
    text = Path(log_path).read_text() if Path(log_path).exists() else ""
    m = re.search(r"passing E-value threshold\s*=\s*\d+\s*\(([\d.]+)%?\)", text)
    if m:
        return float(m.group(1)) / 100.0
    total = re.search(r"Total reads\s*=\s*(\d+)", text)
    passed = re.search(r"passing E-value threshold\s*=\s*(\d+)", text)
    if total and passed and int(total.group(1)) > 0:
        return int(passed.group(1)) / int(total.group(1))
    return 0.0


def run_sortmerna_deplete(reads: list[Path], rrna_db: Path, workdir: Path, paired: bool,
                          threads: int = 8, env: str = "rnaforge-seqqc") -> dict:
    """SortMeRNA --fastx --other ile rRNA'sız okumaları üret. Hatada gürültülü yüksel."""
    workdir = Path(workdir)
    if workdir.exists():
        shutil.rmtree(workdir)
    workdir.mkdir(parents=True)
    aligned = workdir / "out" / "aligned"
    other = workdir / "out" / "other"
    cmd = ["conda", "run", "-n", env, "sortmerna", "--ref", str(rrna_db),
           "--workdir", str(workdir), "--fastx",
           "--aligned", str(aligned), "--other", str(other),
           "--threads", str(threads), "-v"]
    if paired:
        cmd += ["--paired_in", "--out2"]
    for r in reads:
        cmd += ["--reads", str(r)]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError(f"sortmerna deplete failed (exit {res.returncode}):\n{res.stderr[-2000:]}")
    out_dir = workdir / "out"
    others = sorted(out_dir.glob("other*.f*q*"))
    if not others:
        raise RuntimeError(f"sortmerna produced no --other output in {out_dir}")
    log = out_dir / "aligned.log"
    return {"other": others, "depletion_rate": parse_depletion_rate(log), "aligned_log": log}
