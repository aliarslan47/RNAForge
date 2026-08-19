"""NanoCount: ONT transkript-düzeyi (izoform) niceleme — çalıştırır ve çıktısını parse eder.

Uzun-okuma bir genin izoformları arasında belirsiz eşleşir; primer-hizalama sayımı bunu
ayrıştıramaz. NanoCount, transkriptoma çok-hizalanmış (minimap2 `-N`) BAM üzerinde EM ile
transkript başına beklenen okuma sayısını (`est_count`, kesirli) tahmin eder. Parser saftır
(string girer, veri çıkar); bowtie2/minimap2 deseni."""
from __future__ import annotations

import subprocess
from pathlib import Path


class NanoCountRunError(RuntimeError):
    """NanoCount çalıştırılamadı ya da beklenen çıktıyı üretmedi."""


def parse_nanocount(tsv_text: str) -> dict[str, float]:
    """NanoCount çıktısı (transcript_name, raw, est_count, tpm) → {transcript: est_count}.

    est_count = EM ile tahmin edilen (kesirli) okuma sayısı; DE için yuvarlanır (çağıran)."""
    lines = [ln for ln in tsv_text.splitlines() if ln.strip()]
    if not lines:
        return {}
    header = lines[0].split("\t")
    try:
        name_i = header.index("transcript_name")
        est_i = header.index("est_count")
    except ValueError:
        raise NanoCountRunError(
            "NanoCount output missing 'transcript_name'/'est_count' columns; "
            f"got header: {header}"
        ) from None
    out: dict[str, float] = {}
    for line in lines[1:]:
        f = line.split("\t")
        if len(f) <= max(name_i, est_i):
            continue
        out[f[name_i]] = float(f[est_i])
    return out


def run_nanocount(bam_path: Path, out_tsv: Path,
                  env: str = "rnaforge-longread") -> dict[str, float]:
    """NanoCount'u çok-hizalanmış BAM üzerinde çalıştır; {transcript: est_count} döndür.
    Nonzero exit → yüksek sesle hata (sessiz kısmi çıktı yok)."""
    out_tsv = Path(out_tsv)
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    cmd = ["conda", "run", "-n", env, "NanoCount",
           "-i", str(bam_path), "-o", str(out_tsv)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0 or not out_tsv.exists():
        raise NanoCountRunError(
            f"NanoCount failed (exit {r.returncode})\ncmd: {' '.join(cmd)}\n"
            f"stderr: {r.stderr.strip()[-800:]}"
        )
    return parse_nanocount(out_tsv.read_text())
