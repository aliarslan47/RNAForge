"""Hizalama-sonrası QC: samtools stats (insert-size), samtools coverage
(kontig derinliği) ve RSeQC read_distribution parser+runner'ları. Parser'lar
saftır (string girer, dict çıkar) — I/O yok, hızlı ve deterministik test edilir."""
from __future__ import annotations

import subprocess
from pathlib import Path


class AlignQCParseError(ValueError):
    """Araç çıktısı beklenen biçimde değil."""


class AlignQCRunError(RuntimeError):
    """Araç çalıştırılamadı ya da beklenen çıktıyı üretmedi."""


# --- F2: samtools stats (insert-size) ---

def parse_samtools_stats(text: str) -> dict:
    """`samtools stats` çıktısından SN metrikleri + IS histogramı.
    Döner: {"reads_paired": int, "insert_size_average": float,
            "insert_size_sd": float, "histogram": [(insert_size, pair_count), ...]}."""
    sn: dict[str, float] = {}
    hist: list[tuple[int, int]] = []
    for line in text.splitlines():
        if line.startswith("SN\t"):
            parts = line.split("\t")
            if len(parts) >= 3:
                key = parts[1].rstrip(":").strip()
                val = parts[2].strip()
                try:
                    sn[key] = float(val)
                except ValueError:
                    pass
        elif line.startswith("IS\t"):
            parts = line.split("\t")
            # IS <insert_size> <all_pairs> <inward> <outward> <other>
            if len(parts) >= 3:
                try:
                    hist.append((int(parts[1]), int(parts[2])))
                except ValueError:
                    pass
    return {
        "reads_paired": int(sn.get("reads paired", 0)),
        "insert_size_average": sn.get("insert size average", 0.0),
        "insert_size_sd": sn.get("insert size standard deviation", 0.0),
        "histogram": hist,
    }


def run_samtools_stats(bam: Path, env: str = "rnaforge-quant-prok") -> str:
    cmd = ["conda", "run", "-n", env, "samtools", "stats", str(bam)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise AlignQCRunError(
            f"samtools stats failed (exit {r.returncode}) for {bam}\nstderr: {r.stderr.strip()}")
    return r.stdout


# --- F4: samtools coverage (kontig başına derinlik) ---

def parse_samtools_coverage(text: str) -> list[dict]:
    """`samtools coverage` çıktısından kontig başına satırlar.
    Döner: [{"contig": str, "numreads": int, "coverage": float, "meandepth": float}, ...]."""
    rows: list[dict] = []
    header: list[str] | None = None
    for line in text.splitlines():
        if line.startswith("#rname"):
            header = line.lstrip("#").split("\t")
            continue
        if not line.strip() or header is None:
            continue
        fields = line.split("\t")
        rec = dict(zip(header, fields))
        try:
            rows.append({
                "contig": rec["rname"],
                "numreads": int(rec.get("numreads", 0)),
                "coverage": float(rec.get("coverage", 0.0)),
                "meandepth": float(rec.get("meandepth", 0.0)),
            })
        except (KeyError, ValueError):
            continue
    return rows


def run_samtools_coverage(bam: Path, env: str = "rnaforge-quant-prok") -> str:
    cmd = ["conda", "run", "-n", env, "samtools", "coverage", str(bam)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise AlignQCRunError(
            f"samtools coverage failed (exit {r.returncode}) for {bam}\nstderr: {r.stderr.strip()}")
    return r.stdout


# --- F3: RSeQC read_distribution ---

# Prokaryot BED CDS-only olduğundan yalnız çakışmayan gruplar anlamlı; iç içe
# (kümülatif) TSS/TES pencereleri çift sayımı önlemek için % tablosundan dışlanır.
_RD_EXCLUSIVE = ("CDS_Exons", "5'UTR_Exons", "3'UTR_Exons", "Introns")


def parse_read_distribution(text: str) -> dict:
    """RSeQC read_distribution.py çıktısını parse eder. Çakışmayan gruplar için
    (CDS/UTR/Intron) Total Tags'e oranı + 'Intergenic' = atanmamış tag oranı.
    Döner: {"total_tags": int, "assigned_tags": int, "groups": {ad: tag_count},
            "percentages": {ad: yüzde}}."""
    total_tags = 0
    assigned_tags = 0
    groups: dict[str, int] = {}
    for line in text.splitlines():
        s = line.strip()
        if s.startswith("Total Tags"):
            total_tags = int(s.split()[-1])
        elif s.startswith("Total Assigned Tags"):
            assigned_tags = int(s.split()[-1])
        else:
            parts = s.split()
            if len(parts) == 4 and parts[0] in _RD_EXCLUSIVE:
                try:
                    groups[parts[0]] = int(parts[2])   # Group Total_bases Tag_count Tags/Kb
                except ValueError:
                    pass
    if total_tags <= 0:
        raise AlignQCParseError("read_distribution: 'Total Tags' bulunamadı veya sıfır")
    pct = {g: round(100.0 * groups.get(g, 0) / total_tags, 3) for g in _RD_EXCLUSIVE}
    pct["Intergenic"] = round(100.0 * (total_tags - assigned_tags) / total_tags, 3)
    return {"total_tags": total_tags, "assigned_tags": assigned_tags,
            "groups": groups, "percentages": pct}


def run_read_distribution(bam: Path, bed: Path, env: str = "rnaforge-seqqc") -> str:
    cmd = ["conda", "run", "-n", env, "read_distribution.py", "-i", str(bam), "-r", str(bed)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise AlignQCRunError(
            f"read_distribution.py failed (exit {r.returncode}) for {bam}\nstderr: {r.stderr.strip()}")
    return r.stdout


def aggregate_histograms(histograms: list[list[tuple[int, int]]],
                         max_insert: int = 1000) -> tuple[list[str], list[int]]:
    """Birden çok örneğin IS histogramını insert-size'a göre toplar; max_insert
    üstünü keser. Döner: (etiketler, toplam_çiftler)."""
    totals: dict[int, int] = {}
    for hist in histograms:
        for size, count in hist:
            if size <= max_insert:
                totals[size] = totals.get(size, 0) + count
    sizes = sorted(totals)
    return [str(s) for s in sizes], [totals[s] for s in sizes]
