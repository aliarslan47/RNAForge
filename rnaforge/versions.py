"""Kurulu araç sürümlerini conda ortamlarından runtime'da okur (tekrarlanabilirlik).

Rapordaki yazılım tablosu eskiden curated (elle yazılmış) sürümler kullanıyordu; bunlar
`envs/*.yml` kesin-pin'leriyle eşleşiyordu ama koşuda GERÇEKTE kurulu olanı yansıtmıyordu.
Burası `conda list --json` (yapısal çıktı — kırılgan `--version` scraping DEĞİL) ile gerçek
sürümleri okur. **Best-effort:** conda ya da env yoksa boş döner → çağıran curated fallback'e
düşer, rapor asla çökmez ([[feedback_gurultulu_hata]] anlamında: sessiz değil, dürüst fallback)."""
from __future__ import annotations

import json
import shutil
import subprocess

# Araç etiketi (report_html._SOFTWARE ilk sütunuyla EŞLEŞMELİ) → (conda env, conda paket adı).
_VERSION_SOURCES: dict[str, tuple[str, str]] = {
    "Python": ("rnaforge-core", "python"),
    "FastQC": ("rnaforge-qc", "fastqc"),
    "fastp": ("rnaforge-qc", "fastp"),
    "Bowtie2": ("rnaforge-quant-prok", "bowtie2"),
    "NanoPlot": ("rnaforge-longread", "nanoplot"),
    "Pychopper": ("rnaforge-longread", "pychopper"),
    "chopper": ("rnaforge-longread", "chopper"),
    "minimap2": ("rnaforge-longread", "minimap2"),
    "Subread/featureCounts": ("rnaforge-quant-prok", "subread"),
    "DESeq2 (R)": ("rnaforge-de", "bioconductor-deseq2"),
    "ggplot2 (R)": ("rnaforge-de", "r-ggplot2"),
    "fgsea (R)": ("rnaforge-de", "bioconductor-fgsea"),
    "abricate": ("rnaforge-amr", "abricate"),
    "SortMeRNA": ("rnaforge-seqqc", "sortmerna"),
    "RSeQC": ("rnaforge-seqqc", "rseqc"),
    "networkx": ("rnaforge-core", "networkx"),
}


def parse_conda_list_json(text: str) -> dict[str, str]:
    """`conda list --json` çıktısını {paket_adı: sürüm} dict'ine çevirir. Bozuk → {}."""
    try:
        data = json.loads(text)
    except (json.JSONDecodeError, ValueError, TypeError):
        return {}
    if not isinstance(data, list):
        return {}
    out: dict[str, str] = {}
    for pkg in data:
        if isinstance(pkg, dict):
            name, ver = pkg.get("name"), pkg.get("version")
            if name and ver:
                out[str(name)] = str(ver)
    return out


def _conda_list(env: str) -> dict[str, str]:
    """Bir conda env'inin kurulu paketlerini oku. Hata → {} (best-effort)."""
    try:
        r = subprocess.run(
            ["conda", "list", "--json", "-n", env],
            capture_output=True, text=True,
        )
    except OSError:
        return {}
    if r.returncode != 0:
        return {}
    return parse_conda_list_json(r.stdout)


def capture_tool_versions(sources: dict[str, tuple[str, str]] | None = None) -> dict[str, str]:
    """Kurulu araç sürümlerini env'lerden topla (env başına TEK sorgu, cache'li).

    Dönen dict {araç_etiketi: sürüm}; kurulu olmayan/okunamayan araç ATLANIR (çağıran o
    araç için curated fallback kullanır). conda yoksa {} → tamamen curated'a düşer."""
    sources = sources if sources is not None else _VERSION_SOURCES
    if shutil.which("conda") is None:
        return {}
    env_cache: dict[str, dict[str, str]] = {}
    out: dict[str, str] = {}
    for tool, (env, pkg) in sources.items():
        if env not in env_cache:
            env_cache[env] = _conda_list(env)
        ver = env_cache[env].get(pkg)
        if ver:
            out[tool] = ver
    return out
