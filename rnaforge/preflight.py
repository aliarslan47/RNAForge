"""Kurulum ön-uçuş: gerekli conda env'leri mevcut mu?

Araç sarmalayıcıları `conda run -n <env> ...` çağırır; env yoksa hata koşunun
ORTASINDA derin patlar. `rnaforge doctor` bunu ÖNDEN, actionable mesajla yakalar:
hangi env eksik + nasıl kurulur (envs/<env>.yml). Kurulumun tek elle-adımı 9 env
oluşturmak; bu kontrol onları unutmayı ucuza yakalar."""
from __future__ import annotations

import subprocess

# env adı → tek satır amaç (envs/<name>.yml ile bire bir).
ENVIRONMENTS = {
    "rnaforge-core": "orchestration + networkx/scipy (PPI)",
    "rnaforge-qc": "FastQC + fastp (short-read QC/trim)",
    "rnaforge-quant-prok": "Bowtie2 + samtools + subread (prokaryote quant)",
    "rnaforge-quant-euk": "Salmon (eukaryote quant)",
    "rnaforge-longread": "minimap2 + NanoPlot + Pychopper + chopper (long-read)",
    "rnaforge-basecall": "pod5 (+ harici dorado) — m00 basecalling",
    "rnaforge-de": "R/DESeq2 + fgsea + tximport (DE + downstream)",
    "rnaforge-amr": "abricate (AMR/virulence, m13)",
    "rnaforge-seqqc": "SortMeRNA + RSeQC + MultiQC (m16/m18)",
}


def list_conda_envs() -> set[str]:
    """`conda env list` çıktısındaki env ADLARI kümesi. conda yoksa boş küme."""
    try:
        proc = subprocess.run(["conda", "env", "list"], capture_output=True, text=True)
    except FileNotFoundError:
        return set()
    names: set[str] = set()
    for line in proc.stdout.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        names.add(line.split()[0])
    return names


def environment_report(available: set[str]) -> list[tuple[str, str, bool]]:
    """(env_adı, açıklama, mevcut_mu) — sabit ENVIRONMENTS sırasında (saf)."""
    return [(name, desc, name in available) for name, desc in ENVIRONMENTS.items()]
