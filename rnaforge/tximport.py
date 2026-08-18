"""tximport (transkript→gen) — R Bioconductor wrapper (deseq2.py deseni).

countsFromAbundance="lengthScaledTPM" → uzunluk-düzeltilmiş sayım; m06 DESeq2 bunu düz
sayım gibi okur ve doğru olur (m06 organizma-agnostik kalır). Env: rnaforge-de."""
from __future__ import annotations

import subprocess
from dataclasses import dataclass
from pathlib import Path

_SCRIPT = Path(__file__).parent / "scripts" / "tximport.R"


class TximportError(RuntimeError):
    """tximport çalıştırılamadı ya da beklenen çıktıyı üretmedi."""


@dataclass
class TximportResult:
    gene_ids: list[str]
    counts: dict[str, list[float]]


def run_tximport(quant_sfs: dict[str, Path], tx2gene: Path, out_dir: Path,
                 env: str = "rnaforge-de") -> TximportResult:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_tsv = out_dir / "gene_counts.tsv"
    args = [str(tx2gene), str(out_tsv)]
    for sid, sf in quant_sfs.items():
        args += [sid, str(sf)]
    cmd = ["conda", "run", "-n", env, "Rscript", str(_SCRIPT), *args]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0 or not out_tsv.exists():
        raise TximportError(f"tximport failed: {r.stderr[-500:]}")
    lines = out_tsv.read_text().splitlines()
    header = lines[0].split("\t")[1:]     # sample ids in file order
    gene_ids: list[str] = []
    cols: dict[str, list[float]] = {sid: [] for sid in header}
    for line in lines[1:]:
        parts = line.split("\t")
        gene_ids.append(parts[0])
        for sid, val in zip(header, parts[1:]):
            cols[sid].append(float(val))
    return TximportResult(gene_ids=gene_ids, counts=cols)
