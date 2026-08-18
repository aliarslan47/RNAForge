"""Salmon (ökaryot transkriptom niceleme) — saf parser + runner (bowtie2.py deseni).

Decoy-aware selective alignment: genome_fasta verilirse genom decoy olarak indekslenir
(anotasyonsuz/intergenik sahte eşleşmeleri eler; Salmon önerisi). Env: rnaforge-quant-euk."""
from __future__ import annotations

import json
import subprocess
from dataclasses import dataclass
from pathlib import Path


class SalmonError(RuntimeError):
    """Salmon çalıştırılamadı ya da beklenen çıktıyı üretmedi."""


@dataclass
class SalmonQuant:
    quant_sf: Path
    mapping_rate: float


def parse_salmon_meta(meta_info_json: Path) -> float:
    data = json.loads(Path(meta_info_json).read_text())
    return float(data["percent_mapped"]) / 100.0


def _contig_names(fasta: Path) -> list[str]:
    names = []
    with Path(fasta).open() as fh:
        for line in fh:
            if line.startswith(">"):
                names.append(line[1:].split()[0])
    return names


def build_salmon_index(transcriptome_fasta: Path, index_dir: Path,
                       genome_fasta: Path | None = None, threads: int = 8,
                       log=None) -> Path:
    index_dir = Path(index_dir)
    index_dir.parent.mkdir(parents=True, exist_ok=True)
    if genome_fasta is not None:
        decoys = index_dir.parent / "decoys.txt"
        decoys.write_text("\n".join(_contig_names(genome_fasta)) + "\n")
        gentrome = index_dir.parent / "gentrome.fa"
        with gentrome.open("wb") as out:
            for src in (transcriptome_fasta, genome_fasta):
                out.write(Path(src).read_bytes())
        cmd = ["conda", "run", "-n", "rnaforge-quant-euk", "salmon", "index",
               "-t", str(gentrome), "-d", str(decoys), "-i", str(index_dir),
               "-k", "31", "-p", str(threads)]
    else:
        if log:
            log("salmon index: genome_fasta verilmedi → transkriptom-only "
                "(doğruluk için genome_fasta decoy önerilir)")
        cmd = ["conda", "run", "-n", "rnaforge-quant-euk", "salmon", "index",
               "-t", str(transcriptome_fasta), "-i", str(index_dir),
               "-k", "31", "-p", str(threads)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0 or not index_dir.exists():
        raise SalmonError(f"salmon index failed: {r.stderr[-500:]}")
    return index_dir


def run_salmon_quant(index_dir: Path, out_dir: Path, fastq_1: Path,
                     fastq_2: Path | None = None, threads: int = 8,
                     env: str = "rnaforge-quant-euk") -> SalmonQuant:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = ["conda", "run", "-n", env, "salmon", "quant", "-i", str(index_dir),
           "-l", "A", "-p", str(threads), "--validateMappings", "-o", str(out_dir)]
    if fastq_2 is not None:
        cmd += ["-1", str(fastq_1), "-2", str(fastq_2)]
    else:
        cmd += ["-r", str(fastq_1)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    quant_sf = out_dir / "quant.sf"
    meta = out_dir / "aux_info" / "meta_info.json"
    if r.returncode != 0 or not quant_sf.exists() or not meta.exists():
        raise SalmonError(f"salmon quant failed for {fastq_1.name}: {r.stderr[-500:]}")
    return SalmonQuant(quant_sf=quant_sf, mapping_rate=parse_salmon_meta(meta))
