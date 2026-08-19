"""minimap2 long-read hizalama: çalıştırır ve samtools flagstat özetini parse eder.

Parser saftır. bowtie2.py deseni izlenir; fark: minimap2 stderr'de "overall alignment
rate" satırı yazmaz → hizalama oranı `samtools flagstat`'ın primary-mapped/primary
sayımından hesaplanır (yazdırılan %'e değil, sayımlara bakılır — biçim kaymasına dayanıklı)."""
from __future__ import annotations

import re
import subprocess
from dataclasses import dataclass
from pathlib import Path

# "N + M primary" (satır sonu) ve "N + M primary mapped ..." satırları.
_PRIMARY_RE = re.compile(r"^(\d+)\s+\+\s+\d+\s+primary\s*$", re.MULTILINE)
_PRIMARY_MAPPED_RE = re.compile(r"^(\d+)\s+\+\s+\d+\s+primary mapped\b", re.MULTILINE)

_PRESETS = {
    "ont": "map-ont",
    "pacbio_hifi": "map-hifi",
}


class Minimap2ParseError(ValueError):
    """samtools flagstat çıktısı beklenen biçimde değil."""


class Minimap2RunError(RuntimeError):
    """minimap2/samtools çalıştırılamadı ya da beklenen çıktıyı üretmedi."""


@dataclass(frozen=True)
class AlignmentResult:
    bam: Path
    alignment_rate: float


def minimap2_preset(platform: str) -> str:
    """Tespit edilen platform -> minimap2 preset. Kısa okuma buraya ulaşmaz."""
    try:
        return _PRESETS[platform]
    except KeyError:
        raise Minimap2RunError(
            f"no minimap2 preset for platform {platform!r}; "
            f"long-read platforms: {', '.join(_PRESETS)}"
        ) from None


def parse_flagstat_mapped(text: str) -> float:
    """primary-mapped / primary oranı (0..1). primary satırları yoksa ya da
    primary=0 ise (boş girdi) yüksek sesle hata (sessiz 0.0 döndürmez)."""
    primary_m = _PRIMARY_RE.search(text)
    mapped_m = _PRIMARY_MAPPED_RE.search(text)
    if primary_m is None or mapped_m is None:
        raise Minimap2ParseError(
            "samtools flagstat output has no 'primary' / 'primary mapped' lines — "
            "alignment may have failed"
        )
    primary = int(primary_m.group(1))
    mapped = int(mapped_m.group(1))
    if primary == 0:
        raise Minimap2ParseError(
            "samtools flagstat reports 0 primary reads — no alignable input"
        )
    return mapped / primary


def _run(cmd: list[str]) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, capture_output=True, text=True)


def run_minimap2(genome_fasta: Path, out_dir: Path, fastq: Path, preset: str,
                 threads: int = 4, env: str = "rnaforge-longread",
                 secondary_n: int | None = None) -> AlignmentResult:
    genome_fasta = Path(genome_fasta)
    if not genome_fasta.exists():
        raise Minimap2RunError(f"genome FASTA does not exist: {genome_fasta}")
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    sam = out_dir / "aligned.sam"
    bam = out_dir / "aligned.sorted.bam"
    log = out_dir / "minimap2.log"

    mm = ["conda", "run", "-n", env, "minimap2", "-ax", preset]
    # İzoform EM (NanoCount) transkript başına çok-hizalama ister → -N ile ikincilleri sakla.
    # Default None: primary-only davranış (gen-yolu / prok-long değişmez).
    if secondary_n is not None:
        mm += ["-N", str(secondary_n)]
    mm += ["-t", str(threads), str(genome_fasta), str(fastq), "-o", str(sam)]
    r = _run(mm)
    log.write_text(r.stderr)
    if r.returncode != 0 or not sam.exists():
        raise Minimap2RunError(
            f"minimap2 failed (exit {r.returncode})\ncmd: {' '.join(mm)}\n"
            f"stderr: {r.stderr.strip()[-800:]}"
        )

    sort = _run(["conda", "run", "-n", env, "samtools", "sort",
                 "-@", str(threads), "-o", str(bam), str(sam)])
    if sort.returncode != 0 or not bam.exists():
        raise Minimap2RunError(f"samtools sort failed: {sort.stderr.strip()}")
    idx = _run(["conda", "run", "-n", env, "samtools", "index", str(bam)])
    if idx.returncode != 0:
        raise Minimap2RunError(f"samtools index failed: {idx.stderr.strip()}")

    flag = _run(["conda", "run", "-n", env, "samtools", "flagstat", str(bam)])
    if flag.returncode != 0:
        raise Minimap2RunError(f"samtools flagstat failed: {flag.stderr.strip()}")
    rate = parse_flagstat_mapped(flag.stdout)

    sam.unlink(missing_ok=True)
    return AlignmentResult(bam=bam, alignment_rate=rate)


def count_primary_alignments(bam_path: Path, env: str = "rnaforge-longread") -> dict[str, int]:
    """Primer hizalama sayımı hedef (transkript) başına. -F 2308 = unmapped(4)+
    secondary(256)+supplementary(2048) hariç → okuma başına tek satır; sütun 3 = hedef.
    Ökaryot uzun-okuma gen-düzeyi sayımı (tx2gene ile gene toplanır)."""
    r = _run(["conda", "run", "-n", env, "samtools", "view", "-F", "2308", str(bam_path)])
    if r.returncode != 0:
        raise Minimap2RunError(f"samtools view failed: {r.stderr.strip()[-500:]}")
    counts: dict[str, int] = {}
    for line in r.stdout.splitlines():
        if not line:
            continue
        ref = line.split("\t")[2]
        if ref and ref != "*":
            counts[ref] = counts.get(ref, 0) + 1
    return counts
