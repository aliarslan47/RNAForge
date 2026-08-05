"""m16 — Sekans QC: rRNA% (SortMeRNA) + strandedness (RSeQC). Referans/BED üretimi + runner + parser.

rRNA referansı referans genomdan çıkarılır (indirme yok, agnostik). Strandedness BAM'den çıkarılıp
config beyanıyla karşılaştırılır. Kötü girdiyi yakalayıp WARN kapısına dönüştürür ([[feedback_dogruluk_kontrol]]).
"""
from __future__ import annotations

import gzip
import re
import shutil
import subprocess
from pathlib import Path

_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def _read_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, str] = {}
    name = None
    parts: list[str] = []
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                seqs[name] = "".join(parts)
            name = line[1:].split()[0]
            parts = []
        else:
            parts.append(line.strip())
    if name is not None:
        seqs[name] = "".join(parts)
    return seqs


def rrna_fasta_from_reference(genome_fa: Path, gff: Path, out_fasta: Path) -> int:
    """GFF `rRNA` feature koordinatlarından genom rRNA dizilerini çıkar (− strand ters-tümler).
    SortMeRNA --ref girdisi. Yazılan rRNA dizisi sayısını döndürür."""
    genome = _read_fasta(genome_fa)
    n = 0
    with Path(out_fasta).open("w") as f:
        for line in Path(gff).read_text().splitlines():
            if not line or line.startswith("#") or "\trRNA\t" not in line:
                continue
            c = line.split("\t")
            if len(c) < 8:
                continue
            contig, start, end, strand = c[0], int(c[3]), int(c[4]), c[6]
            seq = genome.get(contig, "")[start - 1:end]
            if not seq:
                continue
            if strand == "-":
                seq = seq.translate(_COMP)[::-1]
            attrs = dict(kv.split("=", 1) for kv in c[8].split(";") if "=" in kv)
            f.write(f">{attrs.get('locus_tag', f'rRNA_{n+1}')}\n{seq}\n")
            n += 1
    return n


def subsample_fastq(src_gz: Path, n_reads: int, out_path: Path) -> int:
    """İlk n_reads okumayı (4·n satır) düz FASTQ'a yaz. Yazılan okuma sayısını döndürür."""
    written = 0
    with gzip.open(src_gz, "rt") as fi, Path(out_path).open("w") as fo:
        for i, line in enumerate(fi):
            if i >= n_reads * 4:
                break
            fo.write(line)
            if i % 4 == 0:
                written += 1
    return written


def run_sortmerna(reads: list[Path], ref_fasta: Path, workdir: Path,
                  threads: int = 8, env: str = "rnaforge-seqqc") -> str:
    """SortMeRNA'yı temiz bir workdir'de çalıştır (kvdb/idx yeniden kullanım hatasını önle).
    aligned.log yolunu döndürür; hatada gürültülü yüksel."""
    workdir = Path(workdir)
    if workdir.exists():
        shutil.rmtree(workdir)
    workdir.mkdir(parents=True)
    cmd = ["conda", "run", "-n", env, "sortmerna", "--ref", str(ref_fasta),
           "--workdir", str(workdir), "--num_alignments", "1", "--threads", str(threads), "-v"]
    for r in reads:
        cmd += ["--reads", str(r)]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError(f"sortmerna failed (exit {res.returncode}):\n{res.stderr[-2000:]}")
    return str(workdir / "out" / "aligned.log")


def parse_sortmerna_log(log_path: Path) -> float:
    """aligned.log'tan rRNA fraksiyonu (0-1). 'Total reads passing E-value threshold = N (X%)'."""
    text = Path(log_path).read_text() if Path(log_path).exists() else ""
    m = re.search(r"passing E-value threshold\s*=\s*\d+\s*\(([\d.]+)%?\)", text)
    if m:
        return float(m.group(1)) / 100.0
    total = re.search(r"Total reads\s*=\s*(\d+)", text)
    passed = re.search(r"passing E-value threshold\s*=\s*(\d+)", text)
    if total and passed and int(total.group(1)) > 0:
        return int(passed.group(1)) / int(total.group(1))
    return 0.0


def gff_to_bed(gff: Path, out_bed: Path) -> int:
    """GFF `gene` feature'ları → BED12 (prokaryotta gen = tek blok). RSeQC infer_experiment girdisi."""
    n = 0
    with Path(out_bed).open("w") as f:
        for line in Path(gff).read_text().splitlines():
            if not line or line.startswith("#") or "\tgene\t" not in line:
                continue
            c = line.split("\t")
            if len(c) < 9:
                continue
            attrs = dict(kv.split("=", 1) for kv in c[8].split(";") if "=" in kv)
            lt = attrs.get("locus_tag")
            if not lt:
                continue
            start0, end, strand = int(c[3]) - 1, int(c[4]), c[6]
            size = end - start0
            f.write(f"{c[0]}\t{start0}\t{end}\t{lt}\t0\t{strand}\t{start0}\t{end}\t0\t1\t{size},\t0,\n")
            n += 1
    return n


def run_infer_experiment(bam: Path, bed: Path, env: str = "rnaforge-seqqc") -> str:
    """RSeQC infer_experiment.py → stdout (metin). Hatada gürültülü yüksel."""
    cmd = ["conda", "run", "-n", env, "infer_experiment.py", "-i", str(bam), "-r", str(bed)]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError(f"infer_experiment.py failed (exit {res.returncode}):\n{res.stderr[-2000:]}")
    return res.stdout


def parse_infer_experiment(output: str, threshold: float = 0.8) -> tuple[str, float, float]:
    """RSeQC çıktısı → (strand, fwd_frac, rev_frac). strand ∈ {unstranded, stranded, reverse}.
    '1++,1--,2+-,2-+' = forward (stranded); '1+-,1-+,2++,2--' = reverse. Baskınlık < threshold → unstranded."""
    fwd = rev = 0.0
    for line in output.splitlines():
        m = re.search(r'explained by "([^"]+)":\s*([\d.]+)', line)
        if not m:
            continue
        pattern, frac = m.group(1), float(m.group(2))
        if pattern.startswith("++") or pattern in ("1++,1--,2+-,2-+", "++,--"):
            fwd = frac
        elif pattern.startswith("+-") or pattern in ("1+-,1-+,2++,2--", "+-,-+"):
            rev = frac
    if fwd >= threshold and fwd > rev:
        return "stranded", fwd, rev
    if rev >= threshold and rev > fwd:
        return "reverse", fwd, rev
    return "unstranded", fwd, rev
