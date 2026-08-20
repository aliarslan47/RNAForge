"""Kraken2 taksonomik profilleme ve Bracken zenginleştirmesi: çalıştırır ve parse eder.

Parser saftır. bowtie2.py/minimap2.py deseni izlenir."""
from __future__ import annotations

import subprocess
from pathlib import Path


class Kraken2ParseError(ValueError):
    """Kraken2/Bracken çıktısı beklenen biçimde değil."""


class Kraken2RunError(RuntimeError):
    """Kraken2/Bracken çalıştırılamadı ya da beklenen çıktıyı üretmedi."""


def parse_kraken2_report(path: Path) -> list[dict]:
    """Kraken2 report parse: fraction(%), clade_reads, taxon_reads, rank, taxid, name.

    Returns: list[dict] with keys {rank, taxid, name, reads, fraction}.
    fraction is in range [0, 1].
    reads is taxon_reads (3rd column).
    """
    path = Path(path)
    if not path.exists():
        raise Kraken2ParseError(f"Kraken2 report does not exist: {path}")

    rows = []
    try:
        with path.open() as f:
            for line_no, line in enumerate(f, 1):
                line = line.strip()
                if not line:
                    continue
                parts = line.split("\t")
                if len(parts) < 6:
                    raise Kraken2ParseError(
                        f"Kraken2 report line {line_no} has < 6 fields: {line}"
                    )
                try:
                    fraction_pct = float(parts[0])
                    # parts[1] is clade_reads (unused)
                    reads = int(parts[2])
                    rank = parts[3].strip()
                    taxid = parts[4].strip()
                    name = parts[5].strip()  # Remove indentation and trailing whitespace
                except (ValueError, IndexError) as e:
                    raise Kraken2ParseError(
                        f"Kraken2 report line {line_no} parse error: {e}"
                    ) from e
                rows.append({
                    "fraction": fraction_pct / 100.0,
                    "reads": reads,
                    "rank": rank,
                    "taxid": taxid,
                    "name": name,
                })
    except Kraken2ParseError:
        raise
    except Exception as e:
        raise Kraken2ParseError(f"Error reading Kraken2 report: {e}") from e

    return rows


def parse_bracken(path: Path) -> dict[str, float]:
    """Bracken output parse: name, taxonomy_id, taxonomy_lvl, kraken_assigned_reads,
    added_reads, new_est_reads, fraction_total_reads.

    Returns: dict[taxon_name -> fraction_total_reads].
    """
    path = Path(path)
    if not path.exists():
        raise Kraken2ParseError(f"Bracken file does not exist: {path}")

    result = {}
    try:
        with path.open() as f:
            lines = f.readlines()
            if not lines:
                raise Kraken2ParseError("Bracken file is empty")
            # Skip header
            for line_no, line in enumerate(lines[1:], 2):
                line = line.strip()
                if not line:
                    continue
                parts = line.split("\t")
                if len(parts) < 7:
                    raise Kraken2ParseError(
                        f"Bracken line {line_no} has < 7 fields: {line}"
                    )
                try:
                    name = parts[0].strip()
                    # parts[1] is taxonomy_id (unused for this output)
                    # parts[2] is taxonomy_lvl (unused)
                    # parts[3] is kraken_assigned_reads (unused)
                    # parts[4] is added_reads (unused)
                    # parts[5] is new_est_reads (unused)
                    fraction = float(parts[6])
                except (ValueError, IndexError) as e:
                    raise Kraken2ParseError(
                        f"Bracken line {line_no} parse error: {e}"
                    ) from e
                result[name] = fraction
    except Kraken2ParseError:
        raise
    except Exception as e:
        raise Kraken2ParseError(f"Error reading Bracken file: {e}") from e

    return result


def _run(cmd: list[str]) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, capture_output=True, text=True)


def run_kraken2(reads: list[Path], db: Path, out_prefix: Path, paired: bool = False,
                threads: int = 4, env: str = "rnaforge-meta") -> Path:
    """Run Kraken2 taxonomic profiling.

    Args:
        reads: List of FASTQ file(s). If paired=True, must be [r1, r2].
               If paired=False, must be single-element list [r1].
        db: Kraken2 database directory.
        out_prefix: Output prefix (report will be <out_prefix>.report).
        paired: If True, expects reads to have 2 elements; adds --paired flag.
        threads: Number of threads.
        env: Conda environment name (default: "rnaforge-meta").

    Returns: Path to the generated report file (<out_prefix>.report).

    Raises: Kraken2RunError on failure.
    """
    reads = [Path(r) for r in reads]
    db = Path(db)
    out_prefix = Path(out_prefix)

    expected_count = 2 if paired else 1
    if len(reads) != expected_count:
        raise Kraken2RunError(
            f"expected {expected_count} read file(s), got {len(reads)}"
        )

    for r in reads:
        if not r.exists():
            raise Kraken2RunError(f"input reads file does not exist: {r}")
    if not db.exists():
        raise Kraken2RunError(f"Kraken2 database does not exist: {db}")

    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    report = Path(str(out_prefix) + ".report")

    cmd = ["conda", "run", "-n", env, "kraken2", "--db", str(db),
           "--report", str(report), "--threads", str(threads)]

    if paired:
        cmd += ["--paired"]

    for r in reads:
        cmd += [str(r)]

    r = _run(cmd)
    if r.returncode != 0 or not report.exists():
        raise Kraken2RunError(
            f"Kraken2 failed (exit {r.returncode})\ncmd: {' '.join(cmd)}\n"
            f"stderr: {r.stderr.strip()[-800:]}"
        )
    return report


def run_bracken(kraken_report: Path, db: Path, out_path: Path, read_len: int = 100,
                level: str = "S", env: str = "rnaforge-meta") -> Path:
    """Run Bracken abundance estimation on Kraken2 report.

    Args:
        kraken_report: Path to Kraken2 report file.
        db: Bracken database directory (same as Kraken2 db).
        out_path: Output file path.
        read_len: Read length (default 100).
        level: Taxonomic level (default "S" for species).
        env: Conda environment name (default: "rnaforge-meta").

    Returns: Path to the output file.

    Raises: Kraken2RunError on failure.
    """
    kraken_report = Path(kraken_report)
    db = Path(db)
    out_path = Path(out_path)

    if not kraken_report.exists():
        raise Kraken2RunError(f"Kraken2 report does not exist: {kraken_report}")
    if not db.exists():
        raise Kraken2RunError(f"Bracken database does not exist: {db}")

    out_path.parent.mkdir(parents=True, exist_ok=True)

    cmd = ["conda", "run", "-n", env, "bracken", "-d", str(db),
           "-i", str(kraken_report), "-o", str(out_path),
           "-r", str(read_len), "-l", level]

    r = _run(cmd)
    if r.returncode != 0 or not out_path.exists():
        raise Kraken2RunError(
            f"Bracken failed (exit {r.returncode})\ncmd: {' '.join(cmd)}\n"
            f"stderr: {r.stderr.strip()[-800:]}"
        )
    return out_path
