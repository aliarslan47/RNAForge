"""FastQC çıktısını parse eder ve çalıştırır. Parser saftır: string girer,
FastQCReport çıkar — I/O yok, hızlı ve deterministik test edilir."""
from __future__ import annotations

import subprocess
import zipfile
from dataclasses import dataclass
from pathlib import Path

_STATUSES = ("PASS", "WARN", "FAIL")


class FastQCParseError(ValueError):
    """FastQC çıktısı beklenen biçimde değil."""


class FastQCRunError(RuntimeError):
    """FastQC çalıştırılamadı ya da beklenen çıktıyı üretmedi."""


@dataclass(frozen=True)
class FastQCReport:
    modules: dict[str, str]
    basic_stats: dict[str, str]
    # Per-base A/T/G/C kompozisyonu: [(pozisyon etiketi, {"A":%, "T":%, "G":%, "C":%}), ...]
    per_base_content: tuple = ()
    # 'Total Deduplicated Percentage' (Sequence Duplication Levels); yoksa None.
    deduplication: float | None = None


def parse_fastqc_report(summary_text: str, data_text: str) -> FastQCReport:
    modules: dict[str, str] = {}
    for line in summary_text.splitlines():
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) < 2:
            raise FastQCParseError(f"summary line is not tab-delimited: {line!r}")
        status, name = parts[0].strip(), parts[1].strip()
        if status not in _STATUSES:
            raise FastQCParseError(
                f"unknown FastQC status {status!r} for module {name!r} "
                f"(expected one of {_STATUSES})"
            )
        modules[name] = status

    basic_stats = _parse_basic_stats(data_text)
    return FastQCReport(
        modules=modules, basic_stats=basic_stats,
        per_base_content=tuple(parse_per_base_content(data_text)),
        deduplication=parse_deduplication(data_text),
    )


def parse_per_base_content(data_text: str) -> list[tuple[str, dict[str, float]]]:
    """>>Per base sequence content modülünden pozisyon başına A/T/G/C yüzdeleri.
    Header sütun sırası FastQC sürümüne göre değişir (ör. 'G A T C') — sıra
    header'dan okunur. Modül yoksa boş liste döner (I/O yok, saf)."""
    rows: list[tuple[str, dict[str, float]]] = []
    header: list[str] | None = None
    in_mod = False
    for line in data_text.splitlines():
        if line.startswith(">>Per base sequence content"):
            in_mod = True
            continue
        if not in_mod:
            continue
        if line.startswith(">>END_MODULE"):
            break
        if line.startswith("#"):
            # #Base<TAB>G<TAB>A<TAB>T<TAB>C
            header = [c.strip() for c in line.lstrip("#").split("\t")][1:]
            continue
        if not line.strip() or header is None:
            continue
        fields = line.split("\t")
        pos = fields[0].strip()
        comp: dict[str, float] = {}
        for name, value in zip(header, fields[1:]):
            try:
                comp[name.upper()] = float(value)
            except ValueError:
                continue
        rows.append((pos, comp))
    return rows


def parse_deduplication(data_text: str) -> float | None:
    """Sequence Duplication Levels modülündeki 'Total Deduplicated Percentage'
    (benzersiz okuma yüzdesi). Bulunmazsa None."""
    for line in data_text.splitlines():
        if "Total Deduplicated Percentage" in line:
            _, _, value = line.partition("\t")
            try:
                return float(value.strip())
            except ValueError:
                return None
    return None


def _parse_basic_stats(data_text: str) -> dict[str, str]:
    stats: dict[str, str] = {}
    in_basic = False
    for line in data_text.splitlines():
        if line.startswith(">>Basic Statistics"):
            in_basic = True
            continue
        if in_basic:
            if line.startswith(">>END_MODULE"):
                return stats
            if line.startswith("#") or not line.strip():
                continue
            key, _, value = line.partition("\t")
            stats[key.strip()] = value.strip()
    raise FastQCParseError(
        "FastQC data has no '>>Basic Statistics' module — output is malformed"
    )


def parse_fastqc_zip(zip_path: Path) -> FastQCReport:
    zip_path = Path(zip_path)
    with zipfile.ZipFile(zip_path) as zf:
        summary_name = _find_member(zf, "summary.txt")
        data_name = _find_member(zf, "fastqc_data.txt")
        summary_text = zf.read(summary_name).decode()
        data_text = zf.read(data_name).decode()
    return parse_fastqc_report(summary_text, data_text)


def _find_member(zf: zipfile.ZipFile, suffix: str) -> str:
    for name in zf.namelist():
        if name.endswith(suffix):
            return name
    raise FastQCParseError(f"FastQC zip has no {suffix}: {zf.filename}")


def run_fastqc(fastq: Path, out_dir: Path, env: str = "rnaforge-qc") -> Path:
    fastq = Path(fastq)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = ["conda", "run", "-n", env, "fastqc", str(fastq), "-o", str(out_dir), "-q"]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise FastQCRunError(
            f"FastQC failed (exit {result.returncode}) for {fastq}\n"
            f"cmd: {' '.join(cmd)}\nstderr: {result.stderr.strip()}"
        )
    # FastQC çıktı adı: <fastq stem, .gz/.fastq atılmış>_fastqc.zip
    stem = fastq.name
    for ext in (".gz", ".fastq", ".fq"):
        if stem.endswith(ext):
            stem = stem[: -len(ext)]
    zip_path = out_dir / f"{stem}_fastqc.zip"
    if not zip_path.exists():
        raise FastQCRunError(
            f"FastQC reported success but produced no zip at {zip_path} "
            f"(out_dir contents: {[p.name for p in out_dir.iterdir()]})"
        )
    return zip_path
