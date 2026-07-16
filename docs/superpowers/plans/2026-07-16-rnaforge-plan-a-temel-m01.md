# RNAForge Plan A — Temel + m01 (Doğrulama & Platform Tespiti)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** `rnaforge validate` komutu çalışır hale gelir: config'i doğrular, metadata'yı okur, FASTQ'lardan platformu tespit eder, Illumina dışını net hatayla reddeder ve run dizinine kalıcı durum kaydı yazar.

**Architecture:** Python paketi (`rnaforge/`), `ali-wgs-pipeline` deseni izlenir: küçük odaklı modüller, YAML config, `Module` ABC sözleşmesi, zaman damgalı `runs/<ts>_<id>/` dizinleri. Bu plandaki hiçbir kod harici biyoinformatik araca ihtiyaç duymaz (saf stdlib + pyyaml) — conda ortamları Plan B'de gelir.

**Tech Stack:** Python ≥3.10, pyyaml, pytest. Conda env `rnaforge-core`.

**Kaynak spec:** `docs/superpowers/specs/2026-07-16-rnaforge-mvp-design.md`
**Referans doküman:** `PLAN.md` v1.2

## Global Constraints

Bu kısıtlar her görevin gereksinimlerine dahildir:

- **`organism_type` ZORUNLUDUR, varsayılanı YOKTUR** (`prokaryote` | `eukaryote`). Yanlış varsayım sessiz hataya yol açar. (PLAN §2.1, §11)
- **Tespit etmek ≠ desteklemek** (PLAN Kural 7): desteklenmeyen girdi sessizce işlenmez, net hatayla reddedilir. MVP yalnızca Illumina.
- **Trimming varsayılanı NAZİK** olmalı: `aggressive_quality: false`, `min_length` zorunlu. Gerekçe Williams et al. 2016 (PLAN §4.2). Bu bir testle sabitlenir.
- **Kod, değişken adları, log mesajları İngilizce.** Docstring'ler Türkçe olabilir (WGS deseni). Kullanıcıya giden hata mesajları İngilizce.
- **Müşteri verisi/PII asla commit edilmez.** Testlerde yalnız sentetik veri; gerçek örnek ID'si/hasta bilgisi kullanılmaz. (PLAN Kural 8)
- **Kapatma dayanıklılığı:** 10 sn heartbeat + her modül bitiminde kalıcı durum kaydı → resume. (PLAN §15)
- Python ≥3.10 (`from __future__ import annotations` her modülde).
- Her dosya tek sorumluluk; WGS'de hiçbir çekirdek dosya 152 satırı geçmiyor — bu ölçek korunur.

---

### Task 1: Paket iskeleti + conda env + CLI

**Files:**
- Create: `pyproject.toml`
- Create: `envs/rnaforge-core.yml`
- Create: `rnaforge/__init__.py`
- Create: `rnaforge/cli.py`
- Test: `tests/test_cli.py`

**Interfaces:**
- Consumes: yok (ilk görev)
- Produces: `rnaforge.__version__: str`; `rnaforge.cli.main(argv: list[str] | None = None) -> int` — 0 başarı, 2 kullanım hatası. Konsol betiği `rnaforge`.

- [ ] **Step 1: Conda ortamını oluştur**

`envs/rnaforge-core.yml`:

```yaml
name: rnaforge-core
channels:
  - conda-forge
dependencies:
  - python=3.11
  - pyyaml>=6
  - pytest>=8
  - pip
```

Çalıştır:

```bash
cd /home/ali/rnaforge-pipeline
conda env create -f envs/rnaforge-core.yml
conda run -n rnaforge-core python --version
```

Beklenen: `Python 3.11.x`

- [ ] **Step 2: pyproject.toml yaz**

```toml
[build-system]
requires = ["setuptools>=68"]
build-backend = "setuptools.build_meta"

[project]
name = "rnaforge"
version = "0.1.0"
description = "RNAForge — Bulk RNA-seq analiz pipeline'ı"
requires-python = ">=3.10"
dependencies = ["pyyaml>=6"]

[project.scripts]
rnaforge = "rnaforge.cli:main"

[tool.setuptools.packages.find]
include = ["rnaforge*"]

[tool.pytest.ini_options]
testpaths = ["tests"]
```

- [ ] **Step 3: Başarısız testi yaz**

`tests/test_cli.py`:

```python
from __future__ import annotations

import pytest

from rnaforge.cli import main


def test_version_flag_prints_version(capsys):
    with pytest.raises(SystemExit) as exc:
        main(["--version"])
    assert exc.value.code == 0
    assert "0.1.0" in capsys.readouterr().out


def test_no_command_returns_usage_error():
    assert main([]) == 2
```

- [ ] **Step 4: Testi çalıştır, BAŞARISIZ olduğunu gör**

```bash
conda run -n rnaforge-core python -m pytest tests/test_cli.py -v
```

Beklenen: FAIL — `ModuleNotFoundError: No module named 'rnaforge'`

- [ ] **Step 5: Paketi yaz**

`rnaforge/__init__.py`:

```python
"""RNAForge — Bulk RNA-seq analiz pipeline'ı."""
from __future__ import annotations

__version__ = "0.1.0"
```

`rnaforge/cli.py`:

```python
"""Komut satırı girişi."""
from __future__ import annotations

import argparse

from rnaforge import __version__


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="rnaforge", description="Bulk RNA-seq pipeline")
    parser.add_argument("--version", action="version", version=__version__)
    parser.add_argument("command", nargs="?", help="validate")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command is None:
        print("error: no command given (try: rnaforge validate)")
        return 2
    return 0
```

- [ ] **Step 6: Kur ve testi çalıştır, GEÇTİĞİNİ gör**

```bash
conda run -n rnaforge-core pip install -e .
conda run -n rnaforge-core python -m pytest tests/test_cli.py -v
```

Beklenen: 2 passed

- [ ] **Step 7: Commit**

```bash
git add pyproject.toml envs/rnaforge-core.yml rnaforge/ tests/
git commit -m "feat: paket iskeleti, conda env ve CLI girişi"
```

---

### Task 2: Config yükleme + şema doğrulama

**Files:**
- Create: `rnaforge/config.py`
- Create: `config/config.yaml`
- Test: `tests/test_config.py`

**Interfaces:**
- Consumes: yok
- Produces:
  - `rnaforge.config.ConfigError(ValueError)`
  - Dataclass'lar (hepsi `frozen=True`): `Reference(genome_fasta: Path|None, annotation_gff: Path|None, transcriptome_fasta: Path|None, tx2gene: Path|None)`, `Library(strandedness: str, selection: str)`, `Trimming(min_length: int, aggressive_quality: bool)`, `DE(design: str, fdr_threshold: float, log2fc_threshold: float)`, `Report(language: str)`, `Resources(threads: int, memory_gb: int)`
  - `Config(organism: str, organism_type: str, platform: str, reference: Reference, library: Library, trimming: Trimming, de: DE, report: Report, resources: Resources)`
  - `load_config(path: Path | str) -> Config` — geçersizse `ConfigError`

- [ ] **Step 1: Başarısız testleri yaz**

`tests/test_config.py`:

```python
from __future__ import annotations

import textwrap

import pytest

from rnaforge.config import ConfigError, load_config


def _write(tmp_path, body: str):
    path = tmp_path / "config.yaml"
    path.write_text(textwrap.dedent(body))
    return path


# DİKKAT: girintisiz tanımlı — testler buna satır EKLİYOR (PROK_BODY + '...').
# Girintili olsaydı eklenen girintisiz satır textwrap.dedent'i bozardı ve
# YAML parse hatası alırdık, beklediğimiz ConfigError'ı değil.
PROK_BODY = """organism: "Escherichia coli"
organism_type: "prokaryote"
reference:
  genome_fasta: "ref/genome.fa"
  annotation_gff: "ref/genes.gff"
de:
  design: "~condition"
"""


def test_valid_prokaryote_config_loads(tmp_path):
    cfg = load_config(_write(tmp_path, PROK_BODY))
    assert cfg.organism_type == "prokaryote"
    assert cfg.reference.genome_fasta.name == "genome.fa"


def test_missing_organism_type_raises(tmp_path):
    path = _write(tmp_path, """
        organism: "Escherichia coli"
        reference:
          genome_fasta: "ref/genome.fa"
          annotation_gff: "ref/genes.gff"
    """)
    with pytest.raises(ConfigError, match="organism_type"):
        load_config(path)


def test_invalid_organism_type_raises(tmp_path):
    path = _write(tmp_path, """
        organism: "X"
        organism_type: "virus"
        reference:
          genome_fasta: "a"
          annotation_gff: "b"
    """)
    with pytest.raises(ConfigError, match="prokaryote"):
        load_config(path)


def test_prokaryote_requires_genome_and_annotation(tmp_path):
    path = _write(tmp_path, """
        organism: "Escherichia coli"
        organism_type: "prokaryote"
        reference:
          transcriptome_fasta: "ref/tx.fa"
    """)
    with pytest.raises(ConfigError, match="genome_fasta"):
        load_config(path)


def test_eukaryote_requires_transcriptome_and_tx2gene(tmp_path):
    path = _write(tmp_path, """
        organism: "Homo sapiens"
        organism_type: "eukaryote"
        reference:
          genome_fasta: "ref/genome.fa"
    """)
    with pytest.raises(ConfigError, match="transcriptome_fasta"):
        load_config(path)


def test_trimming_defaults_are_gentle(tmp_path):
    """PLAN §4.2: agresif trimming ekspresyon tahminlerini bozar (Williams 2016).
    Varsayılan NAZİK olmalı; bu test o kararı sabitler."""
    cfg = load_config(_write(tmp_path, PROK_BODY))
    assert cfg.trimming.aggressive_quality is False
    assert cfg.trimming.min_length >= 1


def test_platform_defaults_to_auto(tmp_path):
    assert load_config(_write(tmp_path, PROK_BODY)).platform == "auto"


def test_invalid_platform_raises(tmp_path):
    path = _write(tmp_path, PROK_BODY + '\nplatform: "ont"\n')
    with pytest.raises(ConfigError, match="platform"):
        load_config(path)


def test_invalid_strandedness_raises(tmp_path):
    path = _write(tmp_path, PROK_BODY + '\nlibrary:\n  strandedness: "sideways"\n')
    with pytest.raises(ConfigError, match="strandedness"):
        load_config(path)
```

- [ ] **Step 2: Testleri çalıştır, BAŞARISIZ olduğunu gör**

```bash
conda run -n rnaforge-core python -m pytest tests/test_config.py -v
```

Beklenen: FAIL — `ModuleNotFoundError: No module named 'rnaforge.config'`

- [ ] **Step 3: config.py'yi yaz**

`rnaforge/config.py`:

```python
"""Config yükleme ve şema doğrulama. Geçersiz config m01'den ÖNCE yakalanır."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import yaml

ORGANISM_TYPES = ("prokaryote", "eukaryote")
PLATFORMS = ("auto", "illumina")
STRANDEDNESS = ("unstranded", "stranded", "reverse")
SELECTIONS = ("rrna_depletion", "polya")
REPORT_LANGUAGES = ("tr", "en")

# organism_type -> zorunlu reference alanları
REQUIRED_REFERENCE = {
    "prokaryote": ("genome_fasta", "annotation_gff"),
    "eukaryote": ("transcriptome_fasta", "tx2gene"),
}


class ConfigError(ValueError):
    """Config dosyası geçersiz."""


@dataclass(frozen=True)
class Reference:
    genome_fasta: Path | None = None
    annotation_gff: Path | None = None
    transcriptome_fasta: Path | None = None
    tx2gene: Path | None = None


@dataclass(frozen=True)
class Library:
    strandedness: str = "unstranded"
    selection: str = "rrna_depletion"


@dataclass(frozen=True)
class Trimming:
    # PLAN §4.2 — nazik varsayılan, literatür temelli. Değiştirmeden önce oku.
    min_length: int = 36
    aggressive_quality: bool = False


@dataclass(frozen=True)
class DE:
    design: str = "~condition"
    fdr_threshold: float = 0.05
    log2fc_threshold: float = 1.0


@dataclass(frozen=True)
class Report:
    language: str = "tr"


@dataclass(frozen=True)
class Resources:
    threads: int = 8
    memory_gb: int = 32


@dataclass(frozen=True)
class Config:
    organism: str
    organism_type: str
    platform: str
    reference: Reference
    library: Library
    trimming: Trimming
    de: DE
    report: Report
    resources: Resources


def _one_of(value, allowed, field: str):
    if value not in allowed:
        raise ConfigError(f"{field}: {value!r} is invalid; allowed: {', '.join(allowed)}")
    return value


def _require_organism_type(raw: dict) -> str:
    value = raw.get("organism_type")
    if value is None:
        raise ConfigError(
            "organism_type is required and has no default "
            f"(allowed: {', '.join(ORGANISM_TYPES)}). "
            "Guessing it would silently produce wrong results."
        )
    return _one_of(value, ORGANISM_TYPES, "organism_type")


def _build_reference(raw: dict, organism_type: str) -> Reference:
    missing = [f for f in REQUIRED_REFERENCE[organism_type] if not raw.get(f)]
    if missing:
        raise ConfigError(
            f"organism_type={organism_type} requires reference fields: {', '.join(missing)}"
        )
    to_path = lambda v: Path(v) if v else None  # noqa: E731
    return Reference(
        genome_fasta=to_path(raw.get("genome_fasta")),
        annotation_gff=to_path(raw.get("annotation_gff")),
        transcriptome_fasta=to_path(raw.get("transcriptome_fasta")),
        tx2gene=to_path(raw.get("tx2gene")),
    )


def load_config(path: Path | str) -> Config:
    path = Path(path)
    if not path.exists():
        raise ConfigError(f"config file not found: {path}")
    raw = yaml.safe_load(path.read_text()) or {}
    if not isinstance(raw, dict):
        raise ConfigError(f"config must be a YAML mapping: {path}")

    organism = raw.get("organism")
    if not organism:
        raise ConfigError("organism is required")

    organism_type = _require_organism_type(raw)
    library_raw = raw.get("library") or {}
    trimming_raw = raw.get("trimming") or {}
    de_raw = raw.get("de") or {}
    report_raw = raw.get("report") or {}
    resources_raw = raw.get("resources") or {}

    trimming = Trimming(
        min_length=int(trimming_raw.get("min_length", 36)),
        aggressive_quality=bool(trimming_raw.get("aggressive_quality", False)),
    )
    if trimming.min_length < 1:
        raise ConfigError(f"trimming.min_length must be >= 1, got {trimming.min_length}")

    return Config(
        organism=str(organism),
        organism_type=organism_type,
        platform=_one_of(raw.get("platform", "auto"), PLATFORMS, "platform"),
        reference=_build_reference(raw.get("reference") or {}, organism_type),
        library=Library(
            strandedness=_one_of(
                library_raw.get("strandedness", "unstranded"), STRANDEDNESS, "library.strandedness"
            ),
            selection=_one_of(
                library_raw.get("selection", "rrna_depletion"), SELECTIONS, "library.selection"
            ),
        ),
        trimming=trimming,
        de=DE(
            design=str(de_raw.get("design", "~condition")),
            fdr_threshold=float(de_raw.get("fdr_threshold", 0.05)),
            log2fc_threshold=float(de_raw.get("log2fc_threshold", 1.0)),
        ),
        report=Report(
            language=_one_of(report_raw.get("language", "tr"), REPORT_LANGUAGES, "report.language")
        ),
        resources=Resources(
            threads=int(resources_raw.get("threads", 8)),
            memory_gb=int(resources_raw.get("memory_gb", 32)),
        ),
    )
```

- [ ] **Step 4: Testleri çalıştır, GEÇTİĞİNİ gör**

```bash
conda run -n rnaforge-core python -m pytest tests/test_config.py -v
```

Beklenen: 9 passed

- [ ] **Step 5: Örnek config'i yaz**

`config/config.yaml`:

```yaml
# RNAForge — örnek config. Alan açıklamaları: PLAN.md Bölüm 11.
organism: "Escherichia coli"

# ZORUNLU, varsayılanı yok: prokaryote | eukaryote (PLAN §2.1)
# Akışı dallandırır: prokaryote -> bowtie2+featureCounts, eukaryote -> Salmon+tximport
organism_type: "prokaryote"

# auto: FASTQ'dan tespit et. MVP yalnızca Illumina destekler (PLAN §4.1).
platform: "auto"

reference:
  # prokaryote için:
  genome_fasta: "references/ecoli/genome.fa"
  annotation_gff: "references/ecoli/genes.gff"
  # eukaryote için (organism_type: eukaryote ise bunlar zorunlu olur):
  # transcriptome_fasta: "references/human/tx.fa"
  # tx2gene: "references/human/tx2gene.tsv"

library:
  # FASTQ'dan TESPİT EDİLEMEZ, buradan gelmeli (PLAN §4.1)
  strandedness: "unstranded"   # unstranded | stranded | reverse
  selection: "rrna_depletion"  # rrna_depletion | polya

trimming:
  # PLAN §4.2: agresif kalite trimming ekspresyon tahminlerini bozar
  # (Williams et al. 2016). Değiştirmeden önce gerekçeyi oku.
  min_length: 36
  aggressive_quality: false

de:
  design: "~condition"         # batch varsa: "~batch + condition"
  fdr_threshold: 0.05
  log2fc_threshold: 1.0

report:
  language: "tr"               # tr | en

resources:
  threads: 8
  memory_gb: 32
```

- [ ] **Step 6: Commit**

```bash
git add rnaforge/config.py config/config.yaml tests/test_config.py
git commit -m "feat: config yükleme + şema doğrulama (organism_type zorunlu, nazik trimming varsayılanı)"
```

---

### Task 3: Platform tespiti

**Files:**
- Create: `rnaforge/platform.py`
- Test: `tests/test_platform.py`
- Test: `tests/conftest.py`

**Interfaces:**
- Consumes: yok
- Produces:
  - `rnaforge.platform.PlatformInfo(platform: str, mean_read_length: float, n50: int, mean_quality: float, n_reads_sampled: int)` — `platform` ∈ `"illumina" | "ont" | "pacbio_hifi" | "unknown"`
  - `UnsupportedPlatformError(RuntimeError)`
  - `detect_platform(fastq: Path, short_read_max_len: int = 350, hifi_min_qual: float = 25.0, max_reads: int = 5000) -> PlatformInfo`
  - `require_supported(info: PlatformInfo, fastq: Path) -> None` — Illumina değilse `UnsupportedPlatformError`
  - `SUPPORTED_PLATFORMS = ("illumina",)`

**Not:** Eşik mantığı `ali-wgs-pipeline/ali_wgs/detect.py`'den uyarlanır (aynı problem, kanıtlanmış yaklaşım): kısa read → Illumina; uzun + yüksek kalite + N50≥5000 → PacBio HiFi; aksi uzun → ONT.

- [ ] **Step 1: Test fixture'larını yaz**

`tests/conftest.py`:

```python
"""Sentetik FASTQ fixture'ları. Gerçek/müşteri verisi ASLA kullanılmaz (PLAN Kural 8)."""
from __future__ import annotations

import gzip
import random
from pathlib import Path

import pytest


def _record(name: str, seq_len: int, qual_char: str) -> str:
    seq = "".join(random.choice("ACGT") for _ in range(seq_len))
    return f"@{name}\n{seq}\n+\n{qual_char * seq_len}\n"


def write_fastq(path: Path, n_reads: int, seq_len, qual_char: str, gzipped: bool = False) -> Path:
    """seq_len: int (sabit) veya (min, max) tuple (değişken uzunluk)."""
    def length() -> int:
        return seq_len if isinstance(seq_len, int) else random.randint(*seq_len)

    body = "".join(_record(f"read{i}", length(), qual_char) for i in range(n_reads))
    if gzipped:
        with gzip.open(path, "wt") as fh:
            fh.write(body)
    else:
        path.write_text(body)
    return path


@pytest.fixture(autouse=True)
def _seed():
    random.seed(1337)


@pytest.fixture
def illumina_fastq(tmp_path) -> Path:
    # 150 bp sabit, Q40 ('I')
    return write_fastq(tmp_path / "illumina.fastq", 200, 150, "I")


@pytest.fixture
def ont_fastq(tmp_path) -> Path:
    # uzun ve gürültülü: 1-20 kb, Q10 ('+')
    return write_fastq(tmp_path / "ont.fastq", 200, (1000, 20000), "+")


@pytest.fixture
def pacbio_fastq(tmp_path) -> Path:
    # uzun ve yüksek kaliteli: 8-20 kb, Q40 ('I')
    return write_fastq(tmp_path / "pacbio.fastq", 200, (8000, 20000), "I")
```

- [ ] **Step 2: Başarısız testleri yaz**

`tests/test_platform.py`:

```python
from __future__ import annotations

import pytest

from rnaforge.platform import (
    UnsupportedPlatformError,
    detect_platform,
    require_supported,
)
from tests.conftest import write_fastq


def test_short_reads_detected_as_illumina(illumina_fastq):
    info = detect_platform(illumina_fastq)
    assert info.platform == "illumina"
    assert info.mean_read_length == pytest.approx(150.0)
    assert info.n_reads_sampled == 200


def test_long_noisy_reads_detected_as_ont(ont_fastq):
    assert detect_platform(ont_fastq).platform == "ont"


def test_long_high_quality_reads_detected_as_pacbio(pacbio_fastq):
    assert detect_platform(pacbio_fastq).platform == "pacbio_hifi"


def test_gzipped_fastq_supported(tmp_path):
    path = write_fastq(tmp_path / "reads.fastq.gz", 100, 150, "I", gzipped=True)
    assert detect_platform(path).platform == "illumina"


def test_illumina_is_supported(illumina_fastq):
    require_supported(detect_platform(illumina_fastq), illumina_fastq)  # raise etmemeli


def test_ont_rejected_with_actionable_message(ont_fastq):
    """PLAN Kural 7: tespit etmek != desteklemek. Sessizce yanlis araçla koşulmaz."""
    info = detect_platform(ont_fastq)
    with pytest.raises(UnsupportedPlatformError) as exc:
        require_supported(info, ont_fastq)
    message = str(exc.value)
    assert "ont" in message.lower()
    assert "illumina" in message.lower()
    assert str(ont_fastq) in message


def test_pacbio_rejected(pacbio_fastq):
    with pytest.raises(UnsupportedPlatformError):
        require_supported(detect_platform(pacbio_fastq), pacbio_fastq)


def test_empty_fastq_is_unknown_and_rejected(tmp_path):
    path = tmp_path / "empty.fastq"
    path.write_text("")
    info = detect_platform(path)
    assert info.platform == "unknown"
    with pytest.raises(UnsupportedPlatformError):
        require_supported(info, path)
```

- [ ] **Step 3: Testleri çalıştır, BAŞARISIZ olduğunu gör**

```bash
conda run -n rnaforge-core python -m pytest tests/test_platform.py -v
```

Beklenen: FAIL — `ModuleNotFoundError: No module named 'rnaforge.platform'`

- [ ] **Step 4: platform.py'yi yaz**

`rnaforge/platform.py`:

```python
"""FASTQ'dan platform tespiti (saf stdlib — harici araç gerekmez).

Eşik mantığı ali-wgs-pipeline/ali_wgs/detect.py'den uyarlandı.
DİKKAT: Kütüphane kimyası (stranded / rRNA-polyA) FASTQ'da YOKTUR;
tespit edilemez, config'ten gelir (PLAN §4.1).
"""
from __future__ import annotations

import gzip
from dataclasses import dataclass
from pathlib import Path

SUPPORTED_PLATFORMS = ("illumina",)


class UnsupportedPlatformError(RuntimeError):
    """Girdi tespit edildi ama MVP'de desteklenmiyor (PLAN Kural 7)."""


@dataclass(frozen=True)
class PlatformInfo:
    platform: str  # illumina | ont | pacbio_hifi | unknown
    mean_read_length: float
    n50: int
    mean_quality: float
    n_reads_sampled: int


def _open(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt")


def _sample_fastq(path: Path, max_reads: int) -> tuple[list[int], list[float]]:
    lengths: list[int] = []
    quals: list[float] = []
    with _open(path) as fh:
        for i, line in enumerate(fh):
            position = i % 4
            if position == 1:
                lengths.append(len(line.strip()))
            elif position == 3:
                q = line.strip()
                if q:
                    quals.append(sum(ord(c) - 33 for c in q) / len(q))
                if len(lengths) >= max_reads:
                    break
    return lengths, quals


def _n50(lengths: list[int]) -> int:
    if not lengths:
        return 0
    ordered = sorted(lengths, reverse=True)
    half, acc = sum(ordered) / 2, 0
    for length in ordered:
        acc += length
        if acc >= half:
            return length
    return ordered[-1]


def detect_platform(
    fastq: Path,
    short_read_max_len: int = 350,
    hifi_min_qual: float = 25.0,
    max_reads: int = 5000,
) -> PlatformInfo:
    lengths, quals = _sample_fastq(Path(fastq), max_reads)
    if not lengths:
        return PlatformInfo("unknown", 0.0, 0, 0.0, 0)

    mean_len = sum(lengths) / len(lengths)
    mean_q = sum(quals) / len(quals) if quals else 0.0
    n50 = _n50(lengths)

    if mean_len <= short_read_max_len:
        platform = "illumina"
    elif mean_q >= hifi_min_qual and n50 >= 5000:
        platform = "pacbio_hifi"
    else:
        platform = "ont"

    return PlatformInfo(
        platform=platform,
        mean_read_length=round(mean_len, 1),
        n50=n50,
        mean_quality=round(mean_q, 1),
        n_reads_sampled=len(lengths),
    )


def require_supported(info: PlatformInfo, fastq: Path) -> None:
    """Desteklenmeyen platformu net mesajla reddet. Sessiz devam YOK."""
    if info.platform in SUPPORTED_PLATFORMS:
        return
    raise UnsupportedPlatformError(
        f"detected platform {info.platform!r} is not supported in the MVP "
        f"(supported: {', '.join(SUPPORTED_PLATFORMS)}).\n"
        f"  file: {fastq}\n"
        f"  mean read length: {info.mean_read_length}, N50: {info.n50}, "
        f"mean quality: {info.mean_quality}, reads sampled: {info.n_reads_sampled}\n"
        f"Long-read support (ONT/PacBio) needs a different route (minimap2) "
        f"and is planned for a later phase. Running the Illumina route on this "
        f"input would produce wrong results, so it is refused."
    )
```

- [ ] **Step 5: Testleri çalıştır, GEÇTİĞİNİ gör**

```bash
conda run -n rnaforge-core python -m pytest tests/test_platform.py -v
```

Beklenen: 8 passed

- [ ] **Step 6: Commit**

```bash
git add rnaforge/platform.py tests/test_platform.py tests/conftest.py
git commit -m "feat: FASTQ'dan platform tespiti; ONT/PacBio net hatayla reddedilir"
```

---

### Task 4: Metadata yükleme + design doğrulama

**Files:**
- Create: `rnaforge/metadata.py`
- Test: `tests/test_metadata.py`

**Interfaces:**
- Consumes: yok
- Produces:
  - `rnaforge.metadata.MetadataError(ValueError)`
  - `Sample(sample_id: str, condition: str, fastq_1: Path, fastq_2: Path | None, batch: str | None)` (frozen)
  - `load_metadata(path: Path | str, base_dir: Path | None = None) -> list[Sample]`
  - `design_variables(design: str) -> list[str]` — `"~batch + condition"` → `["batch", "condition"]`
  - `validate_design(samples: list[Sample], design: str) -> None`

Metadata formatı: TSV, zorunlu sütunlar `sample_id`, `condition`, `fastq_1`; opsiyonel `fastq_2`, `batch`. Göreli yollar `base_dir`'e göre çözülür (verilmezse metadata dosyasının dizini).

- [ ] **Step 1: Başarısız testleri yaz**

`tests/test_metadata.py`:

```python
from __future__ import annotations

import pytest

from rnaforge.metadata import (
    MetadataError,
    design_variables,
    load_metadata,
    validate_design,
)


def _make_fastqs(tmp_path, *names):
    for n in names:
        (tmp_path / n).write_text("@r\nACGT\n+\nIIII\n")


def _write_meta(tmp_path, body: str):
    path = tmp_path / "samples.tsv"
    path.write_text(body)
    return path


def test_loads_paired_end_samples(tmp_path):
    _make_fastqs(tmp_path, "a_R1.fastq", "a_R2.fastq", "b_R1.fastq", "b_R2.fastq")
    path = _write_meta(tmp_path, (
        "sample_id\tcondition\tfastq_1\tfastq_2\n"
        "s1\tcontrol\ta_R1.fastq\ta_R2.fastq\n"
        "s2\ttreated\tb_R1.fastq\tb_R2.fastq\n"
    ))
    samples = load_metadata(path)
    assert [s.sample_id for s in samples] == ["s1", "s2"]
    assert samples[0].fastq_2.name == "a_R2.fastq"
    assert samples[0].batch is None


def test_loads_single_end_samples(tmp_path):
    _make_fastqs(tmp_path, "a.fastq", "b.fastq")
    path = _write_meta(tmp_path, (
        "sample_id\tcondition\tfastq_1\n"
        "s1\tcontrol\ta.fastq\n"
        "s2\ttreated\tb.fastq\n"
    ))
    assert load_metadata(path)[0].fastq_2 is None


def test_missing_fastq_file_raises(tmp_path):
    path = _write_meta(tmp_path, (
        "sample_id\tcondition\tfastq_1\n"
        "s1\tcontrol\tyok.fastq\n"
    ))
    with pytest.raises(MetadataError, match="yok.fastq"):
        load_metadata(path)


def test_duplicate_sample_id_raises(tmp_path):
    _make_fastqs(tmp_path, "a.fastq", "b.fastq")
    path = _write_meta(tmp_path, (
        "sample_id\tcondition\tfastq_1\n"
        "s1\tcontrol\ta.fastq\n"
        "s1\ttreated\tb.fastq\n"
    ))
    with pytest.raises(MetadataError, match="s1"):
        load_metadata(path)


def test_missing_required_column_raises(tmp_path):
    path = _write_meta(tmp_path, "sample_id\tfastq_1\ns1\ta.fastq\n")
    with pytest.raises(MetadataError, match="condition"):
        load_metadata(path)


def test_design_variables_parses_formula():
    assert design_variables("~condition") == ["condition"]
    assert design_variables("~batch + condition") == ["batch", "condition"]


def test_design_variable_missing_from_metadata_raises(tmp_path):
    _make_fastqs(tmp_path, "a.fastq", "b.fastq")
    path = _write_meta(tmp_path, (
        "sample_id\tcondition\tfastq_1\n"
        "s1\tcontrol\ta.fastq\n"
        "s2\ttreated\tb.fastq\n"
    ))
    samples = load_metadata(path)
    with pytest.raises(MetadataError, match="batch"):
        validate_design(samples, "~batch + condition")


def test_design_requires_two_condition_levels(tmp_path):
    _make_fastqs(tmp_path, "a.fastq", "b.fastq")
    path = _write_meta(tmp_path, (
        "sample_id\tcondition\tfastq_1\n"
        "s1\tcontrol\ta.fastq\n"
        "s2\tcontrol\tb.fastq\n"
    ))
    with pytest.raises(MetadataError, match="at least 2 levels"):
        validate_design(load_metadata(path), "~condition")


def test_design_requires_replicates(tmp_path):
    """DESeq2 replika olmadan dispersiyon tahmin edemez — erken ve net uyar."""
    _make_fastqs(tmp_path, "a.fastq", "b.fastq")
    path = _write_meta(tmp_path, (
        "sample_id\tcondition\tfastq_1\n"
        "s1\tcontrol\ta.fastq\n"
        "s2\ttreated\tb.fastq\n"
    ))
    with pytest.raises(MetadataError, match="replicate"):
        validate_design(load_metadata(path), "~condition")


def test_valid_design_passes(tmp_path):
    _make_fastqs(tmp_path, "a.fastq", "b.fastq", "c.fastq", "d.fastq")
    path = _write_meta(tmp_path, (
        "sample_id\tcondition\tfastq_1\n"
        "s1\tcontrol\ta.fastq\n"
        "s2\tcontrol\tb.fastq\n"
        "s3\ttreated\tc.fastq\n"
        "s4\ttreated\td.fastq\n"
    ))
    validate_design(load_metadata(path), "~condition")  # raise etmemeli
```

- [ ] **Step 2: Testleri çalıştır, BAŞARISIZ olduğunu gör**

```bash
conda run -n rnaforge-core python -m pytest tests/test_metadata.py -v
```

Beklenen: FAIL — `ModuleNotFoundError: No module named 'rnaforge.metadata'`

- [ ] **Step 3: metadata.py'yi yaz**

`rnaforge/metadata.py`:

```python
"""Örnek metadata (TSV) yükleme ve DESeq2 design formülü doğrulama.

Hatalar burada, pipeline koşmadan ÖNCE yakalanır (PLAN §13).
"""
from __future__ import annotations

import csv
import re
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

REQUIRED_COLUMNS = ("sample_id", "condition", "fastq_1")


class MetadataError(ValueError):
    """Metadata dosyası veya design formülü geçersiz."""


@dataclass(frozen=True)
class Sample:
    sample_id: str
    condition: str
    fastq_1: Path
    fastq_2: Path | None = None
    batch: str | None = None


def _resolve(value: str, base_dir: Path) -> Path:
    path = Path(value)
    return path if path.is_absolute() else base_dir / path


def load_metadata(path: Path | str, base_dir: Path | None = None) -> list[Sample]:
    path = Path(path)
    if not path.exists():
        raise MetadataError(f"metadata file not found: {path}")
    base_dir = base_dir or path.parent

    with path.open(newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows:
        raise MetadataError(f"metadata file has no data rows: {path}")

    missing_cols = [c for c in REQUIRED_COLUMNS if c not in rows[0]]
    if missing_cols:
        raise MetadataError(
            f"metadata is missing required column(s): {', '.join(missing_cols)} "
            f"(required: {', '.join(REQUIRED_COLUMNS)})"
        )

    samples: list[Sample] = []
    for line_no, row in enumerate(rows, start=2):
        sample_id = (row.get("sample_id") or "").strip()
        condition = (row.get("condition") or "").strip()
        if not sample_id or not condition:
            raise MetadataError(f"line {line_no}: sample_id and condition must not be empty")

        fastqs: list[Path] = []
        for column in ("fastq_1", "fastq_2"):
            value = (row.get(column) or "").strip()
            if not value:
                fastqs.append(None)
                continue
            resolved = _resolve(value, base_dir)
            if not resolved.exists():
                raise MetadataError(
                    f"line {line_no}: {column} file does not exist: {resolved}"
                )
            fastqs.append(resolved)

        batch = (row.get("batch") or "").strip() or None
        samples.append(Sample(sample_id, condition, fastqs[0], fastqs[1], batch))

    duplicates = [s for s, n in Counter(x.sample_id for x in samples).items() if n > 1]
    if duplicates:
        raise MetadataError(f"duplicate sample_id(s): {', '.join(sorted(duplicates))}")
    return samples


def design_variables(design: str) -> list[str]:
    """'~batch + condition' -> ['batch', 'condition']"""
    body = design.strip().lstrip("~")
    return [v for v in (part.strip() for part in re.split(r"[+*:]", body)) if v]


def validate_design(samples: list[Sample], design: str) -> None:
    variables = design_variables(design)
    if not variables:
        raise MetadataError(f"design formula has no variables: {design!r}")

    known = {"condition", "batch"}
    unknown = [v for v in variables if v not in known]
    if unknown:
        raise MetadataError(
            f"design formula references unknown variable(s): {', '.join(unknown)}. "
            f"Metadata columns available for design: {', '.join(sorted(known))}"
        )

    if "batch" in variables and any(s.batch is None for s in samples):
        raise MetadataError(
            "design formula uses 'batch' but the metadata has no batch value for every sample. "
            "Add a 'batch' column, or use design '~condition'."
        )

    counts = Counter(s.condition for s in samples)
    if len(counts) < 2:
        raise MetadataError(
            f"condition must have at least 2 levels for differential expression, "
            f"found: {', '.join(sorted(counts))}"
        )
    without_replicates = sorted(c for c, n in counts.items() if n < 2)
    if without_replicates:
        raise MetadataError(
            f"condition level(s) without replicate: {', '.join(without_replicates)}. "
            "DESeq2 cannot estimate dispersion without replicates."
        )
```

- [ ] **Step 4: Testleri çalıştır, GEÇTİĞİNİ gör**

```bash
conda run -n rnaforge-core python -m pytest tests/test_metadata.py -v
```

Beklenen: 10 passed

- [ ] **Step 5: Commit**

```bash
git add rnaforge/metadata.py tests/test_metadata.py
git commit -m "feat: metadata yükleme + design formülü doğrulama (replika/seviye kontrolü)"
```

---

### Task 5: Run context + durum kaydı (resume + heartbeat)

**Files:**
- Create: `rnaforge/state.py`
- Test: `tests/test_state.py`

**Interfaces:**
- Consumes: yok
- Produces:
  - `rnaforge.state.RunState(run_dir: Path)` — `mark_done(module: str, outputs: list[str]) -> None`, `is_done(module: str) -> bool`, `heartbeat() -> None`, `completed_modules() -> list[str]`
  - Durum dosyası: `<run_dir>/state.json`; heartbeat dosyası: `<run_dir>/heartbeat.txt`
  - `new_run_dir(base: Path, run_id: str, now: datetime | None = None) -> Path` — `runs/<YYYYmmdd_HHMMSS>_<run_id>/`

Kapatma dayanıklılığı (PLAN §15): durum her modül bitiminde diske **atomik** yazılır (geçici dosya + `os.replace`), böylece yarıda kapanma bozuk JSON bırakmaz.

- [ ] **Step 1: Başarısız testleri yaz**

`tests/test_state.py`:

```python
from __future__ import annotations

import json
from datetime import datetime

from rnaforge.state import RunState, new_run_dir


def test_new_run_dir_has_timestamp_and_id(tmp_path):
    now = datetime(2026, 7, 16, 14, 30, 22)
    run_dir = new_run_dir(tmp_path, "demo", now=now)
    assert run_dir.name == "20260716_143022_demo"
    assert run_dir.exists()


def test_module_not_done_initially(tmp_path):
    assert RunState(tmp_path).is_done("m01_validate") is False


def test_mark_done_persists_across_instances(tmp_path):
    RunState(tmp_path).mark_done("m01_validate", ["a.json"])
    # yeni instance = süreç yeniden başladı
    assert RunState(tmp_path).is_done("m01_validate") is True


def test_completed_modules_listed_in_order(tmp_path):
    state = RunState(tmp_path)
    state.mark_done("m01_validate", [])
    state.mark_done("m02_qc", [])
    assert state.completed_modules() == ["m01_validate", "m02_qc"]


def test_state_file_is_valid_json(tmp_path):
    RunState(tmp_path).mark_done("m01_validate", ["out/a.json"])
    data = json.loads((tmp_path / "state.json").read_text())
    assert data["modules"]["m01_validate"]["outputs"] == ["out/a.json"]
    assert "completed_at" in data["modules"]["m01_validate"]


def test_heartbeat_writes_valid_timestamp(tmp_path):
    state = RunState(tmp_path)
    state.heartbeat()
    path = tmp_path / "heartbeat.txt"
    assert path.exists()
    # İçerik geçerli bir ISO timestamp olmalı; bozuksa fromisoformat raise eder.
    datetime.fromisoformat(path.read_text().strip())


def test_corrupt_state_file_does_not_crash(tmp_path):
    (tmp_path / "state.json").write_text("{ bozuk json")
    # bozuk durum = hiç ilerleme yok kabul edilir, çökme YOK
    assert RunState(tmp_path).is_done("m01_validate") is False
```

- [ ] **Step 2: Testleri çalıştır, BAŞARISIZ olduğunu gör**

```bash
conda run -n rnaforge-core python -m pytest tests/test_state.py -v
```

Beklenen: FAIL — `ModuleNotFoundError: No module named 'rnaforge.state'`

- [ ] **Step 3: state.py'yi yaz**

`rnaforge/state.py`:

```python
"""Run durumu: resume + heartbeat (PLAN §15 — kapatma dayanıklılığı).

Durum atomik yazılır (geçici dosya + os.replace): yarıda kapanma
bozuk state.json bırakmaz.
"""
from __future__ import annotations

import json
import os
from datetime import datetime
from pathlib import Path

STATE_FILE = "state.json"
HEARTBEAT_FILE = "heartbeat.txt"
HEARTBEAT_INTERVAL_SECONDS = 10


def new_run_dir(base: Path | str, run_id: str, now: datetime | None = None) -> Path:
    stamp = (now or datetime.now()).strftime("%Y%m%d_%H%M%S")
    run_dir = Path(base) / f"{stamp}_{run_id}"
    run_dir.mkdir(parents=True, exist_ok=True)
    return run_dir


class RunState:
    def __init__(self, run_dir: Path | str):
        self.run_dir = Path(run_dir)
        self.run_dir.mkdir(parents=True, exist_ok=True)

    @property
    def _path(self) -> Path:
        return self.run_dir / STATE_FILE

    def _read(self) -> dict:
        if not self._path.exists():
            return {"modules": {}}
        try:
            data = json.loads(self._path.read_text())
        except (json.JSONDecodeError, OSError):
            # Bozuk durum dosyası = ilerleme yok say. Çökmek yerine baştan koş.
            return {"modules": {}}
        data.setdefault("modules", {})
        return data

    def _write(self, data: dict) -> None:
        tmp = self._path.with_suffix(".json.tmp")
        tmp.write_text(json.dumps(data, indent=2))
        os.replace(tmp, self._path)  # atomik

    def mark_done(self, module: str, outputs: list[str]) -> None:
        data = self._read()
        data["modules"][module] = {
            "completed_at": datetime.now().isoformat(timespec="seconds"),
            "outputs": [str(o) for o in outputs],
        }
        self._write(data)

    def is_done(self, module: str) -> bool:
        return module in self._read()["modules"]

    def completed_modules(self) -> list[str]:
        return list(self._read()["modules"].keys())

    def heartbeat(self) -> None:
        (self.run_dir / HEARTBEAT_FILE).write_text(
            datetime.now().isoformat(timespec="seconds") + "\n"
        )
```

- [ ] **Step 4: Testleri çalıştır, GEÇTİĞİNİ gör**

```bash
conda run -n rnaforge-core python -m pytest tests/test_state.py -v
```

Beklenen: 7 passed

- [ ] **Step 5: Commit**

```bash
git add rnaforge/state.py tests/test_state.py
git commit -m "feat: run durumu — atomik kalıcı kayıt, resume ve heartbeat"
```

---

### Task 6: m01 modülü + `rnaforge validate` komutu

**Files:**
- Create: `rnaforge/modules/__init__.py`
- Create: `rnaforge/modules/m01_validate.py`
- Modify: `rnaforge/cli.py` (tamamen değiştirilir, kod aşağıda)
- Test: `tests/test_m01_validate.py`

**Interfaces:**
- Consumes: `load_config` (Task 2), `detect_platform`/`require_supported`/`PlatformInfo`/`UnsupportedPlatformError` (Task 3), `load_metadata`/`validate_design`/`Sample`/`MetadataError` (Task 4), `RunState`/`new_run_dir` (Task 5)
- Produces:
  - `rnaforge.modules.m01_validate.run_validation(config: Config, metadata_path: Path, run_dir: Path) -> dict` — özet sözlük döner ve `<run_dir>/logs/validation.log` + `<run_dir>/statistics/raw_statistics.json` yazar
  - `rnaforge.cli.main(argv)` → `rnaforge validate --config X --metadata Y [--run-id Z] [--runs-dir D]`; çıkış kodları: 0 başarı, 1 doğrulama hatası, 2 kullanım hatası

- [ ] **Step 1: Başarısız testleri yaz**

`tests/test_m01_validate.py`:

```python
from __future__ import annotations

import json
import textwrap

import pytest

from rnaforge.cli import main
from rnaforge.config import load_config
from rnaforge.modules.m01_validate import run_validation
from rnaforge.platform import UnsupportedPlatformError
from tests.conftest import write_fastq


def _setup(tmp_path, fastq_maker) -> tuple:
    """Geçerli config + metadata + FASTQ üretir; (config_path, metadata_path) döner."""
    (tmp_path / "ref").mkdir()
    (tmp_path / "ref" / "genome.fa").write_text(">c1\nACGT\n")
    (tmp_path / "ref" / "genes.gff").write_text("##gff-version 3\n")

    names = ["c1.fastq", "c2.fastq", "t1.fastq", "t2.fastq"]
    for n in names:
        fastq_maker(tmp_path / n)

    config_path = tmp_path / "config.yaml"
    config_path.write_text(textwrap.dedent(f"""
        organism: "Escherichia coli"
        organism_type: "prokaryote"
        reference:
          genome_fasta: "{tmp_path / 'ref' / 'genome.fa'}"
          annotation_gff: "{tmp_path / 'ref' / 'genes.gff'}"
        de:
          design: "~condition"
    """))

    metadata_path = tmp_path / "samples.tsv"
    metadata_path.write_text(
        "sample_id\tcondition\tfastq_1\n"
        "s1\tcontrol\tc1.fastq\n"
        "s2\tcontrol\tc2.fastq\n"
        "s3\ttreated\tt1.fastq\n"
        "s4\ttreated\tt2.fastq\n"
    )
    return config_path, metadata_path


def _illumina(path):
    write_fastq(path, 50, 150, "I")


def _ont(path):
    write_fastq(path, 50, (1000, 20000), "+")


def test_validation_succeeds_on_illumina(tmp_path):
    config_path, metadata_path = _setup(tmp_path, _illumina)
    run_dir = tmp_path / "run"
    summary = run_validation(load_config(config_path), metadata_path, run_dir)

    assert summary["n_samples"] == 4
    assert summary["platform"] == "illumina"
    assert summary["organism_type"] == "prokaryote"
    assert summary["conditions"] == {"control": 2, "treated": 2}


def test_validation_writes_log_and_statistics(tmp_path):
    config_path, metadata_path = _setup(tmp_path, _illumina)
    run_dir = tmp_path / "run"
    run_validation(load_config(config_path), metadata_path, run_dir)

    assert (run_dir / "logs" / "validation.log").exists()
    stats = json.loads((run_dir / "statistics" / "raw_statistics.json").read_text())
    assert stats["n_samples"] == 4
    assert len(stats["samples"]) == 4
    assert stats["samples"][0]["mean_read_length"] == pytest.approx(150.0)


def test_validation_rejects_ont(tmp_path):
    config_path, metadata_path = _setup(tmp_path, _ont)
    with pytest.raises(UnsupportedPlatformError, match="ont"):
        run_validation(load_config(config_path), metadata_path, tmp_path / "run")


def test_validation_marks_module_done_for_resume(tmp_path):
    from rnaforge.state import RunState

    config_path, metadata_path = _setup(tmp_path, _illumina)
    run_dir = tmp_path / "run"
    run_validation(load_config(config_path), metadata_path, run_dir)
    assert RunState(run_dir).is_done("m01_validate") is True


def test_cli_validate_returns_zero_on_success(tmp_path, capsys):
    config_path, metadata_path = _setup(tmp_path, _illumina)
    code = main([
        "validate",
        "--config", str(config_path),
        "--metadata", str(metadata_path),
        "--runs-dir", str(tmp_path / "runs"),
        "--run-id", "demo",
    ])
    assert code == 0
    assert "illumina" in capsys.readouterr().out


def test_cli_validate_returns_one_on_ont(tmp_path, capsys):
    config_path, metadata_path = _setup(tmp_path, _ont)
    code = main([
        "validate",
        "--config", str(config_path),
        "--metadata", str(metadata_path),
        "--runs-dir", str(tmp_path / "runs"),
        "--run-id", "demo",
    ])
    assert code == 1
    assert "not supported" in capsys.readouterr().err
```

- [ ] **Step 2: Testleri çalıştır, BAŞARISIZ olduğunu gör**

```bash
conda run -n rnaforge-core python -m pytest tests/test_m01_validate.py -v
```

Beklenen: FAIL — `ModuleNotFoundError: No module named 'rnaforge.modules'`

- [ ] **Step 3: m01 modülünü yaz**

`rnaforge/modules/__init__.py`:

```python
"""Pipeline modülleri. Sıra ve sözleşmeler: PLAN.md Bölüm 5."""
```

`rnaforge/modules/m01_validate.py`:

```python
"""m01 — Girdi doğrulama + platform tespiti.

Bu modül pipeline'ın kapısıdır: config, metadata ve FASTQ'lar burada
doğrulanır. Hata varsa BURADA durulur, sessiz devam yoktur (PLAN §13, Kural 7).
"""
from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

from rnaforge.config import Config
from rnaforge.metadata import Sample, load_metadata, validate_design
from rnaforge.platform import PlatformInfo, detect_platform, require_supported
from rnaforge.state import RunState

MODULE_NAME = "m01_validate"


def _check_reference(config: Config) -> None:
    fields = {
        "prokaryote": ("genome_fasta", "annotation_gff"),
        "eukaryote": ("transcriptome_fasta", "tx2gene"),
    }[config.organism_type]
    for field in fields:
        path = getattr(config.reference, field)
        if not path.exists():
            raise FileNotFoundError(
                f"reference.{field} does not exist: {path} "
                f"(required for organism_type={config.organism_type})"
            )


def _sample_stats(sample: Sample, info: PlatformInfo) -> dict:
    return {
        "sample_id": sample.sample_id,
        "condition": sample.condition,
        "batch": sample.batch,
        "paired": sample.fastq_2 is not None,
        "platform": info.platform,
        "mean_read_length": info.mean_read_length,
        "mean_quality": info.mean_quality,
        "n_reads_sampled": info.n_reads_sampled,
    }


def run_validation(config: Config, metadata_path: Path, run_dir: Path) -> dict:
    run_dir = Path(run_dir)
    logs_dir = run_dir / "logs"
    stats_dir = run_dir / "statistics"
    logs_dir.mkdir(parents=True, exist_ok=True)
    stats_dir.mkdir(parents=True, exist_ok=True)
    state = RunState(run_dir)
    lines: list[str] = []

    def log(message: str) -> None:
        lines.append(message)

    log(f"organism={config.organism} organism_type={config.organism_type}")
    _check_reference(config)
    log("reference files: OK")

    samples = load_metadata(metadata_path)
    log(f"metadata: {len(samples)} sample(s) loaded from {metadata_path}")

    validate_design(samples, config.de.design)
    log(f"design formula {config.de.design!r}: OK")

    per_sample: list[dict] = []
    platforms: set[str] = set()
    for sample in samples:
        state.heartbeat()
        info = detect_platform(sample.fastq_1)
        require_supported(info, sample.fastq_1)  # desteklenmiyorsa BURADA durur
        platforms.add(info.platform)
        per_sample.append(_sample_stats(sample, info))
        log(f"{sample.sample_id}: platform={info.platform} "
            f"mean_read_length={info.mean_read_length}")

    if len(platforms) > 1:
        raise ValueError(
            f"samples come from mixed platforms: {', '.join(sorted(platforms))}. "
            "A single run must use one platform."
        )
    platform = platforms.pop()

    if config.platform != "auto" and config.platform != platform:
        raise ValueError(
            f"config says platform={config.platform!r} but the FASTQ files look like "
            f"{platform!r}. Fix the config, or set platform: auto."
        )

    conditions = dict(Counter(s.condition for s in samples))
    summary = {
        "organism": config.organism,
        "organism_type": config.organism_type,
        "platform": platform,
        "n_samples": len(samples),
        "conditions": conditions,
        "design": config.de.design,
        "samples": per_sample,
    }

    stats_path = stats_dir / "raw_statistics.json"
    stats_path.write_text(json.dumps(summary, indent=2))
    log(f"raw statistics written: {stats_path}")
    (logs_dir / "validation.log").write_text("\n".join(lines) + "\n")

    state.mark_done(MODULE_NAME, [str(stats_path)])
    return summary
```

- [ ] **Step 4: CLI'yi yaz (Task 1'deki içerik tamamen değişir)**

`rnaforge/cli.py`:

```python
"""Komut satırı girişi."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from rnaforge import __version__
from rnaforge.config import ConfigError, load_config
from rnaforge.metadata import MetadataError
from rnaforge.modules.m01_validate import run_validation
from rnaforge.platform import UnsupportedPlatformError
from rnaforge.state import new_run_dir


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="rnaforge", description="Bulk RNA-seq pipeline")
    parser.add_argument("--version", action="version", version=__version__)
    sub = parser.add_subparsers(dest="command")

    validate = sub.add_parser("validate", help="validate config, metadata and FASTQ inputs")
    validate.add_argument("--config", required=True, type=Path)
    validate.add_argument("--metadata", required=True, type=Path)
    validate.add_argument("--runs-dir", type=Path, default=Path("runs"))
    validate.add_argument("--run-id", default="run")
    return parser


def _cmd_validate(args) -> int:
    config = load_config(args.config)
    run_dir = new_run_dir(args.runs_dir, args.run_id)
    summary = run_validation(config, args.metadata, run_dir)
    print(
        f"validation OK: {summary['n_samples']} sample(s), "
        f"platform={summary['platform']}, "
        f"organism_type={summary['organism_type']}, "
        f"conditions={summary['conditions']}"
    )
    print(f"run directory: {run_dir}")
    return 0


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command is None:
        print("error: no command given (try: rnaforge validate --help)")
        return 2
    try:
        return _cmd_validate(args)
    except (ConfigError, MetadataError, UnsupportedPlatformError,
            FileNotFoundError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
```

- [ ] **Step 5: Testleri çalıştır, GEÇTİĞİNİ gör**

```bash
conda run -n rnaforge-core python -m pytest tests/test_m01_validate.py -v
```

Beklenen: 6 passed

- [ ] **Step 6: TÜM testleri çalıştır**

```bash
conda run -n rnaforge-core python -m pytest -v
```

Beklenen: 42 passed (2 CLI + 9 config + 8 platform + 10 metadata + 7 state + 6 m01)

Not: `tests/test_cli.py::test_no_command_returns_usage_error` hâlâ geçmelidir; `--version` testi de öyle.

- [ ] **Step 7: Commit**

```bash
git add rnaforge/modules/ rnaforge/cli.py tests/test_m01_validate.py
git commit -m "feat: m01 doğrulama modülü + rnaforge validate komutu"
```

---

### Task 7: README (TR/EN) + Plan A kapanışı

**Files:**
- Create: `README.md`
- Create: `README.tr.md`
- Modify: `DURUM.md`

**Interfaces:**
- Consumes: `rnaforge validate` CLI (Task 6)
- Produces: yok (dokümantasyon)

- [ ] **Step 1: README.md (EN) yaz**

````markdown
# RNAForge

Reproducible, modular Bulk RNA-seq analysis pipeline.

Turkish version: [README.tr.md](README.tr.md) · Reference document: [PLAN.md](PLAN.md) (v1.2)

## Status

Early development. Currently implemented: input validation and platform detection (`m01`).

## Install

```bash
conda env create -f envs/rnaforge-core.yml
conda activate rnaforge-core
pip install -e .
```

## Usage

```bash
rnaforge validate --config config/config.yaml --metadata samples.tsv --run-id demo
```

### Metadata format (TSV)

| Column | Required | Description |
|---|---|---|
| `sample_id` | yes | Unique sample identifier |
| `condition` | yes | Experimental group; needs ≥2 levels and ≥2 replicates each |
| `fastq_1` | yes | Path to R1 (or single-end reads) |
| `fastq_2` | no | Path to R2 for paired-end |
| `batch` | no | Batch/covariate; required if the design formula uses `batch` |

## Key design decisions

- **`organism_type` is required and has no default** (`prokaryote` | `eukaryote`). It routes
  quantification: prokaryote uses genome alignment + featureCounts, eukaryote uses
  Salmon + tximport. Both converge on the same gene × sample count matrix.
- **Illumina only (MVP).** ONT/PacBio inputs are detected and refused with a clear error
  rather than silently processed through the wrong route.
- **Trimming is deliberately gentle.** Aggressive quality trimming distorts expression
  estimates ([Williams et al. 2016](https://doi.org/10.1186/s12859-016-0956-2)); a minimum
  length filter is what prevents the distortion.

## Development

```bash
conda run -n rnaforge-core python -m pytest -v
```

## Privacy

Customer data is never committed. `runs/`, `raw/` and `references/` are git-ignored.
````

- [ ] **Step 2: README.tr.md (TR) yaz — içerik eşdeğer**

````markdown
# RNAForge

Yeniden üretilebilir, modüler Bulk RNA-seq analiz pipeline'ı.

İngilizce sürüm: [README.md](README.md) · Referans doküman: [PLAN.md](PLAN.md) (v1.2)

## Durum

Erken geliştirme. Şu an hazır olan: girdi doğrulama ve platform tespiti (`m01`).

## Kurulum

```bash
conda env create -f envs/rnaforge-core.yml
conda activate rnaforge-core
pip install -e .
```

## Kullanım

```bash
rnaforge validate --config config/config.yaml --metadata samples.tsv --run-id demo
```

### Metadata formatı (TSV)

| Sütun | Zorunlu | Açıklama |
|---|---|---|
| `sample_id` | evet | Benzersiz örnek kimliği |
| `condition` | evet | Deney grubu; ≥2 seviye ve her seviyede ≥2 replika gerekir |
| `fastq_1` | evet | R1 yolu (veya single-end okumalar) |
| `fastq_2` | hayır | Paired-end için R2 yolu |
| `batch` | hayır | Batch/kovaryat; design formülü `batch` kullanıyorsa zorunlu |

## Temel tasarım kararları

- **`organism_type` zorunludur, varsayılanı yoktur** (`prokaryote` | `eukaryote`).
  Kantifikasyonu yönlendirir: prokaryotta genom hizalama + featureCounts, ökaryotta
  Salmon + tximport. İkisi de aynı gen × örnek count matrisinde buluşur.
- **Yalnızca Illumina (MVP).** ONT/PacBio girdileri tespit edilir ve sessizce yanlış
  yoldan işlenmek yerine net bir hatayla reddedilir.
- **Trimming bilinçli olarak naziktir.** Agresif kalite trimming ekspresyon tahminlerini
  bozar ([Williams et al. 2016](https://doi.org/10.1186/s12859-016-0956-2)); sapmayı
  engelleyen şey minimum uzunluk filtresidir.

## Geliştirme

```bash
conda run -n rnaforge-core python -m pytest -v
```

## Gizlilik

Müşteri verisi asla commit edilmez. `runs/`, `raw/` ve `references/` git tarafından yok sayılır.
````

- [ ] **Step 3: DURUM.md'yi güncelle**

`DURUM.md` içindeki "Şu an nerede kaldık" bölümünü şununla değiştir:

```markdown
## Şu an nerede kaldık
- **Plan A BİTTİ (2026-07-16).** `rnaforge validate` çalışıyor: config + metadata +
  platform tespiti, ONT/PacBio reddi, resume/heartbeat altyapısı. Tüm testler geçiyor.
- Sıradaki: **Plan B** (m02 FastQC → m03 fastp → m04 quant router → m05 count matrisi).
  Plan B conda araç ortamlarını gerektirir: fastqc, fastp, salmon, bowtie2 KURULU DEĞİL
  (featureCounts ve samtools sistemde var).
```

- [ ] **Step 4: Tüm testleri son kez çalıştır**

```bash
conda run -n rnaforge-core python -m pytest -v
```

Beklenen: 42 passed

- [ ] **Step 5: Commit ve push**

```bash
git add README.md README.tr.md DURUM.md
git commit -m "docs: çift dilli README (TR/EN) + DURUM güncellemesi — Plan A tamam"
git push origin main
```

---

## Plan A tamamlanma ölçütü

- [ ] `rnaforge validate --config config/config.yaml --metadata samples.tsv` sentetik
      Illumina girdisinde 0 döner ve `runs/<ts>_<id>/` altına `logs/validation.log` +
      `statistics/raw_statistics.json` yazar.
- [ ] Aynı komut ONT girdisinde 1 döner ve stderr'de "not supported" içeren, dosya adı
      ve ölçülen değerleri gösteren net bir mesaj basar.
- [ ] `organism_type` eksik config reddedilir.
- [ ] `trimming.aggressive_quality` varsayılanı `False` (testle sabitlenmiş).
- [ ] 42 test geçer.

## Sonraki planlar

- **Plan B — QC→Count:** conda araç ortamları (fastqc/fastp/salmon/bowtie2), m02, m03,
  m04 (yönlendirici), m05. Çıktı: `quantification/counts.tsv` (gen × örnek sözleşmesi).
- **Plan C — DE→Rapor:** m06 (DESeq2 Rscript köprüsü), m07 (PCA/Volcano/Heatmap),
  m08 (HTML rapor), pydeseq2 çapraz kontrolü.
