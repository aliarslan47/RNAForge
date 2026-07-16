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
