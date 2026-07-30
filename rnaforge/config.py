"""Config yükleme ve şema doğrulama. Geçersiz config m01'den ÖNCE yakalanır."""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import yaml

ORGANISM_TYPES = ("prokaryote", "eukaryote")
PLATFORMS = ("auto", "illumina")
STRANDEDNESS = ("unstranded", "stranded", "reverse")
SELECTIONS = ("rrna_depletion", "polya")
REPORT_LANGUAGES = ("tr", "en")

# İzin verilen üst seviye anahtarlar. Bunun DIŞINDA bir anahtar (ör. `design:`
# doğrusu `de.design`, ya da `refernce` yazım hatası) SESSİZCE yutulmamalı:
# kullanıcı config'ini değiştirdiğini sanıp eski varsayılanla koşar → makul
# görünen sahte sonuç. Yeni üst anahtar eklerken buraya da ekle.
KNOWN_TOP_LEVEL_KEYS = frozenset({
    "organism", "organism_type", "platform", "reference", "library",
    "trimming", "de", "report", "resources", "paired", "quality",
    "quantification",
})

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
    reference: str | None = None


@dataclass(frozen=True)
class Quantification:
    # featureCounts (prokaryot) parametreleri. Anotasyon kaynağı değişir → config-driven.
    # Prokaryot GFF3 kullanıcısı tipik CDS/locus_tag ile ezer.
    feature_type: str = "exon"
    attribute: str = "gene_id"


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
    paired: bool | None = None
    quality: dict = field(default_factory=dict)
    quantification: Quantification = field(default_factory=Quantification)


def _one_of(value, allowed, field: str):
    if value not in allowed:
        raise ConfigError(f"{field}: {value!r} is invalid; allowed: {', '.join(allowed)}")
    return value


def _section(raw: dict, field: str) -> dict:
    """Bir config bölümünü mapping olarak al. `library: "foo"` ham AttributeError
    yerine ConfigError vermeli: hata mesajı kullanıcıya ne yapacağını söylemeli."""
    value = raw.get(field)
    if value is None:
        return {}
    if not isinstance(value, dict):
        raise ConfigError(
            f"{field} must be a mapping of settings, got {type(value).__name__} ({value!r})"
        )
    return value


def _as_int(value, field: str) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        raise ConfigError(f"{field}: expected an integer, got {value!r}") from None


def _as_float(value, field: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        raise ConfigError(f"{field}: expected a number, got {value!r}") from None


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

    unknown = [k for k in raw if k not in KNOWN_TOP_LEVEL_KEYS]
    if unknown:
        raise ConfigError(
            f"unknown top-level config key(s): {', '.join(map(repr, sorted(map(str, unknown))))}. "
            f"Allowed: {', '.join(sorted(KNOWN_TOP_LEVEL_KEYS))}. "
            "A misplaced key (e.g. top-level 'design' instead of 'de.design') is "
            "silently ignored otherwise, so the run would use stale defaults."
        )

    organism = raw.get("organism")
    if not organism:
        raise ConfigError("organism is required")

    organism_type = _require_organism_type(raw)
    library_raw = _section(raw, "library")
    trimming_raw = _section(raw, "trimming")
    de_raw = _section(raw, "de")
    report_raw = _section(raw, "report")
    resources_raw = _section(raw, "resources")
    quantification_raw = _section(raw, "quantification")

    trimming = Trimming(
        min_length=_as_int(trimming_raw.get("min_length", 36), "trimming.min_length"),
        aggressive_quality=bool(trimming_raw.get("aggressive_quality", False)),
    )
    if trimming.min_length < 1:
        raise ConfigError(f"trimming.min_length must be >= 1, got {trimming.min_length}")

    return Config(
        organism=str(organism),
        organism_type=organism_type,
        platform=_one_of(raw.get("platform", "auto"), PLATFORMS, "platform"),
        reference=_build_reference(_section(raw, "reference"), organism_type),
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
            fdr_threshold=_as_float(de_raw.get("fdr_threshold", 0.05), "de.fdr_threshold"),
            log2fc_threshold=_as_float(
                de_raw.get("log2fc_threshold", 1.0), "de.log2fc_threshold"
            ),
            reference=(str(de_raw["reference"]) if de_raw.get("reference") else None),
        ),
        report=Report(
            language=_one_of(report_raw.get("language", "tr"), REPORT_LANGUAGES, "report.language")
        ),
        resources=Resources(
            threads=_as_int(resources_raw.get("threads", 8), "resources.threads"),
            memory_gb=_as_int(resources_raw.get("memory_gb", 32), "resources.memory_gb"),
        ),
        paired=None if raw.get("paired") is None else bool(raw.get("paired")),
        quality=_section(raw, "quality"),
        quantification=Quantification(
            feature_type=str(quantification_raw.get("feature_type", "exon")),
            attribute=str(quantification_raw.get("attribute", "gene_id")),
        ),
    )
