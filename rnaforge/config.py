"""Config yükleme ve şema doğrulama. Geçersiz config m01'den ÖNCE yakalanır."""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import yaml

ORGANISM_TYPES = ("prokaryote", "eukaryote")
PLATFORMS = ("auto", "illumina", "ont", "pacbio_hifi")
STRANDEDNESS = ("unstranded", "stranded", "reverse")
SELECTIONS = ("rrna_depletion", "polya")
CHEMISTRY = ("cdna", "direct_rna")
REPORT_LANGUAGES = ("tr", "en")

# İzin verilen üst seviye anahtarlar. Bunun DIŞINDA bir anahtar (ör. `design:`
# doğrusu `de.design`, ya da `refernce` yazım hatası) SESSİZCE yutulmamalı:
# kullanıcı config'ini değiştirdiğini sanıp eski varsayılanla koşar → makul
# görünen sahte sonuç. Yeni üst anahtar eklerken buraya da ekle.
KNOWN_TOP_LEVEL_KEYS = frozenset({
    "organism", "organism_type", "platform", "reference", "library",
    "trimming", "de", "report", "resources", "paired", "quality",
    "quantification", "enrichment", "amr", "operon", "ppi", "basecall",
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
class Basecall:
    # ONT ham-sinyal (FAST5/POD5) → FASTQ (m00). dorado GPU zorunlu (CPU pratik değil).
    dorado_bin: str = "dorado"      # PATH'te değilse mutlak yol verilir
    model: str = "hac"              # dorado kompleks; tam modeli POD5 metadata'sından seçer
    device: str = "cuda:all"
    env: str = "rnaforge-basecall"  # pod5 dönüşümü için
    models_dir: str | None = None   # model önbelleği (yeniden indirmeyi önler)


@dataclass(frozen=True)
class Library:
    strandedness: str = "unstranded"
    selection: str = "rrna_depletion"
    # Long-read only: cDNA needs Pychopper full-length orientation; direct-RNA does not.
    # NOT detectable from FASTQ (spec 2026-08-05). None = unset (fine for short reads).
    chemistry: str | None = None
    # Long-read cDNA only: True (default) = ONT-kit cDNA with SSP/VNP strand-switch
    # primers -> Pychopper orients/trims full-length. False = non-ONT-kit cDNA
    # (random-primed dscDNA, native barcoding) that has NO strand-switch primers ->
    # Pychopper would discard ~all reads, so skip it and run chopper-only. Ignored for
    # direct_rna (already chopper-only) and short reads. NOT detectable from FASTQ.
    full_length_cdna: bool = True


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
    # Çok-seviyeli condition için açık kontrastlar: her biri (test_seviyesi, referans_seviyesi).
    # Boş → DESeq2 varsayılanı (son-vs-ilk seviye), 2-seviyeli koşular için yeterli.
    # 3+ gruplu koşu (kontrol/düşük/yüksek, zaman serisi) tüm karşılaştırmaları buradan ister.
    contrasts: tuple = ()


@dataclass(frozen=True)
class Quantification:
    # featureCounts (prokaryot) parametreleri. Anotasyon kaynağı değişir → config-driven.
    # Prokaryot GFF3 kullanıcısı tipik CDS/locus_tag ile ezer.
    feature_type: str = "exon"
    attribute: str = "gene_id"


@dataclass(frozen=True)
class Enrichment:
    # m09 GO ORA. min_term_size: gürültü filtresi (arka planda < bu kadar genli terim atlanır).
    # top_n: figürde gösterilen zenginleşmiş terim sayısı. obo/gaf: referans yolları (yoksa
    # modül yüksek sesle hata; gaf None ise GAF doldurma atlanır — GFF+propagation yeterli).
    min_term_size: int = 3
    top_n: int = 15
    obo: Path | None = None
    gaf: Path | None = None
    # m10 KEGG: organizma kodu (eco/hsa/mmu…) ZORUNLU (yoksa m10 çalışmaz, net hata).
    # kegg_dir: KEGG REST dosyalarının bulunduğu dizin (yoksa references/kegg/<org>).
    kegg_organism: str | None = None
    kegg_dir: Path | None = None
    # m11 GSEA (fgsea): gen-seti boyut eşikleri (fgsea minSize/maxSize). Eşik = veri/config ilkesi.
    gsea_min_size: int = 15
    gsea_max_size: int = 500
    # m12 REVIGO: semantik benzerlik eşiği (Lin, [0,1]); üstündeki terimler tek temsilcide toplanır.
    revigo_similarity: float = 0.7


@dataclass(frozen=True)
class AMR:
    # m13 AMR/virülans (abricate). amr_db: direnç DB (card/ncbi/resfinder); virulence_db: vfdb.
    # env: izole abricate ortamı; min_identity/coverage: abricate hit eşikleri (%).
    amr_db: str = "card"
    virulence_db: str = "vfdb"
    env: str = "rnaforge-amr"
    min_identity: float = 80.0
    min_coverage: float = 80.0
    # İkinci AMR aracı AMRFinderPlus (yan-yana konkordans). organism verilirse koşar (ör. "Escherichia").
    amrfinder_organism: str | None = None
    amrfinder_env: str = "ali-amrfinder"


@dataclass(frozen=True)
class Operon:
    # m14 operon tahmini: aynı yönde bitişik + intergenik gap ≤ max_gap (bp) genler aynı operon.
    max_gap: int = 50


@dataclass(frozen=True)
class PPI:
    # m15 STRING PPI + community. taxid: STRING organizma (ör. "511145" E.coli K-12); ZORUNLU (yoksa çalışmaz).
    # min_score: STRING combined_score eşiği (700=yüksek); min_community_size: rapora giren en küçük modül.
    taxid: str | None = None
    string_dir: Path | None = None
    min_score: int = 700
    min_community_size: int = 3


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
    enrichment: Enrichment = field(default_factory=Enrichment)
    amr: AMR = field(default_factory=AMR)
    operon: Operon = field(default_factory=Operon)
    ppi: PPI = field(default_factory=PPI)
    basecall: Basecall = field(default_factory=Basecall)


def _one_of(value, allowed, field: str):
    if value not in allowed:
        raise ConfigError(f"{field}: {value!r} is invalid; allowed: {', '.join(allowed)}")
    return value


# Bölüm-içi bilinen anahtarlar. Top-level koruması (KNOWN_TOP_LEVEL_KEYS) gibi:
# bölüm içindeki yazım hatası (ör. de.fdr_treshold) sessizce default kullanmasın.
# quality HARİÇ: gate-override adları serbest-biçimdir (değişken), sabit sete bağlanamaz.
_KNOWN_SECTION_KEYS = {
    "reference": {"genome_fasta", "annotation_gff", "transcriptome_fasta", "tx2gene"},
    "library": {"strandedness", "selection", "chemistry", "full_length_cdna"},
    "trimming": {"min_length", "aggressive_quality"},
    "de": {"design", "fdr_threshold", "log2fc_threshold", "reference", "contrasts"},
    "report": {"language"},
    "resources": {"threads", "memory_gb"},
    "quantification": {"feature_type", "attribute"},
    "enrichment": {"min_term_size", "top_n", "obo", "gaf", "kegg_organism", "kegg_dir",
                   "gsea_min_size", "gsea_max_size", "revigo_similarity"},
    "amr": {"amr_db", "virulence_db", "env", "min_identity", "min_coverage",
            "amrfinder_organism", "amrfinder_env"},
    "operon": {"max_gap"},
    "ppi": {"taxid", "string_dir", "min_score", "min_community_size"},
    "basecall": {"dorado_bin", "model", "device", "env", "models_dir"},
}


def _section(raw: dict, field: str) -> dict:
    """Bir config bölümünü mapping olarak al. `library: "foo"` ham AttributeError
    yerine ConfigError vermeli: hata mesajı kullanıcıya ne yapacağını söylemeli.
    Bilinen-anahtar seti olan bölümlerde yazım hatası (bilinmeyen anahtar) reddedilir."""
    value = raw.get(field)
    if value is None:
        return {}
    if not isinstance(value, dict):
        raise ConfigError(
            f"{field} must be a mapping of settings, got {type(value).__name__} ({value!r})"
        )
    known = _KNOWN_SECTION_KEYS.get(field)
    if known is not None:
        unknown = [k for k in value if k not in known]
        if unknown:
            raise ConfigError(
                f"unknown key(s) in '{field}' section: "
                f"{', '.join(map(repr, sorted(map(str, unknown))))}. "
                f"Allowed: {', '.join(sorted(known))}. "
                "A typo'd key is silently ignored otherwise, so the run would use "
                "stale defaults."
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


def _as_bool(value, field: str) -> bool:
    if isinstance(value, bool):
        return value
    if isinstance(value, str) and value.strip().lower() in {"true", "false"}:
        return value.strip().lower() == "true"
    raise ConfigError(f"{field}: expected a boolean (true/false), got {value!r}")


def _build_contrasts(raw) -> tuple:
    """de.contrasts'ı doğrulayıp ((test, ref), ...) döndürür. Boş/None → ()."""
    if raw is None:
        return ()
    if not isinstance(raw, list):
        raise ConfigError(
            f"de.contrasts must be a list of [test, reference] pairs, got {type(raw).__name__}"
        )
    out: list[tuple[str, str]] = []
    for i, pair in enumerate(raw):
        if not isinstance(pair, (list, tuple)) or len(pair) != 2:
            raise ConfigError(
                f"de.contrasts[{i}] must be a [test, reference] pair, got {pair!r}"
            )
        test, ref = str(pair[0]).strip(), str(pair[1]).strip()
        for lvl in (test, ref):
            if not lvl:
                raise ConfigError(f"de.contrasts[{i}] has an empty condition level")
            if ":" in lvl or ";" in lvl:
                raise ConfigError(
                    f"de.contrasts[{i}] level {lvl!r} must not contain ':' or ';' "
                    "(used as internal delimiters)"
                )
        if test == ref:
            raise ConfigError(
                f"de.contrasts[{i}] test and reference are identical ({test!r}); "
                "a contrast compares two different condition levels"
            )
        out.append((test, ref))
    return tuple(out)


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
    enrichment_raw = _section(raw, "enrichment")
    amr_raw = _section(raw, "amr")
    operon_raw = _section(raw, "operon")
    ppi_raw = _section(raw, "ppi")
    basecall_raw = _section(raw, "basecall")

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
            chemistry=(
                _one_of(library_raw["chemistry"], CHEMISTRY, "library.chemistry")
                if library_raw.get("chemistry") is not None
                else None
            ),
            full_length_cdna=(
                _as_bool(library_raw["full_length_cdna"], "library.full_length_cdna")
                if library_raw.get("full_length_cdna") is not None
                else True
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
            contrasts=_build_contrasts(de_raw.get("contrasts")),
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
        enrichment=Enrichment(
            min_term_size=_as_int(enrichment_raw.get("min_term_size", 3), "enrichment.min_term_size"),
            top_n=_as_int(enrichment_raw.get("top_n", 15), "enrichment.top_n"),
            obo=(Path(enrichment_raw["obo"]) if enrichment_raw.get("obo") else None),
            gaf=(Path(enrichment_raw["gaf"]) if enrichment_raw.get("gaf") else None),
            kegg_organism=(str(enrichment_raw["kegg_organism"])
                           if enrichment_raw.get("kegg_organism") else None),
            kegg_dir=(Path(enrichment_raw["kegg_dir"]) if enrichment_raw.get("kegg_dir") else None),
            gsea_min_size=_as_int(enrichment_raw.get("gsea_min_size", 15), "enrichment.gsea_min_size"),
            gsea_max_size=_as_int(enrichment_raw.get("gsea_max_size", 500), "enrichment.gsea_max_size"),
            revigo_similarity=_as_float(
                enrichment_raw.get("revigo_similarity", 0.7), "enrichment.revigo_similarity"),
        ),
        amr=AMR(
            amr_db=str(amr_raw.get("amr_db", "card")),
            virulence_db=str(amr_raw.get("virulence_db", "vfdb")),
            env=str(amr_raw.get("env", "rnaforge-amr")),
            min_identity=_as_float(amr_raw.get("min_identity", 80.0), "amr.min_identity"),
            min_coverage=_as_float(amr_raw.get("min_coverage", 80.0), "amr.min_coverage"),
            amrfinder_organism=(str(amr_raw["amrfinder_organism"])
                                if amr_raw.get("amrfinder_organism") else None),
            amrfinder_env=str(amr_raw.get("amrfinder_env", "ali-amrfinder")),
        ),
        operon=Operon(max_gap=_as_int(operon_raw.get("max_gap", 50), "operon.max_gap")),
        basecall=Basecall(
            dorado_bin=str(basecall_raw.get("dorado_bin", "dorado")),
            model=str(basecall_raw.get("model", "hac")),
            device=str(basecall_raw.get("device", "cuda:all")),
            env=str(basecall_raw.get("env", "rnaforge-basecall")),
            models_dir=(str(basecall_raw["models_dir"])
                        if basecall_raw.get("models_dir") else None),
        ),
        ppi=PPI(
            taxid=(str(ppi_raw["taxid"]) if ppi_raw.get("taxid") else None),
            string_dir=(Path(ppi_raw["string_dir"]) if ppi_raw.get("string_dir") else None),
            min_score=_as_int(ppi_raw.get("min_score", 700), "ppi.min_score"),
            min_community_size=_as_int(
                ppi_raw.get("min_community_size", 3), "ppi.min_community_size"),
        ),
    )
