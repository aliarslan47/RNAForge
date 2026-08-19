"""Örnek metadata (TSV) yükleme ve DESeq2 design formülü doğrulama.

Hatalar burada, pipeline koşmadan ÖNCE yakalanır (PLAN §13).
"""
from __future__ import annotations

import csv
import re
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Mapping

from rnaforge.gates import FAIL, PASS, GateResult

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
    subject: str | None = None
    # Çekirdek-dışı metadata sütunları (sex, lane, RIN, genotype...). Design'a keyfi
    # adlı kovaryat olarak katılabilir; condition ANA tetkik faktörü olarak kalır (Faz 3).
    covariates: Mapping[str, str] = field(default_factory=dict)


# batch/subject'in gate'de özel işlemesi var; bunlar + kimlik/fastq sütunları
# covariates'e GİRMEZ, geri kalan her sütun kovaryat olur.
_CORE_COLUMNS = frozenset(
    {"sample_id", "condition", "fastq_1", "fastq_2", "batch", "subject"}
)


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
                if column == "fastq_1":
                    raise MetadataError(f"line {line_no}: fastq_1 must not be empty")
                fastqs.append(None)
                continue
            resolved = _resolve(value, base_dir)
            if not resolved.exists():
                raise MetadataError(
                    f"line {line_no}: {column} file does not exist: {resolved}"
                )
            fastqs.append(resolved)

        batch = (row.get("batch") or "").strip() or None
        subject = (row.get("subject") or "").strip() or None
        covariates = {
            col: (val or "").strip()
            for col, val in row.items()
            if col is not None and col not in _CORE_COLUMNS
        }
        samples.append(
            Sample(sample_id, condition, fastqs[0], fastqs[1], batch, subject, covariates)
        )

    duplicates = [s for s, n in Counter(x.sample_id for x in samples).items() if n > 1]
    if duplicates:
        raise MetadataError(f"duplicate sample_id(s): {', '.join(sorted(duplicates))}")
    return samples


def design_variables(design: str) -> list[str]:
    """'~batch + condition' -> ['batch', 'condition']"""
    body = design.strip().lstrip("~")
    return [v for v in (part.strip() for part in re.split(r"[+*:]", body)) if v]


MODULE = "m01"


def validate_design(
    samples: list[Sample], design: str, paired: bool | None = None
) -> list[GateResult]:
    """Tasarım kapılarını döndürür. Bozuk FORMÜL hâlâ MetadataError'dır.

    Kapı döndürmenin sebebi: exception fırlatan kontrol gates.json'a yazılamaz,
    dolayısıyla FAIL anında teşhis raporunun gösterecek verisi olmaz (spec §3.5).
    Zorlama çağırana aittir (raise_if_failed).
    """
    variables = design_variables(design)
    if not variables:
        raise MetadataError(f"design formula has no variables: {design!r}")

    # 'condition' ANA tetkik faktörü (her zaman var). 'batch'/'subject' özel gate
    # işlemeli çekirdek faktörler. Geri kalan HER değişken keyfi bir kovaryat
    # sütunudur ve TÜM örneklerde değeri olmalı (yoksa DESeq2 sessizce/kriptik çöker).
    for variable in variables:
        if variable in ("condition", "batch", "subject"):
            if variable != "condition" and any(
                getattr(s, variable) is None for s in samples
            ):
                raise MetadataError(
                    f"design formula uses {variable!r} but the metadata has no {variable} "
                    f"value for every sample. Add a {variable!r} column, or drop it from "
                    "the design."
                )
            continue
        if any(not (s.covariates.get(variable) or "").strip() for s in samples):
            raise MetadataError(
                f"design formula references {variable!r} but the metadata has no "
                f"{variable!r} column with a value for every sample. Add a {variable!r} "
                "column, or drop it from the design."
            )

    return [
        _rank_gate(samples, variables),
        _replication_gate(samples),
        _paired_gate(samples, variables, paired),
    ]


def validate_read_layout(samples: list[Sample]) -> GateResult:
    """Tüm örnekler aynı düzende (hepsi tek-uçlu VEYA hepsi çift-uçlu) mi?

    Karışık SE/PE tek koşuda desteklenmez: m05 tek global `paired` bayrağını TÜM
    BAM'lere uyguluyordu → prokaryot featureCounts paired-mode'u SE BAM'ler üzerinde
    çalıştırıp sessizce yanlış sayardı. Tahmin/otomatik-onarım yerine yüksek sesle
    FAIL: kullanıcı ayrı koşular çalıştırmalı (Kural 7, gürültülü hata)."""
    paired = [s for s in samples if s.fastq_2 is not None]
    single = [s for s in samples if s.fastq_2 is None]
    if paired and single:
        # Azınlıktaki grup "sorumlu" örnekler olarak adlandırılır (teşhis raporu render eder).
        offending = single if len(single) <= len(paired) else paired
        return GateResult(
            name="read_layout", module=MODULE, status=FAIL,
            message=(
                f"samples mix single-end and paired-end layouts in one run "
                f"({len(paired)} paired, {len(single)} single-end); this is not "
                "supported in a single run and would silently miscount"
            ),
            remedy=(
                "run single-end and paired-end samples as separate runs (separate "
                "metadata sheets), or make the layout consistent"
            ),
            samples=[s.sample_id for s in offending],
        )
    layout = "paired-end" if paired else "single-end"
    return _ok("read_layout", f"all samples share one read layout ({layout})")


def _ok(name: str, message: str) -> GateResult:
    return GateResult(
        name=name, module=MODULE, status=PASS, message=message, remedy="no action needed"
    )


def _factor_value(sample: Sample, factor: str):
    """Bir faktörün örnekteki değeri: çekirdek faktörler attribute, kovaryatlar dict."""
    if factor in ("condition", "batch", "subject"):
        return getattr(sample, factor)
    return sample.covariates.get(factor)


def _rank_gate(samples: list[Sample], variables: list[str]) -> GateResult:
    # Tasarımdaki condition-DIŞI HER kategorik faktör (batch, subject, sex, lane...)
    # ayrı ayrı kontrol edilmeli — yalnızca batch'e bakmak, subject/kovaryat üzerinden
    # doygun/confounded tasarımları sessizce PASS ettirir (bkz. Finding 1: 4 örnek + 4
    # benzersiz subject -> DESeq2 "model matrix is not full rank" ile patlar).
    factors = [v for v in variables if v != "condition"]
    if not factors:
        return _ok(
            "design_rank",
            "the design uses no factor beyond condition, so it is full rank",
        )

    n_samples = len(samples)
    for factor in factors:
        levels = {_factor_value(s, factor) for s in samples}

        # (a) Tek seviye: faktör hiçbir ayrım bilgisi eklemez.
        if len(levels) < 2:
            return GateResult(
                name="design_rank", module=MODULE, status=FAIL,
                message=(
                    f"the design uses {factor!r} but every sample is in the same {factor} "
                    f"({sorted(levels)[0]!r}), which makes the model matrix rank-deficient"
                ),
                remedy=f"use design '~condition', or supply samples from more than one {factor}",
                samples=[s.sample_id for s in samples],
            )

        # (c) Doygun (saturated): seviye sayısı örnek sayısına eşitse, her
        # seviye tam olarak BİR örnekte görülür (pigeonhole ilkesi) -> residual
        # serbestlik derecesi kalmaz. Bu, 'subject' için özellikle sinsi bir
        # durumdur: her hastadan tek ölçüm varsa (tekrarlanan ölçüm yok),
        # tasarım "confounded" değil ama yine de rank-deficient'tır — bu yüzden
        # confounded kontrolünden ÖNCE kontrol ediyoruz ki teşhis mesajı doğru
        # nedeni (doygunluk) adlandırsın.
        if len(levels) == n_samples:
            offending = sorted(str(x) for x in levels)
            return GateResult(
                name="design_rank", module=MODULE, status=FAIL,
                message=(
                    f"the design uses {factor!r} but it has as many levels as there are "
                    f"samples ({len(levels)} levels for {n_samples} samples: "
                    f"{', '.join(offending)}), so the design is saturated and leaves no "
                    "residual degrees of freedom"
                ),
                remedy=(
                    f"drop {factor!r} from the design, or provide repeated measurements so "
                    f"that {factor!r} levels repeat across samples"
                ),
                samples=[s.sample_id for s in samples],
            )

        # (b) Tam confounded: her seviye yalnızca tek bir condition'da
        # görülüyorsa, faktörün etkisi condition'ınkinden ayrıştırılamaz.
        by_factor: dict[str, set[str]] = {}
        for sample in samples:
            by_factor.setdefault(_factor_value(sample, factor), set()).add(sample.condition)
        if all(len(conditions) == 1 for conditions in by_factor.values()):
            mapping = ", ".join(f"{lvl}->{next(iter(c))}" for lvl, c in sorted(by_factor.items()))
            offending_samples = [
                s.sample_id for s in samples if len(by_factor[_factor_value(s, factor)]) == 1
            ]
            return GateResult(
                name="design_rank", module=MODULE, status=FAIL,
                message=(
                    f"{factor} is completely confounded with condition, so their effects "
                    f"cannot be separated ({mapping})"
                ),
                remedy=(
                    f"drop {factor!r} from the design, or use a layout where at least one "
                    f"{factor} contains more than one condition"
                ),
                samples=offending_samples,
            )

    return _ok(
        "design_rank", "design factors are not confounded with condition; the design is full rank"
    )


def _replication_gate(samples: list[Sample]) -> GateResult:
    counts = Counter(s.condition for s in samples)
    if len(counts) < 2:
        return GateResult(
            name="replication", module=MODULE, status=FAIL,
            message=(
                "condition has fewer than 2 levels, so there is nothing to compare "
                f"(found: {', '.join(sorted(counts))})"
            ),
            remedy="provide samples from at least two condition levels",
        )
    without = sorted(c for c, n in counts.items() if n < 2)
    if without:
        return GateResult(
            name="replication", module=MODULE, status=FAIL,
            message=(
                f"condition level(s) without replicate: {', '.join(without)}; "
                "DESeq2 cannot estimate dispersion without replicates"
            ),
            remedy="add at least one more sample for each condition level listed above",
            samples=[s.sample_id for s in samples if s.condition in without],
        )
    return _ok("replication", "every condition level has at least two replicates")


def _paired_gate(
    samples: list[Sample], variables: list[str], paired: bool | None
) -> GateResult:
    if "subject" in variables:
        return _ok("paired_declared", "the design accounts for the paired structure")
    if paired is not None:
        return _ok("paired_declared", f"pairing was declared explicitly (paired={paired})")
    if not looks_paired(samples):
        return _ok("paired_declared", "the data is not paired")

    # Hangi subject(ler) eşleşmiş görünüyor (>1 condition'da)? Kapı mesajı ve
    # samples alanı bunları ADLANDIRMALI — teşhis raporu FAIL veren her kapının
    # sorumlu örneklerini render eder (bkz. Finding 2).
    by_subject: dict[str, set[str]] = {}
    for sample in samples:
        if sample.subject is None:
            continue
        by_subject.setdefault(sample.subject, set()).add(sample.condition)
    paired_subjects = sorted(subj for subj, conditions in by_subject.items() if len(conditions) > 1)
    offending_samples = [s.sample_id for s in samples if s.subject in paired_subjects]

    return GateResult(
        name="paired_declared", module=MODULE, status=FAIL,
        message=(
            "the metadata looks PAIRED (subject(s) appearing in more than one condition: "
            f"{', '.join(paired_subjects)}) but the design does not use 'subject'; an "
            "unpaired analysis leaves subject-to-subject variation in the noise and hides "
            "real differences"
        ),
        remedy=(
            "use design '~subject + condition', or declare 'paired: false' in the config "
            "to run unpaired on purpose"
        ),
        samples=offending_samples,
    )


def looks_paired(samples: list[Sample]) -> bool:
    """En az bir subject birden fazla condition'da görünüyorsa veri eşleşmiş demektir.

    Bu YALNIZCA bir tespittir; design'a karar VERMEZ (spec 2026-07-20). Tahmin etmek
    yanlış modeli sessizce koşturabilirdi — burada amaç kullanıcıya sormak.
    """
    by_subject: dict[str, set[str]] = {}
    for sample in samples:
        if sample.subject is None:
            continue
        by_subject.setdefault(sample.subject, set()).add(sample.condition)
    return any(len(conditions) > 1 for conditions in by_subject.values())
