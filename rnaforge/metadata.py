"""Örnek metadata (TSV) yükleme ve DESeq2 design formülü doğrulama.

Hatalar burada, pipeline koşmadan ÖNCE yakalanır (PLAN §13).
"""
from __future__ import annotations

import csv
import re
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

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
        samples.append(Sample(sample_id, condition, fastqs[0], fastqs[1], batch, subject))

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

    known = {"condition", "batch", "subject"}
    unknown = [v for v in variables if v not in known]
    if unknown:
        raise MetadataError(
            f"design formula references unknown variable(s): {', '.join(unknown)}. "
            f"Metadata columns available for design: {', '.join(sorted(known))}"
        )
    for variable in ("batch", "subject"):
        if variable in variables and any(getattr(s, variable) is None for s in samples):
            raise MetadataError(
                f"design formula uses {variable!r} but the metadata has no {variable} "
                f"value for every sample. Add a {variable!r} column, or drop it from "
                "the design."
            )

    return [
        _rank_gate(samples, variables),
        _replication_gate(samples),
        _paired_gate(samples, variables, paired),
    ]


def _ok(name: str, message: str) -> GateResult:
    return GateResult(
        name=name, module=MODULE, status=PASS, message=message, remedy="no action needed"
    )


def _rank_gate(samples: list[Sample], variables: list[str]) -> GateResult:
    if "batch" not in variables:
        return _ok("design_rank", "the design uses no batch term, so it is full rank")

    batches = {s.batch for s in samples}
    if len(batches) < 2:
        return GateResult(
            name="design_rank", module=MODULE, status=FAIL,
            message=(
                f"the design uses 'batch' but every sample is in the same batch "
                f"({sorted(batches)[0]!r}), which makes the model matrix rank-deficient"
            ),
            remedy="use design '~condition', or supply samples from more than one batch",
        )

    by_batch: dict[str, set[str]] = {}
    for sample in samples:
        by_batch.setdefault(sample.batch, set()).add(sample.condition)
    if all(len(conditions) == 1 for conditions in by_batch.values()):
        mapping = ", ".join(f"{b}->{next(iter(c))}" for b, c in sorted(by_batch.items()))
        return GateResult(
            name="design_rank", module=MODULE, status=FAIL,
            message=(
                "batch is completely confounded with condition, so their effects cannot "
                f"be separated ({mapping})"
            ),
            remedy=(
                "drop 'batch' from the design, or use a layout where at least one batch "
                "contains more than one condition"
            ),
            samples=[s.sample_id for s in samples],
        )
    return _ok("design_rank", "batch and condition are not confounded; the design is full rank")


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
    return GateResult(
        name="paired_declared", module=MODULE, status=FAIL,
        message=(
            "the metadata looks PAIRED (at least one subject appears in more than one "
            "condition) but the design does not use 'subject'; an unpaired analysis "
            "leaves subject-to-subject variation in the noise and hides real differences"
        ),
        remedy=(
            "use design '~subject + condition', or declare 'paired: false' in the config "
            "to run unpaired on purpose"
        ),
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
