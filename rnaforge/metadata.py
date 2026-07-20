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

    if "batch" in variables:
        # Rank-deficient design'lar DESeq2'de kriptik bir matris hatasına dönüşür
        # ("model matrix is not full rank"). Burada yakalayıp NE yapılacağını söylemek
        # çok daha ucuz — ve sessizce yanlış bir modele koşmaktan güvenli.
        batches = {s.batch for s in samples}
        if len(batches) < 2:
            raise MetadataError(
                f"design formula uses 'batch' but every sample is in the same batch "
                f"({batches.pop()!r}). A single-level batch adds no information and makes "
                "the model matrix rank-deficient. Use design '~condition'."
            )
        by_batch: dict[str, set[str]] = {}
        for sample in samples:
            by_batch.setdefault(sample.batch, set()).add(sample.condition)
        if all(len(conds) == 1 for conds in by_batch.values()):
            mapping = ", ".join(
                f"{b}->{next(iter(c))}" for b, c in sorted(by_batch.items())
            )
            raise MetadataError(
                "batch is completely confounded with condition, so their effects cannot be "
                f"separated ({mapping}). Either drop 'batch' from the design, or use a layout "
                "where at least one batch contains more than one condition."
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
