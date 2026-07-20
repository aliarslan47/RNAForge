# Kalite Kapıları Çerçevesi — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Her koşunun sonucuna güvenilip güvenilemeyeceğini ölçen kapı çerçevesini kurmak; m01'in mevcut kontrollerini bu sözleşmeye taşımak ve eşleşmiş tasarım desteği eklemek.

**Architecture:** `gates.py` bir `GateResult` sözleşmesi tanımlar (PASS/WARN/FAIL). Eşikler koda gömülü değil, `profiles/*.yml` içinde **veri** olarak durur; `organism_type` profili seçer, config `quality:` bölümü ezer, ezilen eşik damgalanır. Modüller kapılarını koşup `runs/<id>/quality/gates.json`'a ekler; `FAIL` varsa `GateFailure` fırlatılır ve pipeline durur.

**Tech Stack:** Python 3.11, PyYAML, pytest. Harici biyoinformatik aracı YOK (bu plan saf Python).

## Global Constraints

- Referans doküman `PLAN.md` **v1.3**; yeniden yazılmaz, yalnız revize edilip sürümü yükseltilir.
- Kod, değişken adları ve log mesajları **İngilizce**; yorumlar Türkçe olabilir (mevcut emsal).
- `organism_type` zorunlu, varsayılanı yok — tahmin sessiz hataya yol açar.
- **FAIL = sonuç geçersiz** (biyolojik çıktı üretilmez). **WARN = sonuç şüpheli** (üretilir, damgalanır).
- Ezilen her eşik `overridden=True` alır ve güvence kartında görünür — **sessiz gevşetme yok**.
- Ökaryot profili bilinçli gevşektir ve "geniş toleranslı" damgası taşır.
- Hata mesajı **ne yapılacağını** söylemeli (`remedy` alanı zorunlu, boş bırakılamaz).
- Test env: `conda run -n rnaforge-core python -m pytest -q`
- Her task sonunda tüm suite yeşil olmalı; kırmızı bırakılmaz.

## File Structure

| Dosya | Sorumluluk |
|---|---|
| `rnaforge/gates.py` (yeni) | `GateStatus`, `GateResult`, `GateFailure`, `run_gates`, `write_gate_results` |
| `rnaforge/profiles/prokaryote.yml` (yeni) | Prokaryot eşikleri (veri) |
| `rnaforge/profiles/eukaryote.yml` (yeni) | Ökaryot eşikleri (veri, gevşek + damga) |
| `rnaforge/quality.py` (yeni) | Profil yükleme + config ile ezme + override kaydı |
| `rnaforge/metadata.py` (mevcut) | `subject` sütunu; eşleşmiş tasarım tespiti |
| `rnaforge/modules/m01_validate.py` (mevcut) | Kontrolleri kapı sözleşmesine taşı |
| `rnaforge/config.py` (mevcut) | `quality:` ve `paired:` bölümleri |
| `rnaforge/report/confidence.py` (yeni) | Güvence kartı JSON üretimi |

---

### Task 1: Kapı sözleşmesi (`gates.py`)

**Files:**
- Create: `rnaforge/gates.py`
- Test: `tests/test_gates.py`

**Interfaces:**
- Consumes: yok (çekirdek)
- Produces: `GateStatus` (PASS/WARN/FAIL sabitleri), `GateResult` dataclass, `GateFailure(Exception)` (`.failures: list[GateResult]`), `write_gate_results(run_dir: Path, results: list[GateResult]) -> Path`, `raise_if_failed(results: list[GateResult]) -> None`

- [ ] **Step 1: Write the failing test**

```python
# tests/test_gates.py
from __future__ import annotations

import json

import pytest

from rnaforge.gates import (
    FAIL,
    PASS,
    WARN,
    GateFailure,
    GateResult,
    raise_if_failed,
    write_gate_results,
)


def _result(name="alignment_rate", status=PASS, **kw):
    defaults = dict(
        module="m04", status=status, measured=0.9, threshold=0.7,
        message="alignment rate is the share of reads mapped to the reference",
        remedy="check that the reference genome matches the organism",
        overridden=False,
    )
    defaults.update(kw)
    return GateResult(name=name, **defaults)


def test_remedy_must_not_be_empty():
    """Kapı mesajı NE YAPILACAGINI söylemezse müşteriye faydası yok."""
    with pytest.raises(ValueError, match="remedy"):
        _result(remedy="")


def test_raise_if_failed_is_quiet_when_all_pass():
    raise_if_failed([_result(status=PASS), _result(status=WARN)])


def test_raise_if_failed_raises_on_fail():
    failing = _result(status=FAIL, measured=0.42)
    with pytest.raises(GateFailure) as exc:
        raise_if_failed([_result(status=PASS), failing])
    assert exc.value.failures == [failing]
    assert "alignment_rate" in str(exc.value)


def test_write_gate_results_appends_across_modules(tmp_path):
    """Resume ile uyumlu olmali: m02'nin sonuclari m01'inkileri EZMEMELI."""
    write_gate_results(tmp_path, [_result(name="design_rank", module="m01")])
    write_gate_results(tmp_path, [_result(name="alignment_rate", module="m04")])
    data = json.loads((tmp_path / "quality" / "gates.json").read_text())
    assert [g["name"] for g in data["gates"]] == ["design_rank", "alignment_rate"]


def test_rerunning_a_module_replaces_its_own_results(tmp_path):
    """--force ile yeniden kosulan modul kendi eski sonucunu birakmamali."""
    write_gate_results(tmp_path, [_result(name="design_rank", module="m01", status=FAIL)])
    write_gate_results(tmp_path, [_result(name="design_rank", module="m01", status=PASS)])
    data = json.loads((tmp_path / "quality" / "gates.json").read_text())
    assert len(data["gates"]) == 1
    assert data["gates"][0]["status"] == PASS
```

- [ ] **Step 2: Run test to verify it fails**

Run: `conda run -n rnaforge-core python -m pytest tests/test_gates.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'rnaforge.gates'`

- [ ] **Step 3: Write minimal implementation**

```python
# rnaforge/gates.py
"""Kalite kapıları: bir koşunun sonucuna güvenilip güvenilemeyeceğini ölçer.

Sözleşme (spec 2026-07-20):
  FAIL = sonuç GEÇERSİZ -> pipeline durur, biyolojik çıktı üretilmez
  WARN = sonuç ŞÜPHELİ  -> üretilir ama damgalanır
  PASS = kapı geçildi
"""
from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from datetime import datetime
from pathlib import Path

PASS = "PASS"
WARN = "WARN"
FAIL = "FAIL"
STATUSES = (PASS, WARN, FAIL)

GATES_FILE = "gates.json"
QUALITY_DIR = "quality"


@dataclass(frozen=True)
class GateResult:
    name: str
    module: str
    status: str
    message: str
    remedy: str
    measured: float | None = None
    threshold: float | None = None
    overridden: bool = False
    samples: list[str] = field(default_factory=list)

    def __post_init__(self) -> None:
        if self.status not in STATUSES:
            raise ValueError(f"status must be one of {STATUSES}, got {self.status!r}")
        if not self.message.strip():
            raise ValueError(f"gate {self.name}: message must not be empty")
        if not self.remedy.strip():
            # Ne yapılacağını söylemeyen kapı, müşteriyi çıkmaza sokar.
            raise ValueError(f"gate {self.name}: remedy must not be empty")


class GateFailure(Exception):
    """Bir veya daha fazla kapı FAIL verdi; pipeline burada durur."""

    def __init__(self, failures: list[GateResult]):
        self.failures = failures
        names = ", ".join(f.name for f in failures)
        detail = "\n".join(f"  - {f.name}: {f.message} -> {f.remedy}" for f in failures)
        super().__init__(f"quality gate(s) failed: {names}\n{detail}")


def raise_if_failed(results: list[GateResult]) -> None:
    failures = [r for r in results if r.status == FAIL]
    if failures:
        raise GateFailure(failures)


def write_gate_results(run_dir: Path | str, results: list[GateResult]) -> Path:
    """Sonuçları quality/gates.json'a EKLE. Aynı modülün eski sonuçları değiştirilir
    (--force ile yeniden koşma), diğer modüllerinkine dokunulmaz (resume uyumu)."""
    quality_dir = Path(run_dir) / QUALITY_DIR
    quality_dir.mkdir(parents=True, exist_ok=True)
    path = quality_dir / GATES_FILE

    existing: list[dict] = []
    if path.exists():
        try:
            existing = json.loads(path.read_text()).get("gates", [])
        except (json.JSONDecodeError, OSError):
            existing = []

    modules = {r.module for r in results}
    kept = [g for g in existing if g.get("module") not in modules]
    payload = {
        "updated_at": datetime.now().isoformat(timespec="seconds"),
        "gates": kept + [asdict(r) for r in results],
    }
    path.write_text(json.dumps(payload, indent=2))
    return path
```

- [ ] **Step 4: Run test to verify it passes**

Run: `conda run -n rnaforge-core python -m pytest tests/test_gates.py -v`
Expected: PASS (5 passed)

- [ ] **Step 5: Run full suite**

Run: `conda run -n rnaforge-core python -m pytest -q`
Expected: 57 passed

- [ ] **Step 6: Commit**

```bash
git add rnaforge/gates.py tests/test_gates.py
git commit -m "feat: kalite kapisi sozlesmesi (GateResult/GateFailure/gates.json)"
```

---

### Task 2: Profiller ve eşik ezme (`quality.py`)

**Files:**
- Create: `rnaforge/profiles/prokaryote.yml`, `rnaforge/profiles/eukaryote.yml`, `rnaforge/quality.py`
- Modify: `rnaforge/config.py` (`quality:` bölümü), `pyproject.toml` (paket verisi)
- Test: `tests/test_quality.py`

**Interfaces:**
- Consumes: Task 1 (`GateResult`)
- Produces: `load_profile(organism_type: str, overrides: dict | None = None) -> Profile`; `Profile.threshold(name: str) -> float`; `Profile.is_overridden(name: str) -> bool`; `Profile.name: str`; `Profile.permissive: bool`

- [ ] **Step 1: Write the failing test**

```python
# tests/test_quality.py
from __future__ import annotations

import pytest

from rnaforge.quality import ProfileError, load_profile


def test_prokaryote_profile_is_strict_about_alignment():
    profile = load_profile("prokaryote")
    assert profile.name == "prokaryote"
    assert profile.threshold("alignment_rate") == 0.70
    assert profile.permissive is False


def test_eukaryote_profile_is_marked_permissive():
    """Elde okaryot dogrulamasi YOK; gevsek esikler rapora damgalanmali."""
    profile = load_profile("eukaryote")
    assert profile.permissive is True
    assert profile.threshold("alignment_rate") == 0.50


def test_override_changes_threshold_and_is_recorded():
    profile = load_profile("prokaryote", overrides={"alignment_rate": 0.30})
    assert profile.threshold("alignment_rate") == 0.30
    assert profile.is_overridden("alignment_rate") is True
    assert profile.is_overridden("survival_rate") is False


def test_unknown_override_key_is_rejected():
    """Yazim hatasi sessizce yutulmamali; kullanici esigi ezdigini saniyor olabilir."""
    with pytest.raises(ProfileError, match="alignment_rat"):
        load_profile("prokaryote", overrides={"alignment_rat": 0.3})


def test_unknown_organism_type_is_rejected():
    with pytest.raises(ProfileError, match="no quality profile"):
        load_profile("archaea")
```

- [ ] **Step 2: Run test to verify it fails**

Run: `conda run -n rnaforge-core python -m pytest tests/test_quality.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'rnaforge.quality'`

- [ ] **Step 3: Write profile data files**

```yaml
# rnaforge/profiles/prokaryote.yml
name: prokaryote
permissive: false
description: >
  Bakteriyel RNA-seq. Sıkı eşikler: temiz bir bakteri koşusu tipik olarak
  %90+ hizalama verir. Bu profil doğrulanmış alandır.
thresholds:
  read_depth: 1000000
  base_quality: 30
  survival_rate: 0.50
  alignment_rate: 0.70
  rrna_fraction: 0.20
  replicate_correlation: 0.85
```

```yaml
# rnaforge/profiles/eukaryote.yml
name: eukaryote
permissive: true
description: >
  Ökaryotik RNA-seq. Eşikler BİLİNÇLİ olarak gevşektir: elde doğrulanmış ökaryot
  veri seti yoktur ve uydurma bir eşik kapı sisteminin anlamını bozar. Bu profille
  üretilen rapor "geniş toleranslı" damgası taşır. Ökaryot doğrulaması geldiğinde
  eşikler sıkılacaktır.
thresholds:
  read_depth: 10000000
  base_quality: 30
  survival_rate: 0.50
  alignment_rate: 0.50
  rrna_fraction: 0.30
  replicate_correlation: 0.80
```

- [ ] **Step 4: Write minimal implementation**

```python
# rnaforge/quality.py
"""Kalite eşikleri: koda gömülü DEĞİL, profiles/*.yml içinde veri olarak durur.

Yeni bir profil eklemek kod değişikliği değil, dosya eklemektir (spec 2026-07-20).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import yaml

PROFILE_DIR = Path(__file__).parent / "profiles"


class ProfileError(ValueError):
    """Kalite profili yüklenemedi veya geçersiz eşik ezmesi verildi."""


@dataclass(frozen=True)
class Profile:
    name: str
    permissive: bool
    description: str
    _thresholds: dict[str, float]
    _overridden: set[str] = field(default_factory=set)

    def threshold(self, gate: str) -> float:
        if gate not in self._thresholds:
            raise ProfileError(
                f"profile {self.name!r} has no threshold for gate {gate!r}; "
                f"known gates: {', '.join(sorted(self._thresholds))}"
            )
        return self._thresholds[gate]

    def is_overridden(self, gate: str) -> bool:
        return gate in self._overridden

    def overrides(self) -> dict[str, float]:
        return {g: self._thresholds[g] for g in sorted(self._overridden)}


def load_profile(organism_type: str, overrides: dict | None = None) -> Profile:
    path = PROFILE_DIR / f"{organism_type}.yml"
    if not path.exists():
        available = ", ".join(sorted(p.stem for p in PROFILE_DIR.glob("*.yml")))
        raise ProfileError(
            f"no quality profile for organism_type={organism_type!r} "
            f"(available: {available})"
        )
    raw = yaml.safe_load(path.read_text()) or {}
    thresholds = dict(raw.get("thresholds") or {})

    applied: set[str] = set()
    for gate, value in (overrides or {}).items():
        if gate not in thresholds:
            # Yazım hatasını yutmak, kullanıcıya "eşiği gevşettim" yanılgısı verir.
            raise ProfileError(
                f"quality.{gate}: unknown gate for profile {organism_type!r}; "
                f"known gates: {', '.join(sorted(thresholds))}"
            )
        try:
            thresholds[gate] = float(value)
        except (TypeError, ValueError):
            raise ProfileError(f"quality.{gate}: expected a number, got {value!r}") from None
        applied.add(gate)

    return Profile(
        name=raw.get("name", organism_type),
        permissive=bool(raw.get("permissive", False)),
        description=str(raw.get("description", "")).strip(),
        _thresholds=thresholds,
        _overridden=applied,
    )
```

- [ ] **Step 5: Ship profile YAMLs inside the package**

`pyproject.toml` içinde paket verisi bildirilmeli, aksi halde kurulu pakette `profiles/`
bulunamaz. `[tool.setuptools.package-data]` bölümüne ekle:

```toml
[tool.setuptools.package-data]
rnaforge = ["profiles/*.yml"]
```

- [ ] **Step 6: Run test to verify it passes**

Run: `conda run -n rnaforge-core python -m pytest tests/test_quality.py -v`
Expected: PASS (5 passed)

- [ ] **Step 7: Commit**

```bash
git add rnaforge/quality.py rnaforge/profiles tests/test_quality.py pyproject.toml
git commit -m "feat: kalite profilleri veri olarak (prokaryot siki, okaryot gevsek+damgali)"
```

---

### Task 3: `subject` sütunu ve eşleşmiş tasarım tespiti

**Files:**
- Modify: `rnaforge/metadata.py`
- Test: `tests/test_metadata.py`

**Interfaces:**
- Consumes: yok
- Produces: `Sample.subject: str | None`; `looks_paired(samples: list[Sample]) -> bool`

- [ ] **Step 1: Write the failing test**

```python
# tests/test_metadata.py — dosyanın SONUNA ekle
def test_subject_column_is_loaded(tmp_path):
    _make_fastqs(tmp_path, "a.fastq", "b.fastq")
    path = _write_meta(tmp_path, (
        "sample_id\tcondition\tsubject\tfastq_1\n"
        "s1\tbefore\tp1\ta.fastq\n"
        "s2\tafter\tp1\tb.fastq\n"
    ))
    samples = load_metadata(path)
    assert [s.subject for s in samples] == ["p1", "p1"]


def test_subject_is_none_when_column_absent(tmp_path):
    _make_fastqs(tmp_path, "a.fastq", "b.fastq")
    path = _write_meta(tmp_path, (
        "sample_id\tcondition\tfastq_1\n"
        "s1\tcontrol\ta.fastq\n"
        "s2\ttreated\tb.fastq\n"
    ))
    assert all(s.subject is None for s in load_metadata(path))


def test_looks_paired_detects_repeated_subject_across_conditions(tmp_path):
    _make_fastqs(tmp_path, "a.fastq", "b.fastq", "c.fastq", "d.fastq")
    path = _write_meta(tmp_path, (
        "sample_id\tcondition\tsubject\tfastq_1\n"
        "s1\tbefore\tp1\ta.fastq\n"
        "s2\tafter\tp1\tb.fastq\n"
        "s3\tbefore\tp2\tc.fastq\n"
        "s4\tafter\tp2\td.fastq\n"
    ))
    assert looks_paired(load_metadata(path)) is True


def test_looks_paired_false_when_each_subject_has_one_condition(tmp_path):
    """Her subject tek condition'da -> eslesmis degil, sadece etiketlenmis."""
    _make_fastqs(tmp_path, "a.fastq", "b.fastq", "c.fastq", "d.fastq")
    path = _write_meta(tmp_path, (
        "sample_id\tcondition\tsubject\tfastq_1\n"
        "s1\tcontrol\tp1\ta.fastq\n"
        "s2\tcontrol\tp2\tb.fastq\n"
        "s3\ttreated\tp3\tc.fastq\n"
        "s4\ttreated\tp4\td.fastq\n"
    ))
    assert looks_paired(load_metadata(path)) is False
```

`tests/test_metadata.py` başındaki import bloğuna `looks_paired` ekle:

```python
from rnaforge.metadata import (
    MetadataError,
    design_variables,
    load_metadata,
    looks_paired,
    validate_design,
)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `conda run -n rnaforge-core python -m pytest tests/test_metadata.py -v`
Expected: FAIL — `ImportError: cannot import name 'looks_paired'`

- [ ] **Step 3: Write minimal implementation**

`rnaforge/metadata.py` içinde `Sample` dataclass'ına alan ekle:

```python
@dataclass(frozen=True)
class Sample:
    sample_id: str
    condition: str
    fastq_1: Path
    fastq_2: Path | None = None
    batch: str | None = None
    subject: str | None = None
```

`load_metadata` içinde `batch` satırının hemen ardına ekle ve yapıcıyı güncelle:

```python
        batch = (row.get("batch") or "").strip() or None
        subject = (row.get("subject") or "").strip() or None
        samples.append(Sample(sample_id, condition, fastqs[0], fastqs[1], batch, subject))
```

Modülün sonuna ekle:

```python
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
```

- [ ] **Step 4: Run test to verify it passes**

Run: `conda run -n rnaforge-core python -m pytest tests/test_metadata.py -v`
Expected: PASS

- [ ] **Step 5: Run full suite**

Run: `conda run -n rnaforge-core python -m pytest -q`
Expected: tüm testler yeşil

- [ ] **Step 6: Commit**

```bash
git add rnaforge/metadata.py tests/test_metadata.py
git commit -m "feat: metadata subject sutunu + eslesmis tasarim tespiti"
```

---

### Task 4: `paired_declared` kapısı ve `subject` design değişkeni

**Files:**
- Modify: `rnaforge/metadata.py` (`validate_design`), `rnaforge/config.py` (`paired`)
- Test: `tests/test_metadata.py`, `tests/test_config.py`

**Interfaces:**
- Consumes: Task 3 (`looks_paired`, `Sample.subject`)
- Produces: `Config.paired: bool | None` (None = beyan edilmedi); `validate_design(samples, design, paired: bool | None = None)`

- [ ] **Step 1: Write the failing test**

```python
# tests/test_metadata.py — sona ekle
def test_paired_data_without_subject_in_design_is_rejected(tmp_path):
    """Tespit VAR, karar YOK: dur ve sor (spec 2026-07-20)."""
    _make_fastqs(tmp_path, "a.fastq", "b.fastq", "c.fastq", "d.fastq")
    path = _write_meta(tmp_path, (
        "sample_id\tcondition\tsubject\tfastq_1\n"
        "s1\tbefore\tp1\ta.fastq\n"
        "s2\tafter\tp1\tb.fastq\n"
        "s3\tbefore\tp2\tc.fastq\n"
        "s4\tafter\tp2\td.fastq\n"
    ))
    with pytest.raises(MetadataError, match="paired"):
        validate_design(load_metadata(path), "~condition")


def test_paired_false_declaration_allows_unpaired_design(tmp_path):
    """Bilerek eslesmemis kosmak MUMKUN olmali — ama beyan edilerek."""
    _make_fastqs(tmp_path, "a.fastq", "b.fastq", "c.fastq", "d.fastq")
    path = _write_meta(tmp_path, (
        "sample_id\tcondition\tsubject\tfastq_1\n"
        "s1\tbefore\tp1\ta.fastq\n"
        "s2\tafter\tp1\tb.fastq\n"
        "s3\tbefore\tp2\tc.fastq\n"
        "s4\tafter\tp2\td.fastq\n"
    ))
    validate_design(load_metadata(path), "~condition", paired=False)


def test_subject_in_design_satisfies_the_gate(tmp_path):
    _make_fastqs(tmp_path, "a.fastq", "b.fastq", "c.fastq", "d.fastq")
    path = _write_meta(tmp_path, (
        "sample_id\tcondition\tsubject\tfastq_1\n"
        "s1\tbefore\tp1\ta.fastq\n"
        "s2\tafter\tp1\tb.fastq\n"
        "s3\tbefore\tp2\tc.fastq\n"
        "s4\tafter\tp2\td.fastq\n"
    ))
    validate_design(load_metadata(path), "~subject + condition")


def test_subject_in_design_requires_subject_column(tmp_path):
    _make_fastqs(tmp_path, "a.fastq", "b.fastq", "c.fastq", "d.fastq")
    path = _write_meta(tmp_path, (
        "sample_id\tcondition\tfastq_1\n"
        "s1\tcontrol\ta.fastq\n"
        "s2\tcontrol\tb.fastq\n"
        "s3\ttreated\tc.fastq\n"
        "s4\ttreated\td.fastq\n"
    ))
    with pytest.raises(MetadataError, match="subject"):
        validate_design(load_metadata(path), "~subject + condition")
```

```python
# tests/test_config.py — sona ekle
def test_paired_defaults_to_undeclared(tmp_path):
    assert load_config(_write(tmp_path, PROK_BODY)).paired is None


def test_paired_can_be_declared_false(tmp_path):
    cfg = load_config(_write(tmp_path, PROK_BODY + "\npaired: false\n"))
    assert cfg.paired is False
```

- [ ] **Step 2: Run test to verify it fails**

Run: `conda run -n rnaforge-core python -m pytest tests/test_metadata.py tests/test_config.py -v`
Expected: FAIL — `TypeError: validate_design() got an unexpected keyword argument 'paired'` ve `AttributeError: 'Config' object has no attribute 'paired'`

- [ ] **Step 3: Implement config field**

`rnaforge/config.py` — `Config` dataclass'ına alan ekle (varsayılan `None`, çünkü
"beyan edilmedi" ile "hayır" farklı şeylerdir):

```python
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
```

`from dataclasses import dataclass, field` importunu güncelle. `load_config` içinde
`Config(...)` çağrısına ekle:

```python
        paired=None if raw.get("paired") is None else bool(raw.get("paired")),
        quality=_section(raw, "quality"),
```

- [ ] **Step 4: Implement design validation**

`rnaforge/metadata.py` — `validate_design` imzasını ve bilinen değişkenleri güncelle:

```python
def validate_design(
    samples: list[Sample], design: str, paired: bool | None = None
) -> None:
    variables = design_variables(design)
    if not variables:
        raise MetadataError(f"design formula has no variables: {design!r}")

    known = {"condition", "batch", "subject"}
```

`batch` kontrolünün hemen ardına ekle:

```python
    if "subject" in variables and any(s.subject is None for s in samples):
        raise MetadataError(
            "design formula uses 'subject' but the metadata has no subject value for "
            "every sample. Add a 'subject' column, or drop 'subject' from the design."
        )

    if "subject" not in variables and paired is not True and looks_paired(samples):
        if paired is None:
            raise MetadataError(
                "the metadata looks PAIRED: at least one subject appears in more than one "
                "condition, but the design formula does not use 'subject'. A paired design "
                f"({'~subject + condition'!r}) compares each subject with itself and finds "
                "considerably more true differences. Either use 'subject' in the design, or "
                "declare 'paired: false' in the config to run unpaired on purpose."
            )
```

- [ ] **Step 5: Run test to verify it passes**

Run: `conda run -n rnaforge-core python -m pytest tests/test_metadata.py tests/test_config.py -v`
Expected: PASS

- [ ] **Step 6: Run full suite**

Run: `conda run -n rnaforge-core python -m pytest -q`
Expected: tüm testler yeşil

- [ ] **Step 7: Commit**

```bash
git add rnaforge/metadata.py rnaforge/config.py tests/test_metadata.py tests/test_config.py
git commit -m "feat: paired_declared kapisi — eslesmis veri beyan edilmeden gecemez"
```

---

### Task 5: m01 kontrollerini kapı sözleşmesine taşı

**Files:**
- Modify: `rnaforge/modules/m01_validate.py`, `rnaforge/cli.py`
- Test: `tests/test_m01_validate.py`

**Interfaces:**
- Consumes: Task 1 (`GateResult`, `write_gate_results`, `raise_if_failed`), Task 2 (`load_profile`), Task 4 (`validate_design(..., paired=)`)
- Produces: `run_validation(...)` artık `quality/gates.json` yazar; `summary["quality_profile"]`, `summary["permissive_profile"]`

- [ ] **Step 1: Write the failing test**

`tests/test_m01_validate.py` içinde `_setup(tmp_path, fastq_maker)` yardımcısı zaten
var; `(config_path, metadata_path)` döndürür ve `illumina_fastq` deseninde FASTQ üretir.
Yeni testler onu kullanır — yeni fixture YAZMA.

```python
# tests/test_m01_validate.py — sona ekle
from rnaforge.gates import PASS


def _illumina(path):
    return write_fastq(path, 200, 150, "I")


def test_m01_writes_gate_results(tmp_path):
    """Kapilar gorunur olmali: PASS alan kosuda da neyin kontrol edildigi yazilir."""
    config_path, metadata_path = _setup(tmp_path, _illumina)
    run_dir = tmp_path / "run"
    run_validation(load_config(config_path), metadata_path, run_dir)
    data = json.loads((run_dir / "quality" / "gates.json").read_text())
    names = {g["name"] for g in data["gates"]}
    assert {"design_rank", "replication", "paired_declared"} <= names
    assert all(g["status"] == PASS for g in data["gates"])
    assert all(g["module"] == "m01" for g in data["gates"])


def test_m01_summary_records_quality_profile(tmp_path):
    config_path, metadata_path = _setup(tmp_path, _illumina)
    run_dir = tmp_path / "run"
    summary = run_validation(load_config(config_path), metadata_path, run_dir)
    assert summary["quality_profile"] == "prokaryote"
    assert summary["permissive_profile"] is False
```

- [ ] **Step 2: Run test to verify it fails**

Run: `conda run -n rnaforge-core python -m pytest tests/test_m01_validate.py -v`
Expected: FAIL — `FileNotFoundError` (gates.json yok)

- [ ] **Step 3: Write minimal implementation**

`rnaforge/modules/m01_validate.py` — import bloğuna ekle:

```python
from rnaforge.gates import PASS, GateResult, raise_if_failed, write_gate_results
from rnaforge.quality import load_profile
```

`validate_design` çağrısını `config.paired` geçecek şekilde güncelle ve kapıları kaydet.
`log(f"design formula ...: OK")` satırının yerine:

```python
        profile = load_profile(config.organism_type, config.quality)
        log(f"quality profile: {profile.name} (permissive={profile.permissive})")

        validate_design(samples, config.de.design, paired=config.paired)
        log(f"design formula {config.de.design!r}: OK")

        # Buraya ulaşıldıysa tasarım kontrolleri geçmiştir: validate_design ihlalde
        # MetadataError fırlatır. Kapıları PASS olarak kaydetmek, müşterinin NEYİN
        # kontrol edildiğini görmesini sağlar (spec §3.4 — görünmeyen güvence, güvence değildir).
        gates = [
            GateResult(
                name="design_rank", module="m01", status=PASS,
                message="the design matrix is full rank (batch is not confounded with condition)",
                remedy="no action needed",
            ),
            GateResult(
                name="replication", module="m01", status=PASS,
                message="every condition level has at least two replicates",
                remedy="no action needed",
            ),
            GateResult(
                name="paired_declared", module="m01", status=PASS,
                message="the pairing structure of the data matches the design formula",
                remedy="no action needed",
            ),
        ]
        write_gate_results(run_dir, gates)
        raise_if_failed(gates)
```

`summary` sözlüğüne ekle:

```python
            "quality_profile": profile.name,
            "permissive_profile": profile.permissive,
```

- [ ] **Step 4: Run test to verify it passes**

Run: `conda run -n rnaforge-core python -m pytest tests/test_m01_validate.py -v`
Expected: PASS

- [ ] **Step 5: Run full suite**

Run: `conda run -n rnaforge-core python -m pytest -q`
Expected: tüm testler yeşil

- [ ] **Step 6: Commit**

```bash
git add rnaforge/modules/m01_validate.py tests/test_m01_validate.py
git commit -m "feat: m01 kapilarini gates.json'a yaziyor + kalite profili secimi"
```

---

### Task 6: Güvence kartı

**Files:**
- Create: `rnaforge/report/__init__.py`, `rnaforge/report/confidence.py`
- Modify: `rnaforge/cli.py`
- Test: `tests/test_confidence.py`

**Interfaces:**
- Consumes: Task 1 (`gates.json`), Task 2 (`Profile`)
- Produces: `build_confidence_card(run_dir: Path, profile: Profile) -> dict`; `write_confidence_card(run_dir: Path, profile: Profile) -> Path`

- [ ] **Step 1: Write the failing test**

```python
# tests/test_confidence.py
from __future__ import annotations

import json

from rnaforge.gates import FAIL, PASS, WARN, GateResult, write_gate_results
from rnaforge.quality import load_profile
from rnaforge.report.confidence import build_confidence_card, write_confidence_card


def _gate(name, status, module="m01"):
    return GateResult(
        name=name, module=module, status=status,
        message="check", remedy="do something",
    )


def test_card_summarises_gate_counts(tmp_path):
    write_gate_results(tmp_path, [
        _gate("design_rank", PASS),
        _gate("rrna_fraction", WARN, module="m04"),
    ])
    card = build_confidence_card(tmp_path, load_profile("prokaryote"))
    assert card["counts"] == {"PASS": 1, "WARN": 1, "FAIL": 0}
    assert card["verdict"] == "SUSPECT"


def test_verdict_is_trustworthy_when_all_pass(tmp_path):
    write_gate_results(tmp_path, [_gate("design_rank", PASS)])
    card = build_confidence_card(tmp_path, load_profile("prokaryote"))
    assert card["verdict"] == "TRUSTWORTHY"


def test_verdict_is_invalid_on_any_fail(tmp_path):
    write_gate_results(tmp_path, [_gate("design_rank", PASS), _gate("alignment_rate", FAIL)])
    card = build_confidence_card(tmp_path, load_profile("prokaryote"))
    assert card["verdict"] == "INVALID"


def test_card_records_permissive_profile_and_overrides(tmp_path):
    """Gevsetilen esik ve gevsek profil GORUNMEK zorunda (spec §3.2)."""
    write_gate_results(tmp_path, [_gate("design_rank", PASS)])
    profile = load_profile("eukaryote", overrides={"alignment_rate": 0.2})
    card = build_confidence_card(tmp_path, profile)
    assert card["profile"]["permissive"] is True
    assert card["profile"]["overrides"] == {"alignment_rate": 0.2}


def test_write_confidence_card_creates_file(tmp_path):
    write_gate_results(tmp_path, [_gate("design_rank", PASS)])
    path = write_confidence_card(tmp_path, load_profile("prokaryote"))
    assert path.name == "confidence_card.json"
    assert json.loads(path.read_text())["verdict"] == "TRUSTWORTHY"
```

- [ ] **Step 2: Run test to verify it fails**

Run: `conda run -n rnaforge-core python -m pytest tests/test_confidence.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'rnaforge.report'`

- [ ] **Step 3: Write minimal implementation**

```python
# rnaforge/report/__init__.py
from __future__ import annotations
```

```python
# rnaforge/report/confidence.py
"""Güvence kartı: koşunun sonucuna ne kadar güvenilebileceğinin tek sayfalık özeti.

PASS alan koşuda da üretilir — müşteri NEYİN kontrol edildiğini görmelidir.
Görünmeyen güvence, güvence değildir (spec 2026-07-20 §3.4).
"""
from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path

from rnaforge.gates import FAIL, GATES_FILE, PASS, QUALITY_DIR, WARN
from rnaforge.quality import Profile

CARD_FILE = "confidence_card.json"

TRUSTWORTHY = "TRUSTWORTHY"
SUSPECT = "SUSPECT"
INVALID = "INVALID"


def build_confidence_card(run_dir: Path | str, profile: Profile) -> dict:
    gates_path = Path(run_dir) / QUALITY_DIR / GATES_FILE
    gates = json.loads(gates_path.read_text())["gates"] if gates_path.exists() else []

    counts = {
        PASS: sum(1 for g in gates if g["status"] == PASS),
        WARN: sum(1 for g in gates if g["status"] == WARN),
        FAIL: sum(1 for g in gates if g["status"] == FAIL),
    }
    if counts[FAIL]:
        verdict = INVALID
    elif counts[WARN]:
        verdict = SUSPECT
    else:
        verdict = TRUSTWORTHY

    return {
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "verdict": verdict,
        "counts": counts,
        "profile": {
            "name": profile.name,
            "permissive": profile.permissive,
            "description": profile.description,
            "overrides": profile.overrides(),
        },
        "gates": gates,
    }


def write_confidence_card(run_dir: Path | str, profile: Profile) -> Path:
    path = Path(run_dir) / QUALITY_DIR / CARD_FILE
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(build_confidence_card(run_dir, profile), indent=2))
    return path
```

- [ ] **Step 4: Run test to verify it passes**

Run: `conda run -n rnaforge-core python -m pytest tests/test_confidence.py -v`
Expected: PASS (5 passed)

- [ ] **Step 5: Wire into the CLI**

`rnaforge/cli.py` — `_cmd_validate` içinde, özet basıldıktan sonra kartı üret ve
verdict'i ekrana yaz (sessiz güvence, güvence değildir):

```python
    profile = load_profile(config.organism_type, config.quality)
    card_path = write_confidence_card(run_dir, profile)
    card = json.loads(card_path.read_text())
    print(f"quality verdict: {card['verdict']} "
          f"(PASS={card['counts']['PASS']} WARN={card['counts']['WARN']} "
          f"FAIL={card['counts']['FAIL']}, profile={profile.name})")
    if profile.permissive:
        print("note: this profile is permissive — thresholds are deliberately loose "
              "and results should be read with that in mind.")
    for gate, value in profile.overrides().items():
        print(f"note: threshold {gate} was overridden to {value} in the config.")
```

Gerekli importları ekle: `import json`, `from rnaforge.quality import load_profile`,
`from rnaforge.report.confidence import write_confidence_card`.

`GateFailure` CLI'da yakalanmalı — biyolojik sonuç üretilmediği net söylenmeli:

```python
    except GateFailure as exc:
        print(f"error: {exc}", file=sys.stderr)
        print("no results were produced: the data did not pass the quality gates.",
              file=sys.stderr)
        return 1
```

`from rnaforge.gates import GateFailure` importunu ekle. Bu blok mevcut
`except (ConfigError, ...)` bloğundan ÖNCE gelmeli.

- [ ] **Step 6: Verify live end-to-end**

```bash
cd /tmp && rm -rf gatecheck && mkdir gatecheck && cd gatecheck
python3 -c "
seq='ACGT'*25
for s in ['a','b','c','d']:
    open(s+'.fastq','w').write(''.join('@r%d\n%s\n+\n%s\n'%(i,seq,'I'*100) for i in range(200)))
"
touch genome.fasta ann.gff
printf 'organism: \"E. coli\"\norganism_type: \"prokaryote\"\nreference:\n  genome_fasta: \"genome.fasta\"\n  annotation_gff: \"ann.gff\"\n' > config.yaml
printf 'sample_id\tcondition\tsubject\tfastq_1\ns1\tbefore\tp1\ta.fastq\ns2\tafter\tp1\tb.fastq\ns3\tbefore\tp2\tc.fastq\ns4\tafter\tp2\td.fastq\n' > samples.tsv
conda run -n rnaforge-core rnaforge validate --config config.yaml --metadata samples.tsv --runs-dir runs --run-id p
```

Expected: `error:` ile durur, mesaj eşleşmiş veriyi bildirir ve `paired: false` ya da
`~subject + condition` önerir. Ardından design'ı `~subject + condition` yapıp tekrar
koş: `quality verdict: TRUSTWORTHY` görülmeli.

- [ ] **Step 7: Run full suite**

Run: `conda run -n rnaforge-core python -m pytest -q`
Expected: tüm testler yeşil

- [ ] **Step 8: Commit**

```bash
git add rnaforge/report tests/test_confidence.py rnaforge/cli.py
git commit -m "feat: guvence karti + CLI verdict ciktisi"
```

---

### Task 7: Dokümantasyon ve PLAN güncellemesi

**Files:**
- Modify: `PLAN.md` (v1.3 → v1.4), `README.md`, `README.tr.md`, `config/` altındaki örnek config

**Interfaces:**
- Consumes: Task 1-6
- Produces: yok (dokümantasyon)

- [ ] **Step 1: PLAN.md'ye kalite kapıları bölümü ekle**

Sürümü **v1.4** yap ve changelog'a ekle. Yeni bölüm: kapı sözleşmesi (FAIL/WARN),
profiller, güvence kartı, teşhis kipi. Katalog için spec'e referans ver
(`docs/superpowers/specs/2026-07-20-kalite-kapilari-design.md`) — **kopyalama**, tek
kaynak ilkesi (Kural 3) bozulmasın.

- [ ] **Step 2: README'lere "Quality gates" bölümü ekle**

Her iki dilde: kapıların ne olduğu, FAIL/WARN farkı, `quality:` ile eşik ezmenin
rapora yazıldığı, `paired:` beyanı. Örnek config parçası göster.

- [ ] **Step 3: Doğrula — README iddiaları koda karşı**

README'de yazan her bayrak/alan gerçekten var mı? `rnaforge validate --help` çıktısıyla
karşılaştır. (Task 7 emsali: Plan A'da bu kontrol bulgusuz geçmişti, aynı titizlik.)

- [ ] **Step 4: Commit**

```bash
git add PLAN.md README.md README.tr.md
git commit -m "docs: PLAN v1.4 kalite kapilari + README quality gates bolumu"
```

---

## Sonraki adım (bu planın DIŞINDA)

Veri kapıları kendi modülleriyle gelir: m02 `read_depth`/`base_quality`,
m03 `survival_rate`, m04 `alignment_rate`/`rrna_fraction`, m05 `genes_detected`,
m06 `replicate_correlation`/`sample_swap`/`dispersion_fit`.

**m04 uyarısı:** kurulu salmon **2.3.4**, PLAN 1.x varsayıyordu. Index/CLI davranışı
değişmiş olabilir — m04 yazılmadan önce doğrulanmalı.
