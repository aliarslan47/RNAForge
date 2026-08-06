# Long-Read Arm — Step 1: Detection → Routing + `library.chemistry` Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Introduce the `read_type` (`short`|`long`) dimension — detect it from the platform, stop *rejecting* long reads and *route* them instead, add the `library.chemistry` config field, and record `read_type` as the single source of truth that downstream modules dispatch on (long branches raise a clear "not implemented yet" until later steps land).

**Architecture:** Follow the codebase's existing router convention — `m04_quant` already dispatches on `organism_type` while tool wrappers live in separate modules (`bowtie2.py`). `read_type` is derived from the already-detected platform in `platform.py`, persisted by m01 into `statistics/raw_statistics.json`, and read back by a small `routing.py` helper. m02/m03/m04/m05 gain a one-line `require_short_read()` guard at the router seam where the long-read branches will later attach — so a long-read run gets an honest `NotImplementedError` instead of silently running Illumina tools. This step touches **pure stdlib only**: no new conda env, no external tools, no real ONT data — everything is TDD-able with the synthetic `ont_fastq`/`pacbio_fastq` fixtures already in `tests/conftest.py`.

**Tech Stack:** Python 3.11 (stdlib), pytest. Test runner: `conda run -n rnaforge-core --cwd /home/ali/rnaforge-pipeline python -m pytest -q`.

## Global Constraints

- **Language:** code + logs in English; comments may be Turkish (match surrounding file). Plans/specs Turkish.
- **read_type values:** exactly `short` | `long`. **chemistry values:** exactly `cdna` | `direct_rna`.
- **platform → read_type map:** `illumina → short`; `ont → long`; `pacbio_hifi → long`; `unknown → error`.
- **Rule 7 (detection ≠ support):** unsupported/unidentifiable input is never silently processed. `unknown` platform is still refused; long platforms are now *routed*, and any not-yet-built long stage must fail **loudly** (`NotImplementedError`), never run the wrong tool. (feedback_gurultulu_hata)
- **Quality-gate philosophy unchanged:** FAIL = invalid (stop), WARN = suspect (stamp). This step adds no new data gate; it adds one *validation* error (missing chemistry for ONT-long) which is a config error, not a quality verdict.
- **Config discipline:** `library` is already in `KNOWN_TOP_LEVEL_KEYS`; `chemistry` is a sub-key, so no top-level change. Invalid values raise `ConfigError` with an actionable message.
- **Frozen dataclasses:** `Library` is `@dataclass(frozen=True)` — keep it frozen.
- **Test command (run from repo root via conda):** `conda run -n rnaforge-core --cwd /home/ali/rnaforge-pipeline python -m pytest -q` (running pytest from outside the repo breaks `tests.conftest` import — false alarm).

---

## File Structure

- `rnaforge/platform.py` — **modify**: add `READ_TYPES`, `read_type_for()`; broaden `SUPPORTED_PLATFORMS` to route long reads; rewrite `require_supported` message for the `unknown`-only rejection.
- `rnaforge/config.py` — **modify**: add `CHEMISTRY` enum, `Library.chemistry` field, parse + validate it.
- `rnaforge/modules/m01_validate.py` — **modify**: compute + record `read_type` and `chemistry`; enforce chemistry-required-for-ONT; drop the now-stale "unreachable" comment.
- `rnaforge/routing.py` — **create**: `resolve_read_type(run_dir)` + `require_short_read(run_dir, stage)`.
- `rnaforge/modules/m02_qc.py`, `m03_trim.py`, `m04_quant.py`, `m05_counts.py` — **modify**: one-line `require_short_read()` guard after the m01-done precondition.
- `config/config.yaml` — **modify**: document the new `library.chemistry` field in comments.
- Tests: `tests/test_platform.py`, `tests/test_config.py`, `tests/test_m01_validate.py` (modify) and `tests/test_routing.py` (create).

---

### Task 1: `read_type_for()` derivation in `platform.py`

**Files:**
- Modify: `rnaforge/platform.py`
- Test: `tests/test_platform.py`

**Interfaces:**
- Produces: `READ_TYPES = ("short", "long")`; `read_type_for(platform: str) -> str` — maps `"illumina"→"short"`, `"ont"→"long"`, `"pacbio_hifi"→"long"`; raises `ValueError` for `"unknown"` or any unrecognised platform.

- [ ] **Step 1: Write the failing test**

Add to `tests/test_platform.py`:

```python
def test_read_type_for_maps_platforms():
    from rnaforge.platform import read_type_for
    assert read_type_for("illumina") == "short"
    assert read_type_for("ont") == "long"
    assert read_type_for("pacbio_hifi") == "long"


def test_read_type_for_rejects_unknown():
    from rnaforge.platform import read_type_for
    with pytest.raises(ValueError):
        read_type_for("unknown")
```

- [ ] **Step 2: Run test to verify it fails**

Run: `conda run -n rnaforge-core --cwd /home/ali/rnaforge-pipeline python -m pytest tests/test_platform.py -k read_type -q`
Expected: FAIL — `ImportError: cannot import name 'read_type_for'`.

- [ ] **Step 3: Write minimal implementation**

In `rnaforge/platform.py`, near the top (after `SUPPORTED_PLATFORMS`):

```python
READ_TYPES = ("short", "long")

_PLATFORM_READ_TYPE = {
    "illumina": "short",
    "ont": "long",
    "pacbio_hifi": "long",
}


def read_type_for(platform: str) -> str:
    """Detected platform -> read_type. 'unknown' has no route (Rule 7)."""
    try:
        return _PLATFORM_READ_TYPE[platform]
    except KeyError:
        raise ValueError(
            f"cannot derive read_type for platform {platform!r}; "
            f"known platforms: {', '.join(_PLATFORM_READ_TYPE)}"
        ) from None
```

- [ ] **Step 4: Run test to verify it passes**

Run: `conda run -n rnaforge-core --cwd /home/ali/rnaforge-pipeline python -m pytest tests/test_platform.py -k read_type -q`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add rnaforge/platform.py tests/test_platform.py
git commit -m "feat(platform): derive read_type (short|long) from detected platform"
```

---

### Task 2: Route long reads instead of rejecting them

**Files:**
- Modify: `rnaforge/platform.py`
- Test: `tests/test_platform.py`

**Interfaces:**
- Consumes: `detect_platform`, `read_type_for` (Task 1).
- Produces: `SUPPORTED_PLATFORMS = ("illumina", "ont", "pacbio_hifi")`; `require_supported(info, fastq)` now raises `UnsupportedPlatformError` **only** for `unknown` (or empty) input. The exception type name is unchanged (still imported by `cli.py` and `m01_validate.py`).

- [ ] **Step 1: Rewrite the existing rejection tests to expect routing**

In `tests/test_platform.py`, **replace** `test_ont_rejected_with_actionable_message` and `test_pacbio_rejected` with:

```python
def test_ont_is_now_routed_not_rejected(ont_fastq):
    """Long reads are routed, no longer refused (Step 1 of the long-read arm)."""
    require_supported(detect_platform(ont_fastq), ont_fastq)  # must NOT raise


def test_pacbio_is_now_routed_not_rejected(pacbio_fastq):
    require_supported(detect_platform(pacbio_fastq), pacbio_fastq)  # must NOT raise
```

Keep `test_illumina_is_supported` and `test_empty_fastq_is_unknown_and_rejected` unchanged.

- [ ] **Step 2: Run tests to verify the new expectations fail**

Run: `conda run -n rnaforge-core --cwd /home/ali/rnaforge-pipeline python -m pytest tests/test_platform.py -q`
Expected: FAIL — `test_ont_is_now_routed_not_rejected` / `test_pacbio_is_now_routed_not_rejected` raise `UnsupportedPlatformError` (old behaviour still in place).

- [ ] **Step 3: Broaden support + rewrite the rejection message**

In `rnaforge/platform.py`:

Change the constant:

```python
SUPPORTED_PLATFORMS = ("illumina", "ont", "pacbio_hifi")
```

Replace the body of `require_supported` with an `unknown`-only refusal:

```python
def require_supported(info: PlatformInfo, fastq: Path) -> None:
    """Refuse only input we cannot identify. Long reads (ONT/PacBio) are routed
    by read_type downstream; unidentifiable input has no safe route (Rule 7)."""
    if info.platform in SUPPORTED_PLATFORMS:
        return
    raise UnsupportedPlatformError(
        f"could not identify the sequencing platform for this input "
        f"(detected {info.platform!r}; supported: {', '.join(SUPPORTED_PLATFORMS)}).\n"
        f"  file: {fastq}\n"
        f"  mean read length: {info.mean_read_length}, N50: {info.n50}, "
        f"mean quality: {info.mean_quality}, reads sampled: {info.n_reads_sampled}\n"
        f"Running any route on unidentifiable reads would produce wrong results, "
        f"so it is refused."
    )
```

- [ ] **Step 4: Run the full platform test module**

Run: `conda run -n rnaforge-core --cwd /home/ali/rnaforge-pipeline python -m pytest tests/test_platform.py -q`
Expected: PASS (all, including `test_empty_fastq_is_unknown_and_rejected`).

- [ ] **Step 5: Commit**

```bash
git add rnaforge/platform.py tests/test_platform.py
git commit -m "feat(platform): route ONT/PacBio instead of rejecting; refuse only unknown"
```

---

### Task 3: `library.chemistry` config field

**Files:**
- Modify: `rnaforge/config.py`
- Modify: `config/config.yaml` (document the field)
- Test: `tests/test_config.py`

**Interfaces:**
- Produces: `CHEMISTRY = ("cdna", "direct_rna")`; `Library` dataclass gains `chemistry: str | None = None`. `load_config` parses `library.chemistry`: absent → `None`; present → validated via `_one_of(..., CHEMISTRY, "library.chemistry")`.

- [ ] **Step 1: Write the failing tests**

Add to `tests/test_config.py` (mirror the file's existing `write_text`/`load_config` pattern; keep the required `organism`, `organism_type`, and `reference` fields a valid config in that file already uses):

```python
def test_library_chemistry_parsed(tmp_path):
    from rnaforge.config import load_config
    cfg = tmp_path / "c.yaml"
    cfg.write_text(
        'organism: "E. coli"\n'
        'organism_type: "prokaryote"\n'
        "reference:\n"
        '  genome_fasta: "g.fa"\n'
        '  annotation_gff: "g.gff"\n'
        "library:\n"
        '  chemistry: "direct_rna"\n'
    )
    assert load_config(cfg).library.chemistry == "direct_rna"


def test_library_chemistry_defaults_to_none(tmp_path):
    from rnaforge.config import load_config
    cfg = tmp_path / "c.yaml"
    cfg.write_text(
        'organism: "E. coli"\n'
        'organism_type: "prokaryote"\n'
        "reference:\n"
        '  genome_fasta: "g.fa"\n'
        '  annotation_gff: "g.gff"\n'
    )
    assert load_config(cfg).library.chemistry is None


def test_library_chemistry_invalid_rejected(tmp_path):
    from rnaforge.config import load_config, ConfigError
    cfg = tmp_path / "c.yaml"
    cfg.write_text(
        'organism: "E. coli"\n'
        'organism_type: "prokaryote"\n'
        "reference:\n"
        '  genome_fasta: "g.fa"\n'
        '  annotation_gff: "g.gff"\n'
        "library:\n"
        '  chemistry: "nanopore"\n'
    )
    with pytest.raises(ConfigError):
        load_config(cfg)
```

Ensure `import pytest` is present at the top of `tests/test_config.py` (it is used elsewhere in the file).

- [ ] **Step 2: Run tests to verify they fail**

Run: `conda run -n rnaforge-core --cwd /home/ali/rnaforge-pipeline python -m pytest tests/test_config.py -k chemistry -q`
Expected: FAIL — `TypeError: __init__() got an unexpected keyword argument 'chemistry'` or `AttributeError`.

- [ ] **Step 3: Implement the field**

In `rnaforge/config.py`:

Add the enum next to the other tuples (after `SELECTIONS`):

```python
CHEMISTRY = ("cdna", "direct_rna")
```

Extend the `Library` dataclass:

```python
@dataclass(frozen=True)
class Library:
    strandedness: str = "unstranded"
    selection: str = "rrna_depletion"
    # Long-read only: cDNA needs Pychopper full-length orientation; direct-RNA does not.
    # NOT detectable from FASTQ (spec 2026-08-05). None = unset (fine for short reads).
    chemistry: str | None = None
```

In `load_config`, where the `Library(...)` is built, add chemistry parsing (keep the existing `strandedness`/`selection` lines):

```python
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
        ),
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `conda run -n rnaforge-core --cwd /home/ali/rnaforge-pipeline python -m pytest tests/test_config.py -k chemistry -q`
Expected: PASS.

- [ ] **Step 5: Document the field in the example config**

In `config/config.yaml`, under the `library:` block, add:

```yaml
  # Long-read only (ONT): cdna | direct_rna. FASTQ'dan tespit EDİLEMEZ.
  # cDNA -> Pychopper full-length yönlendirme; direct-RNA -> yok. HiFi platformdan çıkar.
  # chemistry: "cdna"
```

- [ ] **Step 6: Commit**

```bash
git add rnaforge/config.py config/config.yaml tests/test_config.py
git commit -m "feat(config): add library.chemistry (cdna|direct_rna) for long-read arm"
```

---

### Task 4: m01 records `read_type` + `chemistry`, routes long, enforces ONT chemistry

**Files:**
- Modify: `rnaforge/modules/m01_validate.py`
- Test: `tests/test_m01_validate.py`

**Interfaces:**
- Consumes: `read_type_for` (Task 1), `Library.chemistry` (Task 3), `detect_platform`/`require_supported` (Task 2).
- Produces: the validation summary written to `statistics/raw_statistics.json` gains two keys — `read_type` (`"short"`/`"long"`) and `chemistry` (the config value or `null`). New validation rule: `read_type == "long"` **and** `platform == "ont"` **and** `config.library.chemistry is None` → `ValueError` (actionable). Short reads and PacBio-HiFi do not require chemistry.

- [ ] **Step 1: Write the failing tests**

In `tests/test_m01_validate.py` (the file already has `_setup`, `_illumina`, and `_ont` helpers). The existing `_ont` writes only 50 reads with variable length — good enough for detection. Add:

```python
def _setup_ont_with_chemistry(tmp_path, chemistry: str | None):
    """_setup but ONT reads + optional library.chemistry."""
    config_path, metadata_path = _setup(tmp_path, _ont)
    text = config_path.read_text()
    if chemistry is not None:
        text += f'library:\n  chemistry: "{chemistry}"\n'
    config_path.write_text(text)
    return config_path, metadata_path


def test_illumina_records_short_read_type(tmp_path):
    config_path, metadata_path = _setup(tmp_path, _illumina)
    summary = run_validation(load_config(config_path), metadata_path, tmp_path / "run")
    assert summary["read_type"] == "short"
    assert summary["chemistry"] is None


def test_ont_with_chemistry_records_long(tmp_path):
    config_path, metadata_path = _setup_ont_with_chemistry(tmp_path, "cdna")
    summary = run_validation(load_config(config_path), metadata_path, tmp_path / "run")
    assert summary["platform"] == "ont"
    assert summary["read_type"] == "long"
    assert summary["chemistry"] == "cdna"


def test_ont_without_chemistry_is_rejected(tmp_path):
    config_path, metadata_path = _setup_ont_with_chemistry(tmp_path, None)
    with pytest.raises(ValueError) as exc:
        run_validation(load_config(config_path), metadata_path, tmp_path / "run")
    assert "chemistry" in str(exc.value).lower()
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `conda run -n rnaforge-core --cwd /home/ali/rnaforge-pipeline python -m pytest tests/test_m01_validate.py -k "read_type or chemistry or long" -q`
Expected: FAIL — `KeyError: 'read_type'` for the first two; the ONT-no-chemistry test currently does NOT raise (long reads no longer rejected after Task 2), so it fails too.

- [ ] **Step 3: Implement in `m01_validate.py`**

Add the import at the top (with the existing platform import):

```python
from rnaforge.platform import PlatformInfo, detect_platform, read_type_for, require_supported
```

In `run_validation`, after `platform = platforms.pop()` and the existing `config.platform` mismatch check, derive read_type and enforce chemistry. **Delete** the now-stale comment block that says the checks are "MVP'de ulaşılamaz" (they are reachable now). Insert:

```python
        read_type = read_type_for(platform)
        chemistry = config.library.chemistry
        # cDNA vs direct-RNA is undetectable from FASTQ but changes m03 (Pychopper).
        # ONT long reads must declare it; HiFi is inferred, short reads don't care.
        if read_type == "long" and platform == "ont" and chemistry is None:
            raise ValueError(
                "ONT long-read input requires library.chemistry to be set "
                "('cdna' or 'direct_rna'): it cannot be detected from the FASTQ "
                "and it changes preprocessing (cDNA needs Pychopper). "
                "Set library.chemistry in the config and re-run validate."
            )
```

Then add the two keys to the `summary` dict (next to `"platform": platform,`):

```python
            "read_type": read_type,
            "chemistry": chemistry,
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `conda run -n rnaforge-core --cwd /home/ali/rnaforge-pipeline python -m pytest tests/test_m01_validate.py -q`
Expected: PASS (new tests + all pre-existing m01 tests; the old ONT-rejection test in this file, if any, must be reconciled — see Step 5).

- [ ] **Step 5: Reconcile any pre-existing ONT-rejection test in this file**

Search: `grep -n "UnsupportedPlatformError\|_ont" tests/test_m01_validate.py`. If a test asserts m01 raises `UnsupportedPlatformError` on ONT input, it now contradicts routing — update it to instead assert the chemistry requirement (i.e., ONT *without* chemistry raises `ValueError` mentioning "chemistry", ONT *with* chemistry succeeds). Re-run the module.

Run: `conda run -n rnaforge-core --cwd /home/ali/rnaforge-pipeline python -m pytest tests/test_m01_validate.py -q`
Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add rnaforge/modules/m01_validate.py tests/test_m01_validate.py
git commit -m "feat(m01): record read_type + chemistry, route long reads, require ONT chemistry"
```

---

### Task 5: `routing.py` resolver + `require_short_read` guard in m02–m05

**Files:**
- Create: `rnaforge/routing.py`
- Create: `tests/test_routing.py`
- Modify: `rnaforge/modules/m02_qc.py`, `m03_trim.py`, `m04_quant.py`, `m05_counts.py`

**Interfaces:**
- Consumes: the `read_type` key written by m01 into `statistics/raw_statistics.json` (Task 4).
- Produces:
  - `resolve_read_type(run_dir: Path) -> str` — reads `run_dir/statistics/raw_statistics.json`, returns its `read_type`; raises `ValueError` if the file or key is missing (m01 must have run).
  - `require_short_read(run_dir: Path, stage: str) -> None` — no-op for `short`; raises `NotImplementedError` for `long`, naming the stage.

- [ ] **Step 1: Write the failing tests**

Create `tests/test_routing.py`:

```python
from __future__ import annotations

import json

import pytest

from rnaforge.routing import require_short_read, resolve_read_type


def _write_stats(run_dir, read_type):
    stats = run_dir / "statistics"
    stats.mkdir(parents=True)
    (stats / "raw_statistics.json").write_text(json.dumps({"read_type": read_type}))


def test_resolve_read_type_reads_short(tmp_path):
    _write_stats(tmp_path, "short")
    assert resolve_read_type(tmp_path) == "short"


def test_resolve_read_type_missing_file_raises(tmp_path):
    with pytest.raises(ValueError):
        resolve_read_type(tmp_path)


def test_require_short_read_passes_for_short(tmp_path):
    _write_stats(tmp_path, "short")
    require_short_read(tmp_path, "qc")  # must NOT raise


def test_require_short_read_blocks_long(tmp_path):
    _write_stats(tmp_path, "long")
    with pytest.raises(NotImplementedError) as exc:
        require_short_read(tmp_path, "qc")
    assert "qc" in str(exc.value)
    assert "long" in str(exc.value).lower()
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `conda run -n rnaforge-core --cwd /home/ali/rnaforge-pipeline python -m pytest tests/test_routing.py -q`
Expected: FAIL — `ModuleNotFoundError: No module named 'rnaforge.routing'`.

- [ ] **Step 3: Implement `rnaforge/routing.py`**

```python
"""read_type dispatch seam. m01 detects read_type and writes it into the run's
raw_statistics.json (single source of truth); the per-stage routers read it back
here. Long-read stages are not built yet — they must fail loudly, never run the
short-read tool on long reads (Rule 7 / feedback_gurultulu_hata)."""
from __future__ import annotations

import json
from pathlib import Path


def resolve_read_type(run_dir: Path | str) -> str:
    """Return the read_type m01 recorded for this run."""
    stats = Path(run_dir) / "statistics" / "raw_statistics.json"
    if not stats.exists():
        raise ValueError(
            f"cannot resolve read_type: {stats} not found. "
            "Run `rnaforge validate` (m01) with the same --run-id first."
        )
    data = json.loads(stats.read_text())
    read_type = data.get("read_type")
    if read_type is None:
        raise ValueError(
            f"cannot resolve read_type: no 'read_type' key in {stats}. "
            "Re-run `rnaforge validate` (m01) to regenerate it."
        )
    return read_type


def require_short_read(run_dir: Path | str, stage: str) -> None:
    """Guard at a short-read-only stage. No-op for short; loud stop for long."""
    read_type = resolve_read_type(run_dir)
    if read_type == "long":
        raise NotImplementedError(
            f"long-read {stage} is not implemented yet (long-read arm Step 1 only "
            "routes; NanoPlot/Pychopper/minimap2/featureCounts -L come in later "
            f"steps). This run is read_type={read_type!r}."
        )
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `conda run -n rnaforge-core --cwd /home/ali/rnaforge-pipeline python -m pytest tests/test_routing.py -q`
Expected: PASS.

- [ ] **Step 5: Wire the guard into m02–m05**

In each of `m02_qc.py`, `m03_trim.py`, `m04_quant.py`, `m05_counts.py`, add the import:

```python
from rnaforge.routing import require_short_read
```

Then, in each module's `run_*` function, immediately **after** the existing "m01 must have completed" precondition check (and before any tool is invoked), add one line with the stage name:

- `m02_qc.py` (`run_qc`): `require_short_read(run_dir, "qc")`
- `m03_trim.py` (`run_trim`): `require_short_read(run_dir, "trim")`
- `m04_quant.py` (`run_quant`): `require_short_read(run_dir, "quant")`
- `m05_counts.py` (`run_counts`): `require_short_read(run_dir, "counts")`

For any module whose `run_*` does not already assert m01 completion, place the guard right after the resume short-circuit (the `if not force and state.is_done(...)` block) so a resumed short-read run is unaffected. Confirm each module resolves `run_dir` as a `Path` before the call (all four already do).

- [ ] **Step 6: Run the full suite**

Run: `conda run -n rnaforge-core --cwd /home/ali/rnaforge-pipeline python -m pytest -q`
Expected: PASS — all pre-existing tests (they run short-read fixtures, so the guard is a no-op) plus the new routing tests.

- [ ] **Step 7: Commit**

```bash
git add rnaforge/routing.py tests/test_routing.py \
        rnaforge/modules/m02_qc.py rnaforge/modules/m03_trim.py \
        rnaforge/modules/m04_quant.py rnaforge/modules/m05_counts.py
git commit -m "feat(routing): read_type resolver + short-read guard in m02-m05"
```

---

## Post-plan verification (whole step)

- [ ] Full suite green: `conda run -n rnaforge-core --cwd /home/ali/rnaforge-pipeline python -m pytest -q`.
- [ ] Manual smoke — short-read run unaffected: `rnaforge validate` on the existing `raw/GSE300731/` config still reports `read_type=short` in `statistics/raw_statistics.json`, and `qc`/`quant` still run.
- [ ] Manual smoke — long-read honest stop: craft a tiny synthetic ONT config (or reuse a test fixture) → `validate` succeeds and records `read_type=long` (with `library.chemistry` set) → `qc` stops with a clear `NotImplementedError` naming the stage. No Illumina tool runs on long reads.

## Out of scope (subsequent plans — need `rnaforge-longread` env + a selected real ONT dataset)

- Step 2: m02-long NanoPlot · Step 3: m03-long Pychopper + chopper · Step 4: m04-long minimap2 · Step 5: m05-long featureCounts `-L` · Step 6: long-read quality profile/gates · Step 7: report `read_type` badge + end-to-end live smoke on real bacterial ONT data (microbepore *E. coli* candidate — still to be selected; the first Step-2 task).

---

## Self-Review

**Spec coverage (spec §"Uygulama sırası" item 1 = "Tespit+yönlendirme (red→route) + library.chemistry config"):**
- "red→route": Task 2 (broaden support) + Task 4 (m01 no longer rejects long) ✓
- "library.chemistry config": Task 3 ✓
- read_type dimension "m02'den itibaren": Tasks 1 (derive) + 4 (persist) + 5 (dispatch seam in m02–m05) ✓
- Spec detection §: chemistry undetectable → config field (Task 3) + ONT requires it (Task 4) ✓; HiFi chemistry "platformdan çıkar" → not required for pacbio_hifi (Task 4 condition is ONT-only) ✓
- Spec "risk: module naming (m04-prok-long vs dispatcher)": resolved by the router-seam approach (guard in existing modules, matching the m04 `organism_type` convention) — no new module names introduced this step ✓
- Items 2–7 (NanoPlot/Pychopper/minimap2/featureCounts -L/profile/report) explicitly deferred to their own plans (they need the new env + a real ONT dataset that is not yet selected) ✓

**Placeholder scan:** No TBD/TODO; every code step shows concrete code. ✓

**Type consistency:** `read_type_for` (Task 1) returns `short|long`; consumed by m01 (Task 4) and stored under `read_type`; `resolve_read_type`/`require_short_read` (Task 5) read the same key. `Library.chemistry` (Task 3) is read by m01 (Task 4) as `config.library.chemistry`. `require_short_read(run_dir, stage)` signature identical in definition (Task 5 Step 3) and all four call sites (Task 5 Step 5). ✓
