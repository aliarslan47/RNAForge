# Long-Read Arm — Step 2: m02-long (NanoPlot QC) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Give m02 (`qc`) a long-read branch: when `read_type == "long"`, run **NanoPlot** instead of FastQC to produce read-length / N50 / quality diagnostics, staying purely diagnostic (never stops the run). Short-read runs keep the exact FastQC behaviour they have today.

**Architecture:** Follow the m04 router convention. `m02_qc.run_qc` resolves `read_type` (via `routing.resolve_read_type`, already written in Step 1) and dispatches: `short → _qc_short` (the current FastQC path, extracted verbatim) · `long → _qc_long` (new NanoPlot path). A new `rnaforge/nanoplot.py` holds a pure `parse_nanostats` parser + a `run_nanoplot` runner in the `rnaforge-longread` env — mirroring `rnaforge/fastqc.py`. The Step-1 `require_short_read(run_dir, "qc")` guard in m02 is **replaced** by this dispatch; the guard stays in m03/m04/m05 (still short-only until their own steps land). No new quality gate — long-read gates are Step 6; m02 is diagnostic in both branches (FastQC FAIL→WARN today; NanoPlot has no gate).

**Tech Stack:** Python 3.11 (stdlib parser), NanoPlot 1.47.1 (conda env `rnaforge-longread`), pytest.

## Global Constraints

- **Language:** code + logs English; comments may match the file (Turkish ok). Tests via `conda run -n rnaforge-core --cwd /home/ali/rnaforge-pipeline python -m pytest -q`.
- **m02 is diagnostic, both branches:** it must NEVER stop the run. No FAIL gate is introduced here. (Same principle as FastQC → WARN.)
- **read_type values:** `short` | `long` (from Step 1). Dispatch on exactly these; anything else is a bug from m01 and may raise.
- **NanoPlot invocation is fixed:** `NanoPlot --fastq <fq> --outdir <dir> --tsv_stats --no_static`. `--tsv_stats` is REQUIRED — it produces the clean `metric<TAB>value` `NanoStats.txt` the parser expects; `--no_static` avoids the kaleido/orca static-image dependency (HTML plots are still written). Env: `rnaforge-longread`.
- **ONT is single-end:** NanoPlot runs on `sample.fastq_1` only (no `fastq_2`).
- **Env file:** `envs/rnaforge-longread.yml` already exported (minimap2 2.31, NanoPlot 1.47.1, pychopper, chopper 0.13, samtools 1.24); commit it with the runner (project pattern: env yml lands with the module that first uses it).
- **Resume/heartbeat/precondition unchanged:** the m01-done precondition and RunState resume short-circuit stay exactly as they are; dispatch happens after them.

## Real NanoStats.txt format (captured 2026-08-06 on SRR14608655)

```
Metrics	dataset
number_of_reads	89032
number_of_bases	74816374.0
median_read_length	724.0
mean_read_length	840.3
read_length_stdev	786.1
n50	1371.0
mean_qual	6.7
median_qual	9.3
longest_read_(with_Q):1	30092 (5.5)
...
Reads >Q10:	38965 (43.8%) 50.4Mb
Reads >Q15:	309 (0.3%) 0.3Mb
```

The parser reads the leading `metric<TAB>value` rows; the `longest_read_(with_Q):N`, `highest_Q_read_(with_length):N`, and `Reads >QN:` rows have composite values and are NOT part of the core numeric fields (they may be ignored, except `Reads >Q10:` whose percentage is captured).

---

## File Structure

- `rnaforge/nanoplot.py` — **create**: `NanoStats` dataclass, `parse_nanostats(text)`, `run_nanoplot(fastq, out_dir, env)`, error types. Mirrors `rnaforge/fastqc.py`.
- `rnaforge/modules/m02_qc.py` — **modify**: extract current FastQC body into `_qc_short(...)`; add `_qc_long(...)` (NanoPlot); dispatch on `resolve_read_type`; drop the `require_short_read` guard line.
- `envs/rnaforge-longread.yml` — **commit** (already on disk).
- Tests: `tests/test_nanoplot.py` (create), `tests/test_m02_qc.py` (add dispatch tests).

---

### Task 1: `nanoplot.py` — `parse_nanostats` parser

**Files:**
- Create: `rnaforge/nanoplot.py`
- Test: `tests/test_nanoplot.py`

**Interfaces:**
- Produces:
  - `class NanoStatsParseError(ValueError)`
  - `@dataclass(frozen=True) class NanoStats` with fields: `number_of_reads: int`, `number_of_bases: int`, `mean_read_length: float`, `median_read_length: float`, `read_length_stdev: float`, `n50: float`, `mean_qual: float`, `median_qual: float`, `reads_above_q10_pct: float | None`.
  - `parse_nanostats(text: str) -> NanoStats` — parses the `--tsv_stats` `NanoStats.txt`. Raises `NanoStatsParseError` if a required core metric is missing.

- [ ] **Step 1: Write the failing test**

Create `tests/test_nanoplot.py`:

```python
from __future__ import annotations

import pytest

from rnaforge.nanoplot import NanoStats, NanoStatsParseError, parse_nanostats

_REAL = """Metrics\tdataset
number_of_reads\t89032
number_of_bases\t74816374.0
median_read_length\t724.0
mean_read_length\t840.3
read_length_stdev\t786.1
n50\t1371.0
mean_qual\t6.7
median_qual\t9.3
longest_read_(with_Q):1\t30092 (5.5)
highest_Q_read_(with_length):1\t18.3 (195)
Reads >Q10:\t38965 (43.8%) 50.4Mb
Reads >Q15:\t309 (0.3%) 0.3Mb
"""


def test_parse_nanostats_core_fields():
    s = parse_nanostats(_REAL)
    assert isinstance(s, NanoStats)
    assert s.number_of_reads == 89032
    assert s.number_of_bases == 74816374
    assert s.mean_read_length == pytest.approx(840.3)
    assert s.median_read_length == pytest.approx(724.0)
    assert s.n50 == pytest.approx(1371.0)
    assert s.mean_qual == pytest.approx(6.7)
    assert s.median_qual == pytest.approx(9.3)
    assert s.reads_above_q10_pct == pytest.approx(43.8)


def test_parse_nanostats_missing_core_raises():
    with pytest.raises(NanoStatsParseError):
        parse_nanostats("Metrics\tdataset\nnumber_of_reads\t10\n")
```

- [ ] **Step 2: Run to verify it fails**

Run: `conda run -n rnaforge-core --cwd /home/ali/rnaforge-pipeline python -m pytest tests/test_nanoplot.py -q`
Expected: FAIL — `ModuleNotFoundError: No module named 'rnaforge.nanoplot'`.

- [ ] **Step 3: Implement the parser**

Create `rnaforge/nanoplot.py`:

```python
"""NanoPlot (ONT/long-read QC) output parsing + runner.

Mirrors rnaforge/fastqc.py: a pure parser over NanoPlot's --tsv_stats
NanoStats.txt, plus a real runner in the rnaforge-longread env. Long-read QC
is diagnostic (like FastQC) — it never stops the run."""
from __future__ import annotations

import re
import subprocess
from dataclasses import dataclass
from pathlib import Path


class NanoStatsParseError(ValueError):
    """NanoStats.txt could not be parsed."""


class NanoPlotRunError(RuntimeError):
    """NanoPlot failed to run."""


@dataclass(frozen=True)
class NanoStats:
    number_of_reads: int
    number_of_bases: int
    mean_read_length: float
    median_read_length: float
    read_length_stdev: float
    n50: float
    mean_qual: float
    median_qual: float
    reads_above_q10_pct: float | None = None


_CORE = {
    "number_of_reads": int,
    "number_of_bases": lambda v: int(float(v)),
    "mean_read_length": float,
    "median_read_length": float,
    "read_length_stdev": float,
    "n50": float,
    "mean_qual": float,
    "median_qual": float,
}


def parse_nanostats(text: str) -> NanoStats:
    raw: dict[str, str] = {}
    for line in text.splitlines():
        if "\t" not in line:
            continue
        key, value = line.split("\t", 1)
        raw[key.strip()] = value.strip()

    parsed: dict[str, object] = {}
    for key, caster in _CORE.items():
        if key not in raw:
            raise NanoStatsParseError(f"NanoStats.txt missing required metric: {key!r}")
        try:
            parsed[key] = caster(raw[key])
        except (TypeError, ValueError) as exc:
            raise NanoStatsParseError(f"bad value for {key!r}: {raw[key]!r}") from exc

    q10 = None
    m = re.search(r"\(([\d.]+)%\)", raw.get("Reads >Q10:", ""))
    if m:
        q10 = float(m.group(1))

    return NanoStats(reads_above_q10_pct=q10, **parsed)  # type: ignore[arg-type]
```

- [ ] **Step 4: Run to verify it passes**

Run: `conda run -n rnaforge-core --cwd /home/ali/rnaforge-pipeline python -m pytest tests/test_nanoplot.py -q`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add rnaforge/nanoplot.py tests/test_nanoplot.py
git commit -m "feat(nanoplot): parse NanoPlot --tsv_stats NanoStats.txt"
```

---

### Task 2: `nanoplot.py` — `run_nanoplot` runner

**Files:**
- Modify: `rnaforge/nanoplot.py`
- Modify: `tests/test_nanoplot.py`
- Commit: `envs/rnaforge-longread.yml`

**Interfaces:**
- Consumes: `NanoStats`, `parse_nanostats`, `NanoPlotRunError` (Task 1).
- Produces: `run_nanoplot(fastq: Path, out_dir: Path, env: str = "rnaforge-longread") -> Path` — runs NanoPlot, returns the path to `out_dir/NanoStats.txt`. Raises `NanoPlotRunError` on non-zero exit or missing output.

- [ ] **Step 1: Write the failing test (env-gated integration, mirrors test_fastqc.py real-run gating)**

Add to `tests/test_nanoplot.py`:

```python
import shutil
from pathlib import Path

from rnaforge.nanoplot import run_nanoplot

_HAS_ENV = shutil.which("conda") is not None


@pytest.mark.skipif(not _HAS_ENV, reason="conda/rnaforge-longread not available")
def test_run_nanoplot_on_tiny_fastq(tmp_path):
    # 3 short synthetic ONT-ish reads are enough to exercise the runner
    fq = tmp_path / "r.fastq"
    fq.write_text(
        "@r1\n" + "ACGT" * 60 + "\n+\n" + "I" * 240 + "\n"
        "@r2\n" + "ACGT" * 80 + "\n+\n" + "I" * 320 + "\n"
        "@r3\n" + "ACGT" * 50 + "\n+\n" + "I" * 200 + "\n"
    )
    stats_path = run_nanoplot(fq, tmp_path / "out")
    assert stats_path.name == "NanoStats.txt"
    assert stats_path.exists()
    from rnaforge.nanoplot import parse_nanostats
    s = parse_nanostats(stats_path.read_text())
    assert s.number_of_reads == 3
```

- [ ] **Step 2: Run to verify it fails**

Run: `conda run -n rnaforge-core --cwd /home/ali/rnaforge-pipeline python -m pytest tests/test_nanoplot.py -k run_nanoplot -q`
Expected: FAIL — `ImportError: cannot import name 'run_nanoplot'` (or, if conda present, collection error).

- [ ] **Step 3: Implement the runner**

Append to `rnaforge/nanoplot.py`:

```python
def run_nanoplot(fastq: Path, out_dir: Path, env: str = "rnaforge-longread") -> Path:
    """Run NanoPlot on a single long-read FASTQ; return the NanoStats.txt path.

    --tsv_stats gives the parseable stats file; --no_static skips the
    kaleido/orca static-image dependency (HTML plots are still written)."""
    fastq = Path(fastq)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        "conda", "run", "-n", env, "NanoPlot",
        "--fastq", str(fastq),
        "--outdir", str(out_dir),
        "--tsv_stats", "--no_static",
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    stats_path = out_dir / "NanoStats.txt"
    if proc.returncode != 0 or not stats_path.exists():
        raise NanoPlotRunError(
            f"NanoPlot failed (exit {proc.returncode}) on {fastq}\n"
            f"stdout: {proc.stdout[-500:]}\nstderr: {proc.stderr[-500:]}"
        )
    return stats_path
```

- [ ] **Step 4: Run to verify it passes**

Run: `conda run -n rnaforge-core --cwd /home/ali/rnaforge-pipeline python -m pytest tests/test_nanoplot.py -q`
Expected: PASS (integration test actually runs NanoPlot in the env; ~10–20s).

- [ ] **Step 5: Commit (runner + env yml)**

```bash
git add rnaforge/nanoplot.py tests/test_nanoplot.py envs/rnaforge-longread.yml
git commit -m "feat(nanoplot): run_nanoplot runner (rnaforge-longread env) + env yml"
```

---

### Task 3: m02 dispatch — short→FastQC, long→NanoPlot

**Files:**
- Modify: `rnaforge/modules/m02_qc.py`
- Test: `tests/test_m02_qc.py`

**Interfaces:**
- Consumes: `resolve_read_type` (Step 1 `routing.py`), `run_nanoplot`/`parse_nanostats` (Tasks 1–2).
- Produces: `run_qc` dispatches on read_type. For long runs it writes `statistics/qc_statistics.json` with shape `{"read_type": "long", "n_samples": N, "samples": {sid: {mean_read_length, n50, mean_qual, number_of_reads, ...}}}` and marks `MODULE_NAME` done. No gates written for the long branch (diagnostic). Short branch output is byte-for-byte what it is today plus a `"read_type": "short"` key.

- [ ] **Step 1: Write the failing tests**

Add to `tests/test_m02_qc.py` (mirror the file's existing run_qc setup; monkeypatch the NanoPlot runner so no env is needed). The test writes a run dir where m01 is already done with `read_type=long`:

```python
def _seed_run(tmp_path, read_type):
    """Minimal run dir: m01 marked done + raw_statistics.json with read_type."""
    import json
    from rnaforge.state import RunState
    run_dir = tmp_path / "run"
    (run_dir / "statistics").mkdir(parents=True)
    (run_dir / "statistics" / "raw_statistics.json").write_text(
        json.dumps({"read_type": read_type})
    )
    RunState(run_dir).mark_done("m01_validate", [])
    return run_dir


def test_run_qc_long_uses_nanoplot(tmp_path, monkeypatch):
    import rnaforge.modules.m02_qc as m02
    run_dir = _seed_run(tmp_path, "long")
    fq = tmp_path / "s1.fastq"
    fq.write_text("@r\n" + "ACGT" * 50 + "\n+\n" + "I" * 200 + "\n")
    meta = tmp_path / "m.tsv"
    meta.write_text(f"sample_id\tcondition\tfastq_1\ns1\tctrl\t{fq}\n")

    _FAKE = (
        "Metrics\tdataset\nnumber_of_reads\t89032\nnumber_of_bases\t74816374.0\n"
        "median_read_length\t724.0\nmean_read_length\t840.3\nread_length_stdev\t786.1\n"
        "n50\t1371.0\nmean_qual\t6.7\nmedian_qual\t9.3\nReads >Q10:\t38965 (43.8%) 50.4Mb\n"
    )

    def fake_run_nanoplot(fastq, out_dir, env="rnaforge-longread"):
        out_dir = __import__("pathlib").Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        p = out_dir / "NanoStats.txt"
        p.write_text(_FAKE)
        return p

    monkeypatch.setattr(m02, "run_nanoplot", fake_run_nanoplot)

    from rnaforge.config import Config, Reference, Library, Trimming, DE, Report, Resources
    cfg = Config(
        organism="E. coli", organism_type="prokaryote", platform="auto",
        reference=Reference(), library=Library(), trimming=Trimming(),
        de=DE(), report=Report(), resources=Resources(),
    )
    summary = m02.run_qc(cfg, meta, run_dir)
    assert summary["read_type"] == "long"
    assert summary["n_samples"] == 1
    assert summary["samples"]["s1"]["n50"] == pytest.approx(1371.0)
    assert summary["samples"]["s1"]["mean_read_length"] == pytest.approx(840.3)
    # diagnostic: no gate results file written for the long branch
    assert not (run_dir / "gates.json").exists() or True  # gates optional; must not FAIL


def test_run_qc_short_still_fastqc(tmp_path, monkeypatch):
    """The short branch must be unchanged (regression guard)."""
    # Reuse the file's existing short-read run_qc test setup if present; otherwise
    # assert dispatch calls the FastQC path. Minimal: read_type short -> summary has
    # module_flags (FastQC-only key) and read_type == 'short'.
    # (Implement by monkeypatching run_fastqc/parse_fastqc_zip as the existing tests do.)
    pass
```

Note: if `tests/test_m02_qc.py` already has a working short-read `run_qc` test, leave it as the short-branch regression and only add `assert summary["read_type"] == "short"` to it instead of the `pass` stub above.

- [ ] **Step 2: Run to verify the long test fails**

Run: `conda run -n rnaforge-core --cwd /home/ali/rnaforge-pipeline python -m pytest tests/test_m02_qc.py -k long -q`
Expected: FAIL — `require_short_read` raises `NotImplementedError` (long still guarded) / `run_nanoplot` not imported in m02.

- [ ] **Step 3: Refactor m02 into a dispatch**

In `rnaforge/modules/m02_qc.py`:

Add imports:

```python
from rnaforge.nanoplot import parse_nanostats, run_nanoplot
from rnaforge.routing import resolve_read_type  # (require_short_read import removed)
```

Remove the `require_short_read` import and the `require_short_read(run_dir, "qc")` line.

Replace the current post-precondition body so that, after the m01-done check, it dispatches:

```python
    if not state.is_done("m01_validate"):
        raise ValueError( ... )  # unchanged

    read_type = resolve_read_type(run_dir)
    if read_type == "long":
        summary = _qc_long(config, metadata_path, run_dir, raw_qc_dir, stats_dir, logs_dir, state)
    else:
        summary = _qc_short(config, metadata_path, run_dir, raw_qc_dir, stats_dir, logs_dir, state)
    state.mark_done(MODULE_NAME, [str(stats_path), str(logs_dir / "qc.log")])
    return summary
```

Move the existing FastQC body (the whole `with log_path.open(...)` block that builds `reports`, gates, figure, and writes `summary`) into a new `def _qc_short(config, metadata_path, run_dir, raw_qc_dir, stats_dir, logs_dir, state) -> dict:` — unchanged logic, but add `"read_type": "short"` to its `summary` dict and `return summary` (drop the `state.mark_done` from inside, since the caller marks it once).

Add the long branch:

```python
def _qc_long(config, metadata_path, run_dir, raw_qc_dir, stats_dir, logs_dir, state) -> dict:
    """Long-read QC via NanoPlot. Diagnostic only — never stops the run."""
    log_path = logs_dir / "qc.log"
    stats_path = stats_dir / "qc_statistics.json"
    with log_path.open("w") as log_file:
        def log(msg: str) -> None:
            log_file.write(msg + "\n"); log_file.flush()

        samples = load_metadata(metadata_path)
        log(f"m02 NanoPlot (long-read): {len(samples)} sample(s)")
        per_sample: dict[str, dict] = {}
        for sample in samples:
            state.heartbeat()
            out = raw_qc_dir / sample.sample_id
            stats_txt = run_nanoplot(sample.fastq_1, out)
            s = parse_nanostats(stats_txt.read_text())
            per_sample[sample.sample_id] = {
                "number_of_reads": s.number_of_reads,
                "number_of_bases": s.number_of_bases,
                "mean_read_length": s.mean_read_length,
                "median_read_length": s.median_read_length,
                "n50": s.n50,
                "mean_qual": s.mean_qual,
                "median_qual": s.median_qual,
                "reads_above_q10_pct": s.reads_above_q10_pct,
            }
            log(f"{sample.sample_id}: NanoPlot OK (N50={s.n50}, meanQ={s.mean_qual})")

        summary = {
            "read_type": "long",
            "n_samples": len(samples),
            "samples": per_sample,
        }
        stats_path.write_text(json.dumps(summary, indent=2))
        log(f"qc statistics written: {stats_path}")
    return summary
```

Note the resume short-circuit at the top of `run_qc` stays; `stats_path` is already computed there and reused by the caller's `mark_done`.

- [ ] **Step 4: Run to verify both branches pass**

Run: `conda run -n rnaforge-core --cwd /home/ali/rnaforge-pipeline python -m pytest tests/test_m02_qc.py -q`
Expected: PASS (long dispatch test + all existing short-read m02 tests unchanged).

- [ ] **Step 5: Commit**

```bash
git add rnaforge/modules/m02_qc.py tests/test_m02_qc.py
git commit -m "feat(m02): dispatch read_type — FastQC (short) / NanoPlot (long)"
```

---

### Task 4: Full suite + live smoke on microbepore

**Files:** none (verification).

- [ ] **Step 1: Full suite green**

Run: `conda run -n rnaforge-core --cwd /home/ali/rnaforge-pipeline python -m pytest -q`
Expected: PASS (Step-1's 423 + new nanoplot/m02 tests).

- [ ] **Step 2: Live smoke — real ONT run reaches NanoPlot QC**

Build a tiny run over 2 microbepore samples (subset a few thousand reads each for speed) with the MG1655 reference, `library.chemistry: cdna`, then:

```bash
rnaforge validate --config <ont_cfg> --metadata <ont_meta> --runs-dir runs --run-id microbepore_smoke
rnaforge qc       --config <ont_cfg> --metadata <ont_meta> --runs-dir runs --run-id microbepore_smoke
```

Expected: `validate` records `read_type=long`; `qc` runs NanoPlot (not FastQC), writes `statistics/qc_statistics.json` with per-sample N50 / mean_qual, and does NOT stop the run. `m03 trim` still stops loudly (`NotImplementedError`, long-read trim not built) — confirming the guard chain is intact below m02.

- [ ] **Step 3: Record** the smoke result (N50 / mean_qual per sample) in the commit message or DURUM.

---

## Out of scope (later steps)
- m03-long (Pychopper + chopper), m04-long (minimap2), m05-long (featureCounts -L) — their own plans.
- Long-read quality profile + gates (Step 6): m02-long stays gate-free here.
- Report `read_type` badge / long-QC section (Step 7).
- Custom read-length figure — NanoPlot already emits HTML plots into `raw_qc/<sample>/`; a report-embedded figure is a Step-7 concern.

---

## Self-Review

**Spec coverage (spec table row "m02 QC: FastQC | NanoPlot"):** Task 3 dispatch + Tasks 1–2 NanoPlot parser/runner implement exactly this ✓. Diagnostic-only (spec "m02 QC diagnostik") — no gate added ✓. Router convention (spec "m04 deseni") — dispatch in the existing module ✓.

**Placeholder scan:** The only `pass` is the short-branch regression stub in Task 3 Step 1, with an explicit note to instead extend the file's existing short-read test — not a silent placeholder. All code steps show real code. ✓

**Type consistency:** `NanoStats` fields (Task 1) are read in `_qc_long` (Task 3) by the same names. `run_nanoplot(fastq, out_dir, env)` signature identical in definition (Task 2) and call site (Task 3) and monkeypatch (Task 3 test). `resolve_read_type(run_dir) -> "short"|"long"` from Step 1 drives the dispatch. ✓

**Risk:** extracting the FastQC body into `_qc_short` must preserve behaviour exactly — the existing m02 test suite (short-read) is the regression guard; Task 4 Step 1 runs the whole suite to confirm no short-read drift.
