# Long-Read Arm — Step 3: m03-long (Pychopper + chopper) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Give m03 (`trim`) a long-read branch: when `read_type == "long"`, preprocess ONT reads by chemistry — `cdna` → **Pychopper** (orient + trim full-length cDNA) then **chopper** (length/quality filter); `direct_rna` → **chopper** only. Short-read runs keep their exact fastp behaviour. The trimmed output lands at the same `trimmed_reads(run_dir, sample)` path so the (future) m04-long can consume it via the existing contract.

**Architecture:** Follow the m04 router / m02-long convention. `m03_trim.run_trim` resolves `read_type` (via `routing.resolve_read_type`) and dispatches: `short → _trim_short` (current fastp path, extracted verbatim) · `long → _trim_long` (Pychopper/chopper by `config.library.chemistry`). New tool wrappers `rnaforge/chopper.py` and `rnaforge/pychopper.py` mirror `rnaforge/fastp.py` (runner + stats). The Step-1 `require_short_read(run_dir, "trim")` guard in m03 is **replaced** by this dispatch; guards stay in m04/m05.

**Tech Stack:** Python 3.11 (stdlib), Pychopper 2.7.10 + chopper 0.13.0 (conda env `rnaforge-longread`), pytest.

## Global Constraints

- **Language:** code + logs English; comments may match the file (Turkish ok). Tests via `conda run -n rnaforge-core --cwd /home/ali/rnaforge-pipeline python -m pytest -q`.
- **read_type / chemistry:** dispatch on `read_type in {short, long}` (Step 1) and, for long, `config.library.chemistry in {cdna, direct_rna}` (Step 1; m01 already requires it for ONT). A `long` run with `chemistry is None` is an m01 bug, but `_trim_long` must still fail loudly (ValueError) rather than silently pick a branch.
- **Trimmed-output contract:** long writes its final filtered FASTQ to `trimmed_reads(run_dir, sample)[0]` (plain `<stem>.trimmed.fastq`), the SAME single-source path m04 reads. ONT is single-end → `out2` is always None.
- **NO new FAIL gate here.** Short m03 keeps its `survival_rate` FAIL gate (calibrated for Illumina fastp). Long m03 records survival as a **diagnostic** only — the Illumina 0.50 threshold would wrongly fail cDNA (Pychopper discards non-full-length reads by design). The long-read profile/gates are Step 6.
- **Pychopper exits non-zero on a known upstream bug.** Pychopper 2.7.10's end-of-run PDF report crashes under pandas 3 (`_plot_stats` → `float(Series)` TypeError) AFTER writing the oriented FASTQ and the `-S` stats TSV. `run_pychopper` MUST tolerate exactly this: if the output FASTQ + stats TSV exist AND stderr contains `_plot_stats`, log a loud WARNING and proceed; for any other non-zero exit, raise `PychopperRunError` with full stderr. Never swallow a real failure. (feedback_gurultulu_hata)
- **Env:** `rnaforge-longread`. Tools: `pychopper -t <n> -S <stats.tsv> [-k <kit>] <in> <out>` · `chopper -q <q> -l <len> --threads <n>` (stdin→stdout).

## Real tool output (captured 2026-08-06 on microbepore PCR-cDNA subset)

- **chopper** exits 0, stderr: `Kept 1097 reads out of 1109 reads`.
- **pychopper** writes oriented full-length FASTQ (1109/2000 reads) + `-S` stats TSV, then exits 1 on the plotting bug. Stats TSV rows (`Category<TAB>Name<TAB>Value`):
  ```
  ReadStats	PassReads	2000.0
  ReadStats	LenFail	63.0
  ReadStats	QcFail	0.0
  Classification	Primers_found	1166.0
  Classification	Rescue	28.0
  Classification	Unusable	820.0
  ```

---

## File Structure

- `rnaforge/chopper.py` — **create**: `run_chopper` + `parse_kept` + errors.
- `rnaforge/pychopper.py` — **create**: `PychopperStats` dataclass, `parse_pychopper_stats`, `run_pychopper` (crash-tolerant) + errors.
- `rnaforge/modules/m03_trim.py` — **modify**: extract `_trim_short`; add `_trim_long`; dispatch on `resolve_read_type`; drop `require_short_read`.
- Tests: `tests/test_chopper.py`, `tests/test_pychopper.py` (create), `tests/test_m03_trim.py` (add dispatch tests).

---

### Task 1: `chopper.py` — length/quality filter runner

**Files:**
- Create: `rnaforge/chopper.py`
- Test: `tests/test_chopper.py`

**Interfaces:**
- Produces:
  - `class ChopperRunError(RuntimeError)`
  - `parse_kept(stderr: str) -> int | None` — pulls N from `Kept N reads out of M reads`; None if absent.
  - `run_chopper(in_fastq: Path, out_fastq: Path, env: str = "rnaforge-longread", min_qual: int = 7, min_len: int = 50, threads: int = 4) -> int` — decompresses gz input if needed, runs chopper stdin→stdout, returns reads kept (from stderr, else counts output lines/4). Raises `ChopperRunError` on non-zero exit or missing output.

- [ ] **Step 1: Write the failing test**

Create `tests/test_chopper.py`:

```python
from __future__ import annotations

import shutil

import pytest

from rnaforge.chopper import ChopperRunError, parse_kept, run_chopper

_HAS = shutil.which("conda") is not None


def test_parse_kept_reads_number():
    assert parse_kept("Kept 1097 reads out of 1109 reads\n") == 1097


def test_parse_kept_absent_is_none():
    assert parse_kept("no summary here") is None


@pytest.mark.skipif(not _HAS, reason="conda/rnaforge-longread not available")
def test_run_chopper_filters(tmp_path):
    fq = tmp_path / "in.fastq"
    # 2 long clean reads pass; 1 tiny read fails -l 50
    fq.write_text(
        "@a\n" + "ACGT" * 40 + "\n+\n" + "I" * 160 + "\n"
        "@b\n" + "ACGT" * 40 + "\n+\n" + "I" * 160 + "\n"
        "@c\nACGT\n+\nIIII\n"
    )
    out = tmp_path / "out.fastq"
    kept = run_chopper(fq, out, min_qual=7, min_len=50)
    assert kept == 2
    assert out.exists()
```

- [ ] **Step 2: Run to verify it fails**

Run: `conda run -n rnaforge-core --cwd /home/ali/rnaforge-pipeline python -m pytest tests/test_chopper.py -q`
Expected: FAIL — `ModuleNotFoundError: No module named 'rnaforge.chopper'`.

- [ ] **Step 3: Implement**

Create `rnaforge/chopper.py`:

```python
"""chopper (ONT length/quality filter) runner. Reads stdin, writes stdout."""
from __future__ import annotations

import re
import subprocess
from pathlib import Path


class ChopperRunError(RuntimeError):
    """chopper failed to run."""


def parse_kept(stderr: str) -> int | None:
    m = re.search(r"Kept (\d+) reads out of", stderr)
    return int(m.group(1)) if m else None


def run_chopper(in_fastq: Path, out_fastq: Path, env: str = "rnaforge-longread",
                min_qual: int = 7, min_len: int = 50, threads: int = 4) -> int:
    in_fastq = Path(in_fastq)
    out_fastq = Path(out_fastq)
    out_fastq.parent.mkdir(parents=True, exist_ok=True)
    decomp = "zcat" if in_fastq.name.endswith(".gz") else "cat"
    pipe = (
        f"{decomp} {in_fastq!s} | "
        f"chopper -q {min_qual} -l {min_len} --threads {threads} > {out_fastq!s}"
    )
    cmd = ["conda", "run", "-n", env, "bash", "-lc", pipe]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0 or not out_fastq.exists():
        raise ChopperRunError(
            f"chopper failed (exit {proc.returncode}) on {in_fastq}\n"
            f"stderr: {proc.stderr[-500:]}"
        )
    kept = parse_kept(proc.stderr)
    if kept is None:
        kept = sum(1 for _ in out_fastq.open()) // 4
    return kept
```

- [ ] **Step 4: Run to verify it passes**

Run: `conda run -n rnaforge-core --cwd /home/ali/rnaforge-pipeline python -m pytest tests/test_chopper.py -q`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add rnaforge/chopper.py tests/test_chopper.py
git commit -m "feat(chopper): ONT length/quality filter runner"
```

---

### Task 2: `pychopper.py` — full-length cDNA orient/trim runner (crash-tolerant)

**Files:**
- Create: `rnaforge/pychopper.py`
- Test: `tests/test_pychopper.py`

**Interfaces:**
- Produces:
  - `class PychopperParseError(ValueError)`, `class PychopperRunError(RuntimeError)`
  - `@dataclass(frozen=True) class PychopperStats(pass_reads: int, primers_found: int, rescue: int, unusable: int, len_fail: int)`
  - `parse_pychopper_stats(tsv_text: str) -> PychopperStats` — reads the `Category<TAB>Name<TAB>Value` TSV.
  - `run_pychopper(in_fastq: Path, out_fastq: Path, stats_tsv: Path, env: str = "rnaforge-longread", kit: str | None = None, threads: int = 4) -> PychopperStats` — runs pychopper; on the known `_plot_stats` crash (outputs present) logs a warning and proceeds; otherwise raises. gz input is decompressed to a temp file first (pychopper reads plain fastx most reliably). Returns parsed stats.

- [ ] **Step 1: Write the failing tests**

Create `tests/test_pychopper.py`:

```python
from __future__ import annotations

import shutil

import pytest

from rnaforge.pychopper import (
    PychopperParseError, PychopperStats, parse_pychopper_stats, run_pychopper,
)

_HAS = shutil.which("conda") is not None

_STATS = (
    "Category\tName\tValue\n"
    "ReadStats\tPassReads\t2000.0\n"
    "ReadStats\tLenFail\t63.0\n"
    "ReadStats\tQcFail\t0.0\n"
    "Classification\tPrimers_found\t1166.0\n"
    "Classification\tRescue\t28.0\n"
    "Classification\tUnusable\t820.0\n"
    "Strand\t+\t419.0\n"
)


def test_parse_pychopper_stats():
    s = parse_pychopper_stats(_STATS)
    assert isinstance(s, PychopperStats)
    assert s.pass_reads == 2000
    assert s.primers_found == 1166
    assert s.rescue == 28
    assert s.unusable == 820
    assert s.len_fail == 63


def test_parse_pychopper_stats_missing_raises():
    with pytest.raises(PychopperParseError):
        parse_pychopper_stats("Category\tName\tValue\nReadStats\tPassReads\t10.0\n")


@pytest.mark.skipif(not _HAS, reason="conda/rnaforge-longread not available")
def test_run_pychopper_tolerates_plot_crash(tmp_path):
    # ~200 reads is enough for pychopper to classify some full-length reads
    fq = tmp_path / "in.fastq"
    body = "".join(
        f"@r{i}\n" + "ACGT" * 60 + "\n+\n" + "I" * 240 + "\n" for i in range(200)
    )
    fq.write_text(body)
    stats = run_pychopper(fq, tmp_path / "fl.fastq", tmp_path / "stats.tsv")
    assert isinstance(stats, PychopperStats)
    assert stats.pass_reads == 200
    assert (tmp_path / "fl.fastq").exists()
```

- [ ] **Step 2: Run to verify it fails**

Run: `conda run -n rnaforge-core --cwd /home/ali/rnaforge-pipeline python -m pytest tests/test_pychopper.py -q`
Expected: FAIL — `ModuleNotFoundError`.

- [ ] **Step 3: Implement**

Create `rnaforge/pychopper.py`:

```python
"""Pychopper (ONT full-length cDNA orient/trim) runner + stats parser.

Pychopper 2.7.10 writes the oriented FASTQ and the -S stats TSV, then crashes
in its end-of-run PDF report under pandas 3 (_plot_stats -> float(Series)).
The real work is complete before that crash, so run_pychopper tolerates exactly
that signature and raises on anything else (feedback_gurultulu_hata)."""
from __future__ import annotations

import gzip
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path


class PychopperParseError(ValueError):
    """pychopper stats TSV could not be parsed."""


class PychopperRunError(RuntimeError):
    """pychopper failed for a reason other than the known plotting crash."""


@dataclass(frozen=True)
class PychopperStats:
    pass_reads: int
    primers_found: int
    rescue: int
    unusable: int
    len_fail: int


def parse_pychopper_stats(tsv_text: str) -> PychopperStats:
    table: dict[tuple[str, str], float] = {}
    for line in tsv_text.splitlines():
        parts = line.split("\t")
        if len(parts) != 3 or parts[0] == "Category":
            continue
        try:
            table[(parts[0], parts[1])] = float(parts[2])
        except ValueError:
            continue

    def get(cat: str, name: str) -> int:
        if (cat, name) not in table:
            raise PychopperParseError(f"pychopper stats missing {cat}/{name}")
        return int(table[(cat, name)])

    return PychopperStats(
        pass_reads=get("ReadStats", "PassReads"),
        primers_found=get("Classification", "Primers_found"),
        rescue=get("Classification", "Rescue"),
        unusable=get("Classification", "Unusable"),
        len_fail=get("ReadStats", "LenFail"),
    )


def run_pychopper(in_fastq: Path, out_fastq: Path, stats_tsv: Path,
                  env: str = "rnaforge-longread", kit: str | None = None,
                  threads: int = 4) -> PychopperStats:
    in_fastq = Path(in_fastq)
    out_fastq = Path(out_fastq)
    stats_tsv = Path(stats_tsv)
    out_fastq.parent.mkdir(parents=True, exist_ok=True)

    tmp = None
    src = in_fastq
    if in_fastq.name.endswith(".gz"):
        tmp = Path(tempfile.mkstemp(suffix=".fastq")[1])
        with gzip.open(in_fastq, "rt") as fh, tmp.open("w") as out:
            shutil.copyfileobj(fh, out)
        src = tmp
    try:
        cmd = ["conda", "run", "-n", env, "pychopper",
               "-t", str(threads), "-S", str(stats_tsv)]
        if kit:
            cmd += ["-k", kit]
        cmd += [str(src), str(out_fastq)]
        proc = subprocess.run(cmd, capture_output=True, text=True)
    finally:
        if tmp is not None:
            tmp.unlink(missing_ok=True)

    outputs_ok = out_fastq.exists() and stats_tsv.exists()
    known_plot_crash = "_plot_stats" in proc.stderr
    if proc.returncode != 0 and not (outputs_ok and known_plot_crash):
        raise PychopperRunError(
            f"pychopper failed (exit {proc.returncode}) on {in_fastq}\n"
            f"stderr: {proc.stderr[-800:]}"
        )
    if proc.returncode != 0:
        sys.stderr.write(
            f"WARNING: pychopper exited {proc.returncode} on its PDF report "
            f"(known pandas-3 _plot_stats bug); outputs are complete, continuing.\n"
        )
    return parse_pychopper_stats(stats_tsv.read_text())
```

- [ ] **Step 4: Run to verify it passes**

Run: `conda run -n rnaforge-core --cwd /home/ali/rnaforge-pipeline python -m pytest tests/test_pychopper.py -q`
Expected: PASS (integration actually runs pychopper; ~10–30s).

- [ ] **Step 5: Commit**

```bash
git add rnaforge/pychopper.py tests/test_pychopper.py
git commit -m "feat(pychopper): full-length cDNA orient/trim runner (crash-tolerant) + stats"
```

---

### Task 3: m03 dispatch — short→fastp, long→Pychopper/chopper by chemistry

**Files:**
- Modify: `rnaforge/modules/m03_trim.py`
- Test: `tests/test_m03_trim.py`

**Interfaces:**
- Consumes: `resolve_read_type` (Step 1), `run_pychopper`/`PychopperStats` (Task 2), `run_chopper` (Task 1), `trimmed_reads` (existing).
- Produces: `run_trim` dispatches on read_type. Long writes `statistics/trimming_statistics.json` shape `{"read_type": "long", "chemistry": <cdna|direct_rna>, "n_samples": N, "samples": {sid: {reads_before, reads_after, survival_rate}}}`, writes the filtered FASTQ to `trimmed_reads(...)[0]`, and writes NO gates (diagnostic). Short output unchanged plus `"read_type": "short"`.

- [ ] **Step 1: Write the failing tests**

Add to `tests/test_m03_trim.py` (mirror the file's short-read setup; monkeypatch `run_pychopper`/`run_chopper`). Seed a long run dir directly:

```python
def _seed_long(tmp_path, chemistry="cdna"):
    import json
    from rnaforge.state import RunState
    run_dir = tmp_path / "run"
    (run_dir / "statistics").mkdir(parents=True)
    (run_dir / "statistics" / "raw_statistics.json").write_text(
        json.dumps({"read_type": "long", "chemistry": chemistry})
    )
    RunState(run_dir).mark_done("m01_validate", [])
    return run_dir


def _long_cfg(chemistry="cdna"):
    from rnaforge.config import (
        Config, Reference, Library, Trimming, DE, Report, Resources,
    )
    return Config(
        organism="E. coli", organism_type="prokaryote", platform="auto",
        reference=Reference(), library=Library(chemistry=chemistry),
        trimming=Trimming(), de=DE(), report=Report(), resources=Resources(),
    )


def test_run_trim_long_cdna_pychopper_then_chopper(tmp_path, monkeypatch):
    import rnaforge.modules.m03_trim as m03
    from rnaforge.pychopper import PychopperStats
    run_dir = _seed_long(tmp_path, "cdna")
    fq = tmp_path / "s1.fastq"
    fq.write_text("@r\n" + "ACGT" * 50 + "\n+\n" + "I" * 200 + "\n")
    meta = tmp_path / "m.tsv"
    meta.write_text(f"sample_id\tcondition\tfastq_1\ns1\tctrl\t{fq}\n")

    calls = []

    def fake_pychopper(in_fastq, out_fastq, stats_tsv, **k):
        calls.append("pychopper")
        Path(out_fastq).parent.mkdir(parents=True, exist_ok=True)
        Path(out_fastq).write_text("@r\nACGT\n+\nIIII\n")
        return PychopperStats(pass_reads=100, primers_found=60,
                              rescue=5, unusable=35, len_fail=3)

    def fake_chopper(in_fastq, out_fastq, **k):
        calls.append("chopper")
        Path(out_fastq).parent.mkdir(parents=True, exist_ok=True)
        Path(out_fastq).write_text("@r\nACGT\n+\nIIII\n")
        return 57

    monkeypatch.setattr(m03, "run_pychopper", fake_pychopper)
    monkeypatch.setattr(m03, "run_chopper", fake_chopper)

    summary = m03.run_trim(_long_cfg("cdna"), meta, run_dir)
    assert summary["read_type"] == "long"
    assert summary["chemistry"] == "cdna"
    assert calls == ["pychopper", "chopper"]      # order matters
    assert summary["samples"]["s1"]["reads_after"] == 57
    # diagnostic: no FAIL gate written for the long branch
    assert not (run_dir / "quality" / "gates.json").exists()
    # trimmed output lands on the m04 contract path
    from rnaforge.modules.m03_trim import trimmed_reads
    from rnaforge.metadata import load_metadata
    out1, out2 = trimmed_reads(run_dir, load_metadata(meta)[0])
    assert out1.exists() and out2 is None


def test_run_trim_long_direct_rna_chopper_only(tmp_path, monkeypatch):
    import rnaforge.modules.m03_trim as m03
    run_dir = _seed_long(tmp_path, "direct_rna")
    fq = tmp_path / "s1.fastq"
    fq.write_text("@r\n" + "ACGT" * 50 + "\n+\n" + "I" * 200 + "\n")
    meta = tmp_path / "m.tsv"
    meta.write_text(f"sample_id\tcondition\tfastq_1\ns1\tctrl\t{fq}\n")

    calls = []
    monkeypatch.setattr(m03, "run_pychopper",
                        lambda *a, **k: calls.append("pychopper"))

    def fake_chopper(in_fastq, out_fastq, **k):
        calls.append("chopper")
        Path(out_fastq).parent.mkdir(parents=True, exist_ok=True)
        Path(out_fastq).write_text("@r\nACGT\n+\nIIII\n")
        return 90

    monkeypatch.setattr(m03, "run_chopper", fake_chopper)
    summary = m03.run_trim(_long_cfg("direct_rna"), meta, run_dir)
    assert calls == ["chopper"]                   # no pychopper for direct-RNA
    assert summary["samples"]["s1"]["reads_after"] == 90
```

Also add `assert summary["read_type"] == "short"` to the file's existing short-read `run_trim` success test.

- [ ] **Step 2: Run to verify the long tests fail**

Run: `conda run -n rnaforge-core --cwd /home/ali/rnaforge-pipeline python -m pytest tests/test_m03_trim.py -k "long or direct" -q`
Expected: FAIL — `require_short_read` raises `NotImplementedError` / `run_pychopper` not an attribute of m03.

- [ ] **Step 3: Refactor m03 into a dispatch**

In `rnaforge/modules/m03_trim.py`:

Imports — drop `require_short_read`, add:

```python
from rnaforge.chopper import run_chopper
from rnaforge.config import CHEMISTRY, Config
from rnaforge.pychopper import run_pychopper
from rnaforge.routing import resolve_read_type
```

Replace the `require_short_read(...)` line + the fastp body + the trailing `mark_done`/`return` with a dispatch, mirroring m02-long. After the m01-done precondition:

```python
    read_type = resolve_read_type(run_dir)
    if read_type == "long":
        summary = _trim_long(config, metadata_path, run_dir,
                             trimmed_dir, stats_dir, logs_dir, state)
    else:
        summary = _trim_short(config, metadata_path, run_dir,
                             trimmed_dir, stats_dir, logs_dir, state)
    state.mark_done(MODULE_NAME, [str(stats_path), str(logs_dir / "trim.log")])
    return summary
```

Move the existing fastp body (`profile = load_profile(...)` through `raise_if_failed(gates)` and the stats write) into:

```python
def _trim_short(config: Config, metadata_path: Path, run_dir: Path,
                trimmed_dir: Path, stats_dir: Path, logs_dir: Path,
                state: RunState) -> dict:
    """Kısa-okuma trimming (fastp). survival_rate FAIL kapısı korunur."""
    stats_path = stats_dir / "trimming_statistics.json"
    profile = load_profile(config.organism_type, config.quality)
    log_path = logs_dir / "trim.log"
    with log_path.open("w") as log_file:
        ...  # (existing body verbatim)
        summary = {
            "read_type": "short",
            "n_samples": len(samples),
            "samples": per_sample,
            "gate_counts": dict(Counter(g.status for g in gates)),
        }
        stats_path.write_text(json.dumps(summary, indent=2))
        write_gate_results(run_dir, gates)
        for g in gates:
            log(f"gate {g.name}: {g.status} — {g.message}")
        raise_if_failed(gates)
        log(f"trimming statistics written: {stats_path}")
    return summary
```

Add the long branch:

```python
def _trim_long(config: Config, metadata_path: Path, run_dir: Path,
               trimmed_dir: Path, stats_dir: Path, logs_dir: Path,
               state: RunState) -> dict:
    """Uzun-okuma ön-işleme. cdna: Pychopper (yönlendir/kes) + chopper (filtre);
    direct_rna: yalnız chopper. Diagnostik — FAIL kapısı yok (long profil Step 6)."""
    chemistry = config.library.chemistry
    if chemistry not in CHEMISTRY:
        raise ValueError(
            f"long-read trim requires library.chemistry in {CHEMISTRY}, "
            f"got {chemistry!r} (m01 should have enforced this for ONT)."
        )
    stats_path = stats_dir / "trimming_statistics.json"
    log_path = logs_dir / "trim.log"
    min_len = config.trimming.min_length
    with log_path.open("w") as log_file:
        def log(msg: str) -> None:
            log_file.write(msg + "\n"); log_file.flush()

        samples = load_metadata(metadata_path)
        log(f"m03 long-read ({chemistry}): {len(samples)} sample(s)")
        per_sample: dict[str, dict] = {}
        for sample in samples:
            state.heartbeat()
            out1, _ = trimmed_reads(run_dir, sample)
            out1.parent.mkdir(parents=True, exist_ok=True)
            reads_before = _count_fastx(sample.fastq_1)

            if chemistry == "cdna":
                work = out1.parent / "pychopper_full_length.fastq"
                stats_tsv = out1.parent / "pychopper_stats.tsv"
                ps = run_pychopper(sample.fastq_1, work, stats_tsv)
                reads_after = run_chopper(work, out1, min_len=min_len)
                log(f"{sample.sample_id}: pychopper primers={ps.primers_found} "
                    f"rescue={ps.rescue} unusable={ps.unusable}; chopper kept {reads_after}")
            else:  # direct_rna
                reads_after = run_chopper(sample.fastq_1, out1, min_len=min_len)
                log(f"{sample.sample_id}: chopper kept {reads_after}")

            survival = reads_after / reads_before if reads_before else 0.0
            per_sample[sample.sample_id] = {
                "reads_before": reads_before,
                "reads_after": reads_after,
                "survival_rate": round(survival, 4),
            }

        summary = {
            "read_type": "long",
            "chemistry": chemistry,
            "n_samples": len(samples),
            "samples": per_sample,
        }
        stats_path.write_text(json.dumps(summary, indent=2))
        log(f"trimming statistics written: {stats_path}")
    return summary


def _count_fastx(path: Path) -> int:
    import gzip as _gz
    opener = _gz.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        return sum(1 for _ in fh) // 4
```

- [ ] **Step 4: Run to verify all m03 tests pass**

Run: `conda run -n rnaforge-core --cwd /home/ali/rnaforge-pipeline python -m pytest tests/test_m03_trim.py -q`
Expected: PASS (long cdna + direct_rna dispatch + existing short-read regression).

- [ ] **Step 5: Commit**

```bash
git add rnaforge/modules/m03_trim.py tests/test_m03_trim.py
git commit -m "feat(m03): dispatch read_type — fastp (short) / Pychopper+chopper (long)"
```

---

### Task 4: Full suite + live smoke on microbepore cDNA

- [ ] **Step 1: Full suite green**

Run: `conda run -n rnaforge-core --cwd /home/ali/rnaforge-pipeline python -m pytest -q`
Expected: PASS.

- [ ] **Step 2: Live smoke — real ONT cDNA reaches m03 preprocessing**

Reuse the `runs/*_mbp_smoke` run (validate + qc already done there) or a fresh subset run; then:

```bash
rnaforge trim --config <ont_cfg> --metadata <ont_meta> --runs-dir runs --run-id mbp_smoke
```

Expected: `_trim_long` runs Pychopper (loud warning about the PDF crash, but proceeds) then chopper; `trimming_statistics.json` has `read_type=long`, `chemistry=cdna`, per-sample survival; the trimmed FASTQ exists at `trimmed_reads(...)`. `rnaforge quant` still stops loudly (`NotImplementedError`, long-read align not built) — confirming the guard chain below m03.

- [ ] **Step 3: Record** survival numbers in the commit message / DURUM.

---

## Out of scope (later steps)
- m04-long (minimap2), m05-long (featureCounts -L) — own plans.
- Long-read quality profile + gates (Step 6): m03-long stays gate-free.
- Pychopper kit / chopper quality as config fields — defaults for now (kit=None autodetect, min_qual=7, min_len from `config.trimming.min_length`).

---

## Self-Review

**Spec coverage (spec row "m03 ön-işleme: cDNA→Pychopper+chopper; direct-RNA→chopper"):** Task 3 chemistry branch + Tasks 1–2 tools implement exactly this ✓. Diagnostic, no long FAIL gate (deferred to Step 6) ✓. Router convention ✓.

**Placeholder scan:** `_trim_short` body is described as "existing body verbatim" with the summary/gate lines shown — the implementer moves the current fastp block unchanged; not a silent placeholder. All new code is concrete. ✓

**Type consistency:** `PychopperStats` (Task 2) fields used in `_trim_long` (Task 3) logging. `run_chopper(...) -> int` and `run_pychopper(...) -> PychopperStats` signatures match call sites and monkeypatches. `trimmed_reads(run_dir, sample) -> (Path, None)` reused as the output contract. `CHEMISTRY` imported from config (Step 1). `resolve_read_type` from routing. ✓

**Risk:** the crash-tolerant `run_pychopper` is the delicate part — it is unit-tested against the real tool (Task 2 Step 1 integration) and the live smoke (Task 4) confirms it end-to-end; any non-plotting failure still raises.
