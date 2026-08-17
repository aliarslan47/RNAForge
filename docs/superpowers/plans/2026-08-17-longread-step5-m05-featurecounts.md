# Long-Read Arm — Step 5: m05-long (featureCounts `-L`) Implementation Plan

**Goal:** Give m05 (`counts`) a long-read branch: when `read_type == "long"`, count the minimap2
BAMs at gene level with **featureCounts `-L`** (long-read mode), producing the SAME
`quantification/counts.tsv` gene×sample contract that m06 (DESeq2) consumes. Short-read runs keep
their exact behaviour + `assignment_rate` FAIL gate. Long is **diagnostic** (assignment_rate
recorded, NO FAIL gate — mirrors m03-long/m04-long; long profile/gates = Step 6). This removes the
Step-1 `require_short_read(run_dir, "counts")` guard in m05 → a long run now reaches a full count
matrix, and the arm's only remaining stop is the report/profile work.

**Architecture:** mirror the m04 read_type dispatch. `m05_counts.run_counts` resolves `read_type`
and dispatches `short → _counts_short` / `long → _counts_long`. featureCounts is the SAME binary
(subread 2.1.1 in `rnaforge-quant-prok`) — `run_featurecounts` gains a `long_read: bool = False`
param that adds `-L` and forces `paired=False` (long reads are single-molecule; `-L`/`-p`
incompatible). The counts.tsv + tpm/fpkm writing is shared via a helper.

**Tech stack:** Python 3.11, subread/featureCounts 2.1.1 (`rnaforge-quant-prok`), pytest.

## Files
- `rnaforge/featurecounts.py` — **modify**: `run_featurecounts(..., long_read=False)` → `-L`,
  `paired` ignored when `long_read`.
- `rnaforge/modules/m05_counts.py` — **modify**: extract `_counts_short` + shared
  `_write_count_outputs`; add `_counts_long`; dispatch on `resolve_read_type`; drop
  `require_short_read`.
- Tests: `tests/test_featurecounts.py` (add `-L` cmd test), `tests/test_m05_counts.py` (add long
  dispatch tests).

## Tasks (TDD, each: red → green → commit)

### Task 1: `run_featurecounts(long_read=True)` adds `-L`
- Test: monkeypatch `subprocess.run` to capture argv; assert `-L` present and `-p` absent when
  `long_read=True`; assert `-L` absent when `long_read=False` (regression). Reuse the module's
  existing capture style if present, else a fake CompletedProcess writing a minimal counts.txt +
  .summary so parsing succeeds.
- Implement: `if long_read: cmd += ["-L"]` and `paired = paired and not long_read`.
- Commit: `feat(featurecounts): long-read -L mode`.

### Task 2: m05 dispatch — short (gate) / long (diagnostic)
- Refactor `run_counts`: keep resume + m04 precondition outer; `read_type = resolve_read_type(...)`
  replaces `require_short_read`; dispatch.
- `_counts_short`: current body verbatim + `"read_type":"short"` in summary (paired detection,
  featureCounts, empty-matrix raise, `assignment_rate` FAIL gate).
- `_counts_long`: featureCounts `long_read=True` (paired=False); write counts.tsv + tpm/fpkm via the
  shared helper; empty-matrix raise kept; record per-sample `assignment_rate` **diagnostically**;
  **NO gate** (`write_gate_results` not called). Summary `{"read_type":"long", "platform":<p>,
  "n_samples":N,"n_genes":G,"samples":{sid:{assignment_rate}}, "expression_values":[...]}`.
- Shared `_write_count_outputs(result, sample_ids, quant_dir, log)`: positional column→sample_id,
  counts.tsv (`gene\t<sid...>`), tpm/fpkm when lengths present. Both branches call it.
- Tests (`tests/test_m05_counts.py`): seed a long run dir (raw_statistics platform=ont/read_type=long,
  RunState mark_done m01..m04, fake BAMs at `quantification/<sid>/aligned.sorted.bam`), monkeypatch
  `run_featurecounts` capturing `long_read`; assert dispatch passes `long_read=True`, counts.tsv
  written with sample_id headers, summary `read_type=="long"`, **no `quality/gates.json`**. Add
  `assert summary["read_type"]=="short"` to the existing short success test.
- Commit: `feat(m05): dispatch read_type — featureCounts (short gate) / -L (long diagnostic)`.

### Task 3: full suite + live smoke on mbp_smoke
- Full suite green.
- Live: `rnaforge counts` on `runs/*_mbp_smoke` (m04-long BAMs already there) → counts.tsv written,
  `count_statistics.json` read_type=long, per-sample assignment_rate; then `rnaforge de` is the
  next stage (m06 is organism/read-type agnostic — it consumes counts.tsv). Record gene count +
  assignment rates in DURUM.
- Merge to main + push; update DURUM + `reminder_rnaforge_longread` memory (Step 5 done, Step 6 next).

## Out of scope
- Long-read quality profile + FAIL gates (Step 6). m05-long stays gate-free until then.
- Report read_type badge + end-to-end DE smoke (Step 7).
- DE-signal ONT dataset selection (microbepore is single-condition; candidate B glucose-vs-pyruvate).

## Self-Review
- Spec row "m05 sayım: featureCounts `-L`" ✓. Gene-level count-matrix contract preserved → m06+
  unchanged ✓. Diagnostic-only matches m03/m04-long precedent ✓. Guard removed in m05 (was the last
  Step-1 guard) → long arm reaches counts.tsv ✓.
