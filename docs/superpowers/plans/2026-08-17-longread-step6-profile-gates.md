# Long-Read Arm — Step 6: Long-read quality profile + gates

## Context
Steps 1–5 built the long-read arm through the count matrix; m03/m04/m05 long branches are all
**diagnostic** (no gates) because the Illumina thresholds (survival 0.50, alignment 0.70,
assignment 0.50) would wrongly FAIL legitimate ONT data (the mbp_smoke live run showed
alignment 0.71–0.81 and assignment 0.06–0.16). Step 6 adds an ONT-appropriate quality profile
and wires gates into the long branches, preserving the FAIL=invalid / WARN=suspect policy.

**Principle guard (`feedback_dogruluk_kontrol`, "yalancı sonuç asla istemem"):** thresholds must
be defensible, NOT fabricated. ONT is inherently noisier (Q10–15, not Q30); the long profile is
**deliberately permissive and stamped** (like the eukaryote profile), and only a catastrophic
signal (wrong reference) FAILs. Report writes out the profile + thresholds — no silent loosening.

## Design
- **New `rnaforge/profiles/prokaryote_long.yml`** (`permissive: true`, stamped description):
  - `alignment_rate: 0.50` → **FAIL** floor. Wrong reference/contamination is invalid on any
    platform; clean ONT cDNA maps >0.85, so 0.50 is a catastrophic-only floor.
  - `survival_rate: 0.20` → **WARN**. Pychopper discards non-full-length by design → low survival
    is normal; only near-total loss is suspect.
  - `assignment_rate: 0.05` → **WARN**. CDS-only counting on rRNA-replete ONT is legitimately low;
    only near-total mis-assignment (<5%) is suspect. Never FAIL (would false-reject good ONT).
  - Include the remaining gate names (read_depth/base_quality/rrna/dedup/replicate_correlation)
    at permissive values so the profile is complete and the confidence card stamps consistently.
- **`quality.profile_name_for(organism_type, read_type)`** → `f"{organism_type}_long"` for long,
  else `organism_type`. Single source for the read-type→profile mapping (modules + CLI).
- **`gates.build_trim_gates` + `gates`-side `build_count_gates` gain `warn_only: bool = False`**:
  when True and below threshold, status = WARN (not FAIL). Refactor `build_trim_gates` to take
  `survival_rates: dict[str, float]` (mirror `build_count_gates`) so the long branch can pass a
  plain rate dict.
- **Wire long gates:**
  - `m04_quant._quant_long`: `build_alignment_gates(results, long_profile)` (FAIL) →
    write_gate_results → raise_if_failed. Add `gate_counts` to summary.
  - `m03_trim._trim_long`: `build_trim_gates(survivals, long_profile, warn_only=True)` (WARN) →
    write_gate_results, no raise. Add `gate_counts`.
  - `m05_counts._counts_long`: `build_count_gates(assignments, long_profile, warn_only=True)`
    (WARN) → write_gate_results, no raise. Add `gate_counts`.
  - m02-long stays diagnostic (parity with short m02, which is FastQC-derived, no threshold gate).
- **CLI `_load_run_profile(config, run_dir)`**: resolve read_type (fallback short if pre-m01) →
  `load_profile(profile_name_for(...), config.quality)`. Use in `_cmd_trim/_cmd_quant/_cmd_counts`
  so a long run's confidence card stamps `prokaryote_long` (permissive), not `prokaryote`.

## Tasks (TDD)
1. **profile + `profile_name_for`**: test `load_profile("prokaryote_long")` thresholds +
   `permissive`; `profile_name_for("prokaryote","long")=="prokaryote_long"`, short passthrough.
   Commit `feat(profile): prokaryote_long permissive ONT profile + profile_name_for`.
2. **`warn_only` in build_trim_gates/build_count_gates** + trim gates float-based: test below
   threshold with `warn_only=True` → WARN not FAIL; short calls unchanged. Update m03-short call
   site to pass a rate dict. Commit `feat(gates): warn_only + float-based build_trim_gates`.
3. **wire long gates** (m03/m04/m05 long branches): update the Step 4/5 long dispatch tests
   (which asserted "no gates.json") to expect the new gates (m04 alignment PASS at rate>0.50,
   m03 survival PASS, m05 assignment PASS/WARN); add a low-rate case → m04-long FAIL raises,
   m05-long low → WARN (no raise). Commit `feat(m03/m04/m05): long-read gates from prokaryote_long`.
4. **CLI `_load_run_profile`** + full suite + live smoke on mbp_smoke (re-run trim/quant/counts →
   card stamps prokaryote_long; alignment PASS, survival/assignment PASS or honest WARN). Merge +
   push + DURUM/memory. Commit `feat(cli): stamp long-read profile on the confidence card`.

## Out of scope
- Step 7: report read_type badge + end-to-end smoke.
- eukaryote_long profile (euk arm not built).
- m02-long NanoPlot quality WARN (m02 is diagnostic on both arms).
- DE-signal ONT dataset selection (candidate B) — thresholds are calibratable when it lands.

## Self-Review
- FAIL=invalid / WARN=suspect preserved ✓. Only alignment (wrong-reference) FAILs on long ✓.
- Permissive + stamped, thresholds in the card/report — no fabricated strict thresholds ✓.
- mbp_smoke (align 0.71–0.81, survival 0.54–0.63, assignment 0.06–0.16) → PASS/PASS/PASS at the
  chosen floors ✓ (no false FAIL). Profile stamped permissive so the report is honest ✓.
