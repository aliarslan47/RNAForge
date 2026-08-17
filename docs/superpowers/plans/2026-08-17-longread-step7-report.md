# Long-Read Arm — Step 7 (final): report read_type badge + long-read method/tool notes

## Context
Steps 1–6 built the long-read arm through gated counting; a long run flows to DESeq2 and the HTML
report. But the report still described the **short-read** tool chain (FastQC/fastp/Bowtie2) for
every run. Step 7 makes the report honest about read_type: a badge, the long-read software table,
the long-read methods narrative, and long-read citations — so a long run's report cites the tools
that actually ran and none that didn't.

## Design (all in `rnaforge/report_html.py`)
- **read_type badge** in `section_dataset`: `Okuma tipi: uzun` / `Read type: long` next to platform
  (glossed short→kısa/short, long→uzun/long). New `read_type` label in both LABELS.
- **Software table** (`_SOFTWARE` + `section_software`): platform tools gain a `cond` of `"short"` /
  `"long"`; new long entries NanoPlot 1.47.1, Pychopper 2.7.10, chopper 0.13.0, minimap2 2.31.
  featureCounts stays shared with a "-L for long-read" note. `render_report` passes
  `short`/`long` flags from the run's read_type.
- **Methods narrative**: new `_METHODS_TEXT_LONG` (tr/en) describing NanoPlot → Pychopper+chopper →
  minimap2 (-ax map-ont/map-hifi) → featureCounts -L, sharing the DESeq2 tail; notes the permissive
  `prokaryote_long` profile. `section_methods` gains `read_type` and selects the narrative.
- **References**: `_REFERENCES` split into read-type-agnostic base + `_REFERENCES_SHORT`
  (FastQC/fastp/Bowtie2/Williams) + `_REFERENCES_LONG` (NanoPack2 btad311, minimap2 bty191,
  Pychopper repo). `section_references` + `render_report` pick by read_type — no citation of unused
  tools (honesty; `feedback_dogruluk_kontrol`).

## Verification
- Unit (`tests/test_report_longread.py`): dataset badge; software long lists minimap2/NanoPlot/
  Pychopper and hides Bowtie2/FastQC (and vice-versa for short); methods long says minimap2 not
  Bowtie2; references long cite NanoPack/minimap2 not FastQC/Bowtie2.
- Full suite 469 green.
- **Live end-to-end (mbp_smoke, long ONT):** figures → report OK; report.html stamps
  profile=prokaryote_long, shows the read_type badge, lists minimap2/NanoPlot/Pychopper + featureCounts
  -L, cites NanoPack2/minimap2, and contains NO FastQC/Bowtie2 mention. 582 KB.

## Result
**Long-read (ONT/PacBio) arm COMPLETE** (Steps 1–7): detect→route · NanoPlot QC · Pychopper+chopper ·
minimap2 · featureCounts -L · permissive gated profile · honest report. Converges on the same
count-matrix contract → m06+ unchanged.

## Out of scope
- DE-signal ONT dataset (candidate B, *E. coli* glucose-vs-pyruvate) for a biology-meaningful smoke.
- eukaryote long arm; PacBio HiFi live validation.
