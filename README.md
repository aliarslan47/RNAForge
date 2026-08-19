# RNAForge

Reproducible, modular **bulk RNA-seq** analysis pipeline — from raw FASTQ to a single,
self-contained HTML report, with a full functional-analysis layer on top of differential expression.

Turkish version: [README.tr.md](README.tr.md) · Reference document: [PLAN.md](PLAN.md) (v1.4)

[![Pipeline DAG](https://img.shields.io/badge/pipeline-DAG-0d6b8f)](https://aliarslan47.github.io/RNAForge/pipeline_architecture.html)
[![organism](https://img.shields.io/badge/organism-prokaryote%20%C2%B7%20eukaryote-2f8f5b)](https://aliarslan47.github.io/RNAForge/pipeline_architecture.html)
[![reads](https://img.shields.io/badge/reads-short%20%C2%B7%20long-c07211)](https://aliarslan47.github.io/RNAForge/pipeline_architecture.html)

## What it does

A staged pipeline that takes raw reads to biology:

```
validate → qc → trim → quant → counts → de → figures → report
                                          └→ enrich · kegg · gsea · semantic · amr · operon · ppi
```

The full pipeline as an interactive, bilingual node-graph (organism × read-type branching, `m00`–`m18`) —
[**rendered diagram**](https://aliarslan47.github.io/RNAForge/pipeline_architecture.html) · source: `docs/pipeline_architecture.html`.

- **Core**: input/design validation, QC, gentle trimming, alignment/quantification, DESeq2
  differential expression, publication-quality figures, and a bilingual (`tr`/`en`) self-contained
  HTML report. The QC → trim → align tool chain is chosen automatically by **read type** (see below);
  quantification is routed by `organism_type` (prokaryote: Bowtie2/minimap2 + featureCounts ·
  eukaryote: Salmon + tximport). Both dimensions converge on the same gene × sample count matrix.
- **Read types (short / long)**: Illumina reads run the short-read chain (FastQC → fastp → Bowtie2);
  ONT/PacBio reads run the long-read chain (NanoPlot → Pychopper+chopper → minimap2 →
  featureCounts `-L`). The read type is auto-detected in m01; from m05 onward (DESeq2 and all
  functional analysis) the pipeline is read-type-agnostic.
- **Functional analysis** (all optional, organism-agnostic, none produce a new failure gate):
  GO over-representation (ORA), KEGG pathway ORA, GSEA (fgsea), REVIGO-like semantic reduction,
  AMR + virulence overlay (abricate/CARD/VFDB), operon prediction + coordination, and STRING
  protein-interaction modules (Louvain community detection).

### Quality gates (why results are trustworthy)

A correct pipeline still produces a plausible-looking but **fake** result from bad input.
RNAForge enforces gates with a two-tier policy:

- **FAIL** → the result is **invalid**: the run stops, no biological output is produced (exit 1).
- **WARN** → the result is **suspect**: it is produced but stamped as such.

Thresholds are data (`profiles/{prokaryote,eukaryote,prokaryote_long}.yml`); an overridden threshold
is written into the report (no silent loosening). Long-read runs use the `prokaryote_long` profile,
whose thresholds are deliberately permissive and stamped (ONT quality ~Q10–15, not Q30): only a
catastrophic alignment failure (wrong reference) FAILs, while ONT-typical lower survival/assignment
WARN. Every run also writes a confidence card (`UNKNOWN`/`INVALID`/`SUSPECT`/`TRUSTWORTHY`).

## Pipeline modules

| Stage | Subcommand | What it does |
|---|---|---|
| m00 | `basecall` | ONT raw signal (FAST5/POD5) → FASTQ via dorado (GPU); optional, only when input is raw signal |
| m01 | `validate` | Config + metadata + FASTQ validation, platform + read-type detection, design gates |
| m02 | `qc` | short: FastQC · long: NanoPlot (diagnostic; never stops the run) |
| m03 | `trim` | short: fastp (gentle) · long: Pychopper+chopper (cDNA) / chopper (direct-RNA) |
| m04 | `quant` | Alignment (prok short: Bowtie2 · long: minimap2 · euk: Salmon) |
| m05 | `counts` | featureCounts (`-L` for long-read) → gene × sample count matrix |
| m06 | `de` | DESeq2 differential expression |
| m07 | `figures` | PCA, volcano, MA, heatmap, dispersion, … (PNG 300dpi + SVG) |
| m08 | `report` | Single self-contained bilingual HTML report (read-type-aware) |
| m09 | `enrich` | GO over-representation (ORA), hypergeometric + BH |
| m10 | `kegg` | KEGG pathway ORA |
| m11 | `gsea` | GSEA on the ranked gene list (fgsea) |
| m12 | `semantic` | REVIGO-like semantic reduction of GO terms (+ MDS map) |
| m13 | `amr` | AMR (CARD + AMRFinderPlus) + virulence (VFDB) gene overlay onto DE |
| m14 | `operon` | Operon prediction (intergenic distance) + DE coordination |
| m15 | `ppi` | STRING PPI subnetwork + Louvain community modules |
| m16 | `seqqc` | rRNA% (SortMeRNA) + strandedness (RSeQC) — WARN gates |
| m17 | `alignqc` | insert-size + coverage + read-distribution (samtools/RSeQC) |
| m18 | `multiqc` | aggregate MultiQC view across the run (run last) |

The downstream analyses (m09–m15) are organism- and read-type-agnostic and never invalidate a run —
the verdict carries over unchanged from the quality gates. The QC add-ons (m16–m18) are diagnostic.

## Install

One command creates all conda environments and installs the package editable:

```bash
bash install.sh
conda run -n rnaforge-core rnaforge doctor    # verify every required env exists
```

`install.sh` is idempotent (skips envs that already exist). The nine tool environments
it creates (`envs/*.yml`, pinned to exact versions for reproducibility): `rnaforge-core`
(orchestration + networkx/scipy), `rnaforge-qc` (FastQC, fastp), `rnaforge-quant-prok`
(Bowtie2, samtools, featureCounts), `rnaforge-quant-euk` (Salmon), `rnaforge-longread`
(minimap2, NanoPlot, Pychopper, chopper), `rnaforge-basecall` (pod5 + external dorado GPU
binary — raw-signal m00), `rnaforge-de` (R: DESeq2, ggplot2, fgsea, tximport), `rnaforge-amr`
(abricate CARD/VFDB), `rnaforge-seqqc` (SortMeRNA, RSeQC, MultiQC).

> `dorado` (ONT raw-signal basecalling, m00) is a GPU-only binary installed **outside** conda;
> `install.sh` does not fetch it. Install it separately and set `basecall.dorado_bin` only if you
> feed FAST5/POD5 input. FASTQ input skips m00 entirely.

## Usage

The whole pipeline in one command (stop-on-FAIL, resumable — re-run the same command after
a crash and completed stages are skipped):

```bash
rnaforge run --config config/config.yaml --metadata samples.tsv --run-id demo
# add optional stages before the report:
rnaforge run --config config/config.yaml --metadata samples.tsv --run-id demo \
             --include enrich,kegg,gsea,seqqc,alignqc,multiqc
# or run a slice of the core chain:
rnaforge run --config config/config.yaml --metadata samples.tsv --run-id demo --from trim --to counts
```

Or drive each stage by hand (same `--run-id` throughout):

```bash
# core chain
rnaforge validate --config config/config.yaml --metadata samples.tsv --run-id demo
rnaforge qc       --config config/config.yaml --metadata samples.tsv --run-id demo
rnaforge trim     --config config/config.yaml --metadata samples.tsv --run-id demo
rnaforge quant    --config config/config.yaml --metadata samples.tsv --run-id demo
rnaforge counts   --config config/config.yaml --metadata samples.tsv --run-id demo
rnaforge de       --config config/config.yaml --metadata samples.tsv --run-id demo
rnaforge figures  --config config/config.yaml --metadata samples.tsv --run-id demo

# optional QC / diagnostics (require m04; produce diagnostic figures/tables, never FAIL)
rnaforge seqqc    --config config/config.yaml --metadata samples.tsv --run-id demo  # rRNA% + strandedness (m16)
rnaforge alignqc  --config config/config.yaml --metadata samples.tsv --run-id demo  # insert-size + coverage + read-distribution (m17)
rnaforge multiqc  --config config/config.yaml --metadata samples.tsv --run-id demo  # aggregate MultiQC view (m18, run last)

# optional functional analyses (any subset; each needs its reference data — see below)
rnaforge enrich   --config config/config.yaml --metadata samples.tsv --run-id demo
rnaforge kegg     --config config/config.yaml --metadata samples.tsv --run-id demo
rnaforge gsea     --config config/config.yaml --metadata samples.tsv --run-id demo
rnaforge semantic --config config/config.yaml --metadata samples.tsv --run-id demo
rnaforge amr      --config config/config.yaml --metadata samples.tsv --run-id demo
rnaforge operon   --config config/config.yaml --metadata samples.tsv --run-id demo
rnaforge ppi      --config config/config.yaml --metadata samples.tsv --run-id demo

# assemble the report last — it embeds whatever analyses were run
rnaforge report   --config config/config.yaml --metadata samples.tsv --run-id demo
```

> Note: `python -m rnaforge.cli` does not work (no main-guard); use the installed `rnaforge` entry point.

### Metadata format (TSV)

| Column | Required | Description |
|---|---|---|
| `sample_id` | yes | Unique sample identifier |
| `condition` | yes | Experimental group; needs ≥2 levels and ≥2 replicates each |
| `fastq_1` | yes | Path to R1 (or single-end reads) |
| `fastq_2` | no | Path to R2 for paired-end |
| `subject` | no | Paired/subject id; detected and, if it looks paired, must be handled deliberately |
| `batch` | no | Batch/covariate; required if the design formula uses `batch` |

## Reference data (one-time prep, git-ignored)

The functional analyses read local reference files under `references/` (never committed). A
parameterized script fetches them once for your organism (each block is skipped if its argument is
omitted; every download gets a `.sha256` stamp for reproducibility):

```bash
# example for E. coli K-12 (KEGG code eco, STRING taxid 511145)
bash prepare_references.sh \
     --kegg-org eco \
     --string-taxid 511145 \
     --goa-url https://ftp.ebi.ac.uk/pub/databases/GO/goa/proteomes/18.E_coli_MG1655.goa
```

Blocked-source note: QuickGO downloads are blocked on some networks; the script therefore takes the
GO annotation (GAF) from the **EBI-GOA FTP proteome** file via `--goa-url` (the documented fallback),
not QuickGO. AMR/virulence (m13) uses abricate's bundled CARD/VFDB databases (no separate download).

## Key design decisions

- **`organism_type` is required and has no default** (`prokaryote` | `eukaryote`). It routes only
  quantification (m04/m05); both paths converge on the same gene × sample count matrix, so every
  downstream step (m06–m15) is organism-agnostic.
- **Two routing dimensions: `organism_type` × read type.** The read type (short/long) is
  auto-detected from the FASTQ in m01 (Illumina → short; ONT/PacBio → long) and drives m02–m05;
  `organism_type` drives m04/m05. Unidentifiable platforms are still refused with a clear error
  (never silently processed through the wrong route).
- **`library.chemistry` is required for ONT long reads** (`cdna` | `direct_rna`) — it cannot be
  detected from the FASTQ and selects the m03 long-read preprocessing (cDNA → Pychopper+chopper;
  direct-RNA → chopper only). PacBio HiFi does not need it.
- **Raw signal (FAST5/POD5) is supported via m00 `basecall`** — when a sample's `fastq_1` points to a
  POD5/FAST5 file or directory, `rnaforge basecall` runs dorado (GPU, `hac` model) to produce FASTQ,
  then the rest of the pipeline runs unchanged. FASTQ input skips m00 entirely. dorado is a GPU-only
  ONT binary installed separately (set `basecall.dorado_bin`); a GPU is required.
- **Trimming is deliberately gentle.** Aggressive quality trimming distorts expression estimates
  ([Williams et al. 2016](https://doi.org/10.1186/s12859-016-0956-2)); a minimum-length filter
  prevents the distortion.
- **No fabricated results.** Ambiguous annotation joins (by gene symbol) are dropped, not guessed;
  predicted structures (operons, STRING interactions) are stamped as predictions in the report.

## Development

```bash
conda run -n rnaforge-core --cwd "$(pwd)" python -m pytest -q
```

Run from the repository root (the test suite imports `tests.conftest`).

## Privacy

Customer data is never committed. `runs/`, `raw/` and `references/` are git-ignored.
