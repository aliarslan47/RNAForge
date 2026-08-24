# RNAForge

Reproducible, modular **bulk RNA-seq** pipeline — from raw FASTQ to a single, self-contained HTML report, with a functional-analysis layer on top of differential expression.

[![Pipeline DAG](https://img.shields.io/badge/pipeline-DAG-0d6b8f)](https://aliarslan47.github.io/RNAForge/pipeline_architecture.html)
[![organism](https://img.shields.io/badge/organism-prokaryote%20%C2%B7%20eukaryote%20%C2%B7%20metatranscriptome-2f8f5b)](https://aliarslan47.github.io/RNAForge/pipeline_architecture.html)
[![reads](https://img.shields.io/badge/reads-short%20%C2%B7%20long-c07211)](https://aliarslan47.github.io/RNAForge/pipeline_architecture.html)

[Türkçe](README.tr.md) · **English**

## What is it?

RNAForge is the bulk RNA-seq member of the Forge family — same architecture as BacForge (bacteria) and VirusForge (virus/phage), but a separate, isolated installation. It takes raw reads to biology in a single command and ends with a bilingual (`tr`/`en`), self-contained HTML report.

## What it does

A staged pipeline — `validate → qc → trim → quant → counts → de → figures → report` — with an optional functional-analysis layer (GO / KEGG / GSEA / semantic / AMR / operon / PPI).

Two routing dimensions converge on the same gene × sample count matrix, so every step from DE onward is agnostic:

- **Organism** (`organism_type`): prokaryote (Bowtie2/minimap2 + featureCounts) · eukaryote (Salmon + tximport) · metatranscriptome (rRNA depletion + Kraken2/Bracken + gene-catalog).
- **Read type** (auto-detected): short (Illumina: FastQC → fastp → Bowtie2) · long (ONT/PacBio: NanoPlot → Pychopper+chopper → minimap2).

Trustworthy by design: two-tier quality gates (**FAIL** stops the run, **WARN** stamps a suspect result), thresholds are data (`profiles/*.yml`), and every run writes a confidence card. No fabricated results.

Interactive bilingual node-graph: **[rendered diagram](https://aliarslan47.github.io/RNAForge/pipeline_architecture.html)**.

## Installation

```bash
bash install.sh
conda run -n rnaforge-core rnaforge doctor   # verify every required env exists
```

Idempotent; creates nine version-pinned conda environments (`envs/*.yml`). `dorado` (ONT raw-signal basecalling, m00) is a GPU-only binary installed separately — needed only for FAST5/POD5 input.

## Usage

```bash
# whole pipeline (stop-on-FAIL, resumable)
rnaforge run --config config/config.yaml --metadata samples.tsv --run-id demo

# add optional stages before the report
rnaforge run ... --include enrich,kegg,gsea,seqqc,alignqc,multiqc

# or run a slice of the core chain
rnaforge run ... --from trim --to counts
```

Each stage can also be driven by hand with the same `--run-id`. Use the installed `rnaforge` entry point (not `python -m`).

## Modules

| Code | Subcommand | What it does |
|---|---|---|
| m00 | `basecall` | ONT raw signal (FAST5/POD5) → FASTQ via dorado (GPU); optional |
| m01 | `validate` | Config/metadata/FASTQ validation, platform + read-type detection |
| m02 | `qc` | short: FastQC · long: NanoPlot |
| m03 | `trim` | short: fastp (gentle) · long: Pychopper+chopper / chopper |
| m04 | `quant` | Alignment (prok: Bowtie2/minimap2 · euk: Salmon) |
| m05 | `counts` | featureCounts → gene × sample count matrix |
| m06 | `de` | DESeq2 differential expression |
| m07 | `figures` | PCA, volcano, MA, heatmap, dispersion (PNG 300dpi + SVG) |
| m08 | `report` | Single self-contained bilingual HTML report |
| m09 | `enrich` | GO over-representation (ORA) |
| m10 | `kegg` | KEGG pathway ORA |
| m11 | `gsea` | GSEA on the ranked gene list (fgsea) |
| m12 | `semantic` | REVIGO-like semantic reduction of GO terms |
| m13 | `amr` | AMR (CARD + AMRFinderPlus) + virulence (VFDB) overlay |
| m14 | `operon` | Operon prediction + DE coordination |
| m15 | `ppi` | STRING PPI subnetwork + Louvain modules |
| m16 | `seqqc` | rRNA% (SortMeRNA) + strandedness (RSeQC) |
| m17 | `alignqc` | insert-size + coverage + read-distribution |
| m18 | `multiqc` | aggregate MultiQC view (run last) |

Metatranscriptome runs auto-insert `rrna-deplete` and `taxonomy` between `trim` and `quant`. Downstream analyses (m09–m18) are organism- and read-type-agnostic and never invalidate a run. Full design, metadata format and reference-data prep live in `PLAN.md` and `docs/`.

---

Forge family: **RNAForge** (bulk RNA-seq) · [BacForge](https://github.com/aliarslan47/BacForge) (bacteria) · [VirusForge](https://github.com/aliarslan47/VirusForge) (virus/phage) · [MicrobiomeForge](https://github.com/aliarslan47/MicrobiomeForge) (microbiome) · [Vaxforge](https://github.com/aliarslan47/Vaxforge) (reverse vaccinology) · [ImmForge](https://github.com/aliarslan47/ImmForge) (immune simulation) · [PipelineForge](https://github.com/aliarslan47/PipelineForge) (DAG generator). Customer data is never committed (`runs/`, `raw/`, `references/` are git-ignored).
