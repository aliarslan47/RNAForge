# RNAForge

Reproducible, modular **bulk RNA-seq** analysis pipeline — from raw FASTQ to a single,
self-contained HTML report, with a full functional-analysis layer on top of differential expression.

Turkish version: [README.tr.md](README.tr.md) · Reference document: [PLAN.md](PLAN.md) (v1.4)

[![Pipeline DAG](https://img.shields.io/badge/pipeline-DAG-0d6b8f)](https://claude.ai/code/artifact/7d033f10-ade2-4cbe-801d-b468a06b0c5a)
[![organism](https://img.shields.io/badge/organism-prokaryote%20%C2%B7%20eukaryote-2f8f5b)](https://claude.ai/code/artifact/7d033f10-ade2-4cbe-801d-b468a06b0c5a)
[![reads](https://img.shields.io/badge/reads-short%20%C2%B7%20long-c07211)](https://claude.ai/code/artifact/7d033f10-ade2-4cbe-801d-b468a06b0c5a)

## What it does

A staged pipeline that takes raw reads to biology:

```
validate → qc → trim → quant → counts → de → figures → report
                                          └→ enrich · kegg · gsea · semantic · amr · operon · ppi
```

The full pipeline as an interactive, bilingual node-graph (organism × read-type branching, `m00`–`m18`) —
[**rendered diagram**](https://claude.ai/code/artifact/7d033f10-ade2-4cbe-801d-b468a06b0c5a) · source: `docs/pipeline_architecture.html`.

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

```bash
conda env create -f envs/rnaforge-core.yml     # orchestration (Python) + networkx/numpy/scipy
conda activate rnaforge-core
pip install -e .
```

Tool environments (created once, referenced by the modules):

```bash
conda env create -f envs/rnaforge-qc.yml           # FastQC, fastp
conda env create -f envs/rnaforge-quant-prok.yml   # Bowtie2, samtools, featureCounts
conda env create -f envs/rnaforge-quant-euk.yml    # Salmon
conda env create -f envs/rnaforge-longread.yml     # minimap2, NanoPlot, Pychopper, chopper, samtools
conda env create -f envs/rnaforge-basecall.yml     # pod5, samtools (+ dorado binary, GPU) — raw signal m00
conda env create -f envs/rnaforge-de.yml           # R: DESeq2, ggplot2, fgsea
conda env create -f envs/rnaforge-amr.yml          # abricate (CARD/VFDB)
conda env create -f envs/rnaforge-seqqc.yml        # SortMeRNA, RSeQC, MultiQC (m16/m18)
```

## Usage

```bash
# core chain (same --run-id throughout)
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

The functional analyses read local reference files under `references/` (never committed). Download
once for your organism (examples for *E. coli* K-12):

```bash
# GO ontology (m09/m12) + organism GO annotation (EBI-GOA)
curl -L -o references/go/go-basic.obo http://purl.obolibrary.org/obo/go/go-basic.obo
curl -L https://ftp.ebi.ac.uk/pub/databases/GO/goa/proteomes/18.E_coli_MG1655.goa \
     -o references/ecoli_bw25113/ecoli.gaf

# KEGG (m10) — per-organism REST files
curl -s https://rest.kegg.jp/link/pathway/eco > references/kegg/eco/pathway_links.tsv
curl -s https://rest.kegg.jp/list/pathway/eco > references/kegg/eco/pathway_names.tsv
curl -s https://rest.kegg.jp/list/eco        > references/kegg/eco/gene_list.tsv

# STRING (m15) — per-taxon network
curl -s https://stringdb-downloads.org/download/protein.info.v12.0/511145.protein.info.v12.0.txt.gz \
     -o references/string/511145/protein.info.txt.gz
curl -s https://stringdb-downloads.org/download/protein.links.v12.0/511145.protein.links.v12.0.txt.gz \
     -o references/string/511145/protein.links.txt.gz
```

AMR/virulence (m13) uses abricate's bundled CARD/VFDB databases (no separate download).

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
