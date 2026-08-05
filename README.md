# RNAForge

Reproducible, modular **bulk RNA-seq** analysis pipeline — from raw FASTQ to a single,
self-contained HTML report, with a full functional-analysis layer on top of differential expression.

Turkish version: [README.tr.md](README.tr.md) · Reference document: [PLAN.md](PLAN.md) (v1.3)

## What it does

A staged pipeline that takes raw reads to biology:

```
validate → qc → trim → quant → counts → de → figures → report
                                          └→ enrich · kegg · gsea · semantic · amr · operon · ppi
```

- **Core**: input/design validation, FastQC, gentle fastp trimming, alignment/quantification
  (prokaryote: Bowtie2 + featureCounts · eukaryote: Salmon + tximport), DESeq2 differential
  expression, publication-quality figures, and a bilingual (`tr`/`en`) self-contained HTML report.
- **Functional analysis** (all optional, organism-agnostic, none produce a new failure gate):
  GO over-representation (ORA), KEGG pathway ORA, GSEA (fgsea), REVIGO-like semantic reduction,
  AMR + virulence overlay (abricate/CARD/VFDB), operon prediction + coordination, and STRING
  protein-interaction modules (Louvain community detection).

### Quality gates (why results are trustworthy)

A correct pipeline still produces a plausible-looking but **fake** result from bad input.
RNAForge enforces gates with a two-tier policy:

- **FAIL** → the result is **invalid**: the run stops, no biological output is produced (exit 1).
- **WARN** → the result is **suspect**: it is produced but stamped as such.

Thresholds are data (`profiles/{prokaryote,eukaryote}.yml`); an overridden threshold is written into
the report (no silent loosening). Every run also writes a confidence card
(`UNKNOWN`/`INVALID`/`SUSPECT`/`TRUSTWORTHY`).

## Pipeline modules

| Stage | Subcommand | What it does |
|---|---|---|
| m01 | `validate` | Config + metadata + FASTQ validation, platform detection, design gates |
| m02 | `qc` | FastQC (diagnostic; never stops the run) |
| m03 | `trim` | fastp — gentle trimming (adapter + min-length; aggressive quality off) |
| m04 | `quant` | Alignment/quantification (prok: Bowtie2 · euk: Salmon) |
| m05 | `counts` | featureCounts → gene × sample count matrix |
| m06 | `de` | DESeq2 differential expression |
| m07 | `figures` | PCA, volcano, MA, heatmap, dispersion, … (PNG 300dpi + SVG) |
| m08 | `report` | Single self-contained bilingual HTML report |
| m09 | `enrich` | GO over-representation (ORA), hypergeometric + BH |
| m10 | `kegg` | KEGG pathway ORA |
| m11 | `gsea` | GSEA on the ranked gene list (fgsea) |
| m12 | `semantic` | REVIGO-like semantic reduction of GO terms (+ MDS map) |
| m13 | `amr` | AMR (CARD) + virulence (VFDB) gene overlay onto DE (abricate) |
| m14 | `operon` | Operon prediction (intergenic distance) + DE coordination |
| m15 | `ppi` | STRING PPI subnetwork + Louvain community modules |

The downstream analyses (m09–m15) are organism-agnostic and never invalidate a run — the verdict
carries over unchanged from the quality gates.

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
conda env create -f envs/rnaforge-de.yml           # R: DESeq2, ggplot2, fgsea
conda env create -f envs/rnaforge-amr.yml          # abricate (CARD/VFDB)
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
- **Illumina only (MVP).** ONT/PacBio inputs are detected and refused with a clear error rather than
  silently processed through the wrong route.
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
