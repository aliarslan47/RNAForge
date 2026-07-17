# RNAForge

Reproducible, modular Bulk RNA-seq analysis pipeline.

Turkish version: [README.tr.md](README.tr.md) · Reference document: [PLAN.md](PLAN.md) (v1.2)

## Status

Early development. Currently implemented: input validation and platform detection (`m01`).

## Install

```bash
conda env create -f envs/rnaforge-core.yml
conda activate rnaforge-core
pip install -e .
```

## Usage

```bash
rnaforge validate --config config/config.yaml --metadata samples.tsv --run-id demo
```

### Metadata format (TSV)

| Column | Required | Description |
|---|---|---|
| `sample_id` | yes | Unique sample identifier |
| `condition` | yes | Experimental group; needs ≥2 levels and ≥2 replicates each |
| `fastq_1` | yes | Path to R1 (or single-end reads) |
| `fastq_2` | no | Path to R2 for paired-end |
| `batch` | no | Batch/covariate; required if the design formula uses `batch` |

## Key design decisions

- **`organism_type` is required and has no default** (`prokaryote` | `eukaryote`). It routes
  quantification: prokaryote uses genome alignment + featureCounts, eukaryote uses
  Salmon + tximport. Both converge on the same gene × sample count matrix.
- **Illumina only (MVP).** ONT/PacBio inputs are detected and refused with a clear error
  rather than silently processed through the wrong route.
- **Trimming is deliberately gentle.** Aggressive quality trimming distorts expression
  estimates ([Williams et al. 2016](https://doi.org/10.1186/s12859-016-0956-2)); a minimum
  length filter is what prevents the distortion.

## Development

```bash
conda run -n rnaforge-core python -m pytest -v
```

## Privacy

Customer data is never committed. `runs/`, `raw/` and `references/` are git-ignored.
