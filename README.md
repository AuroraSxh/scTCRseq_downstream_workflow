# scTCR-seq Downstream Workflow

Code-only downstream workflow for PBMC-style single-cell `scRNA-seq + scTCR-seq` analysis on 10x multi-omic inputs.

This repository is structured for public sharing of pipeline logic, configuration, and expected results layout without bundling source data, large intermediate objects, or generated outputs.

## Overview

Current included workflow:

- `pbmc_10k_multi`

Main analysis stages:

- gene expression object preparation
- scRNA-seq clustering
- cell-type annotation
- TCR metadata integration
- T-cell clonotype analysis
- report assembly

## Repository Layout

```text
.
├── README.md
├── LICENSE
├── .gitignore
├── TEST_DATASETS.md
├── config/
├── src/
├── workflows/
└── expected_results/
```

## Workflow Entry Point

Run from the repository root:

```bash
./workflows/run_pbmc_10k_pipeline.sh config/pbmc_10k_multi.yaml
```

Main code modules:

- `src/pbmc_10k/00_utils.R`
- `src/pbmc_10k/01_prepare_object.R`
- `src/pbmc_10k/02_cluster_scrna.R`
- `src/pbmc_10k/03_annotate_celltypes.R`
- `src/pbmc_10k/04_integrate_tcr.R`
- `src/pbmc_10k/05_tcell_clonotype_analysis.R`
- `src/pbmc_10k/06_build_report.R`

## Expected Input

The workflow expects upstream-prepared downstream inputs rather than FASTQ files:

- a gene expression H5 matrix
- filtered TCR contig annotations
- clonotype annotations

Paths and parameters are configured in `config/pbmc_10k_multi.yaml`.

Additional dataset guidance is documented in `TEST_DATASETS.md`.

## Repository Policy

Included:

- code
- workflow entry points
- configuration
- dataset guidance
- expected result structure

Excluded:

- raw sequencing data
- public test matrices and contig tables
- generated logs
- generated results objects
- local package environments

## Expected Results

This repository ships only the expected output skeleton:

- `expected_results/pbmc_10k_multi/`

That structure mirrors the staged downstream deliverables for QC, clustering, annotation, TCR integration, clonotype analysis, reports, checkpoint objects, and per-step logs.
