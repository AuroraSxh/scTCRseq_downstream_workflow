# scTCRseq Downstream Workflow

This repository contains a downstream single-cell `scRNA-seq + scTCR-seq` workflow for PBMC 10k style 10x multi-omic inputs.

It intentionally excludes:
- raw sequencing data
- public test matrices and contig tables
- generated logs
- generated results objects

It includes only:
- code
- workflow entry points
- configuration
- dataset guidance
- expected result structure

## Repository Layout

```text
.
├── README.md
├── TEST_DATASETS.md
├── config/
├── src/
├── workflows/
└── expected_results/
```

## Included Workflow

Current workflow:
- `pbmc_10k_multi`

Main runnable entry point:
- `workflows/run_pbmc_10k_pipeline.sh`

Main code modules:
- `src/pbmc_10k/00_utils.R`
- `src/pbmc_10k/01_prepare_object.R`
- `src/pbmc_10k/02_cluster_scrna.R`
- `src/pbmc_10k/03_annotate_celltypes.R`
- `src/pbmc_10k/04_integrate_tcr.R`
- `src/pbmc_10k/05_tcell_clonotype_analysis.R`
- `src/pbmc_10k/06_build_report.R`

## Expected Input

The workflow expects upstream-prepared 10x-style downstream inputs rather than FASTQ:

- a gene expression H5 matrix
- filtered TCR contig annotations
- clonotype annotations

The exact paths and parameters are configured in:
- `config/pbmc_10k_multi.yaml`

More dataset guidance is in:
- `TEST_DATASETS.md`

## How To Run

From repository root:

```bash
./workflows/run_pbmc_10k_pipeline.sh config/pbmc_10k_multi.yaml
```

Each step writes one log file and produces staged outputs under the configured results directory.

## Expected Results

This repository does not include real outputs. Instead it includes the expected output skeleton in:

- `expected_results/pbmc_10k_multi/`

That structure mirrors the staged downstream deliverables:
- QC
- clustering
- annotation
- TCR integration
- T-cell clonotype analysis
- report tables
- checkpoint objects
- per-step logs

See:
- `expected_results/pbmc_10k_multi/README.md`
