# Expected Results: `pbmc_10k_multi`

This directory describes the expected result structure for the downstream PBMC 10k scRNA + scTCR workflow.

```text
expected_results/pbmc_10k_multi/
├── 01_qc/
├── 02_clustering/
├── 03_annotation/
├── 04_tcr/
├── 05_tcell/
├── 06_report/
├── objects/
└── logs/
```

## What Each Directory Represents

- `01_qc/`
  - QC summaries before and after filtering
  - cell-level metadata snapshots
  - QC plots

- `02_clustering/`
  - normalized embeddings
  - cluster labels
  - marker summaries

- `03_annotation/`
  - broad cell-type annotation tables
  - annotation support plots

- `04_tcr/`
  - merged clonotype metadata
  - TRA/TRB summaries
  - clonotype frequency tables

- `05_tcell/`
  - T-cell-only embeddings
  - clonotype expansion summaries in T-cell states
  - state-aware clonotype plots

- `06_report/`
  - final markdown summary
  - compact export tables for interpretation

- `objects/`
  - checkpoint Seurat objects for staged reruns

- `logs/`
  - one log per pipeline step

## Notes

- The repository includes only the structure, not real outputs.
- Real results depend on the chosen dataset and parameters in `config/pbmc_10k_multi.yaml`.
