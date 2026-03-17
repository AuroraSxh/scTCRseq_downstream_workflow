# PBMC 10k Scripts

This folder contains the modular R scripts for the PBMC 10k scRNA + scTCR downstream workflow.

## Steps

- `00_utils.R`: shared helpers for config loading, directory creation, plotting, metadata writing, and clonotype binning.
- `01_prepare_object.R`: reads the 10x H5 matrix, computes QC metrics, applies filtering, and saves the filtered Seurat object.
- `02_cluster_scrna.R`: performs normalization, HVG selection, PCA, neighbors, clustering, UMAP, and cluster-average marker summaries.
- `03_annotate_celltypes.R`: assigns broad cell types from cluster marker patterns and exports annotation tables and plots.
- `04_integrate_tcr.R`: aggregates productive TRA/TRB calls, merges clonotype metadata into Seurat metadata, and saves TCR summary outputs.
- `05_tcell_clonotype_analysis.R`: subsets T cells, reruns embedding, and summarizes clonotype structure within T-cell states.
- `06_build_report.R`: writes the final markdown summary and cell-type count table.

## Rerun strategy

- If QC changes, rerun `01` through `06`.
- If clustering or resolution changes, rerun `02` through `06`.
- If only annotation rules change, rerun `03` through `06`.
- If only report text changes, rerun `06`.
