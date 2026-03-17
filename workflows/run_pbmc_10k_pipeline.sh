#!/usr/bin/env bash
set -euo pipefail

CONFIG_PATH="${1:-config/pbmc_10k_multi.yaml}"
SCRIPTS=(
  "01_prepare_object.R"
  "02_cluster_scrna.R"
  "03_annotate_celltypes.R"
  "04_integrate_tcr.R"
  "05_tcell_clonotype_analysis.R"
  "06_build_report.R"
)

for script in "${SCRIPTS[@]}"; do
  script_name="${script%.R}"
  log_path="logs/pbmc_10k_multi/${script_name}.log"
  echo "[run] ${script}"
  Rscript "src/pbmc_10k/${script}" "${CONFIG_PATH}" >"${log_path}" 2>&1
  echo "[ok] ${script} -> ${log_path}"
done
