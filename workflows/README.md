# Workflows

This folder stores the runnable workflow entry point.

## Files

- `run_pbmc_10k_pipeline.sh`: executes the PBMC 10k workflow in order and writes one log per step under `logs/pbmc_10k_multi/`.

## Use

- Run `./workflows/run_pbmc_10k_pipeline.sh` for a full end-to-end execution.
- For partial reruns, execute the individual scripts in `src/pbmc_10k/` instead.
