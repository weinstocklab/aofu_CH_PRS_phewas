# Workflow portability notes

This directory contains the Snakemake workflow and helper scripts used to build PRS inputs and run targeted/survival PheWAS analyses.

## Environment assumptions

- **AoU/BigQuery dependency is expected** for data extraction scripts (`query_aou_*.py`).
- `workflow/scripts/query_aou_survival_base.py` expects:
  - `WORKSPACE_CDR` (unless `--cdr` is passed)
  - a billing project via `GOOGLE_PROJECT` or `GOOGLE_CLOUD_PROJECT` (unless `--billing-project` is passed)
- You also need authenticated Google Cloud credentials with access to the referenced AoU datasets.

## Non-portable dependencies

### 1) Runtime GitHub fetches for PheWAS resources

`workflow/scripts/build_phecode_resources.R` downloads `.rda` files at runtime from:

- `https://raw.githubusercontent.com/PheWAS/PheWAS/master/data/phecode_map.rda`
- `https://raw.githubusercontent.com/PheWAS/PheWAS/master/data/phecode_rollup_map.rda`
- `https://raw.githubusercontent.com/PheWAS/PheWAS/master/data/pheinfo.rda`

This depends on network access and an upstream `master` branch that may change.

**Recommendation:** pin URLs to a specific commit SHA, or vendor/cache these files in your own storage for reproducible runs.

### 2) Cross-script Python import assumptions

`workflow/scripts/run_somatic_burden_phewas.py` imports helper functions from sibling scripts (for example `run_targeted_phewas` and `run_survival_phewas`). This assumes the scripts directory is importable at runtime.

**Recommendation:** run from repo root (`snakemake` from this repository) or set `PYTHONPATH` so sibling-script imports resolve consistently.

## Input schema expectations

The scripts rely on fixed column names. Common required fields include:

- **Sample IDs**
  - `IID` (primary join key across most TSVs)
  - `ID` (used in age-at-sequencing input, then renamed to `IID`)
- **Score inputs**
  - exactly one standardized score column ending in `*_Z` for `run_phewas.py`
  - default PRS column `PRSFNN_out_final_Z` for `run_survival_phewas.py` and `run_prs_spline_phewas.py` (overridable by CLI)
- **Covariates**
  - `sex` (or configured `--sex-col`)
  - age columns from config (`--age-cols`)
  - principal component columns from config (`--pc-cols`)
- **EHR mapping/evidence tables**
  - `vocabulary_id`, `concept_code`, `code_date`, `phenotype_id`, `map_status`
- **Survival base output**
  - `IID`, `birth_date`, `observation_start_date`, `observation_end_date`, `death_date`

If these names drift, scripts fail with missing-column errors or empty merges.

## Quick portability checklist

- Confirm AoU/BigQuery access and auth are available where the workflow runs.
- Set `WORKSPACE_CDR` and `GOOGLE_PROJECT`/`GOOGLE_CLOUD_PROJECT` (or pass equivalent CLI flags).
- Decide whether to pin or vendor/cache external PheWAS `.rda` resources.
- Run workflow from repo root (or set `PYTHONPATH`) so cross-script imports resolve.
- Validate upstream TSV schemas include the required ID, score, and covariate columns.
