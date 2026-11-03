# Current scripts

## Large benchmark

- `run_large_population.py`: existing experiment coordinator; the PCA is calculated in R.
- `large_population.R`: shared R worker for simulation, PCA and implicit diagnostics.

## Follow-up experiments

- `run_followups.py`: coordinates small exact-reference validation or nested genotyping experiments.
- `followup_worker.R`: exact reference, exact comparisons and fixed-pedigree subset generation.
- `followup_job.sh`: Eddie batch wrapper for those two experiments.

## Requested profiling method: lineprof

- `profile_lineprof.R`: profiles a saved dataset directly in R.
- `lineprof_job.sh`: submits the two saved fixed-pedigree profiling cases to Eddie.
- `view_lineprof.R`: opens the profile locally with `shine()` and remaps source paths.

The earlier `profile_pca.R`/`profile_job.sh` alternative and its documentation were
moved to `archive/cleanup_20261102/`. Use lineprof for the currently requested method.

## Results and reports

- `summarise_eddie.py`: merges and audits downloaded CSVs into comparison tables.
- `export_eddie_tables.R`: exports Excel and figures; defaults to English.
- `render_eddie_report.py`: renders the offline HTML summary; defaults to English.

Older validation/report pipelines were archived together, with a supporting copy
of their R helpers. They are no longer mixed with these active scripts.
