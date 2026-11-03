# randPedGenPCA benchmarks

## Start here

Open **[the English results summary](reports/eddie_comparison/summary_en.html)**
in a browser, or **[the comparison workbook](reports/eddie_comparison/comparison_tables_en.xlsx)**
in Excel. These are local files; no website or server is required.

## Current workflows

| Task | Entry point | Instructions |
|---|---|---|
| Original R simulation/scalability benchmark | `R/run_scalability_benchmark.R` | [R benchmark guide](R/README_scalability_benchmark.md) |
| Large-population experiment used on Eddie | `validation/run_large_population.py` | [Large-population guide](docs/large_population_exercise.md) |
| Exact-reference and fixed-pedigree genotyping experiments | `validation/run_followups.py` | [Follow-up guide](docs/followup_experiments.md) |
| Profiling with the requested Advanced R first-edition method | `validation/profile_lineprof.R` and `validation/lineprof_job.sh` | [lineprof guide](docs/profiling_lineprof.md) |

The scientific computations are in R. The Python launchers used for completed
experiments are retained to reproduce those experiments. The current lineprof
profiling job runs directly in R. Submit compute jobs through `qsub` on Eddie;
do not run large PCA calculations on a login node.

## Folder map

- `R/`: original R benchmark, its dependencies and population generators.
- `validation/`: current experiment, profiling and report scripts. See its [file index](validation/README.md).
- `docs/`: guides for the active workflows.
- `tests/`: numerical smoke checks for the benchmark and follow-up datasets.
- `eddie_results/`: unmodified CSVs/configurations downloaded from Eddie.
- `reports/eddie_comparison/`: English HTML summary and Excel workbook; technical tables in `data/`, figures in `figures/`.
- `archive/cleanup_20261102/`: previous local outputs, older validation workflows, alternative profiler and Spanish reports. Files were archived, not deleted.

New experiment results are created in `output/` when an experiment is launched.
Previous local results were moved to `archive/cleanup_20261102/output/`.
The cleanup did not change any files on Eddie or the neighbouring scientific packages.

## Update the English summary locally

From this repository root, after downloading the latest Eddie CSVs:

```sh
python3 validation/summarise_eddie.py
Rscript validation/export_eddie_tables.R
python3 validation/render_eddie_report.py
```

The report defaults to English. Spanish/bilingual text can be generated with
`summarise_eddie.py --language es` or `--language both`; the R exporter accepts
`reports/eddie_comparison es`, and the HTML renderer accepts `--language es`.

For the next profiling step, read [the lineprof guide](docs/profiling_lineprof.md)
before installing or submitting anything. It explains what is executed and what
the time/memory measurements mean.

The [Hutch++ trace algorithm guide](docs/hutch_trace_algorithm.md) explains the
estimator, its vector budget and its connection to PCA.

The [algorithm alignment note](docs/algorithm_alignment.md) maps the specified
randomised-PCA steps to the implementation and documents the benchmark-only
correction to eigenvalue-scaled score comparisons.

For the reference-version, correctness and before/after measurement workflow,
read [the optimization protocol](docs/optimization_protocol.md). The initial
[change log](reports/optimization/change_log.csv) records historical reference
measurements; baseline commits and any optimization results remain to be recorded.
