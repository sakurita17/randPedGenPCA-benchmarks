# Advanced R first-edition profiling with lineprof

The requested chapter is https://adv-r.had.co.nz/Profiling.html. Its original
site timed out during this session; the lineprof implementation and documentation
were checked directly in https://github.com/hadley/lineprof. The package is
deprecated in favour of profvis; Advanced R's current edition uses profvis at
https://adv-r.hadley.nz/perf-measure.html. This workflow deliberately uses
lineprof to follow the user's request, without silently switching packages.

lineprof uses Rprof for statistical time/memory sampling and provides `shine()`
for interactive inspection. It reports `time`, `alloc`, `release` and `dups`.
Allocation/release are sampled estimates of R-managed memory changes, not peak
process RSS. Duplication counts are the sampled duplicate counters collected by R.
The script uses a 20 ms interval and `torture=FALSE`; forced garbage collection
would heavily perturb this long computation. It does not install packages or
modify the scientific package automatically.

## Install in Eddie (user executes)

```sh
module load R/4.4
export R_LIBS_USER="$HOME/HighlanderLab/jortiz/pca-benchmark/R-library/4.4"
Rscript -e 'if (!requireNamespace("remotes", quietly=TRUE)) install.packages("remotes", repos="https://cloud.r-project.org", lib=Sys.getenv("R_LIBS_USER"))'
Rscript -e 'remotes::install_github("hadley/lineprof", lib=Sys.getenv("R_LIBS_USER"), upgrade="never")'
Rscript -e 'library(lineprof); print(packageVersion("lineprof"))'
```

This legacy package needs compilation. Successful installation on the local Mac
does not guarantee installation or compatibility on Eddie; if it fails, keep the
error log and diagnose it before launching a long run. No automatic fallback to
profvis is configured.

## Copy from the Mac and submit in Eddie

Mac:

```sh
cd ~/Documents/GitHub/randPedGenPCA-benchmarks
rsync -av validation/profile_lineprof.R validation/lineprof_job.sh \
  validation/view_lineprof.R \
  s2745345@eddie.ecdf.ed.ac.uk:HighlanderLab/jortiz/pca-benchmark/benchmarks/validation/
```

Eddie:

```sh
cd ~/HighlanderLab/jortiz/pca-benchmark/benchmarks
ls -lh output/followup_genotyped_20261031_182624_492600/data_g{500,4000}.rds
bash -n validation/lineprof_job.sh
qsub -m ea -M s2745345@ed.ac.uk validation/lineprof_job.sh
```

The job runs directly in R on the two saved nested datasets, H q=4,l=20, trace=50,
seed=123. It is separate from the earlier Rprof-plus-counters job. Use this
lineprof job for the requested method; do not submit both unnecessarily. Results
go to `output/lineprof_JOBID/g500` and `g4000` with source snapshots and session
information. Each directory contains the profile object, flattened CSV, summary,
instrumented runtime and PCA fit. The instrumented runtime includes profile parsing
and must not replace the original benchmark time.

lineprof catches errors internally. The driver records successful completion and
propagates a PCA error after saving the available profile; a profile file alone
does not prove the PCA finished.

## View in R on the Mac after downloading

Download the entire profile directory, including `source/`, not just its CSVs.
Install `lineprof` in the Mac's own R library too. In a local R/RStudio session:

```r
source("validation/view_lineprof.R")
view_lineprof("eddie_profiles/lineprof_JOBID/g500")
```

`shine()` runs a local Shiny application and opens an interactive browser view;
it is not a standalone HTML file. The viewer remaps source paths from Eddie to
the downloaded snapshots. No Shiny server is launched inside the Eddie batch job.
If the legacy viewer fails with current Shiny, report the actual compatibility
error; using profvis instead would be an explicit subsequent choice.

The earlier Rprof profiler was moved to
`archive/cleanup_20261102/validation/profile_pca.R`. lineprof provides a different
display of the same underlying sampling approach.
