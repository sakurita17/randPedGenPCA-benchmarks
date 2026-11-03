# End-to-end scalability benchmark.
# It starts at simulation, builds hybrid no/single/two-metafounder cases,
# calculates gamma where needed, times getLInv/gammaL/rhpca, and compares
# rhpca against prcomp(getH) when the dense H matrix is still feasible.

#script_dir <- file.path(getwd(), "R") #Run scenarios line by line
script_dir <- tryCatch({
  # source() from RStudio or console: look for ofile in active frames
  ofile <- NULL
  for (i in rev(seq_along(sys.frames()))) {
    if (!is.null(sys.frames()[[i]]$ofile)) {
      ofile <- sys.frames()[[i]]$ofile
      break
    }
  }
  if (is.null(ofile)) stop()
  dirname(normalizePath(ofile))
}, error = function(e) tryCatch({
  # Rscript from terminal: parse --file= argument
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- args[grep("^--file=", args)][1]
  dirname(normalizePath(sub("^--file=", "", file_arg)))
}, error = function(e) {
  # Fallback: assume getwd() is the project root
  file.path(normalizePath(getwd()), "R")
}))

source(file.path(script_dir, "00_bench_config.R"))
source(file.path(script_dir, "01_bench_utils.R"))
source(file.path(script_dir, "02_bench_simulation.R"))
source(file.path(script_dir, "03_bench_population.R"))
source(file.path(script_dir, "04_bench_comparison.R"))

load_benchmark_packages(bench_config)
dir.create(bench_config$output_dir, recursive = TRUE, showWarnings = FALSE)

results <- empty_result_table()
run_id <- 0L

for (i in seq_len(nrow(bench_config$simulation_grid))) {
  params <- as.list(bench_config$simulation_grid[i, , drop = FALSE])
  run_id <- run_id + 1L

  message("[", Sys.time(), "] simulation run_id=", run_id,
          " popSize=", params$popSize,
          " nChr=", params$nChr,
          " segSites=", params$segSites)

  sim_t <- timed(simulate_benchmark_population(params, bench_config))
  sim <- sim_t$value

  sim_row_scenario <- list(name = "simulation", nmf = NA_integer_, n_pedigree = length(sim$ped@label), gamma = NULL)
  results <- append_result(results, make_result_row(
    run_id, params, sim_row_scenario, centred = NA_character_, step = "simulation",
    method = "AlphaSimR", x_type = NA_character_, timing = sim_t,
    object = sim, gamma_dim = NA_character_, gamma_values = NA_character_
  ))
  write_incremental_results(results, bench_config)

  for (scenario_name in bench_config$scenarios) {
    message("[", Sys.time(), "] scenario=", scenario_name)
    scenario <- scenario_builders[[scenario_name]](sim, bench_config)
    results <- benchmark_scenario(run_id, params, scenario, bench_config, results)
    write_incremental_results(results, bench_config)
  }
}

utils::capture.output(utils::sessionInfo(), file = file.path(bench_config$output_dir, "session_info.txt"))
message("Done. Results saved in: ", bench_config$output_dir)



# TODO: Bigger populations really big
