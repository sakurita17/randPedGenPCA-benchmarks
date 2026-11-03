# Run from the repository root: Rscript tests/smoke_benchmark.R
for (file in c("00_bench_config.R", "01_bench_utils.R", "02_bench_simulation.R",
               "03_bench_population.R", "04_bench_comparison.R")) {
  source(file.path("R", file))
}
load_benchmark_packages(bench_config)
config <- bench_config
config$rhpca$rho <- 3L
config$rhpca$q <- 2L
config$rhpca$numVectorsTrace <- 5L
config$rhpca$epsilon <- 0.1
config$geno_generations <- "last"
params <- list(design = "smoke", geno_sampling = "fraction_last_generation",
               popSize = 6L, geno_pct = 100, nChr = 1L, segSites = 40L,
               split = 20L, nGenPhase1 = 1L, nGenPhase2 = 1L)
sim <- simulate_benchmark_population(params, config)
config$rhpca$l <- length(sim$ped@label) # Full sketch satisfying l <= n.
for (centredH in c(FALSE, TRUE)) {
  config$rhpca$centredH <- centredH
  for (name in config$scenarios) {
    scenario <- scenario_builders[[name]](sim, config)
    results <- benchmark_scenario(1L, params, scenario, config, empty_result_table())
    rh <- results[results$step == "rhpca", ]
    exact <- results[results$step == "getH_prcomp", ]
    stopifnot(nrow(rh) == 4L, all(rh$pca_rank == 3L),
              all(is.finite(rh$elapsed_sec)), nrow(exact) == 2L,
              all(exact$exact_status == "ok"),
              all(exact$pc1_abs_cor > 0.999), all(exact$pc2_abs_cor > 0.999),
              all(exact$pc1_rmse < 1e-6), all(exact$pc2_rmse < 1e-6))
    message("Passed: ", name, ", centredH = ", centredH)
  }
}
# Load only the stress-test function, without running the large benchmark.
exprs <- parse("R/06_bench_rhpca_limit.R")
for (expr in exprs) {
  if (is.call(expr) && identical(expr[[1L]], as.name("<-")) &&
      identical(expr[[2L]], as.name("run_rhpca_with_timeout"))) eval(expr)
}
rh <- run_rhpca_with_timeout(build_no_metafounder(sim), config, 60)
stopifnot(rh$status == "ok", ncol(rh$object$x) == config$rhpca$rho)
config$exact_pca$enabled <- FALSE
stopifnot(benchmark_exact_pca(scenario, "0.5", config)$status == "disabled")
config$exact_pca$enabled <- TRUE
config$exact_pca$max_n_for_getH <- 1L
stopifnot(benchmark_exact_pca(scenario, "0.5", config)$status == "skipped_n_too_large")
message("Benchmark smoke checks passed.")
