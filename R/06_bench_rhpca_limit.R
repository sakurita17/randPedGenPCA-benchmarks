# Stress-test rhpca until it fails or reaches a configurable upper bound.
# Usage: Rscript R/06_bench_rhpca_limit.R

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

script_dir <- tryCatch({
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
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- args[grep("^--file=", args)][1]
  dirname(normalizePath(sub("^--file=", "", file_arg)))
}, error = function(e) {
  file.path(normalizePath(getwd()), "R")
}))

source(file.path(script_dir, "00_bench_config.R"))
source(file.path(script_dir, "01_bench_utils.R"))
source(file.path(script_dir, "02_bench_simulation.R"))
source(file.path(script_dir, "03_bench_population.R"))

load_benchmark_packages(bench_config)

limit_config <- bench_config
limit_config$output_dir <- file.path(
  "/Users/s2745345/Documents/GitHub/randPedGenPCA-benchmarks",
  "output",
  paste0(format(Sys.time(), "%Y%m%d_%H%M%S"), "_rhpca_limit")
)
dir.create(limit_config$output_dir, recursive = TRUE, showWarnings = FALSE)

# Edit here if you want a deeper or faster stress test.
pop_sizes <- c(20, 50, 100, 200, 500, 1000)
geno_pct <- 80
seg_sites <- 5000       # nChr = 2, so this is 10,000 markers.
max_elapsed_sec <- 600  # Treat longer runs as "timeout".
stop_after_first_failure <- TRUE

empty_limit_results <- function() {
  data.frame(
    popSize = integer(), geno_pct = numeric(), nChr = integer(), segSites = integer(),
    n_pedigree = integer(), n_genotyped = integer(), n_markers = integer(), pca_rank = integer(),
    status = character(), elapsed_sec = numeric(), elapsed_fmt = character(),
    mem_mb = numeric(), object_mb = numeric(), error_message = character(),
    stringsAsFactors = FALSE
  )
}

write_limit_results <- function(results) {
  saveRDS(results, file.path(limit_config$output_dir, "rhpca_limit_results.rds"))
  utils::write.csv(
    results,
    file.path(limit_config$output_dir, "rhpca_limit_results.csv"),
    row.names = FALSE
  )
}

run_rhpca_with_timeout <- function(scenario, config, max_elapsed_sec) {
  setTimeLimit(elapsed = max_elapsed_sec, transient = TRUE)
  on.exit(setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE), add = TRUE)

  tryCatch({
    rh_t <- timed(randPedGenPCA::rhpca(
      X = scenario$ped,
      M = scenario$M,
      rho = config$rhpca$rho,
      q = config$rhpca$q,
      l = config$rhpca$l,
      gamma = scenario$gamma,
      epsilon = config$rhpca$epsilon,
      centredH = config$rhpca$centredH,
      centred = "0.5",
      numVectorsTrace = config$rhpca$numVectorsTrace
    ))
    list(status = "ok", timing = rh_t, object = rh_t$value, error_message = NA_character_)
  }, error = function(e) {
    msg <- conditionMessage(e)
    status <- if (grepl("reached elapsed time limit|time limit", msg, ignore.case = TRUE)) "timeout" else "error"
    elapsed <- if (identical(status, "timeout")) max_elapsed_sec else NA_real_
    list(
      status = status,
      timing = list(elapsed = elapsed, user = NA_real_, system = NA_real_, mem_mb = NA_real_),
      object = NULL,
      error_message = msg
    )
  })
}

results <- empty_limit_results()

for (pop_size in pop_sizes) {
  params <- list(
    design = "rhpca_limit",
    geno_sampling = "fraction_last_generation",
    popSize = pop_size,
    geno_pct = geno_pct,
    nChr = 2L,
    segSites = seg_sites,
    split = 20L,
    nGenPhase1 = generations,
    nGenPhase2 = generations
  )

  message("[", Sys.time(), "] popSize=", pop_size, " segSites=", seg_sites)

  sim_result <- tryCatch({
    sim <- simulate_benchmark_population(params, limit_config)
    scenario <- build_no_metafounder(sim, limit_config)
    rh <- run_rhpca_with_timeout(scenario, limit_config, max_elapsed_sec)

    data.frame(
      popSize = pop_size,
      geno_pct = geno_pct,
      nChr = params$nChr,
      segSites = params$segSites,
      n_pedigree = scenario$n_pedigree,
      n_genotyped = nrow(scenario$M),
      n_markers = ncol(scenario$M),
      pca_rank = limit_config$rhpca$rho,
      status = if (identical(rh$status, "ok") && is.na(rh$timing$elapsed)) "incomplete" else rh$status,
      elapsed_sec = rh$timing$elapsed,
      elapsed_fmt = format_elapsed(rh$timing$elapsed),
      mem_mb = rh$timing$mem_mb,
      object_mb = if (is.null(rh$object)) NA_real_ else object_mb(rh$object),
      error_message = rh$error_message,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    data.frame(
      popSize = pop_size,
      geno_pct = geno_pct,
      nChr = params$nChr,
      segSites = params$segSites,
      n_pedigree = NA_integer_,
      n_genotyped = NA_integer_,
      n_markers = params$nChr * params$segSites,
      pca_rank = limit_config$rhpca$rho,
      status = "simulation_error",
      elapsed_sec = NA_real_,
      elapsed_fmt = NA_character_,
      mem_mb = NA_real_,
      object_mb = NA_real_,
      error_message = conditionMessage(e),
      stringsAsFactors = FALSE
    )
  })

  results <- rbind(results, sim_result)
  write_limit_results(results)

  if (stop_after_first_failure && !identical(sim_result$status, "ok")) {
    break
  }
}

plot_data <- results %>%
  mutate(
    status = factor(status, levels = c("ok", "incomplete", "timeout", "error", "simulation_error")),
    n_pedigree_plot = ifelse(is.na(n_pedigree), popSize * 52, n_pedigree),
    elapsed_plot = ifelse(is.na(elapsed_sec), max_elapsed_sec, elapsed_sec)
  )

p <- ggplot(plot_data, aes(x = n_pedigree_plot, y = elapsed_plot, color = status)) +
  geom_line(data = filter(plot_data, status == "ok"), linewidth = 1.2, color = "#009E73") +
  geom_point(size = 5) +
  geom_text(
    data = filter(plot_data, status != "ok"),
    aes(label = status),
    vjust = -0.9,
    size = 5,
    show.legend = FALSE
  ) +
  scale_x_log10(labels = function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)) +
  scale_y_log10(labels = function(x) signif(x, 2), na.value = NA_real_) +
  scale_color_manual(
    values = c(ok = "#009E73", incomplete = "#999999", timeout = "#D55E00", error = "#CC79A7", simulation_error = "#0072B2"),
    drop = FALSE
  ) +
  labs(
    title = "rhpca stress limit",
    subtitle = paste0(
      "10,000 markers, ", geno_pct, "% genotyped, k = ", limit_config$rhpca$rho,
      ", timeout = ", max_elapsed_sec, "s"
    ),
    x = "Number of individuals in pedigree (log)",
    y = "rhpca execution time (s, log)",
    color = "Status"
  ) +
  theme_bw(base_size = 18) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  )

ggsave(
  file.path(limit_config$output_dir, "plot7_rhpca_limit.png"),
  plot = p,
  width = 11,
  height = 7,
  dpi = 300
)

message("Saved results in: ", limit_config$output_dir)
