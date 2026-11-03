load_benchmark_packages <- function(config) {
  required <- c("AlphaSimR", "pedigreeTools", "randPedGenPCA", "bench")

  if (!is.null(config$package_roots$randPedGenPCA) &&
      dir.exists(config$package_roots$randPedGenPCA) &&
      file.exists(file.path(config$package_roots$randPedGenPCA, "DESCRIPTION")) &&
      requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all(config$package_roots$randPedGenPCA, quiet = TRUE)
  }

  invisible(lapply(required, function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package not installed: ", pkg, call. = FALSE)
    }
  }))
}

as_seconds <- function(proc_time) {
  unname(proc_time[["elapsed"]])
}

timed <- function(expr) {
  gc(reset = TRUE)
  if (requireNamespace("bench", quietly = TRUE)) {
    # substitute() captures the unevaluated expression before bench sees it.
    # eval() inside bench::mark re-evaluates it fresh on every pass (including
    # any memory-measurement dry run), so R never serves a cached promise value.
    expr_q   <- substitute(expr)
    call_env <- parent.frame()
    .env <- environment()
    .env$.result <- NULL
    b <- bench::mark(
      { .env$.result <- eval(expr_q, envir = call_env) },
      iterations = 1L, memory = TRUE, check = FALSE
    )
    list(
      value   = .env$.result,
      elapsed = as.numeric(b$median),
      user    = NA_real_,
      system  = NA_real_,
      mem_mb  = as.numeric(b$mem_alloc) / 1024^2
    )
  } else {
    warning("Package 'bench' not available: mem_mb will be NA.", call. = FALSE)
    start <- proc.time()
    value <- force(expr)
    elapsed <- proc.time() - start
    list(
      value   = value,
      elapsed = as_seconds(elapsed),
      user    = unname(elapsed[["user.self"]]),
      system  = unname(elapsed[["sys.self"]]),
      mem_mb  = NA_real_
    )
  }
}

object_mb <- function(x) {
  as.numeric(utils::object.size(x)) / 1024^2
}

format_elapsed <- function(sec) {
  ifelse(is.na(sec), NA_character_,
    ifelse(sec < 1,    paste0(round(sec * 1000, 1), " ms"),
    ifelse(sec < 60,   paste0(round(sec, 1),        " s"),
    ifelse(sec < 3600, paste0(round(sec / 60, 1),   " min"),
                       paste0(round(sec / 3600, 1), " h")))))
}

report_cols <- c(
  "run_id", "scenario", "step", "method", "x_type", "centred",
  "popSize", "geno_pct", "n_pedigree", "n_genotyped", "n_markers", "pca_rank", "nmf",
  "elapsed_sec", "elapsed_fmt", "mem_mb", "object_mb", "exact_status",
  "pc1_abs_cor", "pc2_abs_cor", "pc1_rmse", "pc2_rmse"
)

write_incremental_results <- function(results, config) {
  dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(results, file.path(config$output_dir, "benchmark_results.rds"))
  keep <- intersect(report_cols, names(results))
  utils::write.csv(results[, keep, drop = FALSE],
                   file.path(config$output_dir, "benchmark_results.csv"),
                   row.names = FALSE)
}

empty_result_table <- function() {
  data.frame(
    run_id = integer(), design = character(), geno_sampling = character(),
    scenario = character(), centred = character(), step = character(),
    method = character(), x_type = character(), popSize = integer(), geno_pct = numeric(),
    nChr = integer(), segSites = integer(), split = integer(), nGenPhase1 = integer(), nGenPhase2 = integer(),
    n_pedigree = integer(), n_genotyped = integer(), n_markers = integer(), pca_rank = integer(), nmf = integer(),
    elapsed_sec = numeric(), elapsed_fmt = character(), user_sec = numeric(), system_sec = numeric(), mem_mb = numeric(), object_mb = numeric(),
    gamma_dim = character(), gamma_values = character(), exact_status = character(),
    pc1_abs_cor = numeric(), pc2_abs_cor = numeric(),
    pc1_rmse = numeric(), pc2_rmse = numeric(),
    stringsAsFactors = FALSE
  )
}

append_result <- function(results, row) {
  missing_cols <- setdiff(names(results), names(row))
  for (nm in missing_cols) row[[nm]] <- NA
  extra_cols <- setdiff(names(row), names(results))
  for (nm in extra_cols) results[[nm]] <- NA
  rbind(results, row[names(results)])
}

# Meention the use of PCA on H matrix for pop correction whenn gwas.
# 
