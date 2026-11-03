# Small regression checks for the algorithm's X = U diag(d) score convention.
# Run from the repository root: Rscript tests/check_score_convention.R PACKAGE_ROOT
args <- commandArgs(TRUE)
stopifnot(length(args) == 1L)
package_root <- normalizePath(args[1], mustWork = TRUE)
worker <- normalizePath("validation/followup_worker.R", mustWork = TRUE)
large <- normalizePath("validation/large_population.R", mustWork = TRUE)
folder <- tempfile("score_convention_")
dir.create(folder)

run_worker <- function(script, ...) {
  arguments <- vapply(list(script, ...), function(x) shQuote(as.character(x)), character(1))
  output <- system2(file.path(R.home("bin"), "Rscript"), arguments, stdout = TRUE, stderr = TRUE)
  if (!is.null(attr(output, "status"))) stop(paste(output, collapse = "\n"))
}

# Exact scores must yield zero aligned error, including when a stored legacy
# reference contains sqrt(d)-scaled scores. This catches the original defect.
K <- diag(c(7, 5, 3, 1))
u <- diag(4)[, 1:2, drop = FALSE]
d <- c(7, 5)
fit <- list(u = u, d = d, x = sweep(u, 2L, d, `*`), varProps = d / sum(diag(K)))
attr(fit, "diagnostics") <- list(trace_estimate = sum(diag(K)))
ref <- list(K = K, u = u, d = d, x = sweep(u, 2L, sqrt(d), `*`), trace = sum(diag(K)))
saveRDS(fit, file.path(folder, "fit.rds"))
saveRDS(ref, file.path(folder, "legacy_reference.rds"))
target <- file.path(folder, "comparison.csv")
run_worker(worker, "compare", file.path(folder, "fit.rds"),
           file.path(folder, "legacy_reference.rds"), target)
metrics <- read.csv(target)
stopifnot(metrics$aligned_scores_relative_error < 1e-12,
          metrics$eigenvalue_max_rel_error < 1e-12,
          metrics$max_relative_eigenpair_residual < 1e-12,
          metrics$reference_score_scaling == "U*d")

# Test the actual saved worker output against a full-size-sketch PCA of H/CHC.
dataset <- file.path(folder, "data.rds")
run_worker(large, "prepare", dataset, 80, 20, 4, 2, 0.5, 203)
for (center in c(FALSE, TRUE)) {
  csv <- file.path(folder, paste0("pca_", center, ".csv"))
  fit_file <- file.path(folder, paste0("pca_", center, ".rds"))
  run_worker(large, "pca", package_root, dataset, csv, fit_file, 0, 80, 3,
             center, 10, 123)
  fit <- readRDS(fit_file)
  diagnostic <- read.csv(csv)
  stopifnot(diagnostic$score_scaling_max_error < 1e-12,
            diagnostic$score_scaling == "U*d",
            isTRUE(all.equal(fit$x, sweep(fit$u, 2L, fit$d, `*`), tolerance = 1e-12)))
  ref_file <- file.path(folder, paste0("reference_", center, ".rds"))
  run_worker(worker, "reference", dataset, center, 3, ref_file,
             file.path(folder, paste0("reference_", center, ".csv")))
  ref <- readRDS(ref_file)
  stopifnot(ref$score_scaling == "U*d",
            isTRUE(all.equal(ref$x, sweep(ref$u, 2L, ref$d, `*`), tolerance = 1e-12)))
  run_worker(worker, "compare", fit_file, ref_file, target)
  metrics <- read.csv(target)
  stopifnot(metrics$eigenvalue_max_rel_error < 1e-8,
            metrics$aligned_scores_relative_error < 1e-7)

  # Independent transcription of the supplied algorithm using a dense K and
  # eigen(B), compared to the implicit implementation using its reduced SVD.
  pkg <- new.env(parent = asNamespace("pedigreeTools"))
  for (file in c("utils.R", "rhpca.R", "hutchH.R"))
    sys.source(file.path(package_root, "R", file), envir = pkg)
  data <- readRDS(dataset)
  for (iterations in c(0L, 2L)) {
    set.seed(321)
    Y <- ref$K %*% matrix(rnorm(80 * 8), 80, 8)
    for (i in seq_len(iterations)) {
      Q <- qr.Q(qr(Y, LAPACK = TRUE))
      Y <- ref$K %*% Q
    }
    Q <- qr.Q(qr(Y, LAPACK = TRUE))
    B <- crossprod(Q, ref$K %*% Q)
    eig <- eigen((B + t(B)) / 2, symmetric = TRUE)
    expected_d <- eig$values[1:3]
    expected_x <- sweep(Q %*% eig$vectors[, 1:3, drop = FALSE], 2L, expected_d, `*`)
    set.seed(321)
    actual <- pkg$rhpca.pedigree(data$ped, data$M, rho = 3, q = iterations, l = 8,
                                 centredH = center, centred = "0.5", epsilon = 0.01)
    stopifnot(max(abs(actual$d - expected_d) / abs(expected_d)) < 1e-8,
              norm(tcrossprod(actual$x) - tcrossprod(expected_x), "F") /
                norm(tcrossprod(expected_x), "F") < 1e-7)
  }
}

# Check the original benchmark's prcomp reference without rescaling its scores.
source("R/04_bench_comparison.R")
timed <- function(expr) {
  value <- force(expr)
  list(value = value, timing = list(elapsed = 0, user = 0, system = 0))
}
data <- readRDS(dataset)
scenario <- list(ped = data$ped, M = data$M, gamma = NULL, n_pedigree = data$n)
config <- list(exact_pca = list(enabled = TRUE, max_n_for_getH = 100L),
               rhpca = list(rho = 3L, epsilon = 0.01))
for (center in c(FALSE, TRUE)) {
  config$rhpca$centredH <- center
  result <- benchmark_exact_pca(scenario, "0.5", config)$object
  ref <- readRDS(file.path(folder, paste0("reference_", center, ".rds")))
  # Score products are invariant to signs and rotations within equal eigenvalues.
  stopifnot(norm(tcrossprod(result$x) - tcrossprod(ref$x), "F") /
              norm(tcrossprod(ref$x), "F") < 1e-8)
}
unlink(folder, recursive = TRUE)
cat("Passed: supplied algorithm vs implicit PCA, U*d worker scores, exact references, legacy-reference comparison and prcomp score scale.\n")
