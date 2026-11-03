# Internal worker for run_followups.py. Each invocation is a fresh R process.
args <- commandArgs(TRUE)
mode <- args[1]

if (mode == "subset") {
  master <- readRDS(args[2])
  target <- args[3]
  count <- as.integer(args[4])
  set.seed(as.integer(args[5]))
  order <- sample.int(nrow(master$M))
  stopifnot(count > 0L, count <= length(order))
  data <- master
  index <- sort(order[seq_len(count)])
  data$M <- master$M[index, , drop = FALSE]
  # All datasets retain exactly the same pedigree and shared marker values.
  stopifnot(identical(data$ped, master$ped),
            identical(data$M, master$M[index, , drop = FALSE]))
  saveRDS(data, target, compress = FALSE)
  write.csv(data.frame(id = rownames(data$M)), paste0(target, "_ids.csv"),
            row.names = FALSE)
  cat("Subset saved:", count, "genotyped of", data$n, "pedigree individuals\n")
} else if (mode == "reference") {
  data <- readRDS(args[2])
  center <- as.logical(args[3])
  rho <- as.integer(args[4])
  fit_file <- args[5]
  csv_file <- args[6]
  getH_sec <- system.time(
    H <- as.matrix(pedigreeTools::getH(data$ped, data$M,
                                     centred = "0.5", epsilon = 0.01))
  )[["elapsed"]]
  K <- H
  if (center) {
    K <- sweep(K, 2L, colMeans(K))
    K <- sweep(K, 1L, rowMeans(K))
  }
  symmetry_error <- norm(K - t(K), "F") / max(norm(K, "F"), .Machine$double.eps)
  K <- (K + t(K)) / 2
  exact_sec <- system.time(exact <- eigen(K, symmetric = TRUE))[["elapsed"]]
  reference <- list(K = K, d = exact$values[seq_len(rho)],
                    u = exact$vectors[, seq_len(rho), drop = FALSE],
                    trace = sum(diag(K)))
  reference$x <- sweep(reference$u, 2L, sqrt(pmax(reference$d, 0)), `*`)
  saveRDS(reference, fit_file, compress = FALSE)
  write.csv(data.frame(
    n = data$n, n_genotyped = nrow(data$M), markers = ncol(data$M),
    centredH = center, rho = rho, getH_sec = getH_sec, exact_eigen_sec = exact_sec,
    exact_trace = reference$trace, exact_min_eigenvalue = min(exact$values),
    reference_symmetry_error = symmetry_error,
    relative_eigengap_at_rho = (exact$values[rho] - exact$values[rho + 1L]) /
      max(abs(exact$values[rho]), .Machine$double.eps)
  ), csv_file, row.names = FALSE)
  cat("Exact reference saved for n =", data$n, "; centered =", center, "\n")
} else if (mode == "compare") {
  fit <- readRDS(args[2])
  reference <- readRDS(args[3])
  target <- args[4]
  rho <- length(fit$d)
  KU <- reference$K %*% fit$u
  residual <- KU - sweep(fit$u, 2L, fit$d, `*`)
  relative_residual <- sqrt(colSums(residual^2)) /
    pmax(sqrt(colSums(KU^2)), .Machine$double.eps)
  overlap <- sum(crossprod(fit$u, reference$u)^2)
  # Orthogonal alignment handles arbitrary signs and rotations in score plots.
  alignment <- svd(crossprod(fit$x, reference$x))
  aligned <- fit$x %*% (alignment$u %*% t(alignment$v))
  diag <- attr(fit, "diagnostics")
  trace_error <- if (is.null(diag$trace_estimate)) NA_real_ else
    as.numeric(diag$trace_estimate) / reference$trace - 1
  prop_error <- if (is.null(fit$varProps)) NA_real_ else
    max(abs(fit$varProps - reference$d / reference$trace))
  write.csv(data.frame(
    eigenvalue_max_rel_error = max(abs(fit$d - reference$d) /
      pmax(abs(reference$d), .Machine$double.eps)),
    normalised_subspace_error = sqrt(max(0, 1 - overlap / rho)),
    aligned_scores_relative_error = norm(aligned - reference$x, "F") /
      max(norm(reference$x, "F"), .Machine$double.eps),
    pc1_abs_cor = abs(cor(fit$x[, 1], reference$x[, 1])),
    pc2_abs_cor = abs(cor(fit$x[, 2], reference$x[, 2])),
    max_relative_eigenpair_residual = max(relative_residual),
    trace_relative_error = trace_error,
    variance_proportion_max_abs_error = prop_error
  ), target, row.names = FALSE)
  cat("Comparison with exact reference saved\n")
} else if (mode == "session") {
  loadNamespace("pedigreeTools")
  writeLines(capture.output(sessionInfo()), args[2])
} else {
  stop("Unknown follow-up worker mode: ", mode)
}
