# Internal worker for run_large_population.py. No dense H or C is constructed.
args <- commandArgs(TRUE)
mode <- args[1]
if (mode == "prepare") {
  target <- args[2]; n <- as.integer(args[3]); markers <- as.integer(args[4])
  generations <- as.integer(args[5]); geno_generations <- as.integer(args[6])
  geno_fraction <- as.numeric(args[7]); seed <- as.integer(args[8])
  set.seed(seed)
  sizes <- rep(n %/% generations, generations)
  sizes[seq_len(n %% generations)] <- sizes[seq_len(n %% generations)] + 1L
  stopifnot(min(sizes) >= 4L, geno_generations <= generations)
  ids <- as.character(seq_len(n))
  sire <- dam <- rep(NA_character_, n)
  frequency <- runif(markers, 0.1, 0.9)
  retained <- list(); retained_ids <- list()
  previous <- NULL; previous_ids <- NULL; offset <- 0L
  tm <- system.time({
    for (generation in seq_len(generations)) {
      size <- sizes[generation]
      current_ids <- ids[offset + seq_len(size)]
      if (generation == 1L) {
        current <- matrix(rbinom(size * markers, 2L, rep(frequency, each = size)), size, markers)
      } else {
        # Disjoint sire/dam pools in the immediately preceding generation.
        split <- nrow(previous) %/% 2L
        fathers <- sample.int(split, size, replace = TRUE)
        mothers <- split + sample.int(nrow(previous) - split, size, replace = TRUE)
        sire[offset + seq_len(size)] <- previous_ids[fathers]
        dam[offset + seq_len(size)] <- previous_ids[mothers]
        paternal <- matrix(rbinom(size * markers, 1L, as.vector(previous[fathers, , drop = FALSE]) / 2), size, markers)
        maternal <- matrix(rbinom(size * markers, 1L, as.vector(previous[mothers, , drop = FALSE]) / 2), size, markers)
        current <- paternal + maternal
        rm(paternal, maternal)
      }
      if (generation > generations - geno_generations) {
        sampled <- sort(sample.int(size, max(1L, ceiling(size * geno_fraction))))
        retained[[length(retained) + 1L]] <- current[sampled, , drop = FALSE]
        retained_ids[[length(retained_ids) + 1L]] <- current_ids[sampled]
      }
      previous <- current; previous_ids <- current_ids; offset <- offset + size
    }
    M <- do.call(rbind, retained)
    rownames(M) <- unlist(retained_ids, use.names = FALSE)
    ped <- pedigreeTools::pedigree(sire = sire, dam = dam, label = ids)
  })[["elapsed"]]
  saveRDS(list(ped = ped, M = M, n = n, generation_sizes = sizes,
               generation_count = generations, geno_fraction = geno_fraction,
               geno_generations = geno_generations, seed = seed,
               simulation_sec = tm), target, compress = FALSE)
  cat("Prepared", n, "individuals;", nrow(M), "genotyped;", ncol(M), "markers\n")
} else if (mode == "pca") {
  package_root <- args[2]; dataset <- args[3]; result_file <- args[4]; fit_file <- args[5]
  q <- as.integer(args[6]); l <- as.integer(args[7]); rho <- as.integer(args[8])
  center <- as.logical(args[9]); trace_vectors <- as.integer(args[10]); seed <- as.integer(args[11])
  pkg <- new.env(parent = asNamespace("pedigreeTools"))
  for (file in c("utils.R", "rhpca.R", "hutchH.R"))
    sys.source(file.path(package_root, "R", file), envir = pkg)
  data <- readRDS(dataset)
  gc()
  set.seed(seed)
  warnings <- character()
  elapsed <- system.time(fit <- withCallingHandlers(
    pkg$rhpca.pedigree(data$ped, data$M, rho = rho, q = q, l = l,
      centredH = center, centred = "0.5", epsilon = 0.01,
      numVectorsTrace = if (trace_vectors == 0L) NULL else trace_vectors),
    warning = function(w) {warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")}
  ))[["elapsed"]]
  # Save PCA results before the separate diagnostic phase.
  saveRDS(fit, fit_file, compress = FALSE)
  diag <- attr(fit, "diagnostics")
  result <- data.frame(n = data$n, n_genotyped = nrow(data$M), markers = ncol(data$M),
    rho = rho, q = q, l = l, centredH = center, seed = seed,
    trace_vectors = trace_vectors, simulation_sec = data$simulation_sec,
    pca_elapsed_sec = elapsed, symmetry_error = diag$symmetry_error,
    sketch_rank = diag$sketch_rank,
    orth_fro_error = norm(crossprod(fit$u) - diag(rho), "F"),
    score_scaling_max_error = max(abs(fit$x - sweep(fit$u, 2L, sqrt(pmax(fit$d, 0)), `*`))),
    centred_scores_max_mean = if (center) max(abs(colMeans(fit$x))) else NA_real_,
    trace_estimate = if (is.null(diag$trace_estimate)) NA_real_ else as.numeric(diag$trace_estimate),
    explained_sum = if (is.null(fit$varProps)) NA_real_ else sum(fit$varProps),
    warnings = paste(warnings, collapse = "; "))
  write.csv(result, result_file, row.names = FALSE)
} else if (mode == "diagnose") {
  package_root <- args[2]; dataset <- args[3]; fit_file <- args[4]
  result_file <- args[5]; reference_file <- args[6]
  pkg <- new.env(parent = asNamespace("pedigreeTools"))
  for (file in c("utils.R", "rhpca.R", "hutchH.R"))
    sys.source(file.path(package_root, "R", file), envir = pkg)
  data <- readRDS(dataset); fit <- readRDS(fit_file)
  elapsed <- system.time({
    Li <- pedigreeTools::getLInv(data$ped)
    op <- pkg$makeHxOperator(pkg$inputHx(Li, data$M), fit$center, "0.5", NULL, 0.01, TRUE, 1000L)
    KU <- op(fit$u)
    residual <- KU - sweep(fit$u, 2L, fit$d, `*`)
    relative <- sqrt(colSums(residual^2)) / pmax(sqrt(colSums(KU^2)), .Machine$double.eps)
  })[["elapsed"]]
  change <- subspace <- NA_real_
  if (reference_file != "none") {
    ref <- readRDS(reference_file)
    change <- max(abs(fit$d - ref$d) / pmax(abs(ref$d), .Machine$double.eps))
    # Normalised projector distance, invariant to signs and within-space rotations.
    overlap <- sum(crossprod(fit$u, ref$u)^2)
    subspace <- sqrt(max(0, 1 - overlap / ncol(fit$u)))
  }
  write.csv(data.frame(diagnostic_elapsed_sec = elapsed,
    max_relative_eigenpair_residual = max(relative),
    eigenvalue_relative_change_vs_base = change,
    normalised_subspace_change_vs_base = subspace), result_file, row.names = FALSE)
} else stop("Unknown worker mode")
