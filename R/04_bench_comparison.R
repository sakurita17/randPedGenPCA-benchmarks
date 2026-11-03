benchmark_scenario <- function(run_id, params, scenario, config, results) {
  gamma_values <- if (is.null(scenario$gamma)) NA_character_ else paste(round(as.numeric(scenario$gamma), 8), collapse = ";")
  gamma_dim <- if (is.null(scenario$gamma)) NA_character_ else paste(dim(scenario$gamma), collapse = "x")

  if (!is.null(scenario$gamma_timing)) {
    results <- append_result(results, make_result_row(
      run_id, params, scenario, centred = NA_character_, step = "gammaL",
      method = "pedigreeTools::getGammaL", x_type = "ped", timing = scenario$gamma_timing,
      object = scenario$gamma, gamma_dim = gamma_dim, gamma_values = gamma_values
    ))
  }

  getlinv_t <- timed(pedigreeTools::getLInv(ped = scenario$ped, gamma = scenario$gamma))
  Li <- getlinv_t$value
  results <- append_result(results, make_result_row(
    run_id, params, scenario, centred = NA_character_, step = "getLInv",
    method = "pedigreeTools::getLInv", x_type = "ped", timing = getlinv_t,
    object = Li, gamma_dim = gamma_dim, gamma_values = gamma_values
  ))

  for (centred in config$centred_values) {
    rh_ped_t <- timed(randPedGenPCA::rhpca(
      X = scenario$ped,
      M = scenario$M,
      rho = config$rhpca$rho,
      q = config$rhpca$q,
      l = config$rhpca$l,
      gamma = scenario$gamma,
      epsilon = config$rhpca$epsilon,
      centredH = config$rhpca$centredH,
      centred = centred,
      numVectorsTrace = config$rhpca$numVectorsTrace
    ))
    results <- append_result(results, make_result_row(
      run_id, params, scenario, centred = centred, step = "rhpca",
      method = "randPedGenPCA::rhpca", x_type = "ped", timing = rh_ped_t,
      object = rh_ped_t$value, gamma_dim = gamma_dim, gamma_values = gamma_values,
      pca_rank = config$rhpca$rho
    ))

    rh_li_t <- timed(randPedGenPCA::rhpca(
      X = Li,
      M = scenario$M,
      rho = config$rhpca$rho,
      q = config$rhpca$q,
      l = config$rhpca$l,
      epsilon = config$rhpca$epsilon,
      centredH = config$rhpca$centredH,
      centred = centred,
      numVectorsTrace = config$rhpca$numVectorsTrace
    ))
    results <- append_result(results, make_result_row(
      run_id, params, scenario, centred = centred, step = "rhpca",
      method = "randPedGenPCA::rhpca", x_type = "Li", timing = rh_li_t,
      object = rh_li_t$value, gamma_dim = gamma_dim, gamma_values = gamma_values,
      pca_rank = config$rhpca$rho
    ))

    exact <- benchmark_exact_pca(scenario, centred, config)
    results <- append_result(results, make_result_row(
      run_id, params, scenario, centred = centred, step = "getH_prcomp",
      method = "pedigreeTools::getH + stats::prcomp", x_type = "ped", timing = exact$timing,
      object = exact$object, gamma_dim = gamma_dim, gamma_values = gamma_values,
      exact_status = exact$status, comparison = compare_pca_scores(rh_ped_t$value, exact$object),
      pca_rank = config$rhpca$rho
    ))
  }

  results
}

benchmark_exact_pca <- function(scenario, centred, config) {
  if (!isTRUE(config$exact_pca$enabled)) {
    return(list(status = "disabled", timing = list(elapsed = NA, user = NA, system = NA, mem_mb = NA_real_), object = NULL))
  }
  if (scenario$n_pedigree > config$exact_pca$max_n_for_getH) {
    return(list(status = "skipped_n_too_large", timing = list(elapsed = NA, user = NA, system = NA, mem_mb = NA_real_), object = NULL))
  }

  exact_t <- timed({
    H <- as.matrix(pedigreeTools::getH(
      ped = scenario$ped, M = scenario$M, gamma = scenario$gamma,
      centred = centred, epsilon = config$rhpca$epsilon
    ))
    if (isTRUE(config$rhpca$centredH)) {
      # Match the implicit C H C operator used by rhpca.
      H <- sweep(H, 2L, colMeans(H), FUN = "-")
      H <- sweep(H, 1L, rowMeans(H), FUN = "-")
    }
    fit <- stats::prcomp(H, center = FALSE, scale. = FALSE, rank. = config$rhpca$rho)
    # For PSD H, prcomp(H) already returns eigenvalue-scaled scores U * d,
    # matching the randomised-PCA algorithm. Do not rescale by sqrt(d).
    fit
  })
  list(status = "ok", timing = exact_t, object = exact_t$value)
}

compare_pca_scores <- function(rhpca_obj, prcomp_obj) {
  out <- list(pc1_abs_cor = NA_real_, pc2_abs_cor = NA_real_, pc1_rmse = NA_real_, pc2_rmse = NA_real_)
  if (is.null(rhpca_obj) || is.null(prcomp_obj) || is.null(rhpca_obj$x) || is.null(prcomp_obj$x)) return(out)

  common <- intersect(rownames(rhpca_obj$x), rownames(prcomp_obj$x))
  if (length(common) < 3) return(out)

  for (pc in 1:2) {
    if (ncol(rhpca_obj$x) < pc || ncol(prcomp_obj$x) < pc) next
    a <- rhpca_obj$x[common, pc]
    b <- prcomp_obj$x[common, pc]
    cors <- suppressWarnings(stats::cor(a, b, use = "complete.obs"))
    sign_align <- ifelse(is.na(cors) || cors >= 0, 1, -1)
    out[[paste0("pc", pc, "_abs_cor")]] <- abs(cors)
    out[[paste0("pc", pc, "_rmse")]] <- sqrt(mean((as.numeric(a) - sign_align * as.numeric(b))^2, na.rm = TRUE))
  }
  out
}

make_result_row <- function(run_id, params, scenario, centred, step, method, x_type,
                            timing, object, gamma_dim, gamma_values, exact_status = NA_character_,
                            comparison = NULL, pca_rank = NA_integer_) {
  if (is.null(comparison)) comparison <- list(pc1_abs_cor = NA_real_, pc2_abs_cor = NA_real_, pc1_rmse = NA_real_, pc2_rmse = NA_real_)

  data.frame(
    run_id = run_id,
    design = if (!is.null(params$design)) params$design else NA_character_,
    geno_sampling = if (!is.null(params$geno_sampling)) params$geno_sampling else NA_character_,
    scenario = scenario$name,
    centred = centred,
    step = step,
    method = method,
    x_type = x_type,
    popSize = params$popSize,
    geno_pct = params$geno_pct,
    nChr = params$nChr,
    segSites = params$segSites,
    split = params$split,
    nGenPhase1 = params$nGenPhase1,
    nGenPhase2 = params$nGenPhase2,
    n_pedigree = scenario$n_pedigree,
    n_genotyped = if (is.null(scenario$M)) NA_integer_ else nrow(scenario$M),
    n_markers = if (is.null(scenario$M)) NA_integer_ else ncol(scenario$M),
    pca_rank = pca_rank,
    nmf = scenario$nmf,
    elapsed_sec = timing$elapsed,
    elapsed_fmt = format_elapsed(timing$elapsed),
    user_sec = timing$user,
    system_sec = timing$system,
    mem_mb = if (is.null(timing$mem_mb)) NA_real_ else timing$mem_mb,
    object_mb = if (is.null(object)) NA_real_ else object_mb(object),
    gamma_dim = gamma_dim,
    gamma_values = gamma_values,
    exact_status = exact_status,
    pc1_abs_cor = comparison$pc1_abs_cor,
    pc2_abs_cor = comparison$pc2_abs_cor,
    pc1_rmse = comparison$pc1_rmse,
    pc2_rmse = comparison$pc2_rmse,
    stringsAsFactors = FALSE
  )
}
