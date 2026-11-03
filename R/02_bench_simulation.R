simulate_benchmark_population <- function(params, config) {
  set.seed(config$seed + params$popSize + params$segSites + params$nGenPhase1 + params$nGenPhase2)

  nFounderPerPop <- params$popSize
  nFounderTotal <- 2L * nFounderPerPop

  founderGenomes <- AlphaSimR::runMacs(
    nInd = nFounderTotal,
    nChr = params$nChr,
    segSites = params$segSites,
    split = params$split,
    species = "GENERIC"
  )

  old_SP_exists <- exists("SP", envir = .GlobalEnv, inherits = FALSE)
  old_SP <- if (old_SP_exists) get("SP", envir = .GlobalEnv, inherits = FALSE) else NULL
  SP <- AlphaSimR::SimParam$new(founderPop = founderGenomes)
  assign("SP", SP, envir = .GlobalEnv)
  on.exit({
    if (old_SP_exists) assign("SP", old_SP, envir = .GlobalEnv) else rm("SP", envir = .GlobalEnv)
  }, add = TRUE)
  SP$setTrackPed(isTrackPed = TRUE)

  varG <- matrix(c(1.0, -0.3, -0.3, 1.0), byrow = TRUE, nrow = 2)
  varE <- matrix(c(2.0, 0.0, 0.0, 2.0), byrow = TRUE, nrow = 2)
  n_qtl <- max(1L, min(20L, floor(params$segSites * 0.20)))
  SP$addTraitA(nQtlPerChr = n_qtl, mean = c(0, 0), var = diag(varG), corA = varG)

  founders <- AlphaSimR::newPop(rawPop = founderGenomes, simParam = SP)
  founders <- AlphaSimR::setPheno(pop = founders, varE = diag(varE), simParam = SP)
  popA <- founders[seq_len(nFounderPerPop)]
  popB <- founders[(nFounderPerPop + 1L):nFounderTotal]

  data <- collect_benchmark_data(pop = popA, data = NULL, population = "A", generation = 0L)
  data <- collect_benchmark_data(pop = popB, data = data, population = "B", generation = 0L)

  parentsA <- popA
  parentsB <- popB
  if (params$nGenPhase1 > 0) {
    for (generation in seq_len(params$nGenPhase1)) {
      parentsA <- AlphaSimR::selectInd(pop = popA, nInd = min(length(popA@id), params$popSize), trait = 1, simParam = SP)
      parentsB <- AlphaSimR::selectInd(pop = popB, nInd = min(length(popB@id), params$popSize), trait = 2, simParam = SP)

      popA <- AlphaSimR::randCross(pop = parentsA, nCrosses = params$popSize, simParam = SP)
      popB <- AlphaSimR::randCross(pop = parentsB, nCrosses = params$popSize, simParam = SP)
      popA <- AlphaSimR::setPheno(pop = popA, varE = diag(varE), simParam = SP)
      popB <- AlphaSimR::setPheno(pop = popB, varE = diag(varE), simParam = SP)

      data <- collect_benchmark_data(pop = popA, data = data, population = "A", generation = generation)
      data <- collect_benchmark_data(pop = popB, data = data, population = "B", generation = generation)
    }
  }

  popAB <- AlphaSimR::randCross(pop = c(parentsA, parentsB), nCrosses = params$popSize, simParam = SP)
  popAB <- AlphaSimR::setPheno(pop = popAB, varE = diag(varE), simParam = SP)
  selIndexWeights <- AlphaSimR::smithHazel(econWt = c(1, 1), varG = varG, varP = varG + varE)

  if (params$nGenPhase2 > 0) {
    for (generation in (params$nGenPhase1 + 1L):(params$nGenPhase1 + params$nGenPhase2)) {
      parentsA <- AlphaSimR::selectInd(pop = popA, nInd = min(length(popA@id), params$popSize), trait = 1, simParam = SP)
      parentsB <- AlphaSimR::selectInd(pop = popB, nInd = min(length(popB@id), params$popSize), trait = 2, simParam = SP)
      parentsAB <- AlphaSimR::selectInd(
        pop = popAB,
        nInd = min(length(popAB@id), params$popSize),
        trait = AlphaSimR::selIndex,
        scale = TRUE,
        b = selIndexWeights,
        simParam = SP
      )

      popA <- AlphaSimR::randCross(pop = parentsA, nCrosses = params$popSize, simParam = SP)
      popB <- AlphaSimR::randCross(pop = parentsB, nCrosses = params$popSize, simParam = SP)
      popAB <- AlphaSimR::randCross(pop = parentsAB, nCrosses = params$popSize, simParam = SP)

      popA <- AlphaSimR::setPheno(pop = popA, varE = diag(varE), simParam = SP)
      popB <- AlphaSimR::setPheno(pop = popB, varE = diag(varE), simParam = SP)
      popAB <- AlphaSimR::setPheno(pop = popAB, varE = diag(varE), simParam = SP)

      data <- collect_benchmark_data(pop = popA, data = data, population = "A", generation = generation)
      data <- collect_benchmark_data(pop = popB, data = data, population = "B", generation = generation)
      data <- collect_benchmark_data(pop = popAB, data = data, population = "AB", generation = generation)
    }
  }

  ped <- data$pedigree
  ped$id <- as.character(ped$id)
  ped$fid <- as.character(ped$fid)
  ped$mid <- as.character(ped$mid)
  ped$fid[ped$fid %in% c("0", "NA", "")] <- NA
  ped$mid[ped$mid %in% c("0", "NA", "")] <- NA
  ped$fid[!is.na(ped$fid) & !(ped$fid %in% ped$id)] <- NA
  ped$mid[!is.na(ped$mid) & !(ped$mid %in% ped$id)] <- NA

  geno_ids <- select_genotyped_ids(ped, params, config)
  M <- data$genoIBS[geno_ids, , drop = FALSE]
  rownames(M) <- as.character(geno_ids)

  list(
    ped_df = ped,
    ped = pedigreeTools::pedigree(sire = ped$fid, dam = ped$mid, label = ped$id),
    M = M
  )
}

collect_benchmark_data <- function(pop, data = NULL, population, generation) {
  remove_first <- FALSE
  if (is.null(data)) {
    remove_first <- TRUE
    data <- vector(mode = "list", length = 2)
    names(data) <- c("pedigree", "genoIBS")
    data$pedigree <- data.frame(id = NA, population = NA, generation = NA, mid = NA, fid = NA)
    data$genoIBS <- matrix(data = NA, ncol = sum(pop@nLoci))
  }

  genoIBS <- AlphaSimR::pullSegSiteGeno(pop = pop)
  if (nrow(genoIBS) == length(pop@id)) rownames(genoIBS) <- pop@id

  data$pedigree <- rbind(
    data$pedigree,
    data.frame(
      id = pop@id,
      population = population,
      generation = generation,
      mid = pop@mother,
      fid = pop@father
    )
  )
  data$genoIBS <- rbind(data$genoIBS, genoIBS)

  if (remove_first) {
    data$pedigree <- data$pedigree[-1, ]
    data$genoIBS <- data$genoIBS[-1, ]
  }
  data
}

select_genotyped_ids <- function(ped, params, config) {
  if (!"generation" %in% names(ped) || !"population" %in% names(ped)) {
    stop("The simulated pedigree must contain generation and population columns.", call. = FALSE)
  }

  max_gen <- max(ped$generation, na.rm = TRUE)
  keep_generation <- if (identical(config$geno_generations, "last")) {
    max_gen
  } else if (is.numeric(config$geno_generations) && length(config$geno_generations) == 1L) {
    seq(max(0L, max_gen - config$geno_generations + 1L), max_gen)
  } else {
    config$geno_generations
  }
  candidates <- ped[ped$generation %in% keep_generation, , drop = FALSE]
  ids <- unlist(tapply(candidates$id, candidates$population, function(x) {
    n <- max(1L, ceiling(length(x) * params$geno_pct / 100))
    sample(x, min(n, length(x)))
  }), use.names = FALSE)
  as.character(ids[!is.na(ids)])
}
