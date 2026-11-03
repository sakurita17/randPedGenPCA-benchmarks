# Generate only the simulated populations, without running PCA benchmarks.

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

load_benchmark_packages(bench_config)

simulate_population_only <- function(params, config) {
  set.seed(config$seed + params$popSize + params$segSites + params$nGenPhase1 + params$nGenPhase2)

  max_gen <- params$nGenPhase1 + params$nGenPhase2
  keep_generation <- if (identical(config$geno_generations, "last")) {
    max_gen
  } else if (is.numeric(config$geno_generations) && length(config$geno_generations) == 1L) {
    seq(max(0L, max_gen - config$geno_generations + 1L), max_gen)
  } else {
    config$geno_generations
  }

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

  ped_chunks <- list()
  geno_chunks <- list()

  append_population <- function(pop, population, generation) {
    ped_chunks[[length(ped_chunks) + 1L]] <<- data.frame(
      id = pop@id,
      population = population,
      generation = generation,
      mid = pop@mother,
      fid = pop@father
    )

    if (generation %in% keep_generation) {
      genoIBS <- AlphaSimR::pullSegSiteGeno(pop = pop)
      if (nrow(genoIBS) == length(pop@id)) rownames(genoIBS) <- pop@id
      geno_chunks[[length(geno_chunks) + 1L]] <<- genoIBS
    }
  }

  founders <- AlphaSimR::newPop(rawPop = founderGenomes, simParam = SP)
  founders <- AlphaSimR::setPheno(pop = founders, varE = diag(varE), simParam = SP)
  popA <- founders[seq_len(nFounderPerPop)]
  popB <- founders[(nFounderPerPop + 1L):nFounderTotal]

  append_population(popA, "A", 0L)
  append_population(popB, "B", 0L)

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

      append_population(popA, "A", generation)
      append_population(popB, "B", generation)
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

      append_population(popA, "A", generation)
      append_population(popB, "B", generation)
      append_population(popAB, "AB", generation)
    }
  }

  ped <- do.call(rbind, ped_chunks)
  ped$id <- as.character(ped$id)
  ped$fid <- as.character(ped$fid)
  ped$mid <- as.character(ped$mid)
  ped$fid[ped$fid %in% c("0", "NA", "")] <- NA
  ped$mid[ped$mid %in% c("0", "NA", "")] <- NA
  ped$fid[!is.na(ped$fid) & !(ped$fid %in% ped$id)] <- NA
  ped$mid[!is.na(ped$mid) & !(ped$mid %in% ped$id)] <- NA

  geno_ids <- select_genotyped_ids(ped, params, config)
  genoIBS <- do.call(rbind, geno_chunks)
  M <- genoIBS[geno_ids, , drop = FALSE]
  rownames(M) <- as.character(geno_ids)

  list(
    ped_df = ped,
    ped = pedigreeTools::pedigree(sire = ped$fid, dam = ped$mid, label = ped$id),
    M = M
  )
}

population_sizes <- data.frame(
  population_size_label = c("small", "medium", "large"),
  popSize = c(5000L, 10000L, 20000L),
  stringsAsFactors = FALSE
)

base_params <- as.list(bench_config$simulation_grid[1, , drop = FALSE])
base_params$geno_pct <- 100
base_params$segSites <- min(bench_config$simulation_grid$segSites)

output_dir <- file.path(
  dirname(bench_config$output_dir),
  paste0(format(Sys.time(), "%Y%m%d_%H%M%S"), "_populations_only")
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

manifest <- data.frame(
  population_size_label = character(),
  popSize = integer(),
  n_pedigree = integer(),
  n_genotyped = integer(),
  n_markers = integer(),
  file = character(),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(population_sizes))) {
  params <- base_params
  params$popSize <- population_sizes$popSize[i]

  label <- population_sizes$population_size_label[i]
  message("[", Sys.time(), "] generating ", label, " population",
          " popSize=", params$popSize,
          " nChr=", params$nChr,
          " segSites=", params$segSites)

  sim <- simulate_population_only(params, bench_config)
  file_name <- paste0("population_", label, "_popSize_", params$popSize, ".rds")
  saveRDS(
    list(label = label, params = params, population = sim),
    file = file.path(output_dir, file_name)
  )

  manifest <- rbind(
    manifest,
    data.frame(
      population_size_label = label,
      popSize = params$popSize,
      n_pedigree = length(sim$ped@label),
      n_genotyped = nrow(sim$M),
      n_markers = ncol(sim$M),
      file = file_name,
      stringsAsFactors = FALSE
    )
  )

  utils::write.csv(
    manifest,
    file.path(output_dir, "population_manifest.csv"),
    row.names = FALSE
  )
}

utils::capture.output(
  utils::sessionInfo(),
  file = file.path(output_dir, "session_info.txt")
)

message("Done. Populations saved in: ", output_dir)
