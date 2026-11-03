# Exact-size, two-source hybrid population, simulated with AlphaSimR.
# A and B are maintained separately; C starts with A x B crosses.
simulate_abc_population <- function(n, markers = 1000L, seed = 12345L) {
  generations <- 10L
  hybrid_generation <- 5L
  n_chr <- 2L
  # Five generations of A/B, followed by five generations of A/B/C.
  group_slots <- 2L * hybrid_generation + 3L * (generations - hybrid_generation)
  if (n %% group_slots != 0L || n / group_slots < 2L) {
    stop("Total size must be divisible by 25, with at least two individuals per group")
  }
  if (markers %% n_chr != 0L) stop("Marker count must be divisible by two")
  group_size <- as.integer(n / group_slots)
  set.seed(seed)
  founders <- AlphaSimR::runMacs(
    nInd = 2L * group_size, nChr = n_chr,
    segSites = as.integer(markers / n_chr), split = 20L,
    species = "GENERIC", nThreads = 1L
  )
  SP <- AlphaSimR::SimParam$new(founders)
  SP$nThreads <- 1L
  SP$setTrackPed(TRUE)
  founder_pop <- AlphaSimR::newPop(founders, simParam = SP)
  A <- founder_pop[seq_len(group_size)]
  B <- founder_pop[group_size + seq_len(group_size)]
  C <- NULL
  pedigrees <- list()
  genotypes <- list()
  collect <- function(pop, population, generation) {
    pedigrees[[length(pedigrees) + 1L]] <<- data.frame(
      id = as.character(pop@id), population = population,
      fid = as.character(pop@father), mid = as.character(pop@mother),
      generation = generation
    )
    M <- AlphaSimR::pullSegSiteGeno(pop, simParam = SP)
    rownames(M) <- as.character(pop@id)
    genotypes[[length(genotypes) + 1L]] <<- M
  }
  collect(A, "A", 0L)
  collect(B, "B", 0L)
  for (generation in seq_len(generations - 1L)) {
    # All groups use parents from the immediately preceding generation.
    next_C <- if (generation == hybrid_generation) {
      AlphaSimR::randCross2(
        females = A, males = B, nCrosses = group_size,
        ignoreSexes = TRUE, simParam = SP
      )
    } else if (generation > hybrid_generation) {
      AlphaSimR::randCross(C, nCrosses = group_size, simParam = SP)
    } else NULL
    A <- AlphaSimR::randCross(A, nCrosses = group_size, simParam = SP)
    B <- AlphaSimR::randCross(B, nCrosses = group_size, simParam = SP)
    C <- next_C
    collect(A, "A", generation)
    collect(B, "B", generation)
    if (!is.null(C)) collect(C, "C", generation)
  }
  ped_df <- do.call(rbind, pedigrees)
  ped_df$fid[ped_df$fid == "0"] <- NA_character_
  ped_df$mid[ped_df$mid == "0"] <- NA_character_
  M <- do.call(rbind, genotypes)
  colnames(M) <- paste0("SNP_", seq_len(ncol(M)))
  list(
    ped = pedigreeTools::pedigree(
      sire = ped_df$fid, dam = ped_df$mid, label = ped_df$id
    ),
    ped_df = ped_df, M = M, n = n, seed = seed,
    generation_count = generations, hybrid_generation = hybrid_generation,
    group_size_per_generation = group_size,
    n_chr = n_chr, seg_sites_per_chr = as.integer(markers / n_chr), split = 20L,
    model = "AlphaSimR: separate A/B; C founded by A x B, followed by C x C"
  )
}
