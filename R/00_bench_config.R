# Scalability benchmark configuration.
# Edit this file to change the starting simulation parameters.

popSize    <- c(20, 50, 100) #seq(10, 1000, 100)
genoPerc   <- c(20, 50, 80)
segSites   <- c(100, 5000)   # 200 or 10,000 markers (nChr=2 * segSites)
generations <- 10



bench_config <- list(
  seed = 12345,

  package_roots = list(
    randPedGenPCA = "/Users/s2745345/Documents/GitHub/RandPedGenPCA/randPedGenPCA"
  ),

  output_dir = file.path(
    "/Users/s2745345/Documents/GitHub/randPedGenPCA-benchmarks",
    "output",
    format(Sys.time(), "%Y%m%d_%H%M%S")
  ),

  simulation_grid = expand.grid(
    design = "baseline_hybrid",
    geno_sampling = "fraction_last_generation",
    popSize = popSize,
    geno_pct = genoPerc,
    nChr = 2,
    segSites = segSites,
    split = 20,
    nGenPhase1 = generations,
    nGenPhase2 = generations,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  ),

  # geno_pct (0-100, see simulation_grid) controls what percentage of the
  # eligible animals get genotyped - it is swept per-run like popSize.
  geno_generations = 5L,

  scenarios = c(
    "hybrid_no_metafounder",
    "hybrid_single_metafounder",
    "hybrid_two_metafounders"
  ),
  
  centred_values = c("0.5", "observed"),

  rhpca = list(
    rho = 5,           # Number of PCs returned.
    q = 10,            # Power iterations.
    l = 50,            # Total range-finder vectors; oversampling = l - rho.
    epsilon = 0.01,
    centredH = FALSE,
    numVectorsTrace = 50
  ),

  gammaL = list(
    epsilon = 0.01
  ),

  exact_pca = list(
    enabled = TRUE,
    max_n_for_getH = Inf
  )
)
