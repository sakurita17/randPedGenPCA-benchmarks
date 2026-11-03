# Scalability benchmark scripts

Entry point:

```r
source("R/run_scalability_benchmark.R")
```

or from shell:

```sh
Rscript R/run_scalability_benchmark.R
```

Edit `R/00_bench_config.R` to change only the simulation starting parameters. The runner then executes the full bottom-to-top path:

1. Simulate population with AlphaSimR.
2. Select genotyped individuals with `geno_pct` (percentage, 0-100), so `n_genotyped` grows with `popSize`.
3. Build three hybrid cases: no metafounder, single metafounder, two metafounders.
4. Estimate gamma with `pedigreeTools::getGammaL()` for the single/two metafounder cases.
5. Time `pedigreeTools::getLInv()`.
6. Time `randPedGenPCA::rhpca()` with `X = ped`.
7. Time `randPedGenPCA::rhpca()` with `X = Li`.
8. When `n_pedigree <= exact_pca$max_n_for_getH`, time `pedigreeTools::getH()` plus `stats::prcomp(rank. = rhpca$rho)` and compare PC1/PC2 against `rhpca`.

Outputs are written to `output/<timestamp>/`:

- `benchmark_results.csv`
- `benchmark_results.rds`
- `session_info.txt`

No plots are created.

The benchmark uses the current `rhpca` API: `rho` (number of PCs), `q`
(power iterations), and `l` (total range-finder vectors, with `l >= rho`).
These replace `k`, `numIter`, and `numVectors`; edit them in
`bench_config$rhpca`. `numVectorsTrace` still controls trace estimation.

The exact reference uses the same `epsilon` and, when `centredH = TRUE`,
double-centers H before `prcomp`. For a positive-semidefinite H, the unmodified
`prcomp(H, center=FALSE, scale.=FALSE)` scores already use `u * d`, matching the
algorithm's eigenvalue-scaled scores. No square-root rescaling is applied.
Older versions of this runner incorrectly rescaled the exact scores to
`u * sqrt(d)`; their coordinate RMSE must not be interpreted as an aligned
comparison with the current algorithm. For a positive-semidefinite H, `rhpca$d` contains
approximate eigenvalues. The `sdev` and variance proportions from `prcomp(H)`
are not directly comparable to those returned by `rhpca`.
The output column `pca_rank` remains unchanged for plotting compatibility.
