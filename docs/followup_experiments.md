# Two follow-up experiments on Eddie

Three files were added without modifying `run_large_population.py` or jobs
already submitted. Both experiments reuse `large_population.R` as their worker.

## Designs

**`exact`:** 200, 500 and 1,000 total individuals, 10 generations, 1,000 markers,
and 50% genotyped in each of the last two generations. For each dataset, H and
CHC are constructed with `centred="0.5"` and `epsilon=0.01`. An exact reference
is computed using `eigen(K, symmetric=TRUE)`. This decomposes the same target
operator; it is not PCA of the genotype matrix. Settings `q=4,l=20` and
`q=6,l=40` are compared, retaining five PCs and using three PCA seeds (123,
124, 125) on each fixed dataset: 36 approximate PCA fits and six references.
There are no replicate population simulations; the three seeds measure PCA
randomness.

Recorded metrics include eigenvalue and subspace errors, aligned-score errors,
absolute per-PC correlations, residuals and trace/variance-proportion errors.
Individual-PC correlations can be less informative for nearly repeated
eigenvalues; inspect subspaces and the recorded reference eigengap. Approximate
variance proportions include Hutch++ error (50 trace vectors). The reference
orders eigenvalues algebraically and records the minimum eigenvalue to flag
potential positive-semidefiniteness problems. Score alignment assumes the
reference score scale, now `U * diag(d)` as required by the algorithm.
Historical reference files used `U * diag(sqrt(d))`. The corrected comparison
worker reconstructs the reference scores from its saved eigenvectors and
eigenvalues, allowing legacy reference RDS files to be reused without rerunning
the eigendecomposition or PCA. Old CSV coordinate-error metrics require this
recomputation; the summary marks them as not aligned rather than reporting them
as valid errors.

**`genotyped`:** a single pedigree of 20,000 individuals with 1,000 markers.
Genotypes for all 4,000 animals in the last two generations are simulated once.
A fixed permutation defines nested subsets of 500, 1,000, 2,000 and 4,000
genotyped animals. Shared animals retain identical genotypes and every dataset
retains exactly the same pedigree. Sampling is random across the two generations
together, not an exact quota per generation. These counts represent 2.5%, 5%,
10% and 20% of the full pedigree.

H and CHC are evaluated with `q=4,l=20`, five PCs, 50 trace vectors and one PCA
seed (123): eight PCA fits and eight implicit diagnostic runs. This is an
initial cost-versus-genotyping experiment without an exact reference.
Changing the genotyped animals changes the operator; PC differences between
subsets are not by themselves approximation errors. For more seeds, add
`--replicates 3`, with the corresponding additional computation.

## Copy from the Mac

Run in the Mac terminal, outside the SSH session:

```sh
cd ~/Documents/GitHub/randPedGenPCA-benchmarks
rsync -av validation/run_followups.py validation/followup_worker.R \
  validation/followup_job.sh \
  s2745345@eddie.ecdf.ed.ac.uk:HighlanderLab/jortiz/pca-benchmark/benchmarks/validation/
```

Do not overwrite Eddie's `run_large_population.py` with the local copy: the
Eddie copy was previously modified to accept `--timeout 0`.

## Check and submit on Eddie

```sh
cd ~/HighlanderLab/jortiz/pca-benchmark/benchmarks
bash -n validation/followup_job.sh
python3 validation/run_followups.py exact --package-root ../randPedGenPCA --dry-run
python3 validation/run_followups.py genotyped --package-root ../randPedGenPCA --dry-run

# Replace this address if it is not your receiving mailbox.
BENCH_EMAIL="s2745345@ed.ac.uk"
qsub -N pca_exact -m ea -M "$BENCH_EMAIL" validation/followup_job.sh exact
qsub -N pca_geno -m ea -M "$BENCH_EMAIL" validation/followup_job.sh genotyped
qstat -u "$USER"
```

Each job requests one core, 32 GB of virtual memory and 48 hours. Python does
not impose an internal timeout; Eddie's scheduler limit still applies. Both
experiments later completed on Eddie; see the English results summary for
measurements. Nothing is submitted automatically: the user executes the
commands above.

Arguments after the experiment name can override the design:

```sh
qsub -N pca_geno validation/followup_job.sh genotyped \
  --pedigree-size 20000 --genotyped 500 1000 2000 4000
```

## Results

Each invocation creates a new `output/followup_exact_<timestamp>/` or
`output/followup_genotyped_<timestamp>/` directory. It contains configuration,
source SHA-256 checksums, session information, datasets, PCA/reference RDS
files, individual CSVs and logs, and `results.csv`, updated atomically after
each completed process.

In `results.csv`, filter `stage == "pca"` for PCA time and RSS.
`pca_elapsed_sec` excludes loading, saving and diagnostics; `process_wall_sec`
includes the entire process. Peak RSS refers to the entire process, not
incremental PCA memory. References and diagnostics run in separate processes.
`stage == "compare_exact"` contains errors against the exact reference; join
rows using n, centredH, q, l and seed. Genotyping-experiment diagnostic rows
have `stage == "diagnose"` and also require n_genotyped in the join key.

Datasets are retained. There is no automatic resume: an interruption preserves
already written stages but not the internal state of a running PCA. Successful
completion is indicated by `Finished. Results:`.

## Small local checks

Use new output directories and set `--package-root` for the local machine:

```sh
python3 validation/run_followups.py exact --sizes 80 --markers 20 \
  --generations 4 --rho 3 --q 2 --sketch 8 --trace-vectors 10 --replicates 1 \
  --package-root ../RandPedGenPCA/randPedGenPCA --output /tmp/followup_exact_smoke
python3 validation/run_followups.py genotyped --pedigree-size 120 \
  --genotyped 8 16 32 40 --markers 20 --generations 6 --rho 3 --q 2 \
  --sketch 8 --trace-vectors 10 --replicates 1 \
  --package-root ../RandPedGenPCA/randPedGenPCA --output /tmp/followup_geno_smoke
Rscript tests/check_followup_outputs.R /tmp/followup_exact_smoke \
  /tmp/followup_geno_smoke ../RandPedGenPCA/randPedGenPCA
```
