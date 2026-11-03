# Large-population hybrid PCA experiment

This experiment measures time and memory as the pedigree grows, without
explicitly constructing H or CHC. It does not modify the randPedGenPCA core.

## Default design

| Total individuals | Generations, including founders | Genotyped animals | Markers |
|---|---|---|---|
| 10,000 | 10 | 1,000 | 1,000 |
| 50,000 | 10 | 5,000 | 1,000 |
| 100,000 | 10 | 10,000 | 1,000 |

A synthetic population is simulated with random mating between sire and dam
pools from the previous generation. Founder allele frequencies are uniform
between 0.1 and 0.9; offspring receive one allele per marker from each parent.
Markers segregate independently. Linkage, selection and multiple populations
are not simulated. This design measures computational cost and numerical
stability; it is not a complete biological model.

To reduce memory use, the simulator retains the previous generation and saves
genotypes only for selected individuals. By default, 50% of each of the last
two generations is genotyped: 10% of the full pedigree. Fractional sample sizes
are rounded upwards when generation sizes are not divisible.

PCA uses `rho=5`, `q=4`, `l=20`, `numVectorsTrace=50`, `epsilon=0.01` and
`centred="0.5"`. H and CHC are evaluated in separate processes. There are no
metafounders. This configuration is not assumed to guarantee accuracy for
every population.

## Run the experiment

From the repository root:

```sh
# Inspect the design without starting calculations.
python3 validation/run_large_population.py --dry-run

# Initial experiment: 10,000 individuals.
python3 validation/run_large_population.py --sizes 10000

# Full size ladder: 10,000, 50,000 and 100,000 individuals.
python3 validation/run_large_population.py

# Compare q=4,l=20 against q=6,l=40 to assess stability.
python3 validation/run_large_population.py --compare

# Three PCA seeds per configuration, keeping each dataset fixed.
python3 validation/run_large_population.py --compare --replicates 3
```

R and the development version of pedigreeTools must be available. The runner
loads randPedGenPCA sources directly from the neighbouring repository. If the
package is elsewhere:

```sh
python3 validation/run_large_population.py --package-root /path/randPedGenPCA
```

Other configurable arguments:

```sh
python3 validation/run_large_population.py \
  --sizes 10000 25000 50000 \
  --markers 500 \
  --generations 10 \
  --geno-generations 2 \
  --geno-fraction 0.5 \
  --rho 5 --q 4 --sketch 20 \
  --trace-vectors 50 \
  --timeout 1800
```

`--trace-vectors 0` disables Hutch++. `--centering true` evaluates only CHC;
`--centering false` evaluates only H. `--seed` controls simulation and PCA.
Population sizes are processed in ascending order, with duplicates removed.

## Results and progress recovery

Each invocation creates a new `output/large_<timestamp>/` directory containing:

- `config.json`: configuration and SHA-256 checksums of the sources used.
- `results.csv`: accumulated results, updated atomically after each process.
  The `stage` column distinguishes preparation, PCA and diagnostics.
- One `.rds` per PCA, containing scores, eigenvectors, eigenvalues and variance
  proportions.
- Individual CSV files and logs with details or errors.

Temporary datasets are removed after all configurations for a size finish.
Use `--keep-datasets` to retain them. PCA results are saved before diagnostics,
so a later diagnostic failure does not remove the fit. `--output` selects a
new output directory; an existing directory is rejected to avoid mixing runs.
Interrupted runs do not resume automatically. To process remaining sizes,
start another invocation with `--sizes`.

The default time limit is **per process**, including startup: 1,800 seconds.
Python enforces it by terminating the process group, even if R is executing
compiled code. An error or timeout is recorded and stops the size ladder.
Ctrl+C terminates the active process and preserves rows already written.
The copy on Eddie was separately modified during the session to accept
`--timeout 0`; the local original runner still requires a positive timeout.
Scheduler wall-clock limits remain independent of this internal limit.

## Measurements and interpretation

**Scalability:** filter `stage == "pca"` and compare `pca_elapsed_sec` and
`peak_process_rss_mib` against n. PCA time includes pedigree factorization and
Hutch++ when enabled, but excludes simulation, data loading, saving and later
diagnostics. `process_wall_sec` includes the entire process. Libraries that
respect the configured thread environment variables use one thread to improve
comparability.

Memory is the peak RSS of the entire process in MiB, measured using
`wait4/ru_maxrss` on macOS or Linux. It includes R, dependencies, input data,
buffers and output saving. It is not incremental PCA-only memory or the sum of
all allocations. Preparation and diagnostics use other processes so they do
not contaminate the PCA process high-water mark. Avoid inferring linear memory
scaling merely because the algorithm does not construct the full H matrix.

**Numerical robustness:** inspect `symmetry_error` before symmetrizing B,
`orth_fro_error` (the norm of U^T U - I), `score_scaling_max_error`, centered
score means and `warnings`. The score-scaling diagnostic currently assumes
`u * sqrt(d)`; check that convention against the scientific package version
before interpreting this column. No numerical error does not guarantee that
the sketch accurately captures the leading components.

**Approximate accuracy without exact PCA:** diagnostics calculate

```text
max_j ||K u_j - d_j u_j||_2 / max(||K u_j||_2, epsilon_mach)
```

This is an eigenpair residual, not error against an exact eigenvalue. A small
residual also does not prove that all leading eigenvalues were found. H is not
constructed for this check: an additional block of implicit products and a
separate pedigree factorization are used. Their cost is recorded under
`stage == "diagnose"`.

With `--compare`, the refined fit reports the largest relative eigenvalue
change against the base fit and the normalized distance between their subspaces:

```text
sqrt(max(0, 1 - ||U_base^T U_refined||F^2 / rho))
```

This comparison is invariant to signs and rotations within the subspace.
Stability across settings provides evidence of convergence, but the refined
fit is not an exact reference. Repeating several seeds helps assess random
variation. Variance percentages also include uncertainty from the Hutch++
trace estimate.

## Pipeline smoke test

```sh
python3 validation/run_large_population.py \
  --sizes 240 480 --markers 40 --generations 6 \
  --trace-vectors 10 --compare --timeout 60
```

This pipeline completed successfully with preparation, eight PCA fits and
eight diagnostic processes. External termination on timeout was also checked.

## Historical initial local run with 10,000 individuals

The initial local run used `--sizes 10000 --compare` and otherwise default
settings: one seed, 1,000 genotyped animals, 1,000 markers and Hutch++.
These are historical Mac measurements, not the later Eddie timings.

| Operator | q | l | PCA time (s) | Peak process RSS (MiB) | Maximum relative residual |
|---|---|---|---|---|---|
| H | 4 | 20 | 11.083 | 527.8 | 0.0981 |
| H | 6 | 40 | 14.488 | 654.5 | 0.0474 |
| CHC | 4 | 20 | 11.214 | 545.7 | 0.0983 |
| CHC | 6 | 40 | 14.331 | 591.8 | 0.0485 |

The residual is not eigenvalue error against exact PCA. No algorithm failures
or warnings occurred in that run. The initial check demonstrated execution at
10,000; sizes 50,000 and 100,000 were not executed as part of that local check.

Results: [`results.csv`](../archive/cleanup_20261102/output/large_20261016_163845_617262/results.csv).

Moving from q=4,l=20 to q=6,l=40 gave a maximum relative eigenvalue change of
4.46% for H and 4.00% for CHC. The normalized subspace distance was about 0.455
in both cases. Residuals decreased, but these changes do not establish
convergence of the base configuration for this pedigree. To study accuracy,
increase the computational budget and repeat seeds. Do not interpret these
changes as errors against exact eigenvalues or as a numerical stability failure.
