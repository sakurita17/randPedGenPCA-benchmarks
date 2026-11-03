# Reproducible optimization protocol

This is the working checklist for Gregor's request: preserve a reference,
profile it, change one implementation detail, verify numerical equivalence,
then measure the change without the profiler.

## 1. Establish the reference version

Freeze all three code components, not just the benchmark driver:

- This benchmark repository.
- The randPedGenPCA package inside the neighbouring RandPedGenPCA repository.
- pedigreeTools, including the version actually installed on Eddie.

As inspected on 3 November 2026, the benchmark repository has only its initial
commit and most current files are untracked. The RandPedGenPCA repository also
has uncommitted package changes. The pedigreeTools working tree was clean at
`3d8914d37fbe56fa82d58fab4dad72fe67799f6e`. These observations are not a baseline
freeze and do not prove which pedigreeTools binary was installed on Eddie.

Before committing, review the relevant changes in each repository. A baseline
commit must reflect the actual package, including intended removals; a commit
that omits untracked R files does not freeze the working version. Avoid adding
installed libraries, compiler objects, personal histories and unrelated work.
Do not use a blanket `git add .` to bypass this review.

After the reviewed code has been committed, record each full commit ID with
`git rev-parse HEAD`. Label the baseline commits, and make improvements on a
separate development branch. Local commits/tags do not publish code; pushing is
a separate action. No commits or tags were created by this protocol document.

Freeze the measurement environment too: R and package versions, RNG kind, source
checksums, BLAS, thread settings, and actual installed pedigreeTools code. The
profile script already records source checksums, session information and selected
dependency function bodies. A session record is documentation, not a package lock
or a complete archived environment.

## 2. Fix the inputs and numerical settings

Use the already saved Eddie datasets:

```text
output/followup_genotyped_20261031_182624_492600/data_g500.rds
output/followup_genotyped_20261031_182624_492600/data_g4000.rds
```

Both have a pedigree of 20,000 individuals, 1,000 markers and nested genotype
subsets. Use H (centredH=FALSE), rho=5, q=4, l=20, epsilon=0.01,
centred="0.5", 50 Hutch++ vectors, PCA seed=123 and one computational thread.
Record hashes of the datasets before running and reuse those exact files.
Store the original fits alongside the immutable reference. Do not regenerate
new inputs between versions when testing an optimization.

The local download currently contains CSVs and settings, not these datasets or
the original fit RDS files. Check and preserve them on Eddie before making code
changes. Do not call an old git revision a verified reference for Eddie until
source/package correspondence has been checked.

## 3. Profile the unchanged reference

Follow [the lineprof guide](profiling_lineprof.md). The prepared job executes R
directly with `lineprof`, interval=0.02 seconds and torture=FALSE. Inspect the
cost of getLInv, Hutch++, repeated getHx calls, implicitMasuda, sparse solves,
orthogonalization and garbage collection.

Repeated reconstruction/solution of the non-genotyped block is a hypothesis,
not an established bottleneck from the small installation smoke test. The
target-size profile needs to be collected and interpreted before optimization.

## 4. Make one change and verify equivalence

Keep data, seed, q, l, centering, trace budget and thread settings unchanged.
Changing q/l or the number of genotyped animals is a different experimental
configuration, not evidence that an implementation was optimized.

Check the following before accepting a performance change:

- Inputs: the same dataset checksum, individual IDs and genotype ordering.
- Operator: compare original and modified H/CHC products on fixed test vectors.
- PCA: compare eigenvalues, residuals and orthogonality.
- Coordinates: allow arbitrary eigenvector signs; use alignment or subspace
  comparisons when signs/rotations can differ.
- Hutch++ and variance proportions: ensure any differences in random-number use
  or arithmetic have been accounted for.

In R, `identical(old, new)` asks for exact equality, including representation
details; `isTRUE(all.equal(old, new, tolerance=...))` allows a stated tolerance.
Choose and record tolerances before inspecting an optimization's outcome. There
is no universal tolerance for all these diagnostics; relative errors, scale,
conditioning and eigenvalue gaps matter. Any intentional change in score scaling
is a separate numerical/API change and must not be hidden as a speed optimization.

## 5. Measure without profiling and keep a change log

Use fresh Rscript processes and identical thread settings for reference and
candidate runs. Initially aim for three repeated measurements with the same PCA
seed and inputs. These repetitions measure runtime variation, not variation in
the simulated population. Keep raw measurements and report the median runtime.
Record peak process RSS using the same convention in both versions; sampled
lineprof allocations are not a substitute for peak RSS.

Keep compute time separate from complete process time, loading, saving and
diagnostics. Do not use lineprof's instrumented time to claim a speedup.

For each accepted change, record:

| Item | Required information |
|---|---|
| Identity | Baseline/candidate commit IDs for the driver and scientific packages |
| Change | One concrete implementation change and its rationale |
| Inputs | Dataset checksum and complete numerical configuration |
| Correctness | Tests, tolerances, largest differences, pass/fail |
| Performance | All repeats, median runtime, peak process RSS and environment |
| Improvement | Speedup = baseline median / candidate median; memory reduction = 100*(1 - candidate RSS / baseline RSS) |

The initial [change log](../reports/optimization/change_log.csv) contains two
historical reference measurements imported from the original Eddie CSV. They
are single measurements and are explicitly marked as not yet tied to a frozen
baseline. No optimization or speedup has been claimed.

## Presenting progress to Gregor

Show a reference version, the profile, one change at a time, numerical validation,
and before/after measurements. An investigation that identifies a bottleneck is
useful progress even before a successful optimization exists.
