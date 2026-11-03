# Alignment with the specified randomised-PCA algorithm

The reference algorithm supplied by the user returns eigenvalue-scaled scores:

\[
X = \bar U_B\operatorname{diag}(d_\rho),\qquad \bar U_B = Q U_B[:,1:\rho].
\]

The reference output is **U*d**, not U*sqrt(d). The latter is a different score
convention and must not be mixed into this algorithm's validation.

## Mapping to the scientific implementation

| Algorithm step | Current implementation |
|---|---|
| Gaussian sketch Omega with l columns | `randRangeFinderHx()` uses `rnorm(n*l)` |
| Initial Y = H Omega | `oracle(Omega)` |
| q subspace iterations | Repeated `orthHx(Y)` and `oracle(Q)` |
| Final orthogonal basis | `orthHx(Y)` |
| Reduced projection B = Q^T H Q | `HQ <- oracle(Q)` and `crossprod(Q,HQ)` |
| Reduced decomposition | `reducedSVDHx(B)` after symmetry and PSD checks |
| Lift eigenvectors and select rho | `Q %*% e$u[,seq_len(rho)]` |
| Eigenvalue-scaled scores | `sweep(u,2L,d,"*")` |

For a symmetric positive-semidefinite reduced matrix, SVD has the same
nonnegative singular values as the eigendecomposition's eigenvalues, and its
left singular subspaces are eigenspaces. Signs and rotations within equal-value
subspaces need not be identical. The implementation checks symmetry and PSD
within numerical tolerances before using the SVD. This is mathematically
consistent with the specified decomposition on the intended PSD domain, not
a claim of identical intermediate floating-point matrices.

The code also trims numerically dependent QR columns. Consequently, the actual
reduced dimension may be r<=l; centered operators have rank at most n-1. The
supplied algorithm assumes a suitable rank. Our comparisons use rho<=l<=n and
sufficient numerical rank. The current public implementation also accepts some
l>n configurations as an extension; this correction does not change that input
policy or the scientific package. CHC runs apply the same algorithm to the
centered operator rather than to uncentered H. Hutch++ is a separate optional
trace calculation and does not define the score scaling.

## Validation corrections after benchmark baseline 1ba536f

The scientific source files were not changed. Three benchmark-side comparisons
were corrected:

1. `large_population.R`: check fit$x against fit$u*fit$d, recording the convention.
2. `followup_worker.R`: construct exact scores as U*d and reconstruct them when
   reading legacy references before computing coordinate errors.
3. `R/04_bench_comparison.R`: keep prcomp(H) scores at their existing U*d scale
   for PSD H, instead of dividing them by sqrt(d).

Historical raw CSVs are preserved. The report excludes historical aligned-score
errors from valid summaries and retains them separately for audit. Eigenvalue
errors, eigenspace comparisons, residuals, absolute per-PC correlations and PCA
timings are unaffected by this scaling correction. Corrected coordinate errors
can be calculated from saved fits/references without rerunning the PCA itself.

Run the focused regression check locally:

```sh
Rscript tests/check_score_convention.R ../RandPedGenPCA/randPedGenPCA
```

It checks the actual worker outputs for H/CHC, legacy-reference handling and the
prcomp score scale. This is a validation fix, not a performance optimization;
keep it in a separate commit from the saved reference and later speed changes.
