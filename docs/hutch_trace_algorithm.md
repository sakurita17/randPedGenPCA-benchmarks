# Estimating the trace of H with `hutchH()`

This explanation describes the local implementation in
`RandPedGenPCA/randPedGenPCA/R/hutchH.R` and its connection to `rhpca.R`.
The pseudocode documents the existing code; it does not change the estimator.

## What it computes

The trace is the sum of diagonal entries and also the sum of eigenvalues:

\[
\operatorname{tr}(H)=\sum_{i=1}^{n}H_{ii}=\sum_{i=1}^{n}\lambda_i.
\]

`hutchH()` approximates this quantity without constructing H. It uses `getHx()`
to multiply H by blocks of vectors. It implements Hutch++, combining a low-rank
approximation with a Hutchinson estimator applied to the residual.
Reference: [Meyer et al., Hutch++](https://arxiv.org/abs/2010.09649).

Define the operator actually being evaluated:

\[
A=\begin{cases}
H,&\texttt{centredH = FALSE},\\
CHC,&\texttt{centredH = TRUE},
\end{cases}
\qquad C=I-\frac{\mathbf1\mathbf1^T}{n}.
\]

Here A denotes the target operator, not the pedigree relationship matrix.
In the second case the estimator targets `tr(CHC)`, which generally differs
from `tr(H)`. C is not constructed: `centerCols(V)` computes CV by subtracting
each column's mean. The code therefore centers both the input and output of
`getHx()`. The `centred` argument separately controls genotype centering.

## Algorithm: Hutch++ for implicit H

**Inputs:** pedigree or Li, genotype matrix M, H construction parameters,
`centredH`, m = `numVectors`, and f = `sketchFrac`.

**Output:** a scalar estimate of the trace of A.

```text
 1. Obtain Li = getLInv(ped, gamma) if the input is a pedigree.
    Align individual IDs using inputHx().
    Define apply(V) = getHx(Li, M, V, ...).
    If centredH: apply(V) = centerCols(getHx(Li, M, centerCols(V), ...)).

 2. s <- round(m * f / 2)                  # sketchQueries
    t <- m - s                            # hutchQueries
    Check that s >= 1 and t >= 1.

 3. Generate S with dimensions n * s, with independent entries -1 or +1,
    each with probability 1/2.

 4. Y <- apply(S)                         # first getHx call
    Q <- orthonormal basis of Y obtained by pivoted QR
         (retain only numerically independent columns)

 5. Generate R with dimensions n * t, independently of S,
    also with entries -1 or +1.

 6. G <- R - Q(Q^T R)                     # remove the projection onto Q

 7. HQ <- apply(Q)                        # second getHx call
    HG <- apply(G)                        # third getHx call

 8. subspace_trace <- tr(Q^T HQ)
    residual_trace <- tr(G^T HG) / t

 9. Return subspace_trace + residual_trace.
```

R distinguishes the original random vectors from G, which the code overwrites
with projected vectors. This G is not the genomic relationship matrix.
Q is an **orthonormal basis**; the projector is QQ^T.

The resulting formula is:

\[
\boxed{\widehat{\operatorname{tr}(A)}
=\operatorname{tr}(Q^TAQ)
+\frac{1}{t}\sum_{j=1}^{t}g_j^TAg_j,
\qquad g_j=(I-QQ^T)r_j.}
\]

The first term is computed exactly for the sampled subspace Q. It is not
necessarily the exact sum of the first s eigenvalues of A. The second term
estimates the contribution outside that subspace.

## Why it works

For a random vector r with independent entries of +/-1,
\(\mathbb E[rr^T]=I\). Therefore, for a fixed matrix B,

\[
\mathbb E[r^TBr]=\operatorname{tr}(B).
\]

Let P = QQ^T. Since P is an orthogonal projector, the trace decomposes as

\[
\operatorname{tr}(A)
=\operatorname{tr}(Q^TAQ)
+\operatorname{tr}((I-P)A(I-P)).
\]

Conditional on Q, the columns of R remain independent random vectors. The
average of \(g_j^TAg_j\) is an unbiased estimate of the second term. The purpose
of Q is to capture much of A's structure and reduce variability in the residual
estimate. One run returns a random estimate, not an error guarantee or a
confidence interval.

## Vector budget and computational cost

With `numVectors = 50` and `sketchFrac = 2/3`:

```text
s = round(50 * (2/3) / 2) = 17
 t = 50 - 17              = 33
```

If Q has s columns, the three calls to `getHx()` process:

| Operation | Columns |
|---|---:|
| A S | 17 |
| A Q | 17 |
| A G | 33 |
| Total | 67 |

Thus m counts the random vectors generated in S and R; it is not a strict
matrix-vector product budget. The equivalent cost is `2s + t = m + s` when Q
has s columns, or `s + ncol(Q) + t` in general. Likewise, `sketchFrac = 2/3`
does not mean that two thirds of the generated vectors go into S: the code
includes an additional division by 2.

## Connection to `rhpca()`

Both functions construct a subspace Q using implicit H products, but use it
for different outputs:

| Stage | `rhpca()` | `hutchH()` |
|---|---|---|
| Initial vectors | Gaussian, l columns | Rademacher (+/-1), s columns |
| Subspace | QR with q power iterations | QR of A S, without power iterations |
| Reduced matrix | B = Q^T A Q | B = Q^T A Q |
| Calculation on B | Reduced SVD with symmetry and positive-semidefiniteness checks, to approximate principal components | Sum of diagonal entries |
| Complement of Q | Its trace is not estimated by the PCA decomposition | Hutchinson estimation |
| Output | Approximate components and eigenvalues | Approximate total trace |

When `numVectorsTrace` is supplied, `rhpca()` calls `hutchH()` and uses

\[
\texttt{varProps}_i=\frac{\texttt{d}_i}{\widehat{\operatorname{tr}(A)}}.
\]

For a positive-semidefinite target, `d` contains approximate eigenvalues of A.
The denominator includes directions outside the returned `rho` components;
it must not be replaced by `sum(d)`. The implementation constructs the
`hutchH()` and `rhpca()` bases separately: Q and B are not shared. Variance
proportions inherit the random error of trace estimation.

## Usage example

```r
set.seed(123)
trace_est <- randPedGenPCA::hutchH(
  X = scenario$ped,
  M = scenario$M,
  gamma = scenario$gamma,
  centred = "0.5",
  centredH = FALSE,
  epsilon = 0.01,
  numVectors = 50
)
```

For a small example, compare against
`sum(diag(as.matrix(pedigreeTools::getH(...))))`, using exactly the same H
construction arguments. If `centredH = TRUE`, double-center the explicit H
first. `epsilon` regularizes the genomic matrix; it is not a Hutch++ error
tolerance.

## Numerical checks shared with PCA

The implementation shares the implicit operator through `makeHxOperator()`.
Pivoted QR removes dependent columns, including the constant direction in the
centered case. The estimator's Rademacher distribution and original vector
allocation are retained.

`rhpca()` checks symmetry and positive semidefiniteness of the reduced matrix.
The currently inspected source uses reduced SVD and scores `u * d`. Older
benchmark diagnostics and documentation assumed `u * sqrt(d)`; verify the
score convention before interpreting score-scaling errors or coordinate
comparisons. This documentation update does not change either implementation.
If the trace is non-finite, numerically non-positive or below the retained
spectral sum, `rhpca()` emits a warning and returns `varProps = NA`, without
changing the estimate produced by `hutchH()`.
