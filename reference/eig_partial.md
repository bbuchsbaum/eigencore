# Compute a partial eigendecomposition.

Compute a partial eigendecomposition.

## Usage

``` r
eig_partial(
  A,
  k,
  target = largest(),
  B = NULL,
  method = auto(),
  tol = 1e-08,
  maxit = NULL,
  vectors = TRUE,
  seed = NULL,
  certify = TRUE,
  allow_dense_fallback = c("auto", "never", "always"),
  initial_subspace = NULL
)
```

## Arguments

- A:

  Matrix or eigencore operator.

- k:

  Number of eigenpairs to compute.

- target:

  Eigencore eigenvalue target descriptor.

- B:

  Optional metric matrix or operator for generalized problems.

- method:

  Solver method descriptor.

- tol:

  Convergence and certification tolerance.

- maxit:

  Optional iteration limit.

- vectors:

  Whether to compute vectors.

- seed:

  Optional random seed for stochastic solver components.

- certify:

  Whether to compute certification diagnostics.

- allow_dense_fallback:

  Dense fallback policy.

- initial_subspace:

  Optional numeric matrix of starting directions (a warm start).
  Supported on standard real Hermitian Lanczos paths: the native paths
  for explicit dense double or `dgCMatrix` operators, the native
  matrix-free callback path selected by `lanczos(block > 1)`, and the
  scalar matrix-free reference path selected by `lanczos(block = 1)`;
  supplying it on any other planned path (generalized, shift-invert,
  dense fallback) is an error. Pass `method = lanczos()` to guarantee a
  Lanczos route: with the default `method = auto()`, sparse or
  [`nearest()`](https://bbuchsbaum.github.io/eigencore/reference/nearest.md)
  problems may be planned as shift-invert, which does not consume a
  start and will reject the argument. The subspace is only a starting
  hint: projected quantities, residuals, orthogonality, convergence, and
  the certificate are recomputed for the current operator on every
  solve. The columns are orthonormalized at the solver boundary and
  fitted to the method's start block — when the accepted rank exceeds
  the block width the block is a seeded random rotation of the full
  accepted basis, so every supplied direction contributes. Because a
  residual certificate proves eigenpair accuracy but not target
  identity, a fully supplied subspace that is already invariant at `tol`
  is discarded in favor of a cold start; provenance records that guard
  decision. Diagnostics distinguish operator block calls, operator
  columns, and certification columns. `NULL` (the default) preserves the
  cold random start exactly.

## Value

An `eigencore_eigen_result` containing computed values, optional
vectors, certificate diagnostics, method/plan metadata, and convergence
diagnostics.

## Examples

``` r
A <- diag(c(5, 4, 3, 2, 1))
A[1, 2] <- A[2, 1] <- 0.1
fit <- eig_partial(A, k = 2, target = largest())
values(fit)
#> [1] 5.009902 3.990098
certificate(fit)$passed
#> [1] TRUE

# Generalized SPD problem A x = lambda B x
B <- diag(c(2, 1, 1, 1, 1))
gfit <- eig_partial(A, B = B, k = 2, target = smallest())
values(gfit)
#> [1] 1 2
```
