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
  initial_subspace = NULL,
  left_vectors = c("auto", "none", "compute")
)
```

## Arguments

- A:

  Matrix or eigencore operator.

- k:

  Number of eigenpairs to compute.

- target:

  Eigencore eigenvalue target descriptor. With
  `method = shift_invert(sigma)` it defaults to `nearest(sigma)`; other
  targets (except
  [`smallest_magnitude()`](https://bbuchsbaum.github.io/eigencore/reference/smallest_magnitude.md)
  with `sigma = 0`) are an error.

- B:

  Optional metric matrix or operator for generalized problems.

- method:

  Solver method descriptor.

- tol:

  Convergence and certification tolerance.

- maxit:

  Optional iteration limit (`NULL` uses each route's default). It bounds
  outer iterations, never the Krylov subspace size: thick-restart cycles
  for (block, generalized, and shift-invert callback) Lanczos,
  Krylov-Schur restarts for nonsymmetric Arnoldi (including shift-invert
  Arnoldi), restart cycles for the reference Arnoldi, LOBPCG iterations,
  and Lanczos steps for the unrestarted reference and native-kernel
  shift-invert Lanczos routes. Dense direct routes ignore it. The
  resolved limit is recorded in `plan$controls$iteration_limit` and
  `plan$controls$iteration_limit_kind`. To set the subspace size (the
  ARPACK `ncv`), use the method descriptor's `max_subspace`
  ([`lanczos()`](https://bbuchsbaum.github.io/eigencore/reference/lanczos.md),
  [`auto()`](https://bbuchsbaum.github.io/eigencore/reference/auto.md),
  [`shift_invert()`](https://bbuchsbaum.github.io/eigencore/reference/shift_invert.md)).
  A `lanczos(max_restarts =)` or `lobpcg(maxit =)` that disagrees with
  `maxit` is an error.

- vectors:

  Whether to compute vectors.

- seed:

  Optional random seed for stochastic solver components. The global
  random number stream is restored on exit.

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

- left_vectors:

  Left-eigenvector policy for nonsymmetric problems. `"auto"` (default)
  computes and certifies left eigenvectors (and the biorthogonality) on
  routes that support them, such as Krylov-Schur Arnoldi, and records
  the reason when they are unavailable; `"none"` skips the left solve
  and its certificate entirely (roughly halving the cost of a
  nonsymmetric Arnoldi solve); `"compute"` is like `"auto"` but an error
  when the route cannot return left eigenvectors. Hermitian problems are
  unaffected.

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
