# Solve a planned eigenproblem.

S3 method that runs the planned solver for an eigenproblem built by
[`eigen_problem()`](https://bbuchsbaum.github.io/eigencore/reference/eigen_problem.md).
Most users call
[`eig_partial()`](https://bbuchsbaum.github.io/eigencore/reference/eig_partial.md),
which constructs the problem and dispatches here; call
[`solve()`](https://rdrr.io/pkg/Matrix/man/solve-methods.html) directly
when you want to build a problem once and reuse or inspect it. Returns a
certified partial eigendecomposition.

## Usage

``` r
# S3 method for class 'eigencore_eigen_problem'
solve(
  a,
  b,
  k,
  method = auto(),
  tol = 1e-08,
  maxit = NULL,
  vectors = TRUE,
  certify = TRUE,
  allow_dense_fallback = c("auto", "never", "always"),
  initial_subspace = NULL,
  ...
)
```

## Arguments

- a:

  Eigencore eigen problem object.

- b:

  Unused second argument reserved by the base
  [`solve()`](https://rdrr.io/pkg/Matrix/man/solve-methods.html)
  generic.

- k:

  Number of eigenpairs to compute.

- method:

  Solver method descriptor.

- tol:

  Convergence and certification tolerance.

- maxit:

  Optional iteration limit.

- vectors:

  Whether to compute vectors.

- certify:

  Whether to compute certification diagnostics.

- allow_dense_fallback:

  Dense fallback policy.

- initial_subspace:

  Optional numeric matrix of starting directions (a warm start).
  Supported on standard real Hermitian Lanczos paths — native dense
  double / `dgCMatrix`, native matrix-free callbacks for
  `lanczos(block > 1)`, and the scalar matrix-free reference path for
  `lanczos(block = 1)`; supplying it on any other planned path is an
  error. The subspace is only a starting hint, never a source of reused
  convergence: every solve recomputes projected quantities, residuals,
  orthogonality, convergence, and a fresh current-operator certificate.
  An already-invariant supplied subspace is discarded to a cold start
  because residual certification alone cannot establish that it contains
  the requested extremal eigenpairs. `NULL` (the default) preserves the
  cold random start exactly.

- ...:

  Reserved for future solver options.

## Value

An `eigencore_eigen_result`.
