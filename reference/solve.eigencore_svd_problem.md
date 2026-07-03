# Solve a planned SVD problem.

S3 method that runs the planned solver for an SVD problem built by
[`svd_problem()`](https://bbuchsbaum.github.io/eigencore/reference/svd_problem.md).
Most users call
[`svd_partial()`](https://bbuchsbaum.github.io/eigencore/reference/svd_partial.md),
which constructs the problem and dispatches here; call
[`solve()`](https://rdrr.io/pkg/Matrix/man/solve-methods.html) directly
when you want to build a problem once and reuse or inspect it. Returns a
certified partial singular-value decomposition.

## Usage

``` r
# S3 method for class 'eigencore_svd_problem'
solve(
  a,
  b,
  rank,
  method = auto(),
  tol = 1e-08,
  vectors = c("both", "left", "right", "none"),
  certify = TRUE,
  allow_dense_fallback = c("auto", "never", "always"),
  ...
)
```

## Arguments

- a:

  Eigencore SVD problem object.

- b:

  Unused second argument reserved by the base
  [`solve()`](https://rdrr.io/pkg/Matrix/man/solve-methods.html)
  generic.

- rank:

  Number of singular values to compute.

- method:

  Solver method descriptor.

- tol:

  Convergence and certification tolerance.

- vectors:

  Which singular-vector sides to compute.

- certify:

  Whether to compute certification diagnostics.

- allow_dense_fallback:

  Dense fallback policy.

- ...:

  Reserved for future solver options.

## Value

An `eigencore_svd_result`.
