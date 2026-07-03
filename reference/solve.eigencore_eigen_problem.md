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

- ...:

  Reserved for future solver options.

## Value

An `eigencore_eigen_result`.
