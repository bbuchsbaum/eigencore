# Compute a partial singular-value decomposition.

Compute a partial singular-value decomposition.

## Usage

``` r
svd_partial(
  A,
  rank,
  target = largest(),
  method = auto(),
  tol = 1e-08,
  vectors = c("both", "left", "right", "none"),
  seed = NULL,
  certify = TRUE,
  allow_dense_fallback = c("auto", "never", "always")
)
```

## Arguments

- A:

  Matrix or eigencore operator.

- rank:

  Number of singular values to compute.

- target:

  Eigencore singular-value target descriptor.

- method:

  Solver method descriptor.

- tol:

  Convergence and certification tolerance.

- vectors:

  Which singular-vector sides to compute.

- seed:

  Optional random seed for stochastic solver components.

- certify:

  Whether to compute certification diagnostics.

- allow_dense_fallback:

  Dense fallback policy.

## Value

An `eigencore_svd_result` containing singular values, optional left and
right singular vectors, certificate diagnostics, method/plan metadata,
and convergence diagnostics.

## Examples

``` r
set.seed(1)
X <- matrix(rnorm(60), 10, 6)
fit <- svd_partial(X, rank = 3)
values(fit)
#> [1] 4.728358 3.042304 2.415933
certificate(fit)$passed
#> [1] TRUE
```
