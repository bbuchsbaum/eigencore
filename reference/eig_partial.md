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
  allow_dense_fallback = c("auto", "never", "always")
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
