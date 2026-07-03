# Compute a dense generalized singular value decomposition

`generalized_svd()` is eigencore's dense GSVD compatibility surface for
matrix pairs with the same number of columns. The current native path
uses LAPACK `dggsvd` through R's native LAPACK interface, so it requires
a linked LAPACK that still provides that deprecated routine.
Sparse/operator inputs are not silently densified, and complex GSVD
remains explicit future scope until a complex GSVD driver is available
through the eigencore native layer.

## Usage

``` r
generalized_svd(A, B, tol = 1e-08, ...)
```

## Arguments

- A:

  Base dense real matrix with `m` rows and `n` columns.

- B:

  Base dense real matrix with `p` rows and `n` columns.

- tol:

  Finite non-negative reconstruction and orthogonality certification
  tolerance.

- ...:

  Reserved for future options.

## Value

An `eigencore_gsvd_result` with fields `alpha`, `beta`, `values`,
`classification`, `U`, `V`, `Q`, `D1`, `D2`, `R`, `zero_R`, `A_factor`,
`B_factor`, `k`, `l`, `rank`, `method`, `plan`, and `certificate`.

## Examples

``` r
A <- matrix(c(1, 2, 3, 3, 2, 1), nrow = 2, byrow = TRUE)
B <- matrix(1:9, nrow = 3)
fit <- generalized_svd(A, B)
alpha_beta(fit)$values
#> [1] 0.2086137 2.6092781        NA
certificate(fit)$passed
#> [1] TRUE
```
