# Migrating from geigen

``` r

library(eigencore)
library(Matrix)
```

eigencore is not a namespace clone of `geigen`. The migration target is
the problem class: dense SPD/Hermitian metrics, dense general pencils,
dense QZ, or sparse SPD partial solves. Results carry planner labels and
certificates so the numerical path is inspectable.

| geigen call | eigencore call | Status |
|----|----|----|
| `geigen::geigen(A, B, symmetric = TRUE)` | `eig_full(A, B = ...)` or `eig_partial(A, B = ..., k = ...)` | Supported for SPD/Hermitian `B` |
| `geigen::geigen(A, B, symmetric = FALSE)` | `eig_full(A, B = ..., structure = general())` | Supported for dense general pencils |
| `geigen::gqz(A, B)` | `generalized_schur(A, B)` | Supported for dense QZ |
| `geigen::gevalues(qz)` | `values(qz)` and `alpha_beta(qz)` | Supported |
| `geigen::gsvd(A, B)` | `generalized_svd(A, B, ...)` | Real dense only |

## Dense SPD/Hermitian Problems

For a full dense generalized SPD/Hermitian problem, pass `B` directly to
[`eig_full()`](https://bbuchsbaum.github.io/eigencore/reference/eig_full.md).

``` r

A <- diag(c(2, 8, 18))
B <- diag(c(1, 2, 3))

fit <- eig_full(A, B = B)
values(fit)
#> [1] 2 4 6
fit$method
#> [1] "native dense generalized SPD/Hermitian LAPACK full"
certificate(fit)$passed
#> [1] TRUE
```

For partial spectra, use
[`eig_partial()`](https://bbuchsbaum.github.io/eigencore/reference/eig_partial.md).

``` r

part <- eig_partial(A, B = B, k = 2, target = smallest())
values(part)
#> [1] 2 4
part$method
#> [1] "native dense generalized SPD LAPACK fallback"
certificate(part)$passed
#> [1] TRUE
```

## Dense General Pencils

Use `structure = general()` when `B` is indefinite, singular,
nonsymmetric, or when homogeneous `alpha`/`beta` diagnostics matter.

``` r

A <- matrix(c(1, 4, 2, 3), 2, 2)
B <- matrix(c(2, 1, 0, -1), 2, 2)

pencil <- eig_full(A, B = B, structure = general())
values(pencil)
#> [1] -0.75+1.391941i -0.75-1.391941i
alpha_beta(pencil)$classification
#> [1] "finite" "finite"
pencil$method
#> [1] "native dense general pencil LAPACK full"
```

Singular pencils expose finite, infinite, and undefined classifications
rather than hiding beta-zero cases.

``` r

singular <- eig_full(
  diag(c(2, 3, 0)),
  B = diag(c(1, 0, 0)),
  structure = general()
)

values(singular)
#> [1]   2+0i Inf+0i     NA
alpha_beta(singular)$classification
#> [1] "finite"    "infinite"  "undefined"
certificate(singular)$failed_indices
#> [1] 2 3
```

## QZ and gevalues

`geigen::gqz()` maps to
[`generalized_schur()`](https://bbuchsbaum.github.io/eigencore/reference/generalized_schur.md).
Use
[`values()`](https://bbuchsbaum.github.io/eigencore/reference/values.md)
for `alpha / beta`, and
[`alpha_beta()`](https://bbuchsbaum.github.io/eigencore/reference/alpha_beta.md)
when the homogeneous coordinates or classification are needed.

``` r

qz <- generalized_schur(A, B)
values(qz)
#> [1] -0.75+1.391941i -0.75-1.391941i
coords <- alpha_beta(qz)
coords$alpha
#> [1] -1.1726039+2.176261i -0.9594032-1.780577i
coords$beta
#> [1] 1.563472+0i 1.279204+0i
coords$classification
#> [1] "finite" "finite"
qz$method
#> [1] "native dense generalized Schur QZ LAPACK full"
```

The current QZ sorting contract is eigencore-specific: use no sorting,
finite-first sorting, or infinite-first sorting. Other geigen sort modes
are not silently accepted.

``` r

beta_zero <- generalized_schur(
  diag(c(2, 3, 0)),
  diag(c(1, 0, 0)),
  sort = "infinite"
)
alpha_beta(beta_zero)$classification
#> [1] "infinite"  "finite"    "undefined"
```

## Sparse Users

Do not migrate sparse problems by relying on implicit densification.
Full dense surfaces require base dense matrices and reject
sparse/operator inputs.

``` r

eig_full(Matrix::Diagonal(3), B = Matrix::Diagonal(3))
#> Error:
#> ! A must be a base dense matrix for eig_full(); sparse/operator full decompositions are not silently densified
```

Sparse SPD/Hermitian generalized partial problems use
[`eig_partial()`](https://bbuchsbaum.github.io/eigencore/reference/eig_partial.md)
with an explicit no-densification policy when that boundary matters.

``` r

A_sparse <- Matrix::Diagonal(x = c(1, 4, 9, 16, 25, 36))
B_sparse <- Matrix::Diagonal(x = c(1, 2, 3, 4, 5, 6))

sparse_fit <- eig_partial(
  A_sparse,
  B = B_sparse,
  k = 3,
  target = smallest(),
  method = lanczos(max_subspace = 6),
  allow_dense_fallback = "never"
)

values(sparse_fit)
#> [1] 1 2 3
sparse_fit$method
#> [1] "native transformed generalized SPD B-orthogonal Lanczos"
sparse_fit$restart$metric_solve
#> [1] "diagonal scaling similarity transform for B"
certificate(sparse_fit)$passed
#> [1] TRUE
```

## GSVD

Use `generalized_svd(A, B)` for real dense GSVD workloads. The result
exposes homogeneous generalized singular values through
[`alpha_beta()`](https://bbuchsbaum.github.io/eigencore/reference/alpha_beta.md),
orthogonal factors `U`, `V`, and `Q`, and reconstructable `D1`, `D2`,
`R`, and `zero_R` factors.

``` r

A <- matrix(c(1, 2, 3, 3, 2, 1), nrow = 2, byrow = TRUE)
B <- matrix(1:9, nrow = 3)

gfit <- generalized_svd(A, B)
alpha_beta(gfit)$values
#> [1] 0.2086137 2.6092781        NA
certificate(gfit)$passed
#> [1] TRUE
```

Complex and sparse GSVD inputs remain future scope and fail explicitly
rather than routing through a dense or external-package fallback.
