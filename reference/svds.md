# RSpectra-compatible SVD shim.

RSpectra-compatible SVD shim.

## Usage

``` r
svds(A, k, nu = k, nv = k, opts = list(), ...)
```

## Arguments

- A:

  Matrix or eigencore operator.

- k:

  Number of singular values to compute.

- nu:

  Number of left singular vectors requested.

- nv:

  Number of right singular vectors requested.

- opts:

  Compatibility options list; currently accepted for API compatibility
  and not interpreted directly.

- ...:

  Additional arguments passed to
  [`svd_partial()`](https://bbuchsbaum.github.io/eigencore/reference/svd_partial.md).

## Value

A list compatible with
[`RSpectra::svds()`](https://rdrr.io/pkg/RSpectra/man/svds.html),
including `d`, optional `u` and `v`, convergence counts, operation
counts, certificate diagnostics, and eigencore diagnostics.

## Examples

``` r
set.seed(1)
X <- matrix(rnorm(60), 10, 6)
res <- svds(X, k = 2)
res$d
#> [1] 4.728358 3.042304
```
