# RSpectra-compatible eigen shim.

RSpectra-compatible eigen shim.

## Usage

``` r
eigs(A, k, which = "LM", opts = list(), ...)
```

## Arguments

- A:

  Matrix or eigencore operator.

- k:

  Number of eigenpairs to compute.

- which:

  RSpectra-style target selector.

- opts:

  Compatibility options list; currently accepted for API compatibility
  and not interpreted directly.

- ...:

  Additional arguments passed to
  [`eig_partial()`](https://bbuchsbaum.github.io/eigencore/reference/eig_partial.md).

## Value

A list compatible with
[`RSpectra::eigs()`](https://rdrr.io/pkg/RSpectra/man/eigs.html),
including `values`, `vectors`, convergence counts, operation counts,
certificate diagnostics, and left/right vector fields when available.

## Examples

``` r
A <- diag(c(5, 4, 3, 2, 1))
A[1, 2] <- 0.5
res <- eigs(A, k = 2, which = "LM")
res$values
#> [1] 5 4
```
