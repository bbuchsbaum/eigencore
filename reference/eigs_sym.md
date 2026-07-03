# RSpectra-compatible symmetric eigen shim.

RSpectra-compatible symmetric eigen shim.

## Usage

``` r
eigs_sym(A, k, which = "LA", opts = list(), ...)
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
  [`solve.eigencore_eigen_problem()`](https://bbuchsbaum.github.io/eigencore/reference/solve.eigencore_eigen_problem.md).

## Value

A list compatible with
[`RSpectra::eigs_sym()`](https://rdrr.io/pkg/RSpectra/man/eigs.html),
including `values`, `vectors`, convergence counts, operation counts,
certificate diagnostics, and eigencore diagnostics.

## Examples

``` r
A <- diag(c(5, 4, 3, 2, 1))
res <- eigs_sym(A, k = 2, which = "LA")
res$values
#> [1] 5 4
```
