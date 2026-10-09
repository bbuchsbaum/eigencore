# RSpectra-compatible symmetric eigen shim.

Mirrors
[`RSpectra::eigs_sym()`](https://rdrr.io/pkg/RSpectra/man/eigs.html):
only the `lower` (or upper) triangle of a dense or `Matrix` input is
read, values are returned in decreasing order, and `sigma` requests the
eigenvalues nearest the shift.

## Usage

``` r
eigs_sym(
  A,
  k,
  which = "LM",
  sigma = NULL,
  opts = list(),
  lower = TRUE,
  ...,
  n = NULL,
  args = NULL
)
```

## Arguments

- A:

  Matrix, eigencore operator, or a function `f(x, args)` returning
  `A %*% x` (then `n` is required).

- k:

  Number of eigenpairs to compute.

- which:

  RSpectra-style target selector (`"LM"`, `"SM"`, `"LA"`, `"SA"`,
  `"BE"`).

- sigma:

  Optional shift; eigenvalues nearest `sigma` are returned.

- opts:

  RSpectra options list. `tol`, `ncv` (Krylov subspace size, passed as
  `auto(max_subspace = ncv)`, or `lanczos(max_subspace = ncv)` when
  `initvec` forces a Lanczos route), `maxitr` (passed as the iteration
  limit `maxit`), `retvec` and `initvec` are honoured; other keys raise
  a warning.

- lower:

  Whether to read the lower (`TRUE`) or upper (`FALSE`) triangle of a
  matrix input.

- ...:

  Additional arguments passed to
  [`eig_partial()`](https://bbuchsbaum.github.io/eigencore/reference/eig_partial.md).

- n:

  Dimension of the operator when `A` is a function.

- args:

  Extra argument passed to a function `A`.

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
