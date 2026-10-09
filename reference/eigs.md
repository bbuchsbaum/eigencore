# RSpectra-compatible eigen shim.

Mirrors
[`RSpectra::eigs()`](https://rdrr.io/pkg/RSpectra/man/eigs.html).
`which` uses ARPACK codes (`"LM"`, `"SM"`, `"LR"`, `"SR"`, `"LI"`,
`"SI"`). Unlike ARPACK, `"LI"`/`"SI"` rank by the signed imaginary part
(see
[`largest_imaginary()`](https://bbuchsbaum.github.io/eigencore/reference/largest_imaginary.md)),
not by its magnitude. A non-`NULL` real `sigma` with `which = "LM"`
requests the eigenvalues nearest `sigma`, computed by shift-invert
Krylov-Schur Arnoldi on a factorised `A - sigma I` (dense QR or sparse
LU).

## Usage

``` r
eigs(
  A,
  k,
  which = "LM",
  sigma = NULL,
  opts = list(),
  ...,
  n = NULL,
  args = NULL,
  left = FALSE
)
```

## Arguments

- A:

  Matrix, eigencore operator, or a function `f(x, args)` returning
  `A %*% x` (then `n` is required).

- k:

  Number of eigenpairs to compute.

- which:

  RSpectra-style target selector.

- sigma:

  Optional real shift; eigenvalues nearest `sigma` are returned.

- opts:

  RSpectra options list. `tol`, `ncv` (Krylov subspace size, passed as
  `auto(max_subspace = ncv)` so the planner still chooses the route),
  `maxitr` (passed as the iteration limit `maxit`) and `retvec` are
  honoured; other keys raise a warning.

- ...:

  Additional arguments passed to
  [`eig_partial()`](https://bbuchsbaum.github.io/eigencore/reference/eig_partial.md).

- n:

  Dimension of the operator when `A` is a function.

- args:

  Extra argument passed to a function `A`.

- left:

  Whether to compute and certify left eigenvectors
  (`left_vectors = "auto"` in
  [`eig_partial()`](https://bbuchsbaum.github.io/eigencore/reference/eig_partial.md)).
  `FALSE` (default) skips the adjoint solve and its certificate
  entirely.

## Value

A list compatible with
[`RSpectra::eigs()`](https://rdrr.io/pkg/RSpectra/man/eigs.html),
including `values`, `vectors`, convergence counts, operation counts,
certificate diagnostics, and, with `left = TRUE`, `left_vectors`,
`right_vectors`, `left_certificate` and `biorthogonality`.

## Details

Like [`RSpectra::eigs()`](https://rdrr.io/pkg/RSpectra/man/eigs.html),
only right eigenvectors are computed by default; `left = TRUE` also
computes and certifies left eigenvectors and their biorthogonality with
the right ones (the
[`eig_partial()`](https://bbuchsbaum.github.io/eigencore/reference/eig_partial.md)
default).

## Examples

``` r
A <- diag(c(5, 4, 3, 2, 1))
A[1, 2] <- 0.5
res <- eigs(A, k = 2, which = "LM")
res$values
#> [1] 5 4
```
