# RSpectra-compatible SVD shim.

RSpectra-compatible SVD shim.

## Usage

``` r
svds(
  A,
  k,
  nu = k,
  nv = k,
  opts = list(),
  ...,
  Atrans = NULL,
  dim = NULL,
  args = NULL
)
```

## Arguments

- A:

  Matrix, eigencore operator, or a function `f(x, args)` returning
  `A %*% x` (then `Atrans` and `dim` are required).

- k:

  Number of singular values to compute.

- nu:

  Number of left singular vectors returned.

- nv:

  Number of right singular vectors returned.

- opts:

  RSpectra options list. `tol`, `center` and `scale` are honoured
  (`center`/`scale` are applied as operators, without densifying); `ncv`
  is passed as `auto(max_subspace = ncv)` (used by the Golub-Kahan and
  implicit-Gram Lanczos routes; the explicit Gram route has no Krylov
  subspace). `maxitr` is accepted but not used, since
  [`svd_partial()`](https://bbuchsbaum.github.io/eigencore/reference/svd_partial.md)
  has no iteration limit yet; other keys raise a warning.

- ...:

  Additional arguments passed to
  [`svd_partial()`](https://bbuchsbaum.github.io/eigencore/reference/svd_partial.md).

- Atrans:

  Function `f(x, args)` returning `t(A) %*% x` when `A` is a function.

- dim:

  Dimensions of `A` when `A` is a function.

- args:

  Extra argument passed to function inputs.

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
