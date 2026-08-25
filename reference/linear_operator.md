# Create a block-native linear operator.

Create a block-native linear operator.

## Usage

``` r
linear_operator(
  dim,
  apply,
  apply_adjoint = NULL,
  dtype = "double",
  structure = general(),
  name = NULL,
  metadata = list(),
  operator_id = NULL,
  revision = NULL,
  portable = FALSE
)
```

## Arguments

- dim:

  Integer vector of length two giving row and column dimensions.

- apply:

  Function implementing block multiplication by the operator.

- apply_adjoint:

  Optional function implementing block multiplication by the adjoint
  operator.

- dtype:

  Scalar character type label, currently `"double"` or `"complex"`.

- structure:

  Eigencore structure descriptor such as
  [`general()`](https://bbuchsbaum.github.io/eigencore/reference/general.md)
  or
  [`hermitian()`](https://bbuchsbaum.github.io/eigencore/reference/hermitian.md).

- name:

  Optional operator label used in plans and diagnostics.

- metadata:

  Optional list of implementation metadata.

- operator_id:

  Optional non-empty character identity for a logical callback-operator
  lineage. Supply together with `revision`.

- revision:

  Optional non-empty character revision for the values and behavior
  implemented by the callback. Supply together with `operator_id`.

- portable:

  Whether an explicitly identified callback operator may be restored and
  validated in another R session. Built-in matrix-backed operators
  derive portable deterministic identity automatically.

## Value

An `eigencore_operator` list containing dimensions, apply callbacks,
scalar type, structure metadata, a display name, and implementation
metadata.

## Examples

``` r
A <- diag(c(3, 2, 1))
op <- linear_operator(
  dim = dim(A),
  apply = function(X, alpha = 1, beta = 0, Y = NULL) {
    Z <- alpha * (A %*% X)
    if (is.null(Y) || beta == 0) Z else Z + beta * Y
  },
  structure = hermitian(),
  metadata = list(frobenius_norm = sqrt(sum(A^2)))
)
fit <- eig_partial(op, k = 1, target = largest())
values(fit)
#> [1] 3
```
