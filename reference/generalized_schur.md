# Compute a dense generalized Schur decomposition

`generalized_schur()` is eigencore's dense QZ surface for general matrix
pencils `A x = lambda B x`. It computes the generalized Schur pair `S`,
`T` and, when requested, left/right Schur vectors `Q`, `Z` such that
`A = Q S Z*` and `B = Q T Z*` for complex inputs, with transpose
replacing conjugate-transpose for real inputs. Sparse and operator
inputs are not silently densified.

## Usage

``` r
generalized_schur(A, B, sort = NULL, vectors = TRUE, ...)
```

## Arguments

- A:

  Base dense square matrix.

- B:

  Base dense square matrix with the same dimension as `A`.

- sort:

  Optional LAPACK sorting class. Use `NULL` or `"none"` for no sorting,
  `"finite"` to move finite generalized eigenvalues first, `"infinite"`
  to move beta-zero nonzero-alpha eigenvalues first. Custom predicates
  and undefined alpha-zero/beta-zero sorting are not part of the public
  contract.

- vectors:

  Whether to compute Schur vectors `Q` and `Z`.

- ...:

  Reserved for future options.

## Value

A classed generalized Schur result with fields `S`, `T`, `Q`, `Z`,
`alpha`, `beta`, `values`, `classification`, `sdim`, `method`, `plan`,
and `certificate`.

## Examples

``` r
A <- matrix(c(0, -1, 1, 0), 2, 2)
B <- diag(2)
qz <- generalized_schur(A, B)
values(qz)
#> [1] 0+1i 0-1i

pencil <- generalized_schur(diag(c(2, 3, 0)), diag(c(1, 0, 0)),
                            sort = "infinite")
pencil$classification
#> [1] "infinite"  "finite"    "undefined"
```
