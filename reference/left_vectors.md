# Extract left singular vectors or left eigenvectors.

For SVD results this returns the left singular vectors `U`. For
nonsymmetric and dense general-pencil eigen results this returns left
eigenvectors when the solver computed them (for example, the dense
general-pencil
[`eig_full()`](https://bbuchsbaum.github.io/eigencore/reference/eig_full.md)
path, which computes left generalized eigenvectors satisfying
`w^H A = lambda w^H B`).

## Usage

``` r
left_vectors(x, ...)
```

## Arguments

- x:

  An eigencore SVD or eigen result object.

- ...:

  Reserved for future methods.

## Value

A matrix of left singular vectors or left eigenvectors, or `NULL` when
the result does not contain a left-vector field.
