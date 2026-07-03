# Extract eigenvectors.

Extract eigenvectors.

## Usage

``` r
vectors(x, ...)
```

## Arguments

- x:

  An eigencore eigen result object.

- ...:

  Reserved for future methods.

## Value

A matrix whose columns are computed right eigenvectors, or `NULL` when
vectors were not requested or are unavailable.

## Examples

``` r
fit <- eig_partial(diag(c(3, 2, 1)), k = 2, target = largest())
dim(vectors(fit))
#> [1] 3 2
```
