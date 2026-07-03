# Extract computed values.

Extract computed values.

## Usage

``` r
values(x, ...)
```

## Arguments

- x:

  An eigencore result object.

- ...:

  Reserved for future methods.

## Value

A numeric or complex vector of computed eigenvalues, singular values, or
generalized singular values.

## Examples

``` r
fit <- eig_partial(diag(c(3, 2, 1)), k = 2, target = largest())
values(fit)
#> [1] 3 2
```
