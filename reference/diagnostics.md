# Extract diagnostics.

Extract diagnostics.

## Usage

``` r
diagnostics(x, ...)
```

## Arguments

- x:

  An eigencore result object.

- ...:

  Reserved for future methods.

## Value

A named list of diagnostic fields, including residuals, backward errors,
orthogonality diagnostics, iteration counts, method/plan metadata,
warnings, and any available left-eigenvector diagnostics.

## Examples

``` r
fit <- eig_partial(diag(c(3, 2, 1)), k = 1, target = largest())
d <- diagnostics(fit)
d$nconv
#> [1] 1
d$method
#> [1] "native dense Hermitian LAPACK fallback"
```
