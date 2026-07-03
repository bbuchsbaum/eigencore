# Extract a result certificate.

Extract a result certificate.

## Usage

``` r
certificate(x, ...)
```

## Arguments

- x:

  An eigencore result object.

- ...:

  Reserved for future methods.

## Value

The `eigencore_certificate` object stored on `x`, or `NULL` if the
result does not carry a certificate field.

## Examples

``` r
fit <- eig_partial(diag(c(3, 2, 1)), k = 1, target = largest())
cert <- certificate(fit)
cert$passed
#> [1] TRUE
cert$max_residual
#> [1] 0
```
