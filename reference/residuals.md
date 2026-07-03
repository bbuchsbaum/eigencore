# Extract residual diagnostics.

Methods for the
[`stats::residuals()`](https://rdrr.io/r/stats/residuals.html) generic:
return the per-pair (or per-triplet) residual norms stored in an
eigencore result or certificate.

## Usage

``` r
# S3 method for class 'eigencore_eigen_result'
residuals(object, ...)

# S3 method for class 'eigencore_svd_result'
residuals(object, ...)

# S3 method for class 'eigencore_certificate'
residuals(object, ...)
```

## Arguments

- object:

  An eigencore result or certificate object.

- ...:

  Reserved for future methods.
