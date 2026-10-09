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

## Details

Backward errors use the standard normwise 2-norm definition:
`||A x - lambda B x|| / ((||A||_2 + |lambda| ||B||_2) ||x||)` for
eigenpairs and
`sqrt(||A v - sigma u||^2 + ||A^H u - sigma v||^2) / ||A||_2` for
singular triplets. The norms in the denominator are exact or LOWER
bounds (never estimates), so the reported backward error is never
smaller than the true one and `passed` is sound. Fields describing the
scale:

- `norm_bound_type`:

  `"two_norm_exact"` or `"two_norm_lower_bound"` (eigen certificates
  report `A+B` parts; `B = I` is `"identity_exact"`).

- `norm_source`:

  Where each value came from: `"diagonal"`, `"full_spectrum"`,
  `"metadata"` (an operator's `metadata$two_norm`), `"column_norms"`,
  `"applied_vectors"` (`||A x|| / ||x||` for the certified vectors),
  `"ritz"` (residual-corrected Ritz values), `"frobenius_rank_bound"`,
  or `"lanczos"` (a short deterministic Krylov estimate run only when it
  could change the verdict).

- `norm_values`:

  The values used, named `A` (and `B`).

- `scale_is_estimate`:

  Always `FALSE` for built-in certificates.

## Examples

``` r
fit <- eig_partial(diag(c(3, 2, 1)), k = 1, target = largest())
cert <- certificate(fit)
cert$passed
#> [1] TRUE
cert$max_residual
#> [1] 0
```
