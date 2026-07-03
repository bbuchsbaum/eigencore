# Golub-Kahan bidiagonalization method descriptor.

Golub-Kahan bidiagonalization method descriptor.

## Usage

``` r
golub_kahan(max_subspace = NULL, reorthogonalize = TRUE)
```

## Arguments

- max_subspace:

  Optional maximum Krylov subspace size.

- reorthogonalize:

  Whether to apply full two-sided reorthogonalization. `FALSE` selects
  the native one-sided small-side policy where supported, with final
  acceptance still controlled by the exact two-sided certificate.

## Value

An `eigencore_method` descriptor selecting Golub-Kahan
bidiagonalization.
