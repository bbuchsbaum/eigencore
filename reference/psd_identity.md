# Construct a certified identity PSD factor

Construct a certified identity PSD factor

## Usage

``` r
psd_identity(dim, policy = psd_policy())
```

## Arguments

- dim:

  Non-negative whole dimension.

- policy:

  A PSD policy from
  [`psd_policy()`](https://bbuchsbaum.github.io/eigencore/reference/psd_policy.md).

## Value

An `eigencore_psd_factor`.
