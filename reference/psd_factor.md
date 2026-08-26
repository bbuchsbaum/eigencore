# Construct a certified real-double PSD factor

A double vector is interpreted as a diagonal. Base double matrices and
the admitted dense/diagonal Matrix classes are validated as complete
sources.

## Usage

``` r
psd_factor(x, policy = psd_policy(), source = c("snapshot", "live"))
```

## Arguments

- x:

  A supported real-double diagonal or square matrix source.

- policy:

  A PSD policy from
  [`psd_policy()`](https://bbuchsbaum.github.io/eigencore/reference/psd_policy.md).

- source:

  Source ownership, currently immutable snapshots for complete
  identity/diagonal/dense paths.

## Value

An `eigencore_psd_factor`.
