# Expose a certified PSD action as an eigencore operator

Expose a certified PSD action as an eigencore operator

## Usage

``` r
psd_operator(
  x,
  action = c("form", "sqrt", "inverse_sqrt", "pseudoinverse", "image_projector",
    "null_projector")
)
```

## Arguments

- x:

  A certified PSD factor.

- action:

  Certified square action.

## Value

An `eigencore_operator`.
