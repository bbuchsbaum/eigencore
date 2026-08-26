# Apply a certified PSD action to a vector or block

Apply a certified PSD action to a vector or block

## Usage

``` r
psd_apply(
  x,
  X,
  action = c("form", "sqrt", "inverse_sqrt", "pseudoinverse", "image_projector",
    "null_projector")
)
```

## Arguments

- x:

  A certified PSD factor.

- X:

  A finite real-double vector or matrix.

- action:

  Certified action to apply.

## Value

A vector or matrix matching the input block shape.
