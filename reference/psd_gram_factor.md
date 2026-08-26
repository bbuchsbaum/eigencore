# Construct a supplied Gram PSD factor

Dense base and `dgeMatrix` inputs receive a complete compact-SVD
certificate and canonical spectral actions. A `dgCMatrix` remains sparse
and exposes only the form and Gram actions justified by `K = L L^T` or
`K = L^T L`; sparse storage alone does not certify rank or a principal
square root.

## Usage

``` r
psd_gram_factor(
  x,
  orientation = c("columns", "rows"),
  policy = psd_policy(),
  source = c("snapshot", "live")
)
```

## Arguments

- x:

  A finite real-double base matrix, `dgeMatrix`, or `dgCMatrix`.

- orientation:

  Whether columns define `K = x x^T` or rows define `K = x^T x`.

- policy:

  A PSD policy.

- source:

  Source ownership.

## Value

A certified PSD factor when the representation is admitted.
