# Compute a PSD Gram or cross-Gram block

Compute a PSD Gram or cross-Gram block

## Usage

``` r
psd_gram(x, X, Y = NULL)
```

## Arguments

- x:

  A certified PSD factor.

- X:

  A finite real-double vector or matrix.

- Y:

  Optional second block.

## Value

`crossprod(X, K Y)` as a dense matrix.
