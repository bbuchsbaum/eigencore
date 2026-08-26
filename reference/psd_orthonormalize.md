# PSD-orthonormalize a block modulo the certified null space

PSD-orthonormalize a block modulo the certified null space

## Usage

``` r
psd_orthonormalize(x, X, required_rank = NULL, tolerance = NULL)
```

## Arguments

- x:

  A certified PSD factor.

- X:

  A finite real-double vector or matrix.

- required_rank:

  Optional exact output width.

- tolerance:

  Optional typed Gram-rank tolerance.

## Value

An `eigencore_psd_block_result`.
