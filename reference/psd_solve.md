# Strictly solve a compatible singular PSD system

Strictly solve a compatible singular PSD system

## Usage

``` r
psd_solve(x, B, tolerance = NULL)
```

## Arguments

- x:

  A certified PSD factor.

- B:

  A right-hand-side vector or matrix.

- tolerance:

  Optional typed RHS compatibility override.

## Value

An `eigencore_psd_solve_result`.
