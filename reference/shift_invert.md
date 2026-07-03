# Shift-invert method descriptor.

Shift-invert method descriptor.

## Usage

``` r
shift_invert(sigma, solve = NULL, factorization = NULL)
```

## Arguments

- sigma:

  Shift value `sigma`.

- solve:

  Optional user-supplied solve operator for `(A - sigma B)`.

- factorization:

  Optional precomputed factorization handle.

## Value

An `eigencore_method` descriptor selecting shift-invert.
