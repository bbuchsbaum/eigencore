# Shift-invert method descriptor.

Shift-invert method descriptor.

## Usage

``` r
shift_invert(sigma, solve = NULL, factorization = NULL, max_subspace = NULL)
```

## Arguments

- sigma:

  Shift value `sigma`.

- solve:

  Optional user-supplied solve operator for `(A - sigma B)`.

- factorization:

  Optional precomputed factorization handle.

- max_subspace:

  Optional maximum Krylov subspace size for the Lanczos (Hermitian) or
  Krylov-Schur Arnoldi (nonsymmetric) iteration on the inverted
  operator. `NULL` uses the route default.

## Value

An `eigencore_method` descriptor selecting shift-invert.
