# Shifted tridiagonal preconditioner.

Shifted tridiagonal preconditioner.

## Usage

``` r
shifted_tridiagonal_preconditioner(A, shift = 0)
```

## Arguments

- A:

  Real symmetric tridiagonal matrix.

- shift:

  Non-negative diagonal shift added to the tridiagonal system.

## Value

A typed preconditioner function mapping residual blocks to
preconditioned blocks.
