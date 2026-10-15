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
  `A + shift * I` must be nonsingular: it is factored once at
  construction, and a singular system is an error naming this
  preconditioner.

## Value

A typed preconditioner function mapping residual blocks to
preconditioned blocks.
