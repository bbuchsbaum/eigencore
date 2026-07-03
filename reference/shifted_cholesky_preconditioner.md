# Shifted Cholesky preconditioner.

Shifted Cholesky preconditioner.

## Usage

``` r
shifted_cholesky_preconditioner(A, shift = 0)
```

## Arguments

- A:

  Symmetric positive semidefinite or positive definite matrix.

- shift:

  Non-negative diagonal shift added before factorization.

## Value

A typed preconditioner function mapping residual blocks to
preconditioned blocks.
