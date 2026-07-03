# Shifted diagonal preconditioner.

Shifted diagonal preconditioner.

## Usage

``` r
shifted_diagonal_preconditioner(A, shift = 0)
```

## Arguments

- A:

  Diagonal matrix or numeric vector containing the diagonal entries.

- shift:

  Non-negative diagonal shift added before inversion.

## Value

A typed native-backed preconditioner function mapping residual blocks to
preconditioned blocks.
