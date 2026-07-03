# Create A^\* A as an operator.

Create A^\* A as an operator.

## Usage

``` r
crossprod_operator(A, name = NULL)
```

## Arguments

- A:

  Operator-like object with an adjoint implementation.

- name:

  Optional label for the cross-product operator.

## Value

A Hermitian `eigencore_operator` representing `A^* A`.
