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

## Details

For a built-in explicit `A`, `A^* A` is formed once and wrapped as a
native explicit operator (`metadata$materialized_crossprod`) only when
that is cheap and memory-safe: always when it has at most 65,536
entries; for a larger dense `A` only when `A` has no more columns than
rows and forming it costs at most 2^30 multiply-adds; for a larger
sparse `A` only when a structural bound on its nonzeros stays within
four times `nnz(A)` and below a quarter of its entries (a sparse
cross-product is never densified). Otherwise the result is lazy and
applies `A` then its adjoint.
