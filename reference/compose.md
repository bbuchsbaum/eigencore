# Compose two operators.

Compose two operators.

## Usage

``` r
compose(A, B, name = NULL)
```

## Arguments

- A:

  Left operator-like object.

- B:

  Right operator-like object.

- name:

  Optional label for the composed operator.

## Value

An `eigencore_operator` representing the composition `A %*% B`.

## Details

When both operands are built-in explicit matrices (dense double or
`Matrix` storage), the product is formed once and wrapped as a native
explicit operator (`metadata$fused == "compose"`) only when that is
cheap and memory-safe: always for a small product (at most 65,536
entries); for a larger dense product only when it holds no more entries
than the two factors together and costs at most 2^30 multiply-adds; for
a larger sparse product only when a structural bound on its nonzeros
stays within four times the factors' nonzeros and below a quarter of its
entries, so a sparse product is never densified. Otherwise (and for
callback or mixed operands) the result is a lazy composition that
applies `B` then `A` and stores no product.
