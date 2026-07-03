# Scale operator rows.

Scale operator rows.

## Usage

``` r
scale_rows(A, weights, name = NULL)
```

## Arguments

- A:

  Operator-like object.

- weights:

  Numeric vector of row weights.

- name:

  Optional label for the scaled operator.

## Value

An `eigencore_operator` representing row-wise left scaling of `A`.
