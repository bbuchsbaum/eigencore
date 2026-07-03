# Scale operator columns.

Scale operator columns.

## Usage

``` r
scale_cols(A, weights, name = NULL)
```

## Arguments

- A:

  Operator-like object.

- weights:

  Numeric vector of column weights.

- name:

  Optional label for the scaled operator.

## Value

An `eigencore_operator` representing column-wise right scaling of `A`.
