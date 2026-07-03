# Center an operator by rows or columns.

Center an operator by rows or columns.

## Usage

``` r
center(
  A,
  rows = FALSE,
  columns = TRUE,
  row_means = NULL,
  col_means = NULL,
  name = NULL
)
```

## Arguments

- A:

  Operator-like object.

- rows:

  Whether to subtract row means.

- columns:

  Whether to subtract column means.

- row_means:

  Optional row means. Required for matrix-free row centering when they
  cannot be derived without densifying.

- col_means:

  Optional column means. Required for matrix-free column centering when
  they cannot be derived without densifying.

- name:

  Optional label for the centered operator.

## Value

An `eigencore_operator` representing the centered linear map.
