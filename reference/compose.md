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
