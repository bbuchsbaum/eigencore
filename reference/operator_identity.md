# Return operator identity and revision provenance.

Return operator identity and revision provenance.

## Usage

``` r
operator_identity(x)
```

## Arguments

- x:

  An eigencore operator, problem, executable plan, restart state, or
  result.

## Value

An `eigencore_operator_identity` for an operator, or a named list with
entries `A` and, when present, `B` for larger workflow objects.
