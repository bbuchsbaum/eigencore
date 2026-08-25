# Return retained-memory accounting.

Return retained-memory accounting.

## Usage

``` r
retained_bytes(x)
```

## Arguments

- x:

  An eigencore plan, restart state, result, or `eigencore_memory`
  record.

## Value

The known retained byte count. The full `eigencore_memory` record is
attached as the `memory` attribute; when `complete` is false the value
is a documented lower bound.
