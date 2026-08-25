# Cold-start provenance record.

Returned when no `initial_subspace` is supplied so downstream
diagnostics carry a uniform provenance schema for cold and warm solves
alike.

## Usage

``` r
warm_start_cold_provenance()
```

## Value

A named list with `start_source = "cold"` and zeroed
supplied/accepted/rejected/augmented/rank counts.
