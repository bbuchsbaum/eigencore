# Create a scale-aware PSD tolerance

Create a scale-aware PSD tolerance

## Usage

``` r
psd_tolerance(abs = 0, rel = sqrt(.Machine$double.eps))
```

## Arguments

- abs:

  Non-negative finite absolute tolerance.

- rel:

  Non-negative finite relative tolerance.

## Value

An immutable `eigencore_psd_tolerance` record.
