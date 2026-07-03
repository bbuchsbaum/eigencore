# Convert an object to an eigencore operator.

Convert an object to an eigencore operator.

## Usage

``` r
as_operator(x, ...)
```

## Arguments

- x:

  Object to convert.

- ...:

  Additional arguments passed to methods.

## Value

An `eigencore_operator` representation of `x`.

## Examples

``` r
op <- as_operator(diag(c(3, 2, 1)))
op$dim
#> [1] 3 3
op$structure$kind
#> [1] "hermitian"
```
