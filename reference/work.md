# Extract typed logical work diagnostics.

Returns shared work counters and phase timings without redefining the
route-specific legacy `matvecs` field. Operator, adjoint, metric,
preconditioner, and certification work are reported separately. Unknown
legacy counters are `NA` and make `complete` false.

## Usage

``` r
work(x, ...)
```

## Arguments

- x:

  An eigencore result or an existing `eigencore_work` record.

- ...:

  Reserved for future methods.

## Value

An `eigencore_work` schema-version-1 record.

## Details

`schema_version` identifies the frozen record layout and `complete` is
true only when every logical counter is known. Each `*_block_calls`
field counts logical applications, while its paired `*_columns` field
counts the vector columns processed by those applications. Certification
has distinct forward-operator and adjoint ledgers. `iterations` and
`restarts` retain route-level progress, the four `*_seconds` fields
report measured phase or total elapsed time when available, and
`legacy_matvecs` preserves the result's historical route-specific
`matvecs` value without treating it as a common unit across solver
families.

## Examples

``` r
fit <- eig_partial(diag(c(5, 3, 1)), k = 2)
work(fit)$operator_columns
#> [1] 0
work(fit)$legacy_matvecs
#> [1] 0
```
