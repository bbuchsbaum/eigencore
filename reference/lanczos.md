# Hermitian Lanczos method descriptor.

Hermitian Lanczos method descriptor.

## Usage

``` r
lanczos(
  max_subspace = NULL,
  max_restarts = NULL,
  block = 1L,
  reorthogonalize = TRUE
)
```

## Arguments

- max_subspace:

  Optional maximum active Krylov subspace size `m`. Must be at least
  `k + 1`. The native thick-restart path keeps the active basis bounded
  by this value across restart cycles.

- max_restarts:

  Optional non-negative integer giving the maximum number of
  thick-restart cycles allowed before stopping with whatever has
  converged. Default `100L`.

- block:

  Native block size. `1L` selects the scalar thick-restart path; values
  greater than one select the native block Krylov prototype where
  supported.

- reorthogonalize:

  Whether to apply full reorthogonalization. The native path always
  reorthogonalizes (DGKS x2) and ignores this flag; it is preserved for
  the R reference solver's public API.

## Value

An `eigencore_method` descriptor selecting Lanczos iteration.
