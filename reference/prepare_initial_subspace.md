# Validate, orthonormalize, and fit a supplied subspace to a start block.

Produces an `n x width` start block for the Lanczos paths together with
provenance counts. Orthonormalization happens here at the solver
boundary for rank detection and explicit reporting; the native kernel
re-orthonormalizes and rank-guards the block again before iterating.

## Usage

``` r
prepare_initial_subspace(
  initial_subspace,
  n,
  width,
  tol = sqrt(.Machine$double.eps)
)
```

## Arguments

- initial_subspace:

  User-supplied numeric matrix (or vector) of start directions.

- n:

  Operator domain dimension.

- width:

  Number of start-block columns the chosen method consumes (the Lanczos
  block size).

- tol:

  Numerical-rank tolerance for boundary orthonormalization.

## Value

A list with `start` (the `n x width` block) and provenance fields
`start_source`, `supplied`, `accepted`, `rejected`, `augmented`, `rank`,
and `compressed`. The internal `accepted_basis` field retains the full
accepted basis for the invariant-subspace safety guard.

## Details

When the accepted numerical rank exceeds `width` (e.g. a k-column
continuation subspace handed to a scalar method), the start block is a
seeded random rotation of the full accepted basis rather than a
truncation: every supplied direction contributes generic weight to every
start column, so a width-1 start still overlaps all k target directions
instead of only the first. Provenance records this as
`compressed = TRUE`.

Seed policy: augmented directions are drawn from the active RNG stream
(controlled by the `seed` argument of
[`eig_partial()`](https://bbuchsbaum.github.io/eigencore/reference/eig_partial.md)
/ [`solve()`](https://rdrr.io/pkg/Matrix/man/solve-methods.html)) and
orthonormalized against the accepted directions, matching how the cold
random start is generated.
