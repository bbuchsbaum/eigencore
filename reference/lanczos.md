# Hermitian Lanczos method descriptor.

Hermitian Lanczos method descriptor.

## Usage

``` r
lanczos(
  max_subspace = NULL,
  max_restarts = NULL,
  block = 1L,
  check_stride = 0L,
  reorthogonalize = TRUE,
  completeness = NULL
)
```

## Arguments

- max_subspace:

  Optional maximum active Krylov subspace size `m` (the ARPACK `ncv`).
  Must be at least `k + 1` (`k + block` for block Lanczos). The native
  thick-restart path keeps the active basis bounded by this value across
  restart cycles; unrestarted reference paths build at most this many
  Lanczos vectors. This is the only subspace-size control: the solve's
  `maxit` argument is an iteration (restart) limit.

- max_restarts:

  Optional non-negative integer giving the maximum number of
  thick-restart cycles allowed before stopping with whatever has
  converged. Default `100L`. Equivalent to the solve's `maxit` argument
  on thick-restart routes; supplying both with different values is an
  error.

- block:

  Native block size. `1L` selects the scalar path; for a matrix-free
  operator this remains the reference Hermitian Lanczos boundary. Values
  greater than one select the native block Krylov path where supported,
  including real Hermitian matrix-free callbacks.

- check_stride:

  Native block thick-restart mid-sweep convergence stride. `0L`
  (default) evaluates convergence once per full sweep (legacy). A
  positive `N` evaluates convergence every `N` block iterations within a
  sweep, letting a warm start that converges after a few blocks stop
  early instead of paying a full cold-sized sweep. Mid-sweep checks
  never consume extra operator applications and never change results at
  `check_stride = 0L`. This control applies only to native block paths.

- reorthogonalize:

  Whether to apply full reorthogonalization. The native path always
  reorthogonalizes (DGKS x2) and ignores this flag; it is preserved for
  the R reference solver's public API.

- completeness:

  Target-completeness check run after a certified Hermitian solve with a
  [`largest()`](https://bbuchsbaum.github.io/eigencore/reference/largest.md),
  [`smallest()`](https://bbuchsbaum.github.io/eigencore/reference/smallest.md),
  [`largest_magnitude()`](https://bbuchsbaum.github.io/eigencore/reference/largest_magnitude.md)
  (or, for the inertia check,
  [`nearest()`](https://bbuchsbaum.github.io/eigencore/reference/nearest.md))
  target. `"inertia"` proves completeness deterministically by counting
  eigenvalues with an \\LDL^T\\ factorisation of `A - t B` (see
  [`eigen_count()`](https://bbuchsbaum.github.io/eigencore/reference/eigen_count.md));
  it needs an explicit dense or sparse matrix and reports
  `"inertia_verified"`, `"inertia_failed"` (repaired when possible,
  otherwise `passed = FALSE`) or `"inertia_inconclusive"` (an eigenvalue
  cluster straddles the target edge within the residual bound).
  `"probe"` runs a short block Lanczos process on the operator deflated
  against the returned eigenvectors (from a fixed-seed start; the global
  random stream is not touched) and, if it finds a more-preferred
  eigenvalue outside the returned set (for example a missed copy of a
  repeated eigenvalue), repairs the result with a deflated complement
  solve; it is probabilistic (it can prove a set incomplete but not
  complete) and is the check for matrix-free operators. `"auto"` uses
  the inertia check for edge targets (not
  [`nearest()`](https://bbuchsbaum.github.io/eigencore/reference/nearest.md))
  when the matrix is explicit and its predicted factorisation time is at
  most
  `max(getOption("eigencore.completeness_inertia_seconds", 0.5), getOption("eigencore.completeness_inertia_ratio", 1) * solve time)`,
  and the probe otherwise. `"none"` skips the check. `NULL` (default)
  uses `getOption("eigencore.target_completeness", "auto")`. The outcome
  is recorded in `certificate(fit)$target_completeness`; see the
  "Certificates" vignette.

## Value

An `eigencore_method` descriptor selecting Lanczos iteration.
