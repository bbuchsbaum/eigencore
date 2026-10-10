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

  Target-completeness check run after a certified Hermitian solve (every
  target except
  [`interval()`](https://bbuchsbaum.github.io/eigencore/reference/interval.md),
  which verifies its own set;
  [`smallest_magnitude()`](https://bbuchsbaum.github.io/eigencore/reference/smallest_magnitude.md)
  is checked as `nearest(0)` and
  [`both_ends()`](https://bbuchsbaum.github.io/eigencore/reference/both_ends.md)
  one end at a time). Unless the check verifies the set,
  `certificate(fit)$passed` is `FALSE` (with `residual_passed = TRUE`).
  `"inertia"` proves completeness deterministically by counting
  eigenvalues with an \\LDL^T\\ factorisation of `A - t B` (see
  [`eigen_count()`](https://bbuchsbaum.github.io/eigencore/reference/eigen_count.md));
  it needs an explicit dense or sparse matrix (a matrix-free operator
  with n at most
  `getOption("eigencore.completeness_materialize_limit", 2000)` is
  materialised by n applies) and reports `"inertia_verified"`,
  `"inertia_failed"` (repaired when possible, otherwise
  `passed = FALSE`) or `"inertia_inconclusive"` (an eigenvalue cluster
  straddles the target edge within the residual bound and the counts
  cannot show that a correct choice was returned; a repeated eigenvalue
  at the edge is verified as a tie, recorded in
  `certificate(fit)$completeness$tie`). `"probe"` runs a short block
  Lanczos process on the operator deflated against the returned
  eigenvectors (from a fixed-seed start; the global random stream is not
  touched) and, if it finds a more-preferred eigenvalue outside the
  returned set (for example a missed copy of a repeated eigenvalue),
  repairs the result with a deflated complement solve; it is
  probabilistic (it can prove a set incomplete but not complete) and is
  the check for large matrix-free operators; for
  [`nearest()`](https://bbuchsbaum.github.io/eigencore/reference/nearest.md)
  targets it runs on the route's own \\(A - \sigma I)^{-1}\\ after
  shift-invert, otherwise on \\(A - \sigma I)^2\\. Generalized problems
  are probed in the standard space \\R^{-T} A R^{-1}\\, \\B = R^T R\\.
  `"auto"` uses the inertia check when the matrix is explicit (or
  materialisable) and its predicted factorisation time is at most
  `max(getOption("eigencore.completeness_inertia_seconds", 0.5), getOption("eigencore.completeness_inertia_ratio", 1) * solve time)`,
  and the probe otherwise. `"none"` skips the check. `NULL` (default)
  uses `getOption("eigencore.target_completeness", "auto")`. The outcome
  is recorded in `certificate(fit)$target_completeness`; see the
  "Certificates" vignette.

## Value

An `eigencore_method` descriptor selecting Lanczos iteration.
