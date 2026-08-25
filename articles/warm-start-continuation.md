# Reuse an eigenspace across a parameter sweep

When you solve a sequence of nearby Hermitian eigenproblems, the
eigenvectors from one solve are often good starting directions for the
next. Passing those directions through `initial_subspace` can reduce
operator work while preserving the ordinary result contract: eigencore
still recomputes the Ritz values, residuals, orthogonality checks, and
certificate for the current operator.

This vignette follows the smallest eigenpairs of `A - rho B` as `rho`
changes. At the end you will have a certified result at every step,
provenance showing whether the supplied start was used, and
operator-column counts for comparing the warm and cold solves.

``` r

library(eigencore)
```

## What problem are we solving?

Here `B` is a perturbation term in the standard eigenproblem, not the
optional generalized-eigenproblem metric argument of
[`eig_partial()`](https://bbuchsbaum.github.io/eigencore/reference/eig_partial.md).
Both `A` and `B` are symmetric, so every member of the family
`A - rho B` is Hermitian.

``` r

n <- 180L
k <- 5L
A <- diag(2, n)
A[row(A) == col(A) + 1L] <- -1
A[row(A) + 1L == col(A)] <- -1
x <- seq_len(n) / (n + 1)
B <- diag(0.5 + x + 0.2 * sin(4 * pi * x))
```

Use an explicit Lanczos method for a warm start. This prevents
[`auto()`](https://bbuchsbaum.github.io/eigencore/reference/auto.md)
from routing a structured or nearest-target problem to shift-invert,
which does not consume `initial_subspace`.

``` r

method <- lanczos(
  block = k,
  max_subspace = 8L * k,
  max_restarts = 200L
)
first <- eig_partial(
  A, k = k, target = smallest(), method = method,
  tol = 1e-8, seed = 2026
)
certificate(first)$passed
#> [1] TRUE
```

## How do you reuse the previous eigenspace?

Pass the returned eigenvectors to the next solve. The values alone are
not enough: `initial_subspace` expects directions with one row per
operator dimension.

``` r

rho <- 0.02
warm <- eig_partial(
  A - rho * B, k = k, target = smallest(), method = method,
  tol = 1e-8, seed = 2026, initial_subspace = vectors(first)
)
certificate(warm)$passed
#> [1] TRUE
```

The same pattern extends to a sweep. This example also runs a cold solve
at each step so the work and answers can be compared on equal terms.

``` r

rhos <- c(0.02, 0.04, 0.06)
start <- vectors(first)
rows <- vector("list", length(rhos))

for (i in seq_along(rhos)) {
  op <- A - rhos[[i]] * B
  cold <- eig_partial(
    op, k, smallest(), method = method, tol = 1e-8, seed = 2026
  )
  warm <- eig_partial(
    op, k, smallest(), method = method, tol = 1e-8, seed = 2026,
    initial_subspace = start
  )
  rows[[i]] <- data.frame(
    rho = rhos[[i]],
    cold_columns = cold$operator_columns,
    warm_columns = warm$operator_columns,
    max_value_difference = max(abs(sort(values(cold)) - sort(values(warm))))
  )
  start <- vectors(warm)
}
comparison <- do.call(rbind, rows)
```

| rho | cold operator columns | warm operator columns | max \|cold value - warm value\| |
|---:|---:|---:|:---|
| 0.02 | 1160 | 405 | 1.75e-13 |
| 0.04 | 1040 | 365 | 1.37e-13 |
| 0.06 | 960 | 565 | 1.70e-13 |

On this reproducible family, every continuation step uses fewer operator
columns than its cold counterpart and agrees on the requested values
well within the solve tolerance. These counts include the work used to
guard and certify the supplied subspace; `operator_block_calls` is a
different metric because one block call may apply the operator to
several columns.

Do not assume the same speedup for every sweep. When the target
eigenspace changes abruptly or the supplied directions have little
overlap with it, a warm solve can cost about as much as a cold one. The
reproducible benchmark `inst/benchmarks/bench-warm-start-continuation.R`
includes both high-overlap continuation and a deliberate overlap-loss
jump.

## How can you tell what happened to the start?

[`diagnostics()`](https://bbuchsbaum.github.io/eigencore/reference/diagnostics.md)
exposes the start source and the boundary-processing counts.

``` r

d <- diagnostics(warm)
data.frame(
  start_source = d$start_source,
  supplied = d$initial_subspace$supplied,
  accepted = d$initial_subspace$accepted,
  rejected = d$initial_subspace$rejected,
  augmented = d$initial_subspace$augmented,
  rank = d$initial_subspace$rank
)
#>    start_source supplied accepted rejected augmented rank
#> 1 user_supplied        5        5        0         0    5
```

The supplied matrix must be numeric, finite, and have `n` rows.
Eigencore orthonormalizes its columns and detects numerical rank at the
solver boundary. If too few independent directions survive, the start
block is augmented; if more directions survive than the method’s block
width, a seeded rotation lets all accepted directions contribute to the
fitted block. The `seed` therefore makes augmentation and compression
reproducible as well as controlling a cold random start.

Passing `initial_subspace = NULL` is exactly the cold-start contract. It
leaves the previous random-start sequence and result unchanged.

## Which solver plans accept a warm start?

The supported surface is deliberately narrow:

| Problem and plan | Warm-start status |
|----|----|
| Standard real Hermitian Lanczos, explicit dense double matrix | Native and supported |
| Standard real Hermitian Lanczos, `dgCMatrix` | Native and supported |
| Standard real Hermitian Lanczos, matrix-free operator, `block > 1` | Native callback path and supported |
| Standard real Hermitian Lanczos, matrix-free operator, `block = 1` | Reference-labelled path and supported |
| Generalized problem using the `B` argument | Rejected |
| Shift-invert or [`nearest()`](https://bbuchsbaum.github.io/eigencore/reference/nearest.md) plan | Rejected |
| Dense fallback, nonsymmetric, or non-Lanczos plan | Rejected |

Use `method = lanczos()` when the start is required. With
`method = auto()`, inspect
[`plan_solver()`](https://bbuchsbaum.github.io/eigencore/reference/plan_solver.md)
first: sparse algebraic-edge targets and
[`nearest()`](https://bbuchsbaum.github.io/eigencore/reference/nearest.md)
may select a factorized shift-invert plan. Supplying a start to an
unsupported plan raises an error instead of silently ignoring it.

For native block paths, `lanczos(check_stride = N)` with a positive
integer checks convergence every `N` block iterations and can stop a
strong warm start before a full cold-sized sweep completes. The default
`check_stride = 0L` preserves full-sweep behavior. Mid-sweep checks do
not add operator applications; enable them when workload-specific
measurements show a benefit.

## What does the fresh certificate prove?

A starting subspace is only a hint. Each solve certifies the returned
eigenpairs against the current `A - rho B`; it never reuses the previous
certificate.

There is one important distinction: a small residual proves that a
returned pair is an eigenpair, but not that it belongs to the requested
end of the spectrum. A fully supplied subspace that is already invariant
at the requested tolerance is therefore discarded in favor of a cold
start, and `start_source` records the invariant-guard decision. For
reliable continuation, supply directions that overlap the requested
target eigenspace and always check `certificate(fit)$passed`.

## Where should you go next?

- [`?eig_partial`](https://bbuchsbaum.github.io/eigencore/reference/eig_partial.md)
  documents the complete `initial_subspace` contract.
- [`vignette("certificates")`](https://bbuchsbaum.github.io/eigencore/articles/certificates.md)
  explains residuals, backward error, and orthogonality checks.
- `diagnostics(fit)` reports start provenance and separates operator
  block calls, operator columns, and certification columns.
- Run `Rscript inst/benchmarks/bench-warm-start-continuation.R --strict`
  for the non-quick continuation evidence on your own machine.
