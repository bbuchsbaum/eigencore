# Method Selection And Workflows

Date: 2026-05-17

This guide gives the current user-facing path through eigencore. It covers the
main workflows required for V1 documentation: partial eigen, partial SVD,
generalized SPD, shift-invert, operator algebra, certificates/diagnostics, and
method selection. It is intentionally conservative: when a path is still a
prototype, reference, or dense oracle, the text says so.

## Working Rule

Start with the high-level function, then inspect the plan and certificate:

```r
fit <- eigencore::eig_partial(A, k = 5, target = eigencore::largest())

fit$method
fit$plan
fit$certificate
eigencore::diagnostics(fit)
```

The method string tells you what actually ran. The certificate tells you
whether the returned values and vectors passed residual, backward-error, and
orthogonality checks.

## Standard Hermitian Eigenproblems

Use `eig_partial()` when the input is symmetric/Hermitian and the eigenproblem
is standard:

```r
fit <- eigencore::eig_partial(
  A,
  k = 10,
  target = eigencore::largest()
)
```

Useful targets:

| Target | Use when |
|---|---|
| `largest()` | You want largest algebraic eigenvalues of a Hermitian problem. |
| `smallest()` | You want smallest algebraic eigenvalues. |
| `largest_magnitude()` | You want eigenvalues with largest modulus. |
| `smallest_magnitude()` | You want eigenvalues with smallest modulus. |
| `nearest(sigma)` | You want values nearest a shift and will use shift-invert where supported. |
| `both_ends(k_low, k_high)` | You want both algebraic ends. |

Method guidance:

| Request | Current planner behavior |
|---|---|
| `method = auto()` | Uses promoted native Hermitian paths where the benchmark-backed planner permits them. Sparse block Lanczos is currently diagnostic opt-in after red G1 benchmark evidence, so sparse `auto()` stays scalar by default. |
| `method = lanczos()` | Requests Lanczos explicitly; supported dense/CSC Hermitian operators use native paths, unsupported targets use honest reference labels. |
| Dense near-full request | May use a native dense LAPACK fallback rather than an iterative solver. |

For sparse memory tests, pass `allow_dense_fallback = "never"`.

## General Eigenproblems

For RSpectra-shaped general eigen calls, use `eigs()`:

```r
res <- eigencore::eigs(A, k = 4, which = "LR")
```

For explicit eigencore calls, construct a general problem:

```r
problem <- eigencore::eigen_problem(
  A,
  structure = eigencore::general(),
  target = eigencore::largest_real()
)
fit <- solve(problem, k = 4)
```

Current limitation: dense and sparse CSC nonsymmetric problems with supported
real/imaginary/magnitude targets can use a native Arnoldi cycle with native
projected Ritz extraction, right-residual certification, a planner-wired
restart budget, and best-attempt retention across restart attempts. Matrix-free
general problems remain reference-labelled. These are compatibility and
correctness bridges, not the final V1 nonsymmetric production solver.

## Partial SVD

Use `svd_partial()` for rectangular matrices or operators:

```r
fit <- eigencore::svd_partial(
  X,
  rank = 5,
  target = eigencore::largest()
)

fit$d
fit$u
fit$v
fit$certificate
```

Method guidance:

| Request | Current planner behavior |
|---|---|
| `method = auto()` | Uses the certified tiny sparse Gram special case where bounded and certified, native Golub-Kahan prototypes for supported explicit operators, or dense LAPACK fallback for dense/small cases. |
| `method = golub_kahan()` | Requests the Golub-Kahan SVD path; explicit dense/CSC operators can use a native prototype. |
| `method = randomized()` | Uses the reference randomized SVD prototype with residual certification/refinement policy. |

Current limitation: production thick-restart SVD for general sparse and
matrix-free workloads remains open. Do not treat `native prototype Golub-Kahan`
or `reference randomized SVD prototype` as completed V1 production paths.

## Generalized SPD Eigenproblems

Use `B = ...` with `eig_partial()` or build an `eigen_problem()` with
`metric = ...`:

```r
fit <- eigencore::eig_partial(
  A,
  k = 5,
  B = B,
  target = eigencore::largest(),
  method = eigencore::lobpcg(maxit = 200),
  allow_dense_fallback = "never"
)
```

The certificate is computed on the original generalized residual
`A v - lambda B v`, and orthogonality is measured in the `B` inner product
where appropriate.

Current limitation: dense generalized `auto()` uses the native dense LAPACK
fallback until iterative gates pass. Native generalized SPD LOBPCG slices
exist for sparse/structured metrics, explicit SPD matrix-free metrics,
constraints, and typed shifted-diagonal/shifted-tridiagonal preconditioners.
The full non-quick strict release gate is still red, so this family is not
fully V1-promoted.

## Shift-Invert

Use `method = shift_invert(sigma)` when you need interior or smallest-near-shift
eigenvalues:

```r
fit <- eigencore::eig_partial(
  A,
  k = 4,
  target = eigencore::nearest(0),
  method = eigencore::shift_invert(sigma = 0)
)
```

Current state:

| Problem class | Current path |
|---|---|
| Dense standard Hermitian | Native factorized Lanczos; certifies original residuals. |
| Dense generalized SPD | Native factorized Lanczos; certifies `A x - lambda B x`. |
| Sparse standard | Honest reference-labelled sparse LU path. |
| Sparse or diagonal-metric generalized SPD | Native only for tridiagonal `A` with diagonal `B`; otherwise honest reference-labelled or rejected depending on support. |
| User-supplied solve | Reference boundary path. |

Sparse native factorization ownership remains a V1 blocker.

## Operator Algebra

Use operators when you need to preserve structure or avoid materializing an
intermediate matrix:

```r
op <- eigencore::as_operator(A)
op_scaled <- eigencore::scale_cols(op, weights)
op_centered <- eigencore::center(op_scaled, columns = TRUE)
gram <- eigencore::crossprod_operator(op_centered)

fit <- eigencore::eig_partial(
  gram,
  k = 5,
  target = eigencore::largest(),
  allow_dense_fallback = "never"
)
```

Native-backed explicit operator transforms currently cover dense/CSC/diagonal
adjoints, scaling, dense/CSC composition, dense/CSC crossproducts, dense
centering, and CSC centering with low-rank corrections. Matrix-free operators
are supported through callbacks, but callback-driven paths are not the same as
built-in native kernels.

For a custom operator, provide both `apply` and `apply_adjoint` if the workload
may need SVD:

```r
op <- eigencore::linear_operator(
  dim = c(n, p),
  apply = function(X, alpha = 1, beta = 0, Y = NULL) {
    out <- A %*% X
    if (is.null(Y)) alpha * out else alpha * out + beta * Y
  },
  apply_adjoint = function(X, alpha = 1, beta = 0, Y = NULL) {
    out <- crossprod(A, X)
    if (is.null(Y)) alpha * out else alpha * out + beta * Y
  },
  structure = eigencore::general(),
  name = "custom_operator"
)
```

## Certificates And Diagnostics

Every solver result should be read in this order:

```r
fit$method
fit$certificate$passed
fit$certificate$max_backward_error
fit$certificate$norm_bound_type
fit$certificate$scale_is_estimate
eigencore::diagnostics(fit)
```

Interpretation:

| Field | Meaning |
|---|---|
| `passed` | Whether all required certificate checks passed. |
| `max_residual` | Worst raw residual over returned pairs/triplets. |
| `max_backward_error` | Worst residual scaled by the certificate norm bound. |
| `max_orthogonality_loss` | Worst basis orthogonality defect. |
| `norm_bound_type` | Provenance of the scale used in backward-error checks. |
| `scale_is_estimate` | Whether the scale was estimated rather than exact. |

A stochastic or estimated scale is not the same as an exact certificate scale.
V1 requires those cases to remain explicitly labelled.

## Method Selection Summary

| Workload | Start with | Use explicitly when |
|---|---|---|
| Standard Hermitian eigen | `eig_partial(A, k)` | `lanczos()` to force Lanczos controls. |
| RSpectra-compatible Hermitian | `eigs_sym(A, k, which)` | Migrating existing RSpectra call sites. |
| General eigen | `eigs(A, k, which)` or `eigen_problem(..., general())` | You need `LR/SR/LI/SI` target semantics; expect dense oracle until Arnoldi lands. |
| Partial SVD | `svd_partial(X, rank)` | `golub_kahan()` for explicit GK, `randomized()` for the reference randomized prototype. |
| Generalized SPD | `eig_partial(A, k, B = B, method = lobpcg())` | You need original generalized residual certification. |
| Interior/shifted eigen | `eig_partial(..., method = shift_invert(sigma))` | Dense standard/generalized, sparse diagonal/symmetric-tridiagonal standard, and tridiagonal generalized cases with diagonal `B` can use native shift-invert; general sparse remains open. |
| Structured transforms | `as_operator()`, `center()`, `scale_cols()`, `crossprod_operator()` | You need planner-visible operator provenance without silent sparse densification. |

See `docs/rspectra-migration.md`, `docs/known-limitations.md`, and
`docs/v1-readiness-audit.md` for migration notes, current limitations, and the
release gate checklist.
