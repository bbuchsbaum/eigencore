# Method Selection And Workflows

Date: 2026-06-07

This guide gives the current user-facing path through eigencore. It covers the
main workflows required for the V2 CRAN release: partial eigen, partial SVD,
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
| `nearest(sigma)` | You want values nearest a shift; `auto()` routes supported Hermitian problems through `shift_invert(sigma)`. |
| `both_ends(k_low, k_high)` | You want both algebraic ends. |

Method guidance:

| Request | Current planner behavior |
|---|---|
| `method = auto()` | Uses promoted native Hermitian paths where the benchmark-backed planner permits them. Symmetric tridiagonal sparse/diagonal Hermitian sources use the native selected tridiagonal solver; general sparse block Lanczos remains diagnostic opt-in, so non-tridiagonal sparse `auto()` stays scalar by default. |
| `method = lanczos()` | Requests Lanczos explicitly; supported dense/CSC Hermitian operators use native paths, unsupported targets use honest reference labels. |
| Dense near-full request | May use a native dense LAPACK fallback rather than an iterative solver. |

For sparse memory tests, pass `allow_dense_fallback = "never"`.

Base complex dense matrices use native dense complex LAPACK labels for
Hermitian eigen, general eigen, and SVD calls, with exact residual/backward-error
certificates. Base complex dense operators also use native `zgemm` block apply,
with `native_operator_kernel = "dense_complex_zgemm"` metadata. This is a base
dense complex promotion, not a complex sparse `Matrix` or matrix-free solver
promotion. Complex-valued `Matrix`/sparse inputs remain rejected so imaginary
components are not silently discarded, and complex matrix-free eigen/SVD
callbacks fail with explicit future-scope messages.

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
real/imaginary/magnitude targets use a native Arnoldi cycle with native refined
Ritz extraction, right-residual certification, a planner-wired restart budget,
and best-attempt retention across restart attempts. Real matrix-free callback
operators with supported targets keep the native callback Arnoldi cycle with
native projected Ritz extraction and the same certification/restart boundary.
This is the scoped compatibility surface; adjoint-capable rows also expose left
vectors with left-residual and biorthogonality diagnostics. Full Krylov-Schur or
harmonic/interior extraction, matrix-free refined extraction, and native
complex-valued input operators remain future scope. Base complex dense
nonsymmetric matrices use the native dense complex general LAPACK label with
exact right-residual certification. These `which = "LR"`, `"SR"`, `"LI"`, and
`"SI"` compatibility targets also apply to real-valued input matrices that may
produce complex eigenpairs; they are not a promoted complex sparse/operator
promise.

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
| `method = auto()` | Uses the certified tiny sparse Gram special case where bounded and certified, native Golub-Kahan prototypes for supported explicit operators, the native callback boundary for real matrix-free operators with adjoints, or dense LAPACK fallback for dense/small cases. |
| `method = golub_kahan()` | Requests the Golub-Kahan SVD path; explicit dense/CSC operators can use a native prototype, and real matrix-free operators with adjoints use the native callback boundary. |
| `method = randomized()` | Dense double and sparse CSC QR-normalized requests for largest singular values use native randomized controllers with exact residual certification diagnostics and native q=0 early stop. Matrix-free, LU/none-normalized, and other unpromoted randomized regimes use the reference-control randomized SVD prototype with native sketch/projection/projected-core kernels where available. |

Current limitation: production thick-restart SVD for general sparse and
matrix-free workloads remains open. Do not treat `native prototype Golub-Kahan`,
the native matrix-free Golub-Kahan callback boundary, or unpromoted
`reference randomized SVD prototype` randomized regimes as completed V2 CRAN
production paths.
Base complex dense SVD uses the native dense complex LAPACK SVD label with an
exact two-sided certificate; native complex sparse/matrix-free SVD remains
future scope.

Complex operator contract:

- complex operators use `dtype = "complex"`; base complex dense operators use
  native `zgemm` block apply, while native `ScalarType::C128` remains reserved
  for broader sparse/matrix-free C++ operator kernels;
- `apply_adjoint()` must implement the conjugate transpose action
  `A^* X = Conj(t(A)) X`;
- native dense complex decomposition labels and native base dense complex
  block-apply metadata are allowed for base complex dense matrices, while
  complex sparse/matrix-free operator labels still require their own
  block-apply and certificate paths to preserve imaginary components;
- explicit dense complex sources use exact Frobenius-scale certificates with
  conjugate-transpose Gram matrices;
- complex matrix-free certificates can pass only when they carry non-estimated
  norm provenance, such as explicit Frobenius norm metadata; estimated scales
  keep `passed = FALSE`.

Smallest and interior SVD targets have an explicit policy boundary:

| Target | Current behavior |
|---|---|
| `smallest()` / `smallest_magnitude()` | Dense inputs use exact dense fallback certificates. Sparse CSC inputs use the native certified smallest Golub-Kahan label, and matrix-free callbacks use the native certified smallest callback label only when explicit norm metadata keeps the certificate scale non-estimated. |
| `nearest(sigma)` on dense input | Uses exact dense fallback selection and two-sided residual certification. |
| `nearest(sigma)` on sparse CSC input | Uses the native full-subspace Golub-Kahan boundary without densifying the original sparse operator. |
| `nearest(sigma)` on matrix-free callbacks | Uses the native full-subspace callback boundary only when explicit norm metadata keeps the certificate scale non-estimated; otherwise it fails loudly rather than returning an estimated-scale certificate. |

Planner controls expose `svd_target_family`, `svd_target_boundary`, and
`svd_target_certificate_policy` so downstream code can distinguish exact dense
fallback, native certified smallest selection, native full-subspace interior
selection, reference/prototype interior selection, and unsupported native
interior routes.

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

### Graph Laplacians And Fiedler Modes

For graph spectral work, keep symmetric problems on the Hermitian/generalized
SPD route:

```r
# Standard Laplacian modes: L x = lambda x
fit <- eigencore::eig_partial(
  L,
  k = 2,
  target = eigencore::smallest(),
  allow_dense_fallback = "never"
)

# Normalized or weighted modes: L x = lambda D x
norm_fit <- eigencore::eig_partial(
  L,
  B = D,
  k = 2,
  target = eigencore::smallest(),
  allow_dense_fallback = "never"
)
```

Connected graph Laplacians have a zero eigenvalue with the constant vector as
the nullspace. If you want the Fiedler mode, do not treat `smallest()` as a
single "skip the zero mode" command. Either request enough pairs and inspect
the second one, or deflate the known nullspace:

```r
nullspace <- matrix(1, nrow = nrow(L), ncol = 1)
fiedler <- eigencore::eig_partial(
  L,
  B = D,
  k = 1,
  target = eigencore::smallest(),
  method = eigencore::lobpcg(maxit = 600, constraints = nullspace),
  allow_dense_fallback = "never"
)

stopifnot(eigencore::certificate(fiedler)$passed)
```

Near-null Fiedler modes can need more LOBPCG iterations than ordinary edge
queries. When a factorized path is available, `shift_invert(sigma)` near zero
is often the more reliable way to certify the mode:

```r
near_zero <- eigencore::eig_partial(
  L,
  B = D,
  k = 2,
  target = eigencore::smallest(),
  method = eigencore::shift_invert(1e-4),
  allow_dense_fallback = "never"
)

stopifnot(eigencore::certificate(near_zero)$passed)
```

Do not force graph Laplacians through the nonsymmetric sparse general-pencil
route. If you truly need a general right-hand pencil, build it explicitly with
`eigen_problem(A, metric = B, structure = general())` and call `solve()`; the
convenience wrapper `eig_partial()` does not take a `structure = general()`
argument.

Current limitation: dense generalized `auto()` uses the native dense LAPACK
fallback. The promoted iterative generalized surface is sparse
shifted-tridiagonal generalized SPD LOBPCG for largest/smallest targets, with
explicit SPD matrix-free metrics, constraints, generalized-Lanczos reference
rows, and adversarial B cases covered by the strict generalized gate.
Generalized Lanczos now distinguishes sparse tridiagonal CSC metrics, which use
a native Thomas metric solve inside the reference-labelled Lanczos refinement,
from general sparse CSC metrics, which retain `Matrix::Cholesky` reference
provenance. Block B-orthogonal Lanczos is covered only inside the native
dense/diagonal transformed generalized-Lanczos boundary; sparse-CSC block
promotion is not claimed.

For nonsymmetric sparse pencils with diagonal nonsingular `B`, build an
explicit general problem with `eigen_problem(A, metric = B, structure =
general())` and solve it. That route uses native Arnoldi on the sparse
transformed operator `B^{-1} A` and certifies the original generalized residual
`A v - lambda B v`. This is not sparse QZ, it is not the graph-Laplacian
Fiedler route, and smallest near-null spectra should be treated as ungated
unless a smallest-specific general-pencil benchmark proves the case.

## Shift-Invert

Use `target = nearest(sigma)` or `method = shift_invert(sigma)` when you need
interior or smallest-near-shift eigenvalues:

```r
fit <- eigencore::eig_partial(
  A,
  k = 4,
  target = eigencore::nearest(0),
  method = eigencore::shift_invert(sigma = 0)
)
```

With `method = auto()`, `target = nearest(sigma)` is planned as an implicit
`shift_invert(sigma)` transform for supported dense and sparse factorized
Hermitian regimes. Matrix-free or otherwise unfactorized requests fail loudly
unless the user supplies an explicit solve.

Current state:

| Problem class | Current path |
|---|---|
| Dense standard Hermitian | Native factorized Lanczos; certifies original residuals. |
| Dense generalized SPD | Native factorized Lanczos; certifies `A x - lambda B x`. |
| Sparse standard | Honest reference-labelled sparse LU path. |
| Sparse or diagonal-metric generalized SPD | Native only for tridiagonal `A` with diagonal `B`; otherwise honest reference-labelled or rejected depending on support. |
| User-supplied solve | Reference boundary path. |

Every shift-invert result exposes a `shift_invert_factorization_contract_v1`
record in `fit$transform$factorization_cache$contract`. Native labels require
`provider = "eigencore_native_factorization"` and
`promotion_status = "promoted_native"`. General sparse rows use
`provider = "Matrix::lu_reference_factorization"` with sparse pivot diagnostics
and no dense rcond. User-supplied solve functions use
`provider = "user_supplied_solve"` with `external_cache = TRUE`. Native general
sparse LU and native ownership of user-supplied solve functions are explicit
PRD non-goals unless a future PRD reopens them.

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
built-in native kernels. Matrix-free SVD with an adjoint now uses a native
callback-cycle label and gate, but the broader sparse/matrix-free SVD
performance moat remains V3 work.

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
The V2 CRAN contract requires those cases to remain explicitly labelled.

## Method Selection Summary

| Workload | Start with | Use explicitly when |
|---|---|---|
| Standard Hermitian eigen | `eig_partial(A, k)` | `lanczos()` to force Lanczos controls. |
| RSpectra-compatible Hermitian | `eigs_sym(A, k, which)` | Migrating existing RSpectra call sites. |
| General eigen | `eigs(A, k, which)` or `eigen_problem(..., general())` | You need `LR/SR/LI/SI` target semantics for real-valued inputs; dense, sparse CSC, and real matrix-free callback operators use native Arnoldi compatibility for supported targets. |
| Partial SVD | `svd_partial(X, rank)` | `golub_kahan()` for explicit GK, `randomized()` for the reference randomized prototype. |
| Generalized SPD | `eig_partial(A, k, B = B, method = lobpcg())` | You need original generalized residual certification. |
| Sparse general pencil | `eig_partial(A, k, B = B)` with sparse general `A` and diagonal nonsingular `B` | You need a few generalized eigenpairs from `A x = lambda B x` without sparse QZ or dense fallback. |
| Interior/shifted eigen | `eig_partial(..., method = shift_invert(sigma))` | Dense standard/generalized, sparse diagonal/symmetric-tridiagonal standard, and tridiagonal generalized cases with diagonal `B` use native shift-invert; general sparse remains reference-labelled with cache provenance. |
| Structured transforms | `as_operator()`, `center()`, `scale_cols()`, `crossprod_operator()` | You need planner-visible operator provenance without silent sparse densification. |

See `docs/rspectra-migration.md`, `docs/known-limitations.md`,
`docs/contribution-methods-artifact.md`, and `docs/v1-readiness-audit.md` for
migration notes, current limitations, the contribution-facing methods summary,
and the release gate checklist.
