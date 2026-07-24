# RSpectra Migration Notes

Date: 2026-06-07

This guide describes the current eigencore migration surface for code that
already calls `RSpectra::eigs()`, `RSpectra::eigs_sym()`, or
`RSpectra::svds()`. It is a compatibility guide for the scoped V2 CRAN release,
not a promise that every future solver family is already production-native. The
planner label and certificate on each result remain the source of truth for
which path actually ran.

## First Rule

Inspect the returned method, plan, certificate, and diagnostics during
migration:

```r
fit <- eigencore::eigs_sym(A, k = 5, which = "LA")

fit$values
fit$certificate
fit$diagnostics$method
fit$diagnostics$plan
```

The shims preserve the familiar list shape and add trust metadata. A result
that used a reference or oracle fallback is still usable for correctness
checks, but it is not evidence that the production native path is complete.

## `eigs_sym()`

Use `eigencore::eigs_sym(A, k, which = ...)` for Hermitian or symmetric
standard eigenproblems.

```r
res <- eigencore::eigs_sym(A, k = 5, which = "LA")
```

Returned fields:

| Field | Meaning |
|---|---|
| `values` | Eigenvalues selected by the requested ARPACK-style target. |
| `vectors` | Eigenvectors when requested by the underlying solve. |
| `nconv` | Number of converged eigenpairs returned. |
| `niter` | Iteration count or backend iteration proxy. |
| `nops` | Operator-application count when available. |
| `certificate` | Residual, backward-error, orthogonality, scale, and pass/fail evidence. |
| `diagnostics` | Method, plan, warnings, and backend-specific metadata. |

Current target mapping:

| RSpectra `which` | eigencore target |
|---|---|
| `"LA"` | `largest()` / largest algebraic |
| `"SA"` | `smallest()` / smallest algebraic |
| `"LM"` | `largest_magnitude()` |
| `"SM"` | `smallest_magnitude()` |
| `"BE"` | `both_ends()` |

For Hermitian regimes, eigencore routes symmetric tridiagonal sparse/diagonal
sources to the promoted native selected tridiagonal solver when the target is
largest or smallest algebraic. Other sparse requests use native scalar Lanczos
by default. Block Lanczos is available for explicit block requests, dense
full-subspace cases, and diagnostic sparse opt-in; it is not promoted for
general sparse `auto()`. Unsupported Hermitian targets or operator classes keep
honest reference/oracle labels rather than pretending to be production native
code.

## `eigs()`

Use `eigencore::eigs(A, k, which = ...)` for the RSpectra-shaped general
eigen shim:

```r
res <- eigencore::eigs(A, k = 4, which = "LR")
```

Additional target mappings:

| RSpectra `which` | eigencore target |
|---|---|
| `"LR"` | `largest_real()` |
| `"SR"` | `smallest_real()` |
| `"LI"` | `largest_imaginary()` |
| `"SI"` | `smallest_imaginary()` |

Current limitation: dense and sparse CSC nonsymmetric matrices with supported
real/imaginary/magnitude targets use a native Arnoldi cycle with native refined
Ritz extraction, right-residual certification, a wired restart budget, and
best-attempt retention. Real matrix-free callback operators with supported
targets keep the native callback Arnoldi cycle with native projected Ritz
extraction. Adjoint-capable dense, sparse CSC, and matrix-free rows also return
left vectors with separate left-residual and biorthogonality diagnostics through
the `eigs()` shim. Treat this as the scoped compatibility surface; full
Krylov-Schur or harmonic/interior extraction, matrix-free refined extraction,
and native complex-valued sparse/operator paths remain future scope. Base
complex dense matrices use native dense complex LAPACK labels with exact
certificates.

## `svds()`

Use `eigencore::svds(A, k, nu = k, nv = k)` for the RSpectra-shaped partial
SVD shim:

```r
res <- eigencore::svds(A, k = 5, nu = 5, nv = 5)

res$d
res$u
res$v
res$certificate
```

Returned fields:

| Field | Meaning |
|---|---|
| `d` | Singular values. |
| `u` | Left singular vectors when `nu > 0`. |
| `v` | Right singular vectors when `nv > 0`. |
| `nconv` | Number of converged singular triplets. |
| `niter` | Iteration count or backend iteration proxy. |
| `nops` | Operator-application count when available. |
| `certificate` | Two-sided SVD residual and orthogonality evidence. |
| `diagnostics` | Planner and backend metadata. |

The production V2 CRAN SVD promise is scoped. The tiny sparse Gram special
case is promoted and certified in original coordinates; dense inputs may use a
native dense LAPACK fallback; explicit dense/CSC Golub-Kahan exists as a native
prototype; general sparse and matrix-free thick-restart SVD remain open
performance work. Randomized SVD has a scoped exact-low-rank release gate with
native fused sketch/projection kernels and native projected-core solves, but
public randomized control remains reference-labelled and broader sparse,
slow-decay, and native-controller promotion remains future work.
Base complex dense SVD is accepted through the native dense complex LAPACK SVD
label with an exact two-sided certificate. Complex sparse or matrix-free SVD
remains future scope and fails explicitly rather than entering real-only
kernels.

For smallest or interior singular-value targets, eigencore is deliberately
stricter than a silent normal-equation fallback. Dense smallest and
`nearest(sigma)` requests are exactly certified through dense fallback.
Sparse CSC `nearest(sigma)` requests use a native full-subspace Golub-Kahan
boundary without densifying the original operator; matrix-free
`nearest(sigma)` requests use the same full-subspace callback boundary only
when explicit norm metadata keeps the certificate scale non-estimated.
Diagonal/reference prototype rows remain labelled as prototypes even when
their two-sided certificates pass.

## Dense Fallback Policy

Use `allow_dense_fallback = "never"` when migration tests must prove that a
sparse or operator workload does not densify:

```r
res <- eigencore::eig_partial(
  A_sparse,
  k = 5,
  target = eigencore::smallest(),
  allow_dense_fallback = "never"
)
```

The accepted values are:

| Value | Meaning |
|---|---|
| `"auto"` | Permit bounded dense fallback only when the planner policy allows it. |
| `"never"` | Error instead of materializing a dense fallback. |
| `"always"` | Explicitly opt into dense fallback where no sparse/reference path owns the case. |

Sparse Hermitian and sparse SVD paths should not silently densify. If migration
code depends on sparse memory behavior, keep `"never"` in the regression tests.

## Generalized SPD Problems

The RSpectra shims do not expose the full generalized-SPD surface. Use
`eig_partial(A, k, B = B, method = lobpcg())` or build an `eigen_problem()`
explicitly:

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

Current state: the promoted iterative V2 CRAN surface is sparse
shifted-tridiagonal generalized SPD LOBPCG for largest/smallest targets. Dense
generalized `auto()` uses the native dense generalized SPD LAPACK fallback, and
explicit generalized-SPD `lanczos()` requests use native transformed Lanczos
for dense and diagonal SPD metrics, including block requests inside that
transformed boundary. Sparse CSC SPD metric solves remain an
honest reference B-orthogonal refinement, but tridiagonal sparse CSC metrics now
use an eigencore-owned native Thomas metric solve instead of the reference
Cholesky boundary. General sparse CSC metrics remain reference Cholesky-labelled.
Matrix-free-B generalized LOBPCG, native shifted-diagonal/tridiagonal
preconditioners, constraints, and adversarial B contract rows are covered by the
current generalized strict gate. Arbitrary sparse-CSC metric factorization and
sparse-CSC block generalized Lanczos promotion are not claimed.

## Shift-Invert

Use `target = nearest(sigma)` for automatic shift-invert planning, or
`method = shift_invert(sigma)` for explicit shift-invert requests:

```r
fit <- eigencore::eig_partial(
  A,
  k = 4,
  target = eigencore::nearest(0),
  method = eigencore::shift_invert(sigma = 0)
)
```

Current state: dense standard Hermitian and dense generalized-SPD shift-invert
use native factorized Lanczos hot loops and certify in original coordinates.
Sparse diagonal/symmetric-tridiagonal standard paths and tridiagonal
generalized paths with diagonal `B` are also native. General sparse standard
shift-invert and general sparse or diagonal-metric generalized shift-invert
remain honest reference-labelled paths with cache provenance rather than native
production claims. `auto()` plus `nearest(sigma)` uses an implicit
`shift_invert(sigma)` transform for supported factorized Hermitian regimes;
matrix-free or otherwise unfactorized requests fail loudly unless a solve is
supplied. User-supplied solve functions are treated as reference boundary code
with external-cache provenance. This is the scoped V2 CRAN boundary; native general
sparse LU ownership and native ownership of user-supplied solve functions are
explicit PRD non-goals unless a future PRD reopens them.

## Planner Labels To Watch

During migration, treat these labels as production-native or promoted only in
their documented regimes:

| Label pattern | Interpretation |
|---|---|
| `native scalar thick-restart Hermitian Lanczos` | Native Hermitian path. |
| `native tridiagonal Hermitian LAPACK selected eigensolver` | Promoted native Hermitian default for symmetric tridiagonal sparse/diagonal sources with largest/smallest algebraic targets. |
| `native block Hermitian Lanczos (thick restart, locking)` | Block Hermitian path for explicit block requests, dense full-subspace cases, and diagnostic sparse opt-in; not promoted for general sparse `auto()`. |
| `native block Hermitian Lanczos (matrix-free callback, thick restart, locking)` | Explicit `lanczos(block > 1)` path for real Hermitian matrix-free callbacks; `block = 1` retains the reference Hermitian label. |
| `native certified Gram SVD special case` | Bounded Gram SVD special case, certified in original coordinates. |
| `native dense ... fallback` | Native dense LAPACK fallback, not an iterative sparse solver. |
| `native prototype ...` | Native prototype or staging path; do not count as V2 CRAN completion by itself. |
| `reference ...` | R reference/prototype/oracle path. |
| `dense ... oracle` | Dense oracle fallback for compatibility and certification. |

The exact string matters because the planner label is the compatibility and
audit contract.

## Migration Checklist

1. Replace RSpectra calls with the eigencore shim or the explicit eigencore
   API.
2. Assert the requested `which` code maps to the intended eigencore target.
3. Add regression checks for `nconv`, sorted values, and certificate pass/fail.
4. Check `diagnostics$method` or `diagnostics$plan$method` in tests for any
   workload that depends on a native path.
5. Add `allow_dense_fallback = "never"` to sparse memory-safety tests.
6. Keep reference/oracle-labelled paths out of production performance claims.
7. Revisit the migration after final release-hardening decisions if public
   labels or recommended methods change.

See `docs/known-limitations.md` and `docs/v1-readiness-audit.md` for the
current release blockers. See `docs/method-selection-and-workflows.md` for the
broader eigencore API path beyond RSpectra-compatible shims.
