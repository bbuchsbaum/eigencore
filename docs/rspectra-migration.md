# RSpectra Migration Notes

Date: 2026-05-17

This guide describes the current eigencore migration surface for code that
already calls `RSpectra::eigs()`, `RSpectra::eigs_sym()`, or
`RSpectra::svds()`. It is a compatibility guide, not a V1 release claim. The
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
real/imaginary/magnitude targets can use a native Arnoldi cycle with native
projected Ritz extraction, right-residual certification, a wired restart
budget, and best-attempt retention. Matrix-free general operators remain
reference-labelled. Treat this as the scoped V1 compatibility surface; fully
restarted matrix-free native Arnoldi is future scope.

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

The production V1 SVD promise is not complete. The tiny sparse Gram special
case is promoted and certified in original coordinates; dense inputs may use a
native dense LAPACK fallback; explicit dense/CSC Golub-Kahan exists as a native
prototype; general sparse and matrix-free thick-restart SVD remain open
performance work. Randomized SVD has a scoped exact-low-rank V1 gate with
native fused sketch/projection kernels, but public randomized control remains
reference-labelled and broader sparse/slow-decay/native-controller promotion
remains future work.

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

Current state: native generalized SPD LOBPCG slices exist for supported built-in
metrics and preconditioners, but the non-quick strict release gate is still red
on convergence and speed. Dense generalized `auto()` uses the native dense
generalized SPD LAPACK fallback until iterative gates pass. Broader generalized
production promotion remains open.

## Shift-Invert

Use `method = shift_invert(sigma)` for explicit shift-invert requests:

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
production claims. User-supplied solve functions are treated as reference
boundary code with external-cache provenance. This is the scoped V1 boundary;
native general sparse LU ownership is future scope.

## Planner Labels To Watch

During migration, treat these labels as production-native or promoted only in
their documented regimes:

| Label pattern | Interpretation |
|---|---|
| `native scalar thick-restart Hermitian Lanczos` | Native Hermitian path. |
| `native tridiagonal Hermitian LAPACK selected eigensolver` | Promoted native Hermitian default for symmetric tridiagonal sparse/diagonal sources with largest/smallest algebraic targets. |
| `native block Hermitian Lanczos (thick restart, locking)` | Block Hermitian path for explicit block requests, dense full-subspace cases, and diagnostic sparse opt-in; not promoted for general sparse `auto()`. |
| `native certified Gram SVD special case` | Bounded Gram SVD special case, certified in original coordinates. |
| `native dense ... fallback` | Native dense LAPACK fallback, not an iterative sparse solver. |
| `native prototype ...` | Native prototype or staging path; do not count as V1 completion by itself. |
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
