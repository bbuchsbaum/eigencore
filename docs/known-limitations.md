# Known Limitations Before V1

Date: 2026-05-17

This file records current user-facing limitations. It should shrink as V1
milestones close. The authoritative release gate map remains
`docs/v1-readiness-audit.md`; this page is the shorter migration-facing view.

## Not V1 Release-Ready

eigencore is still experimental. It has a credible architecture, certificates,
native operator kernels, and certified Hermitian paths, but it has not
satisfied all PRD release gates. Do not treat current successful local checks
as a release signoff.

## Solver Limitations

| Area | Current limitation |
|---|---|
| Hermitian eigen | Native scalar Hermitian Lanczos remains the default sparse path. Native block Hermitian Lanczos exists for explicit requests and diagnostic sparse opt-in, but fresh installed `path_laplacian:1000` strict evidence is red. The native projected-matrix update hotspot is much smaller now, yet scalar and block-candidate rows still fail speed/parity gates. |
| SVD | Tiny sparse Gram rows certify and pass memory on fresh installed probes, but are speed-red; general sparse and matrix-free production thick-restart SVD remains open. |
| Randomized SVD | Reference prototype with certified refinement exists; native randomized sketch/projection engine and broader benchmark gates remain open. |
| Generalized SPD | Dense generalized `auto()` now uses the native dense LAPACK fallback until iterative gates pass. Native generalized LOBPCG slices exist for sparse/structured and explicit requests, but the full non-quick strict gate is still red. |
| B-orthogonal Lanczos | Explicit generalized-SPD `lanczos()` requests now have an honest reference scalar refinement path for dense, diagonal, and CSC SPD metric solves. Native/block production Lanczos and benchmark promotion remain open. |
| Shift-invert | Dense standard, dense generalized SPD, diagonal standard, sparse symmetric-tridiagonal standard, and sparse/diagonal tridiagonal generalized paths are native; general sparse standard and general sparse/diagonal generalized remain reference-labelled. |
| Nonsymmetric eigen | Dense LAPACK oracle with right-residual certification exists and is benchmarked as a compatibility fallback. Sparse CSC real- and imaginary-target cases now use a native Arnoldi cycle with native projected Ritz extraction, a wired restart budget, and best-attempt retention; matrix-free real-spectrum cases remain reference-labelled. Production-grade fully native restarted Arnoldi is still not implemented. |
| Matrix-free operators | Supported through callback boundaries, but callback-driven paths are not the same as built-in native kernels. |

## Dense Fallback And Sparse Memory

Sparse workloads should be tested with `allow_dense_fallback = "never"` when
memory behavior matters. The planner is designed not to silently densify sparse
solver paths, but explicit dense fallbacks still exist for oracle and
compatibility cases. A dense-oracle method label is a warning that the workload
materialized dense state.

## Certificates

Certificates are central to the package contract and are stronger than a raw
eigenvalue or singular-value return. Current limitations:

- Stochastic norm-scale estimates must not produce an unqualified passed
  certificate.
- Any new solver promotion still needs dense/operator certificate agreement
  tests.
- Generalized and shift-invert certificates must continue to certify the
  original problem, not only the transformed operator.

## Planner Labels

Trust the planner label over assumptions about the input type. Labels containing
`reference`, `prototype`, or `oracle` identify paths that are useful for
correctness and migration but are not final V1 production performance claims.
For shift-invert, dense standard/generalized, sparse diagonal or
symmetric-tridiagonal standard cases, and tridiagonal `A` with diagonal `B`
have native labels; general sparse and general sparse diagonal-metric
generalized cases remain reference-labelled.

## Release-Hardening Gaps

The release-hardening issue remains open until:

- CRAN-like checks are rerun from a fresh tarball.
- Full tests are rerun from the current tree.
- Strict benchmark reports are regenerated for promoted solver families; quick
  certification-only runs are not release signoff.
- Sanitizer or valgrind-style checks are complete. Current evidence is partial:
  UBSan install plus `inst/validation/native-smoke.R` passed locally, while
  ASan is blocked by macOS/R dynamic loading and valgrind is absent.
- Migration and user docs pass a final scope audit after the remaining solver
  blockers close.
- Mote issues for completed gates are closed or handed off with remaining scope
  recorded.
