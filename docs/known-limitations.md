# Known V1 Release Boundaries

Date: 2026-06-05

This file records current user-facing limitations and future-scope boundaries
for the scoped V1 surface. The authoritative release gate map remains
`docs/v1-readiness-audit.md`; this page is the shorter migration-facing view.

## Scoped V1 Surface

eigencore is still experimental, and broader solver families remain future
scope, but the promoted V1 surfaces have fresh installed-package evidence and
honest planner labels. Treat this as a scoped V1 surface, not a claim that
every PRD-adjacent future solver is complete.

## Solver Limitations

| Area | Current limitation |
|---|---|
| Hermitian eigen | Symmetric tridiagonal sparse/diagonal Hermitian sources now use a promoted native selected tridiagonal solver and pass the installed `path_laplacian:1000` G1 gate. Native scalar Hermitian Lanczos remains the default for non-tridiagonal sparse paths. Native block Hermitian Lanczos exists for explicit requests and diagnostic sparse opt-in, but is not promoted for general sparse `auto()`. |
| SVD | Promoted non-quick tall/wide sparse rows pass the installed H gate through bounded native Gram special cases. Quick `600x90` probes are diagnostic only, and general sparse/matrix-free production thick-restart SVD remains future scope. |
| Randomized SVD | Scoped V1 release gate is green for the large exact-low-rank planner row with native fused sketch/projection kernels. Public control remains reference-labelled, and quick small-size, sparse-size, slow-decay, and fully native-controller promotion remain future scope. |
| Generalized SPD | The promoted iterative surface is scoped to sparse shifted-tridiagonal generalized SPD LOBPCG for largest/smallest targets. Dense generalized `auto()` uses the native dense LAPACK fallback, and broader generalized preconditioner/block variants remain future scope. |
| B-orthogonal Lanczos | Explicit generalized-SPD `lanczos()` requests have an honest reference scalar refinement path for dense, diagonal, and CSC SPD metric solves, with focused installed contract evidence. Native/block production Lanczos is future scope. |
| Shift-invert | Dense standard, dense generalized SPD, diagonal standard, sparse symmetric-tridiagonal standard, and sparse/diagonal tridiagonal generalized paths are native and pass the scoped V1 installed gate; general sparse standard and general sparse/diagonal generalized remain reference-labelled with cache provenance. Native general sparse LU is future scope. |
| Nonsymmetric eigen | Dense and sparse CSC real/imaginary/magnitude-target cases use a native Arnoldi cycle with native projected Ritz extraction, a wired restart budget, best-attempt retention, and right-residual certification. Matrix-free real-spectrum cases remain reference-labelled. Fully restarted matrix-free native Arnoldi is future scope. |
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

## Release-Hardening Evidence

The scoped release-hardening surface is:

- CRAN-like checks were rerun from a fresh tarball in the final local release
  slice and passed with one expected CRAN incoming/version NOTE.
- Full tests were rerun from the current tree in the final local release slice
  and passed with four expected CRAN skips.
- Strict installed benchmark reports are recorded for promoted solver families;
  quick certification-only runs are not release signoff.
- UBSan install plus `inst/validation/native-smoke.R` passed locally. ASan is
  blocked by macOS/R dynamic loading and valgrind is absent, so those remain
  environment boundaries rather than claimed coverage.
- Migration and user docs have a final scope audit for the current promoted
  surface.
- Mote issues for completed gates are closed or handed off with remaining scope
  recorded during final release handoff.
