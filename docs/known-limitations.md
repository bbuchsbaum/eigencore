# Known V2 CRAN Release Boundaries

Date: 2026-06-07

This file records current user-facing limitations and future-scope boundaries
for the V2 CRAN release surface. The authoritative release gate map remains
`docs/v1-readiness-audit.md`, now used as the preserved evidence inventory for
the V2 CRAN candidate; this page is the shorter migration-facing view.
`prd.json` carries the matching machine-readable status map for scoped current
surfaces, future-scope boundaries, and V3 deferrals.

## V2 CRAN Surface

eigencore is still experimental, and broader solver families remain future
scope, but the promoted CRAN release surfaces have installed-package evidence
and honest planner labels. Treat V2 as the CRAN release boundary, not a claim
that every PRD-adjacent future solver is complete.

## Solver Limitations

| Area | Current limitation |
|---|---|
| Hermitian eigen | Symmetric tridiagonal sparse/diagonal Hermitian sources now use a promoted native selected tridiagonal solver and pass the installed `path_laplacian:1000` G1 gate. Explicitly tagged internal 2D path-grid Laplacian operators have a diagnostic separable prototype boundary for smallest eigenpairs; installed 2026-06-07 evidence on `100x100`, `k = 10` certifies in original coordinates and beats RSpectra, but this is not automatic recognition for arbitrary sparse Hermitian matrices. Native scalar Hermitian Lanczos remains the default for non-tridiagonal sparse paths. Native block Hermitian Lanczos exists for explicit requests and diagnostic sparse opt-in, but is not promoted for general sparse `auto()`. |
| SVD | Promoted non-quick tall/wide sparse largest-SVD rows pass the installed H gate through bounded native Gram special cases. The default Gram-SVD materialization gate remains `eigencore.gram_svd_max_dimension = 512`; opt-in larger runs must also pass `eigencore.gram_svd_memory_mb`, rank-fraction, aspect-ratio, work-budget, and original-coordinate certificate gates. Fresh installed 2026-06-07 cutoff evidence certifies 600/768/1024 small-side tall and wide sparse rows, but tall rows are slower than RSpectra while wide rows pass speed, so V2 does not claim a broad cutoff raise. Sparse CSC smallest-SVD targets use a native certified Gram special case when the smaller Gram dimension is inside the explicit materialization gate, and otherwise retain the native certified Golub-Kahan smallest boundary; exact-norm matrix-free callbacks keep their native smallest callback label. Nearest/interior SVD has a native full-subspace Golub-Kahan boundary for sparse CSC inputs and exact-norm matrix-free callbacks, while dense smallest and nearest-sigma requests use exact dense fallback certificates and diagonal/reference prototype interior rows can certify without densification. Scalable refined/shift-invert interior SVD remains future scope. Real matrix-free operators with adjoints have a native callback-cycle Golub-Kahan label and sidecar gate, but quick `600x90` probes are diagnostic only, and general sparse/matrix-free production thick-restart SVD remains future scope. |
| Randomized SVD | Scoped V1 release gate is green for the large exact-low-rank dense planner row. Dense double `randomized()` requests and sparse CSC `randomized()` requests with QR normalization and largest-value targets now use native randomized controllers for sketch generation, subspace iteration, projected-core SVD, residual certification diagnostics, and q=0 early stop. The sparse CSC controller has a green native-controller certificate/provenance contract, but sparse speed rows remain diagnostic rather than promoted 2x release claims. Matrix-free, LU/none-normalized, slow-decay, and broader adaptive planner regimes remain reference-control or diagnostic surfaces until their own certified native-controller gates pass. |
| Generalized SPD | The promoted iterative surface is scoped to sparse shifted-tridiagonal generalized SPD LOBPCG for largest/smallest targets. Dense generalized `auto()` uses the native dense LAPACK fallback. Broader generalized block variants remain future scope, and standard Hermitian LOBPCG/preconditioner promotion is closed as a diagnostic/prototype-only no-promotion decision because it loses to the current eigencore default. |
| B-orthogonal Lanczos | Explicit generalized-SPD `lanczos()` requests use a native transformed Lanczos path for dense and diagonal SPD metrics, including block requests inside the dense/diagonal similarity-transform boundary. Sparse tridiagonal CSC SPD metrics use a native Thomas metric-solve boundary inside the honest reference Lanczos refinement, with no sparse densification. General sparse CSC metrics remain on `Matrix::Cholesky` reference provenance, and sparse-CSC block generalized Lanczos promotion is not claimed. |
| Shift-invert | Dense standard, dense generalized SPD, diagonal standard, sparse symmetric-tridiagonal standard, and sparse/diagonal tridiagonal generalized paths are native and pass the scoped V1 installed gate. `auto()` with `target = nearest(sigma)` routes through an implicit `shift_invert(sigma)` transform for supported factorized Hermitian regimes and fails loudly for matrix-free/unfactorized inputs without a solve. All shift-invert paths carry `shift_invert_factorization_contract_v1`: native labels require eigencore-owned factorized apply, general sparse uses an honest `Matrix::lu` reference contract, and user solves use an external-cache contract. Native general sparse LU and native ownership of user solve functions are explicit PRD non-goals unless a future PRD reopens them. |
| Nonsymmetric eigen | Dense and sparse CSC cases with supported real/imaginary/magnitude targets use a native Arnoldi cycle with native refined Ritz extraction, a wired restart budget, best-attempt retention, right-residual certification, and adjoint-capable left eigenvectors with biorthogonality diagnostics. Real matrix-free callback cases keep the native callback Arnoldi cycle with native projected Ritz extraction and the same certification/restart boundary. Full Krylov-Schur or harmonic/interior extraction, matrix-free refined extraction, and native complex sparse/operator paths remain future scope; matrix-free operators without `apply_adjoint` remain right-only with an explicit warning. |
| Complex-valued inputs | Base complex dense matrices use native dense complex LAPACK labels for Hermitian eigen, general eigen, and SVD calls, and base complex dense operators use native `zgemm` block apply with exact residual/backward-error certificates. Complex matrix-free solver operators fail with actionable future-scope messages, but exact-norm complex callbacks can be certified directly. Complex-valued `Matrix`/sparse inputs and broader complex sparse/matrix-free operator kernels remain future scope and are rejected rather than silently coerced. Real-valued matrices may still return complex eigenpairs through the `eigs()` compatibility path. |
| Matrix-free operators | Supported through callback boundaries. The real nonsymmetric Arnoldi callback path and real matrix-free SVD callback path have native callback labels and sidecar gates; callback-driven paths remain separate from built-in native kernels unless a gate explicitly promotes them. |

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
- Any new V3 solver promotion still needs dense/operator certificate agreement
  tests.
- Generalized and shift-invert certificates must continue to certify the
  original problem, not only the transformed operator.

## Complex ABI And Certificate Contract

Dense complex Hermitian eigen, general eigen, and SVD calls now use native
LAPACK kernels. Base complex dense operator wrappers also use native `zgemm`
block apply with `storage = "complex_dense_matrix"`, `native = TRUE`, and
`native_operator_kernel = "dense_complex_zgemm"`. Native `ScalarType::C128`
is still reserved for broader sparse/matrix-free C++ operator implementations;
the current R operator contract represents complex operators with
`dtype = "complex"` and requires `apply_adjoint()` to implement the conjugate
transpose action `A^* X = Conj(t(A)) X`.

Complex certificates use the same residual formulas as the real paths, but
with conjugate inner products: eigen/SVD orthogonality is computed as `V^* V`,
and SVD residuals are the exact two-sided pair `A v - sigma u` and
`A^* u - sigma v`. Explicit dense complex sources use exact Frobenius scales;
complex matrix-free callbacks with explicit Frobenius norm metadata can be
certified directly with `norm_bound_type = "frobenius_metadata"` and
`scale_is_estimate = FALSE`. Estimated-norm complex callbacks cannot produce a
passed certificate.

## Planner Labels

Trust the planner label over assumptions about the input type. Labels containing
`reference`, `prototype`, or `oracle` identify paths that are useful for
correctness and migration but are not V2 CRAN production performance claims.
For shift-invert, dense standard/generalized, sparse diagonal or
symmetric-tridiagonal standard cases, and tridiagonal `A` with diagonal `B`
have native labels; general sparse and general sparse diagonal-metric
generalized cases remain reference-labelled.

## Release-Hardening Evidence

The V2 CRAN release-hardening surface is:

- CRAN-like checks must be rerun from a fresh tarball in the final release
  candidate gate and may close only with zero failures/warnings and no
  unaccepted notes.
- Full tests must be rerun from the current tree before the final release gate
  closes.
- Strict installed benchmark reports are recorded for promoted solver families;
  quick certification-only runs are not release signoff.
- UBSan install plus `inst/validation/native-smoke.R` passed locally. ASan is
  blocked by macOS/R dynamic loading and valgrind is absent, so those remain
  environment boundaries rather than claimed coverage.
- Migration and user docs have a final scope audit for the current promoted
  surface.
- Mote issues for completed V2_CRAN gates are closed or handed off with any
  remaining scope recorded during final release handoff.
