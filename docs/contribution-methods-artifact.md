# V2 CRAN Contribution Methods Artifact

Date: 2026-06-07

This document is the contribution-facing package of the V2 CRAN release
evidence. It does not replace `prd.json`, `docs/v1-readiness-audit.md`, or
`docs/v1-benchmark-manifest.md`; it turns the preserved V1/V2 evidence into a
methods story with explicit losses, non-promoted boundaries, and V3 deferrals.

## Contribution Thesis

eigencore's credible contribution is not "another eigensolver wrapper." The
current contribution is a certified, native-first spectral-computation surface
for R with planner labels that tell users what actually ran. The package is
valuable where four pieces line up:

1. A native or explicitly reference-labelled solver path.
2. A target taxonomy that preserves ARPACK/RSpectra-style intent.
3. Residual/backward-error certificates with visible scale provenance.
4. Benchmark gates run from installed-package artifacts, not source-only smoke.

The unqualified product claim remains post-V2 ambition. The scoped CRAN-facing
claim is: promoted V2 release surfaces provide certified native-first answers with
explicit planner labels, no silent sparse densification, and documented
reference/prototype boundaries.

## Method Design

### Planner Labels

Method labels are part of the public trust contract. Labels containing
`reference`, `prototype`, `oracle`, or `diagnostic` identify paths that may be
useful for migration, validation, or future development but are not promoted
performance claims. A solver surface becomes promoted only after its installed
gate passes and the label, certificate provenance, and benchmark report all
agree.

### Certificates

Result certificates report residuals, backward error, orthogonality, the norm
scale, and whether the scale was exact or estimated. Estimated-scale paths must
not silently become passed certificates. Generalized and shift-invert paths
certify the original residual, not only the transformed operator.

### Operator Boundary

Built-in dense, CSC, diagonal, adjoint, centered, scaled, summed, composed, and
crossproduct operators have native-backed provenance where the current V2 CRAN
surface claims it. Matrix-free callbacks are supported, but callback-driven
paths are separate from promoted built-in native kernels. Sparse workloads
should be tested with `allow_dense_fallback = "never"` when memory behavior is
part of the claim.

### API Compatibility

The RSpectra-shaped shims preserve familiar `eigs()`, `eigs_sym()`, and `svds()`
entry points while returning certificate and diagnostics metadata. This is not
unlimited compatibility. Base complex dense matrices use native dense complex
LAPACK labels for Hermitian eigen, general eigen, and SVD calls, and base
complex dense operators use a native `zgemm` block-apply kernel; complex sparse
`Matrix` and promoted complex sparse/matrix-free operator kernels remain future
scope.
Real-valued matrices may still return complex eigenpairs through the supported
general-eigen target taxonomy.

## Benchmark Report

Use `docs/v1-benchmark-manifest.md` for runnable commands and exact artifact
names. The table below is the contribution-level interpretation.

| Surface | Current claim | Evidence | Losses and boundaries |
|---|---|---|---|
| Hermitian eigen | Promoted for the structured-tridiagonal sparse/diagonal default plus scoped dense/CSC native surfaces. | `bench-native-hermitian-gate.R`; installed `path_laplacian:1000` evidence certifies `20/20` with `6.15x` speed versus RSpectra and PRIMME parity. | General sparse block-Lanczos `auto()` has a documented no-promotion decision as of 2026-06-06; explicit block and diagnostic opt-in paths remain available, but new sparse block promotion is V3 scope. |
| SVD | Promoted only for the non-quick tall/wide sparse bounded native Gram special cases. | `bench-svd-surface.R`; installed tall/wide sparse evidence certifies `20/20`, with about `3.85x` and `9.88x` speed over certified references. | Quick `600x90` probes were rejected as signoff; retained IRLBA/LBD/BPRO and general sparse or matrix-free thick-restart SVD are closed as diagnostic/no-promotion in the current PRD under `bd-01KTE8J9SF16Y1832D8HQQ9KEC`. |
| Smallest/interior SVD targets | Smallest-SVD targets are promoted where exact certification is available: sparse CSC inputs use a native certified Golub-Kahan label, and matrix-free callbacks require explicit non-estimated norm metadata. Nearest/interior targets now have a native full-subspace Golub-Kahan boundary for sparse CSC inputs and exact-norm matrix-free callbacks. | `test-svd-adversarial.R`; dense smallest and nearest-sigma rows have exact two-sided certificates, sparse CSC and exact-norm matrix-free smallest/interior rows carry native labels, diagonal/reference prototype interior rows certify without densification, and matrix-free interior rows without exact norm metadata fail loudly. `bench-svd-surface.R` includes benchmark-visible `smallest_sparse` and `interior_sparse` rows. | Scalable refined/shift-invert sparse and matrix-free interior SVD is closed as current-PRD no-promotion under `bd-01KTEH6862GB19JJWX2M3FQP6T`; future promotion needs a revised PRD and strict benchmark rows proving more than a full-subspace exact boundary. |
| Randomized SVD | Scoped V1 row is green for large exact-low-rank dense workloads. Dense double and sparse CSC QR-normalized largest-target requests now use native randomized controllers with native sketch, subspace iteration, projected-core SVD, residual diagnostics, and q=0 early stop; sparse speed/non-QR/matrix-free regimes remain diagnostic or reference-control. | `bench-randomized-rsvd.R`; fresh strict evidence certifies the dense exact-low-rank row and passes at about `2.55x` versus certified `rsvd`, with native controller diagnostics. The same run certifies `low_rank_sparse:2000x500` through `native_csc_randomized_controller` with a green controller contract. | Slow-decay, sparse speed promotion, matrix-free, and LU/none native randomized controllers are not promoted in the current PRD after the closed no-promotion decision under `bd-01KTE8JFKPA90ZJTXK496SBMK4`. |
| Generalized SPD LOBPCG | Promoted for sparse shifted-tridiagonal generalized largest/smallest targets, native matrix-free-B metrics, native constrained rows, and adversarial SPD metric contracts; dense generalized `auto()` remains a native dense LAPACK fallback rather than an iterative promotion. Standard Hermitian LOBPCG remains diagnostic/prototype-only after the 2026-06-06 no-promotion decision. | `bench-generalized-lobpcg.R`; current strict `--iterations=1` evidence certifies the full generalized gate, including sparse performance rows at `10/10`, native shifted-tridiagonal preconditioner provenance, adversarial matrix-free-B / ill-conditioned diagonal-B / sparse-CSC-B contracts, constraints, and generalized-Lanczos native/reference/block contracts. `bench-lobpcg-preconditioned.R` documents the standard Hermitian no-promotion evidence. | Sparse-CSC block B-orthogonal Lanczos and arbitrary sparse-CSC native factorization are not claimed; standard Hermitian LOBPCG/preconditioner promotion is closed under `bd-01KTEH4G1QPR4RT14B4G78PF1M`. |
| Generalized eigen (replacement surface) | Dense full generalized pencils (`eig_full()`), QZ (`generalized_schur()`), and GSVD (`generalized_svd()`) use native dense LAPACK with original-coordinate certificates, left eigenvectors, conditioning diagnostics where available, and `pencil_norm_scaled` alpha/beta classification. Sparse partial general pencils with nonsingular diagonal `B` use native transformed Arnoldi without densification; singular/general sparse `B` and sparse QZ are explicit unsupported boundaries. | `bench-generalized-eigen.R`; quick strict gate certifies dense SPD/general real/complex rows, sparse diagonal-`B` partial Arnoldi, unsupported sparse boundaries, and dense QZ/GSVD reconstruction rows with certificate/oracle agreement and honest planner labels (`bd-01KVWKVFC8QFZSTR6KFQ3MM3NH`). | Sparse QZ, general sparse `B`, and shift-invert/user-solve sparse general-pencil paths are not promoted; speed ratios on dense full solves are diagnostic rather than release performance claims. |
| B-orthogonal Lanczos | Dense and diagonal SPD metrics now use a native transformed generalized Lanczos path, including explicit block requests inside that transformed boundary; sparse CSC metrics remain an honest reference refinement. | Focused generalized-Lanczos rows in `bench-generalized-lobpcg.R`; native dense/diagonal contract rows, the block transformed contract row, and sparse-CSC reference rows certify in original generalized coordinates. | Sparse-CSC block B-orthogonal Lanczos and arbitrary sparse-CSC native factorization are not claimed. |
| Shift-invert | Promoted for dense standard/generalized, diagonal standard, sparse symmetric-tridiagonal standard, and tridiagonal generalized cases. `auto()` with `target = nearest(sigma)` now routes through an implicit `shift_invert(sigma)` transform where factorized shift-invert is supported. General sparse and user-solve paths are explicitly owned as reference/external-cache contracts, not native claims. | `bench-shift-invert.R`; installed contract rows certify native and reference boundaries in original coordinates, current rows expose `shift_invert_factorization_contract_v1`, and `test-shift-invert.R` covers auto nearest planning/dispatch plus matrix-free failure. | Native general sparse LU and stronger user-solve ownership are closed as an explicit PRD non-goal under `bd-01KTEH4ZPPBESDDPR0Z5MRSQCG`; unsupported matrix-free/unfactorized nearest requests fail loudly rather than remapping targets. |
| Nonsymmetric eigen | Compatibility surface is green for dense, sparse CSC, and real matrix-free callback operators with supported real/imaginary/magnitude targets. Dense and sparse CSC explicit matrices include the native refined Ritz release tranche. | `bench-nonsymmetric.R` and the operator sidecar gate; current rows certify native Arnoldi cycle plus native Ritz extraction and right-residual contracts, including refined-Ritz provenance for dense/sparse CSC and the focused 2026-06-06 `matrix_free_nonnormal:30` projected-Ritz callback row. Adjoint-capable rows return left vectors with a separate left-residual/biorthogonality certificate. | `bd-01KTF6H41S9XDN286TR3V184P4` completes native refined Ritz extraction for dense/sparse CSC only. Full Krylov-Schur or harmonic/interior extraction, matrix-free refined extraction, and native complex sparse/operator paths remain V3 scope under `bd-01KTEH5JM64A4CBZG7ECBWT9WB` and `bd-01KTEH60X91VZRSW7NGV65FBDR`. |
| Complex-valued inputs | Dense complex Hermitian eigen, general eigen, and SVD calls now use native dense complex LAPACK labels with exact certificates, and base complex dense operators now use native `zgemm` block apply with `native_operator_kernel = "dense_complex_zgemm"`. Complex callbacks with explicit Frobenius norm metadata can be certified directly without treating the scale as estimated. | `test-operator.R`, `test-target-taxonomy.R`, `test-solvers.R`, `test-result-contracts.R`, `test-svd-adversarial.R`, and `test-validation.R` cover base complex dense native labels, exact certificates, native complex dense block apply, complex adjoint metadata, direct metadata-scale callback certification, actionable matrix-free future-scope errors, and real-input complex eigenpair compatibility. `bench-nonsymmetric.R` and `bench-svd-surface.R` include dense complex benchmark-visible rows. | Complex sparse `Matrix` inputs and promoted complex matrix-free solver operators are closed as current-PRD non-promoted scope under `bd-01KTEH60X91VZRSW7NGV65FBDR`; broader complex sparse/operator kernels remain future-scoped beyond the base dense `zgemm` tranche in `bd-01KTF24G7XC2PXV80XSJSECZN4`. |
| Operator algebra and matrix-free boundaries | Native provenance is covered for dense/CSC/diagonal scaling, sums, composition, crossprod, centering, adjoints, and selected callback solver boundaries. | `test-operator-algebra.R` and `bench-post-v1-operator-sidecars.R`; current quick strict sidecar evidence passes matrix-free SVD, matrix-free nonsymmetric Arnoldi, and matrix-free generalized-B callback rows. | DelayedArray/out-of-core support and automatic adapters for every R matrix ecosystem remain explicit PRD non-goals/deferred scope unless a future PRD reopens them. |

### Closed No-Promotion Decisions

- `bd-01KTEH48HH4X8G9Q69HHAJ983B`: general sparse block Hermitian Lanczos
  `auto()` promotion is closed as a documented no-promotion decision. Current
  non-quick `path_laplacian:1000`, `k = 20` evidence certifies the explicit
  block candidate but shows about `0.458s` time-to-certified-answer versus
  about `0.0446s` for certified RSpectra and about `0.0073s` for the current
  promoted eigencore path. A `path_laplacian:10000` probe was stopped after
  roughly two minutes without a completed row, so it is not promotion evidence.
- `bd-01KTEH4G1QPR4RT14B4G78PF1M`: standard Hermitian
  LOBPCG/preconditioner promotion is closed as a documented no-promotion
  decision. The native shifted-tridiagonal LOBPCG row certifies and beats
  certified external references on `path_laplacian` `n = 200/1000/2000`,
  `k = 5`, but loses to the current eigencore scalar/tridiagonal default with
  speed ratios of about `0.65`, `0.32`, and `0.14`. Keep it diagnostic unless a
  future redesign beats the current default.
- `bd-01KTE8JFKPA90ZJTXK496SBMK4`: randomized SVD broader promotion is closed
  for the current PRD after the native-controller boundary was completed. Dense
  exact-low-rank QR requests pass the strict release gate at about `2.55x`
  versus certified `rsvd`; sparse CSC QR requests have a green native-controller
  certificate/provenance contract, but sparse speed remains diagnostic and
  matrix-free, LU/none, and slow-decay regimes remain unpromoted.
- `bd-01KTE8JVEPGA1EEQYERZS1V7S1`: arbitrary sparse-CSC metric factorization
  and sparse-CSC block generalized B-orthogonal Lanczos promotion are closed as
  current-PRD no-promotion boundaries. Dense/diagonal transformed generalized
  Lanczos, including the focused block row, is covered; sparse tridiagonal B
  exposes native Thomas metric-solve provenance inside the honest reference
  loop; general sparse CSC B stays `Matrix::Cholesky` reference-labelled.
- `bd-01KTE8J9SF16Y1832D8HQQ9KEC`: hard sparse/matrix-free thick-restarted SVD
  promotion is closed as diagnostic-only for the current PRD. The retained
  native workspace and matrix-free callback-boundary slices are implemented, but
  hard SVD candidate evidence remains speed or memory red against certified
  RSpectra/irlba references. Keep retained IRLBA/LBD/BPRO and broad sparse or
  matrix-free thick-restart rows diagnostic unless a future PRD reopens them
  with installed evidence.
- `bd-01KTEH6862GB19JJWX2M3FQP6T`: scalable refined/shift-invert sparse or
  matrix-free smallest/interior SVD promotion is closed for the current PRD.
  The current claim is exact certified boundaries: sparse CSC smallest
  Gram/Golub-Kahan, sparse CSC interior full-subspace Golub-Kahan, exact-norm
  matrix-free callback boundaries, and loud failures for unsupported regimes.

### Closed Non-Goal Decisions

- `bd-01KTEH4ZPPBESDDPR0Z5MRSQCG`: native general sparse shift-invert
  LU/factorized apply and native ownership of user-supplied solve functions are
  closed as an explicit PRD non-goal. The current contribution is the ownership
  contract: general sparse rows use `Matrix::lu_reference_factorization` with
  sparse pivot diagnostics, cache provenance, no dense rcond, and estimated
  scale honesty; user-solve rows use `user_supplied_solve` external-cache
  provenance; native labels remain reserved for `eigencore_native_factorization`.
- `bd-01KTEH5JM64A4CBZG7ECBWT9WB`: full Krylov-Schur-style,
  harmonic/interior, matrix-free refined, and native complex sparse/operator
  nonsymmetric extraction remain future scope, not current-PRD blockers. The
  current contribution is the native Arnoldi compatibility surface: dense and
  sparse CSC rows now carry native refined Ritz labels, real matrix-free
  callback rows carry native projected Ritz callback labels, and all promoted
  rows keep restart diagnostics, best-attempt retention, exact right-residual
  certification, and adjoint-capable left-vector diagnostics where an adjoint is
  available.
- `bd-01KTF6H41S9XDN286TR3V184P4`: native refined Ritz extraction is complete
  for dense and sparse CSC nonsymmetric Arnoldi paths. Planner controls and
  restart diagnostics expose `arnoldi_extraction = "refined_ritz"`,
  `refined_extraction_native = TRUE`, and `krylov_schur = FALSE`; matrix-free
  callback Arnoldi remains projected-Ritz scope.
- `bd-01KTEH60X91VZRSW7NGV65FBDR`: complex sparse Matrix and promoted complex
  matrix-free solver operators are current-PRD non-goals. Dense complex LAPACK
  eigen/SVD paths, base complex dense `zgemm` block apply, and exact
  certificates are implemented; complex callbacks can carry exact metadata-scale
  certificates; unsupported complex sparse/operator solver paths fail with
  actionable future-scope messages instead of silent real coercion.

### Completed Post-V1 Surfaces

- `bd-01KTEH4RA4NJBQSKF9JQV3V4BS`: broader generalized SPD LOBPCG is complete
  against the current gate. The strict generalized benchmark passes sparse
  largest/smallest performance rows, native shifted-tridiagonal preconditioner
  diagnostics, dense-fallback policy rows, constraints, adversarial B metrics,
  and native/reference generalized-Lanczos contract rows.
- `bd-01KTEH57J1RR3SP1SJ27YRC0ZE`: nearest-sigma planner coverage is complete
  for the supported shift-invert surface. `auto()` plus
  `target = nearest(sigma)` now plans and dispatches through an implicit
  `shift_invert(sigma)` transform for dense and sparse factorized regimes,
  keeps sparse general rows reference-labelled, certifies original-coordinate
  residuals, and rejects matrix-free inputs without a solve rather than
  remapping to largest or smallest.
- `bd-01KTEH6HNN33M15YNGW7T35RQR`: current operator algebra and matrix-free
  sidecar boundaries are complete for the scoped contribution surface. Built-in
  dense/CSC/diagonal transforms preserve native provenance and no-densification
  guarantees, while matrix-free SVD, nonsymmetric Arnoldi, and generalized-B
  callback rows pass the strict sidecar gate. Broader ecosystem adapters remain
  deliberately deferred in `prd.json`.
- `bd-01KTEH5SRWDHXBZXK5CPHBT6G2`: nonsymmetric left eigenvectors and
  biorthogonal certificates are complete for adjoint-capable Arnoldi paths.
  Dense, sparse CSC, and matrix-free operators with `apply_adjoint` now expose
  `left_vectors`, `right_vectors`, a separate
  `left_residual_biorthogonal_backward_error` certificate, and
  biorthogonality diagnostics. Matrix-free operators without `apply_adjoint`
  remain right-only with an explicit warning.

## Migration Package

The migration-facing artifacts are:

- `docs/rspectra-migration.md` for `RSpectra::eigs()`, `eigs_sym()`, and
  `svds()` call-site migration.
- `docs/method-selection-and-workflows.md` for choosing targets, methods,
  shift-invert, generalized SPD, and operator workflows.
- `docs/known-limitations.md` for the short current-boundaries table.
- `vignettes/eigencore.Rmd` and `vignettes/certificates.Rmd` for package-level
  usage and certificate interpretation.
- `docs/post-v1-benchmark-gates.md` for broadening claims after scoped V1.

## V3 Deferral Program

V2 does not add solver families beyond the CRAN release surface. Jacobi-Davidson,
Davidson, full nonsymmetric Krylov-Schur or harmonic/interior workflows,
scalable sparse/matrix-free interior SVD, native general sparse LU ownership,
complex sparse/operator kernels, GraphBLAS/GPU/distributed/SLEPc/PRIMME plugins,
and broad matrix ecosystem adapters are V3 work. The post-V1 gates remain useful
as umbrella-owned regression evidence, but broader claims now require a revised
PRD, a new owner issue, and installed strict or long evidence before planner
labels move. Gate owner ids stay current through
`inst/benchmarks/post-v1-gate-manifest.R`, and the benchmark smoke test fails if
a post-V1 gate points at a stale closed owner.
