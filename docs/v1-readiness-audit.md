# eigencore V2 CRAN Readiness Audit

Date: 2026-07-03

This audit maps the V2 CRAN release objective to concrete artifacts and gates.
The promoted solver surfaces retain installed-package evidence; remaining
future work is documented as V3 deferral rather than hidden release debt.
`V2 CRAN` is the internal release-boundary name used by the tracker and audit
files; the package version for this release candidate is `eigencore` 1.0.0.

The final stop-rule checklist lives in `docs/v1-completion-audit.md`. Use that
file before any V2 CRAN completion claim; this audit provides the evidence inventory
that feeds it.

## Objective Restated

V2 CRAN readiness means:

1. Package checks stay clean while the native engine grows.
2. Public solver paths either run native engine code or carry honest
   reference/oracle planner labels.
3. Native operator foundation covers built-in dense, CSC, diagonal, adjoint,
   centered, scaled, summed, composed, and crossprod transforms without silent
   sparse densification.
4. Production Hermitian, SVD, randomized SVD, generalized SPD, shift-invert,
   and nonsymmetric surfaces meet their PRD gates before promotion.
5. Certificates, target taxonomy, result fields, and RSpectra-compatible shims
   remain consistent with the PRD.
6. Release hardening includes full tests, CRAN-like checks, benchmark reports,
   docs, migration guidance, known-limitations documentation, and mote handoff.

## Prompt-to-Artifact Checklist

| Requirement | Evidence artifact | Verification command or check | Current status |
|---|---|---|---|
| PRD scope reconciliation | `prd.json`, `docs/known-limitations.md` | `python -m json.tool prd.json`; search for stale V2 algorithm-expansion language; inspect `release_strategy.v2_cran_release`, `v2_scope`, `v3_scope`, `post_v1_tracker_map`, and milestone `current_status` fields | Green when `prd.json` defines V2 as the CRAN release boundary, moves nonessential solver expansion to V3, qualifies the central product claim, and preserves current-status text for milestones 5-9. |
| Package checks clean | `eigencore_1.0.0.tar.gz` tarball result | `LC_ALL=C LANG=C R CMD build /Users/bbuchsbaum/code/eigencore`; `LC_ALL=C LANG=C R CMD check --as-cran --no-manual eigencore_1.0.0.tar.gz` | Green on 2026-07-03 via `rcmdcheck::rcmdcheck(path = ".", args = c("--as-cran", "--no-manual"), build_args = "--no-manual")` after the generalized-eigen release gate and exported-return-value docs: `0 ERROR`, `0 WARNING`, `1 NOTE` for CRAN new submission. |
| Full test suite clean | `tests/testthat/` | `Rscript -e 'pkgload::load_all("."); testthat::test_dir("tests/testthat", reporter="summary")'`; installed `tests/testthat.R` inside the source-tarball check | Green on 2026-07-03: the full source suite was green after the generalized-eigen gate, and the installed-package `testthat.R` pass in the fresh 1.0.0 CRAN-like check was OK. |
| Diff hygiene | whole repo | `git diff --check` | Green on 2026-07-03 after release metadata, generalized-eigen docs, benchmark, roxygen, and exported return-value documentation updates. |
| CRAN remote checks | GitHub Actions, R-hub, win-builder | GitHub Actions status for the release candidate; R-hub Ubuntu release and Windows; win-builder R-release and R-devel email results | GitHub Actions are green for commit 088d3bd: R-CMD-check, lint, coverage, and pkgdown all completed successfully. R-hub run 28686224852 completed successfully for Ubuntu release and Windows. win-builder R-release and R-devel submissions from 088d3bd were accepted on 2026-07-03, and maintainer email results remain the final external cross-platform evidence gate. |
| Planner honesty | `R/problem.R`, `R/solve.R`, solver result builders, tests | Search method labels; run target/solver tests | Mostly green; nonsymmetric dense/native-CSC Arnoldi compatibility, native matrix-free callback Arnoldi, and remaining reference labels are distinguished. |
| No silent sparse densification | `R/solve.R`, `R/operator_algebra.R`, shift-invert tests, SVD tests | Dense fallback tests; `allow_dense_fallback = "never"` adversarial tests | Mostly green for current public paths. |
| Native operator foundation | `R/operator_algebra.R`, `R/operator.R`, `src/native_operators.cpp`, `src/native_operators.h`, `tests/testthat/test-operator-algebra.R` | `test-operator-algebra.R` | Done for explicit built-ins: dense/CSC/diagonal adjoint, scaling, sum, compose, crossprod, dense centering, CSC centering. Matrix-free centering remains callback-boundary policy. |
| Hermitian G1 default | `R/problem.R`, `R/solve.R`, `R/solve_eigen.R`, `R/certification.R`, `src/small_dense.cpp`, `src/certificates.cpp`, `inst/benchmarks/bench-native-hermitian-gate.R`, `benchmarks/RELEASES.md` | strict Hermitian benchmark script; planner/result/certificate tests | Green for the promoted structured-tridiagonal default. Installed strict `path_laplacian:1000` evidence from 2026-06-05 certifies `20/20` through the native selected tridiagonal solver and tridiagonal residual certificate, with `0.007990162s` median, `579696` bytes, `6.152318x` speed versus RSpectra, `1.073942x` memory versus the best certified reference, and `22.58426x` PRIMME parity. Native block Hermitian Lanczos remains explicit/diagnostic rather than promoted for general sparse `auto()`. The scripts support `--cases=<case[:n]>` and progress messages for further diagnosis. |
| H production SVD | `R/reference_golub_kahan.R`, `R/solve_svd.R`, `R/validation.R`, `inst/benchmarks/bench-svd-surface.R`, `tests/testthat/test-svd-adversarial.R`, `benchmarks/RELEASES.md` | SVD surface gates and target contracts against certified references | Green for the promoted non-quick tall/wide sparse H surface plus sparse smallest/interior target contracts. The SVD surface script supports stable `--cases=` ids, progress messages, `--h-candidate`, and `svd_target_contract` rows. Fresh installed 2026-06-05 strict evidence from `--iterations=1 --h-candidate --methods=eigencore,RSpectra,PRIMME --cases=tall_sparse,wide_sparse --subject=eigencore --strict --save` certifies `20/20` triplets on `tall_sparse` (`100000 x 500`) and `wide_sparse` (`500 x 100000`), uses the bounded native Gram special case with `native_gram_eigensolver = "lapack_dsyevr"` and `materialized_gram = TRUE`, and passes speed/memory (`3.849574x` / `2.049985x` tall, `9.882995x` / `2.049985x` wide). Fresh 2026-06-06 source strict target evidence from `smallest_sparse:5000x500,interior_sparse:5000x500` certifies both target rows and passes contract speed/memory versus dense `base_*` rows: smallest uses `native certified Gram SVD special case`, interior uses `native full-subspace interior Golub-Kahan SVD`. The quick `600x90` fixture is diagnostic only; a final strict rerun was too noisy for signoff and is not the release surface. `complex_dense` is benchmark-visible for the native dense complex SVD label but is excluded from the sparse H speed/memory gate. Retained IRLBA/LBD residual augmentation, BPRO diagnostics, guarded BPRO modes, exact-final lock diagnostics, and cached `A Q` certificate reuse remain benchmark-visible diagnostics, not promoted beyond this H surface. |
| I randomized SVD | `R/reference_golub_kahan.R`, `src/native_operators.cpp`, `inst/benchmarks/bench-randomized-rsvd.R`, `docs/hegelsvd_svd_acceleration.md`, `benchmarks/RELEASES.md` | scoped randomized/rsvd accuracy and time-to-certified-answer release gate plus native controller contracts | Scoped dense release gate and sparse controller contract are green. The randomized-rsvd script supports stable `--cases=` ids, progress messages, `release_gate_required` tags, and `randomized_controller_contract` rows. Fresh 2026-06-06 strict source evidence certifies `exact_low_rank_dense:2000x500` against certified `rsvd` and passes speed/accuracy at about `2.55x`; dense QR randomized requests expose `randomized_controller_native = TRUE`, `randomized_dense_native_controller = TRUE`, `randomized_native_certificate_diagnostics = TRUE`, `randomized_adaptive_stop_used = TRUE`, `randomized_sketch_kind = "native_fused_a_omega"`, `randomized_projection_kind = "native_direct_qt_a"`, and `randomized_core_solver = "native_dense_projected_svd"`. The same strict run certifies `low_rank_sparse:2000x500` through `native_csc_randomized_controller` with `randomized_sparse_native_controller = TRUE` and a green controller contract, while its speed ratio remains diagnostic. Nearly-low-rank, slow-decay, and quick small-size rows remain printed diagnostics; rows with uncertified `rsvd` baselines still record `baseline_certified = FALSE`, `speed_gate = FALSE`, and `passed = FALSE`. Matrix-free/LU/none randomized controllers and slow-decay sparse release candidates remain future scope, not promoted release gates. |
| J generalized SPD LOBPCG | `R/reference_lobpcg.R`, `inst/benchmarks/bench-generalized-lobpcg.R`, `tests/testthat/test-generalized-lobpcg.R` | installed full strict gate plus focused generalized tests | Green for the promoted sparse shifted-tridiagonal surface. Dense generalized `auto()` remains the native dense LAPACK fallback, not an iterative LOBPCG promotion. Fresh installed saved evidence on 2026-06-05 from `bench-generalized-lobpcg.R --iterations=1 --strict --save` certifies `10/10` for both `sparse_generalized_path_smallest:1000` and `sparse_generalized_path_largest:1000`, exposes native generalized/B-orthogonal/shifted-tridiagonal provenance, and passes speed/memory versus dense base at about `1.56x` / `8.8x` for smallest and `1.67x` / `25.0x` for largest. The same run keeps dense fallback, generalized-Lanczos reference, constrained, matrix-free-B, and adversarial-B checks green as scoped contract diagnostics. Standard Hermitian LOBPCG/preconditioner promotion is closed as diagnostic/prototype-only after the 2026-06-06 path-Laplacian probe certified but failed the scalar-speed gate against the current eigencore default. |
| K B-orthogonal Lanczos refinement | `R/reference_generalized_lanczos.R`, `R/problem.R`, `R/solve_eigen.R`, `inst/benchmarks/bench-generalized-lobpcg.R`, `tests/testthat/test-generalized-lobpcg.R` | focused generalized Lanczos native/reference/block contracts; focused benchmark rows; full suite | Green for native dense/diagonal transformed generalized Lanczos, explicit block transformed generalized Lanczos, and sparse CSC metric-boundary diagnostics. Explicit generalized-SPD `method = lanczos()` requests use `native transformed generalized SPD B-orthogonal Lanczos` when `B` is dense or diagonal SPD and the transformed standard operator can run through native Lanczos; the focused strict probe also certifies `diagonal_generalized_block_lanczos_native_smallest` with `block_size = 2`. Sparse CSC SPD metrics keep the honest `reference generalized SPD B-orthogonal Lanczos refinement` solver label; tridiagonal sparse metrics use `native_sparse_tridiagonal_thomas` / `tridiagonal_thomas` metric-solve provenance, while general sparse CSC metrics remain `Matrix::Cholesky` reference provenance. Focused tests cover diagonal, dense, block diagonal, tridiagonal sparse CSC, and general sparse CSC SPD metrics, B-orthogonality, original generalized residual certificates, and no sparse densification. Sparse-CSC block production Lanczos and arbitrary sparse-CSC native factorization are not claimed. |
| L shift-invert | `R/transform_shift_invert.R`, `src/scalar_krylov.cpp`, `src/native_operators.cpp`, `tests/testthat/test-shift-invert.R`, `inst/benchmarks/bench-shift-invert.R` | shift-invert tests and installed benchmark evidence | Green for the scoped V1 surface. Dense standard, dense generalized SPD, diagonal standard, sparse symmetric-tridiagonal standard, and sparse/diagonal tridiagonal generalized shift-invert are native and certify original-coordinate residuals. `auto()` plus `target = nearest(sigma)` now plans and dispatches through an implicit `shift_invert(sigma)` transform for supported factorized Hermitian regimes, and rejects matrix-free/unfactorized inputs without a solve. Fresh installed 2026-06-05 `bench-shift-invert.R --iterations=1 --strict --save` evidence passes all eight native/reference contract rows: the five native rows certify `6/6`; general sparse standard and sparse/diagonal generalized remain honest `Matrix::lu` reference boundaries with estimated-scale certificate honesty; the matrix-free user-solve row records `factorization = user_solve`, `external_cache = TRUE`, and exact original-coordinate certification. General sparse native LU and native ownership of user-supplied solve functions are explicit PRD non-goals unless a future PRD reopens them. |
| Nonsymmetric eigen | `R/reference_arnoldi.R`, `R/solve_eigen.R`, `R/certification.R`, `tests/testthat/test-target-taxonomy.R`, `tests/testthat/test-adversarial.R`, `inst/benchmarks/bench-nonsymmetric.R` | nonsymmetric real/complex residual tests; RSpectra shim tests; nonsymmetric benchmark contract | Green for the scoped V1 compatibility surface, with a current dense/sparse CSC refined-Ritz tranche. Dense and sparse CSC nonsymmetric auto paths run a native Arnoldi cycle with native refined Ritz extraction, exact right-residual certification, planner-wired restart controls, best-attempt retention, and measured native-cycle/native-Ritz stage diagnostics. Real matrix-free callback paths keep native projected Ritz extraction under the native callback Arnoldi label. Base complex dense nonsymmetric paths use the native dense complex general LAPACK label with exact right-residual certification; sparse CSC and matrix-free callback paths stay on bounded `min(n, max(k + 8, 9k))` policy. Fresh installed strict benchmark evidence on 2026-06-05 certifies the original V1 dense/sparse native rows, and a focused installed strict run on 2026-06-06 certifies `matrix_free_nonnormal:30` with `native_matrix_free_arnoldi_label = TRUE`, `matrix_free_native = TRUE`, `ritz_extraction_native = TRUE`, and exact right-residual contracts. Current source tests and benchmark smoke require dense/sparse refined-Ritz provenance. Adjoint-capable rows also return left vectors with separate left-residual and biorthogonality diagnostics; full Krylov-Schur or harmonic/interior extraction, matrix-free refined extraction, and native complex sparse/operator paths remain future scope. |
| Result object fields | `R/solve_result_builders.R`, `R/certification.R`, `tests/testthat/test-result-contracts.R` | result-contract tests plus full suite | Covered for current public paths. The contract test checks eigen and SVD accessor/diagnostics shape across dense Hermitian, native sparse Lanczos, generalized LOBPCG, dense shift-invert, nonsymmetric native Arnoldi, matrix-free Arnoldi, randomized SVD, and RSpectra-compatible shims. Re-audit after new solver promotions. |
| No stochastic estimate marked passed | `R/certification.R`, operator norm metadata, `tests/testthat/test-validation.R` | `test-validation.R`; certificate code inspection | Covered for current certificate entry points. `new_certificate()` with `scale_is_estimate = TRUE`, matrix-free eigen certificates, matrix-free SVD certificates, and residual-backed generalized certificates all withhold `passed` even when all backward-error checks converge. Re-audit only if new certificate constructors or estimated-scale norm paths are added. |
| RSpectra shims | `R/compatibility.R`, target taxonomy tests | `test-target-taxonomy.R` | Current ARPACK target mapping, including `SM`, `LR/SR/LI/SI`, is covered. |
| Benchmark reports | `docs/v1-benchmark-manifest.md`, `benchmarks/RELEASES.md`, `inst/benchmarks/*` | run strict benchmark scripts and save summaries | Green for promoted solver surfaces. The V1 benchmark manifest maps release surfaces to runnable installed-package commands, saved artifacts, and current gate status; the 2026-06-07 all-surface quick-smoke sweep passed every manifest gate as a contract/provenance regression check, while installed strict rows remain the release evidence for speed claims. |
| Documentation scope audit | `docs/v1-doc-scope-audit.md` | manual doc inventory plus `test-bench-smoke.R` source-doc checks | Green for scoped V1 docs after the README/vignette refresh. Future solver promotions must re-open this row. |
| User docs and migration | `README.md`, `README.Rmd`, `vignettes/*.Rmd`, `docs/method-selection-and-workflows.md`, `docs/rspectra-migration.md`, `docs/known-limitations.md`, `docs/v1-doc-scope-audit.md` | doc build / R CMD check plus manual scope audit | Green for the scoped public surface. RSpectra migration, method/workflow selection, known-boundaries, documentation-scope, README, and vignette pages all describe current promoted paths without claiming future solver families. |
| Sanitizer / valgrind-style checks | `inst/validation/native-smoke.R`; temporary UBSan install and smoke run; `benchmarks/RELEASES.md` notes | UBSan `R CMD INSTALL` with sanitizer flags plus `Rscript -e 'source(system.file("validation/native-smoke.R", package = "eigencore"))'`; ASan and valgrind attempts | Scoped green locally. Final 2026-06-05 UBSan-only install and installed native smoke passed, and the smoke is a reusable installed-package artifact. ASan remains blocked locally by macOS/R `dlopen` interceptor ordering, and valgrind is not installed; those are documented environment boundaries rather than hidden V1 solver blockers. |
| Mote handoff | `.mote/` issue state | `mote board`, `mote ls`, issue histories | Final H and release-hardening close operations are recorded in mote history after validation. |

## V2 CRAN Boundaries And Final Checks

- H: the promoted non-quick tall/wide sparse production SVD gate is green on
  fresh installed evidence, and sparse smallest/interior target contracts are
  green on 2026-06-06 source strict evidence. Smallest sparse CSC uses the
  native certified Gram special case inside the explicit Gram-dimension gate;
  interior sparse CSC uses the native full-subspace Golub-Kahan boundary.
  Matrix-free SVD has a native callback-cycle sidecar boundary, but broader
  sparse/matrix-free retained restart work remains documented future scope.
- I: randomized SVD has a scoped V1 release gate for the validated large
  exact-low-rank dense regime, with a native dense QR controller, native fused
  sketch/projection kernels, native projected-core SVD, native certificate
  diagnostics, and a green sparse CSC native-controller contract. Sparse speed,
  quick small-size, slow-decay, and broader matrix-free/LU/none controller work
  remain future scope.
- J: generalized SPD LOBPCG is production-promoted for the sparse
  shifted-tridiagonal largest/smallest surface plus the current broader
  generalized contract bank: native matrix-free-B metrics, shifted native
  preconditioner diagnostics, constraints, adversarial B rows, and honest
  generalized-Lanczos native/reference boundaries. Dense auto still uses the
  native dense generalized LAPACK fallback; arbitrary sparse-CSC metric
  factorization and sparse-CSC block generalized Lanczos variants are not
  claimed by the current gate.
  Standard Hermitian
  LOBPCG/preconditioner promotion is closed as diagnostic/prototype-only after
  it certified but failed the scalar-speed gate against the current eigencore
  default.
- K: B-orthogonal Lanczos now has a native transformed path for dense/diagonal
  SPD metrics, including an explicit block transformed contract, and an
  explicit reference scalar refinement for sparse-CSC SPD metrics. Tridiagonal
  sparse CSC metrics use a native Thomas metric solve inside that
  reference-labelled loop, while general sparse CSC metrics retain
  `Matrix::Cholesky` reference provenance. Sparse-CSC block production Lanczos
  and arbitrary sparse-CSC native factorization are not claimed.
- L: shift-invert is green for the scoped V1 surface. Dense, diagonal, and
  tridiagonal native rows plus general sparse/user-solve reference boundaries
  have installed strict benchmark contracts with explicit cache/factorization
  provenance. `auto()` nearest-sigma planning now uses an implicit
  `shift_invert(sigma)` transform where supported and fails loudly for
  matrix-free/unfactorized inputs without a solve. General sparse native LU and
  native ownership of user-supplied solve functions are explicit PRD non-goals
  unless a future PRD reopens them.
- Nonsymmetric: dense and sparse CSC paths have a native Arnoldi cycle with
  native refined Ritz extraction, while real matrix-free callback paths keep
  native projected Ritz extraction with wired restart controls. The 2026-06-05
  dense/sparse installed strict evidence covers the original V1 native rows,
  and current source tests cover the refined-Ritz tranche. The focused
  installed `matrix_free_nonnormal:30` row requires native callback labels and
  exact right-residual certification. Adjoint-capable rows expose left vectors
  with left-residual and biorthogonality diagnostics; full Krylov-Schur or
  harmonic/interior extraction, matrix-free refined extraction, and native
  complex sparse/operator paths remain future scope.
- N: release docs now include RSpectra migration, method/workflow selection,
  known-boundaries, documentation-scope audit, V1 benchmark manifest,
  README/vignette refresh, and curated release evidence. Solver gates have
  fresh installed evidence; final aggregate signoff includes full testthat,
  installed native smoke, UBSan native smoke, tarball build, CRAN-like check,
  and mote closure.
- Native sanitizer evidence is scoped: UBSan install plus
  `inst/validation/native-smoke.R` are green, while ASan and valgrind-equivalent
  checks are unavailable in this local macOS/R environment.
- G1 is green for the promoted structured-tridiagonal default on
  `path_laplacian:1000`; sparse block `auto()` promotion remains diagnostic
  rather than a production default.

## Stop Rule

Mark the V2 CRAN thread goal complete only after:

1. Every row above is either green with fresh evidence or explicitly scoped out
   by a revised PRD.
2. The strict benchmark gates for promoted solver families have been rerun from
   a clean installed package.
3. `R CMD check --no-manual` is green from a fresh tarball.
4. win-builder R-release and R-devel maintainer emails are confirmed clean for
   the same release candidate as the fresh tarball, GitHub Actions, and R-hub.
5. Mote issues for completed gates are closed or handed off with remaining
   scope recorded.
