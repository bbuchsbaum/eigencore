# eigencore V1 Readiness Audit

Date: 2026-06-05

This audit maps the V1 objective to concrete artifacts and gates. It is not a
release signoff. The current state is still **not V1 ready** because several
PRD-required solver and release-hardening gates remain open.

The final stop-rule checklist lives in `docs/v1-completion-audit.md`. Use that
file before any V1 completion claim; this audit provides the evidence inventory
that feeds it.

## Objective Restated

V1 readiness means:

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
| Package checks clean | `R CMD check --no-manual` tarball result | `R CMD build /Users/bbuchsbaum/code/eigencore`; `env LC_ALL=C LANG=C R CMD check --no-manual eigencore_0.0.0.9000.tar.gz` | Green in latest local runs, but must be rerun before release handoff. |
| Full test suite clean | `tests/testthat/` | `Rscript -e 'pkgload::load_all("."); testthat::test_dir("tests/testthat", reporter="summary")'` | Green with four expected CRAN skips in latest local runs. |
| Diff hygiene | whole repo | `git diff --check` | Green in latest local runs. |
| Planner honesty | `R/problem.R`, `R/solve.R`, solver result builders, tests | Search method labels; run target/solver tests | Mostly green; nonsymmetric dense/native-CSC Arnoldi compatibility and matrix-free reference Arnoldi labels are distinguished. |
| No silent sparse densification | `R/solve.R`, `R/operator_algebra.R`, shift-invert tests, SVD tests | Dense fallback tests; `allow_dense_fallback = "never"` adversarial tests | Mostly green for current public paths. |
| Native operator foundation | `R/operator_algebra.R`, `R/operator.R`, `src/native_operators.cpp`, `src/native_operators.h`, `tests/testthat/test-operator-algebra.R` | `test-operator-algebra.R` | Done for explicit built-ins: dense/CSC/diagonal adjoint, scaling, sum, compose, crossprod, dense centering, CSC centering. Matrix-free centering remains callback-boundary policy. |
| Hermitian G1 default | `R/problem.R`, `R/solve.R`, `R/solve_eigen.R`, `R/certification.R`, `src/small_dense.cpp`, `src/certificates.cpp`, `inst/benchmarks/bench-native-hermitian-gate.R`, `benchmarks/RELEASES.md` | strict Hermitian benchmark script; planner/result/certificate tests | Green for the promoted structured-tridiagonal default. Installed strict `path_laplacian:1000` evidence from 2026-06-05 certifies `20/20` through the native selected tridiagonal solver and tridiagonal residual certificate, with `0.007990162s` median, `579696` bytes, `6.152318x` speed versus RSpectra, `1.073942x` memory versus the best certified reference, and `22.58426x` PRIMME parity. Native block Hermitian Lanczos remains explicit/diagnostic rather than promoted for general sparse `auto()`. The scripts support `--cases=<case[:n]>` and progress messages for further diagnosis. |
| H production SVD | `R/reference_golub_kahan.R`, `R/solve_svd.R`, `R/validation.R`, `inst/benchmarks/bench-svd-surface.R`, `tests/testthat/test-svd-adversarial.R`, `benchmarks/RELEASES.md` | SVD surface gates against certified references | Green for the promoted tall/wide sparse H surface. The SVD surface script supports stable `--cases=` ids and progress messages, and the executable `--h-candidate` gate targets the promoted `eigencore` SVD path. The tall-sparse production row uses native right-normal `implicit_normal_lanczos` without materializing the Gram and returns exact original-coordinate certificate diagnostics; the wide row remains on the certified native left-Gram special case. Warning-free installed 2026-06-05 3-iteration quick probes on `tall_sparse:600x90` and `wide_sparse:90x600` certify all five requested triplets and pass speed/memory (`1.194x` / `2.392x` tall, `1.223x` / `2.392x` wide). Retained IRLBA/LBD residual augmentation, BPRO diagnostics, guarded BPRO modes, exact-final lock diagnostics, and cached `A Q` certificate reuse remain benchmark-visible diagnostics, not promoted beyond this H surface. |
| I randomized SVD | `R/reference_golub_kahan.R`, `src/native_operators.cpp`, `inst/benchmarks/bench-randomized-rsvd.R`, `docs/hegelsvd_svd_acceleration.md`, `benchmarks/RELEASES.md` | scoped randomized/rsvd accuracy and time-to-certified-answer release gate | Scoped V1 gate green. The randomized-rsvd script supports stable `--cases=` ids, progress messages, and `release_gate_required` tags; `--strict` now enforces only release/planner candidate rows. Fresh installed 2026-06-05 strict non-quick evidence certifies `exact_low_rank_dense:2000x500` against certified `rsvd` and passes speed/accuracy at about `3.1x`; native dense/CSC fused sketch and projection diagnostics are benchmark-visible as `randomized_native_sketch = TRUE`, `randomized_sketch_kind = "native_fused_a_omega"`, and `randomized_projection_kind = "native_direct_qt_a"`. Nearly-low-rank, slow-decay, quick exact-low-rank, and quick sparse rows remain printed diagnostics with `release_gate_required = FALSE`; rows with uncertified `rsvd` baselines still record `baseline_certified = FALSE`, `speed_gate = FALSE`, and `passed = FALSE`. Broader fully native randomized controller and slow-decay/sparse release candidates remain future scope, not promoted V1 gates. |
| J generalized SPD LOBPCG | `R/reference_lobpcg.R`, `inst/benchmarks/bench-generalized-lobpcg.R`, `tests/testthat/test-generalized-lobpcg.R` | installed full strict gate plus focused generalized tests | Green for the promoted sparse shifted-tridiagonal surface. Dense generalized `auto()` remains the native dense LAPACK fallback, not an iterative LOBPCG promotion. Fresh installed saved evidence on 2026-06-05 from `bench-generalized-lobpcg.R --iterations=1 --strict --save` certifies `10/10` for both `sparse_generalized_path_smallest:1000` and `sparse_generalized_path_largest:1000`, exposes native generalized/B-orthogonal/shifted-tridiagonal provenance, and passes speed/memory versus dense base at about `1.56x` / `8.8x` for smallest and `1.67x` / `25.0x` for largest. The same run keeps dense fallback, generalized-Lanczos reference, constrained, matrix-free-B, and adversarial-B checks green as scoped contract diagnostics. |
| K B-orthogonal Lanczos refinement | `R/reference_generalized_lanczos.R`, `R/problem.R`, `R/solve_eigen.R`, `inst/benchmarks/bench-generalized-lobpcg.R`, `tests/testthat/test-generalized-lobpcg.R` | focused generalized Lanczos adversarial/agreement tests; installed focused benchmark rows; full suite | Partial reference. Explicit generalized-SPD `method = lanczos()` requests now use an honest `reference generalized SPD B-orthogonal Lanczos refinement` label when `B` has a dense, diagonal, or CSC SPD solve. Focused tests cover diagonal, dense, and sparse CSC SPD metrics, B-orthogonality, original generalized residual certificates, and agreement with LOBPCG certificates. Fresh installed focused quick strict benchmark rows on 2026-05-17 certify `diagonal_generalized_lanczos_ref_smallest` and `sparse_csc_generalized_lanczos_ref_smallest`, pass the generalized-Lanczos reference contract, and record diagonal / sparse-Cholesky metric-solve provenance. Native/block production promotion remains open. |
| L shift-invert | `R/transform_shift_invert.R`, `src/scalar_krylov.cpp`, `src/native_operators.cpp`, `tests/testthat/test-shift-invert.R`, `inst/benchmarks/bench-shift-invert.R` | shift-invert tests and installed benchmark evidence | Green for the scoped V1 surface. Dense standard, dense generalized SPD, diagonal standard, sparse symmetric-tridiagonal standard, and sparse/diagonal tridiagonal generalized shift-invert are native and certify original-coordinate residuals. Fresh installed 2026-06-05 `bench-shift-invert.R --iterations=1 --strict --save` evidence passes all eight native/reference contract rows: the five native rows certify `6/6`; general sparse standard and sparse/diagonal generalized remain honest `Matrix::lu` reference boundaries with estimated-scale certificate honesty; the matrix-free user-solve row records `factorization = user_solve`, `external_cache = TRUE`, and exact original-coordinate certification. General sparse native LU is future scope rather than a scoped V1 blocker. |
| Nonsymmetric eigen | `R/reference_arnoldi.R`, `R/solve_eigen.R`, `R/certification.R`, `tests/testthat/test-target-taxonomy.R`, `tests/testthat/test-adversarial.R`, `inst/benchmarks/bench-nonsymmetric.R` | nonsymmetric real/complex residual tests; RSpectra shim tests; nonsymmetric benchmark contract | Partial. Dense and sparse CSC real- and imaginary-target nonsymmetric auto paths now run a native Arnoldi cycle with native projected Ritz extraction, exact right-residual certification, planner-wired restart controls, best-attempt retention, and measured native-cycle/native-Ritz stage diagnostics. Dense explicit compatibility uses the full dense subspace by default so it does not regress from the old oracle surface; sparse CSC stays on the bounded `min(n, max(k + 8, 9k))` policy. Matrix-free real-spectrum paths remain on the honest `reference Arnoldi (prototype/oracle fallback)`. Installed non-quick strict benchmark evidence on 2026-05-17 certifies `dense_native_arnoldi_lm`, `dense_native_arnoldi_li`, `dense_eigs_native_arnoldi_li`, `sparse_native_arnoldi_lr`, and `sparse_native_arnoldi_li` with `native_arnoldi_label = TRUE`, `ritz_extraction_native = TRUE`, and exact right-residual certificates. Production-grade fully native restarted Arnoldi remains open because restart policy, matrix-free native callback support, and final performance gates remain compatibility-grade. |
| Result object fields | `R/solve_result_builders.R`, `R/certification.R`, `tests/testthat/test-result-contracts.R` | result-contract tests plus full suite | Covered for current public paths. The contract test checks eigen and SVD accessor/diagnostics shape across dense Hermitian, native sparse Lanczos, generalized LOBPCG, dense shift-invert, nonsymmetric native Arnoldi, matrix-free Arnoldi, randomized SVD, and RSpectra-compatible shims. Re-audit after new solver promotions. |
| No stochastic estimate marked passed | `R/certification.R`, operator norm metadata, `tests/testthat/test-validation.R` | `test-validation.R`; certificate code inspection | Covered for current certificate entry points. `new_certificate()` with `scale_is_estimate = TRUE`, matrix-free eigen certificates, matrix-free SVD certificates, and residual-backed generalized certificates all withhold `passed` even when all backward-error checks converge. Re-audit only if new certificate constructors or estimated-scale norm paths are added. |
| RSpectra shims | `R/compatibility.R`, target taxonomy tests | `test-target-taxonomy.R` | Current ARPACK target mapping, including `SM`, `LR/SR/LI/SI`, is covered. |
| Benchmark reports | `docs/v1-benchmark-manifest.md`, `benchmarks/RELEASES.md`, `inst/benchmarks/*` | run strict benchmark scripts and save summaries | Partial but auditable. The V1 benchmark manifest now maps release surfaces to runnable installed-package commands, saved artifacts, and current gate status. The release report is not final because H, I, J, K, sparse-native L, and Arnoldi remain open or red. |
| Documentation scope audit | `docs/v1-doc-scope-audit.md` | manual doc inventory plus `test-bench-smoke.R` source-doc checks | Partial. Current documentation coverage is mapped; final docs/signoff remain blocked by solver gates, final strict benchmark artifacts, and sanitizer/valgrind-style evidence. |
| User docs and migration | `README.md`, `README.Rmd`, `vignettes/*.Rmd`, `docs/method-selection-and-workflows.md`, `docs/rspectra-migration.md`, `docs/known-limitations.md`, `docs/v1-doc-scope-audit.md` | doc build / R CMD check plus manual scope audit | Partial but broader. RSpectra migration, method/workflow selection, known-limitations, and documentation-scope audit pages now exist; final release still needs README/vignette refresh after solver gates close. |
| Sanitizer / valgrind-style checks | `inst/validation/native-smoke.R`; temporary UBSan install and smoke run; `benchmarks/RELEASES.md` notes | UBSan `R CMD INSTALL` with sanitizer flags plus `Rscript -e 'source(system.file("validation/native-smoke.R", package = "eigencore"))'`; ASan and valgrind attempts | Partial. UBSan-only install and native smoke passed on 2026-05-17, and the smoke is now a reusable installed-package artifact. ASan install is blocked locally by macOS/R `dlopen` interceptor ordering, and valgrind is not installed, so final sanitizer/valgrind-style release coverage remains incomplete. |
| Mote handoff | `.mote/` issue state | `mote board`, `mote ls`, issue histories | In progress. Completed slices are noted; several V1 blockers remain open. |

## Current Blocking Gaps

- H: the promoted tall/wide sparse production SVD gate is green on fresh
  installed evidence. General sparse and matrix-free SVD beyond that surface
  still depend on staging paths and benchmark-negative retained restart work,
  and should remain documented limitations until a broader native retained
  restart path is ready.
- I: randomized SVD has a scoped V1 release gate for the validated large
  exact-low-rank dense regime, with native fused sketch/projection kernels and
  honest reference-control public labels. Quick small-size, sparse-size,
  slow-decay, and fully native-controller work remain future scope.
- J: generalized SPD LOBPCG is production-promoted for the sparse
  shifted-tridiagonal largest/smallest surface only. Dense auto still uses the
  native dense generalized LAPACK fallback, and broader generalized
  preconditioner/block variants remain future scope, but the installed full
  strict gate is green for the promoted sparse surface and its adversarial B
  contract bank.
- K: B-orthogonal Lanczos now has an explicit reference scalar refinement path
  for dense/diagonal/CSC SPD metric solves plus focused installed diagonal and
  sparse-CSC benchmark evidence, but native/block production promotion remains
  open.
- L: shift-invert is green for the scoped V1 surface. Dense, diagonal, and
  tridiagonal native rows plus general sparse/user-solve reference boundaries
  have installed strict benchmark contracts with explicit cache/factorization
  provenance. General sparse native LU is future scope.
- Nonsymmetric: sparse CSC real- and imaginary-target paths have a native
  Arnoldi cycle and native projected Ritz extraction with wired restart controls, but
  production-grade fully native restarted Arnoldi is not implemented; dense
  LAPACK oracle and matrix-free reference Arnoldi fallbacks remain honest and
  certified where their current contracts apply.
- N: release docs now include RSpectra migration, method/workflow selection,
  known-limitations, documentation-scope audit, and V1 benchmark manifest
  pages, but sanitizer/valgrind-style runs, final benchmark reports, and final
  README/vignette refresh after solver gates close are incomplete. The
  installed-package Hermitian gate has fresh non-quick evidence, but that
  evidence is red rather than release signoff.
- Native sanitizer evidence is partial: UBSan install plus
  `inst/validation/native-smoke.R` are green, but ASan and
  valgrind-equivalent checks are not yet green in this environment.
- G1 is green for the promoted structured-tridiagonal default on
  `path_laplacian:1000`; sparse block `auto()` promotion remains diagnostic
  rather than a production default.

## Stop Rule

Do not mark the thread goal complete until:

1. Every row above is either green with fresh evidence or explicitly scoped out
   by a revised PRD.
2. The strict benchmark gates for promoted solver families have been rerun from
   a clean installed package.
3. `R CMD check --no-manual` is green from a fresh tarball.
4. Mote issues for completed gates are closed or handed off with remaining
   scope recorded.
