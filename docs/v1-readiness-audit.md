# eigencore V1 Readiness Audit

Date: 2026-05-17

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
| Hermitian block Lanczos G1 | `R/problem.R`, `src/block_lanczos.cpp`, `src/scalar_krylov.cpp`, `inst/benchmarks/bench-native-hermitian-gate.R`, `inst/benchmarks/bench-hermitian-sparse.R`, `benchmarks/RELEASES.md` | strict Hermitian benchmark scripts with dense diagnostics; planner tests | Red on fresh installed evidence. A 2026-05-17 quick strict run certified path-Laplacian and dense rows but failed quick speed/memory diagnostics. A filtered non-quick installed run for `path_laplacian:1000` under the old sparse block auto policy certified and passed memory but failed speed (`0.0064x` vs RSpectra) and PRIMME parity (`0.025x`). Sparse `auto()` promotion is now disabled by default and remains available only via diagnostic opt-in. A post-demotion installed scalar rerun was initially red at `5.41s`; after the native append-only projected-matrix update it improves to `0.327s` and certifies all 20 pairs, but still fails speed (`0.138x`), memory (`0.636x`), and PRIMME parity (`0.513x`). The post-fix block candidate is also red (`0.488s`, memory green, speed `0.093x`, PRIMME parity `0.356x`). Full G1 release signoff is not green. The scripts support `--cases=<case[:n]>` and progress messages for further diagnosis. |
| H production SVD | `R/reference_golub_kahan.R`, `R/solve_svd.R`, `R/validation.R`, `inst/benchmarks/bench-svd-surface.R`, `tests/testthat/test-svd-adversarial.R`, `benchmarks/RELEASES.md` | SVD surface gates against certified references | Not complete. The SVD surface script supports stable `--cases=` ids and progress messages. Native CSC Gram special cases now have fast-result construction on both left- and right-Gram sides, so the tall CSC path avoids R result/certificate assembly. This is useful overhead reduction, not promotion: fresh installed 2026-05-17 3-iteration quick probes still miss the speed gate (`0.40x` on `tall_sparse:600x90`, `0.63x` on `wide_sparse:90x600`) while certifying and passing memory. The retained block-GK H candidate also certifies `tall_sparse:600x90` but fails speed (`0.0188x`) and memory (`0.0619x`). The retained IRLBA/LBD diagnostic now reports the native retained state as `ritz_subspace_only` with `irlba_lbd_recurrence_available = FALSE`; a source-loaded wide probe certifies only after adaptive fallback (`24` scout + `24` retained + `90` fallback matvecs), so true augmented recurrence remains open. A source-loaded normal-scout IRLBA diagnostic runs 8/12/16/20 matrix-free normal scouts and trusts only the final two-sided SVD certificate; the H-shaped wide row certifies but is much slower and higher-allocation than direct one-sided GK and RSpectra, so it is explicit rejection evidence rather than a promotion path. H remains unpromoted. |
| I randomized SVD | `R/reference_golub_kahan.R`, `inst/benchmarks/bench-randomized-rsvd.R`, `docs/hegelsvd_svd_acceleration.md`, `benchmarks/RELEASES.md` | randomized/rsvd accuracy and time-to-certified-answer gates | Not complete. The randomized-rsvd script supports stable `--cases=` ids and progress messages, and its gate now requires the `rsvd` baseline itself to satisfy eigencore certification before parity can pass. Fresh installed evidence separates one green regime from the broader open gate: non-quick `exact_low_rank_dense:2000x500` certifies against certified `rsvd` and passes speed/accuracy at `2.97x`; quick `exact_low_rank_dense:120x80` certifies but reaches only `1.85x`, below the `2x` gate; quick `low_rank_sparse:140x90` certifies and beats certified `rsvd` but is still red at `1.87x`; `slow_decay_dense:140x90` certifies eigencore while the `rsvd` baseline fails eigencore certification and therefore records `baseline_certified = FALSE`, `speed_gate = FALSE`, and `passed = FALSE`. Broader slow-decay/native sketch-projection gates remain open. |
| J generalized SPD LOBPCG | `R/reference_lobpcg.R`, `inst/benchmarks/bench-generalized-lobpcg.R`, `tests/testthat/test-generalized-lobpcg.R` | focused installed/source generalized benchmark rows plus full strict gate before promotion | Not complete. Dense generalized `auto()` is demoted to the native dense LAPACK fallback until iterative gates pass. Fresh installed focused evidence on 2026-05-17 shows the sparse-smallest shifted-tridiagonal native LOBPCG row certifies `10/10` and passes speed/memory versus dense base (`1.04x` speed, `2.79x` memory), while dense auto fallback rows certify as non-performance-gated boundary checks. Current installed sparse-largest evidence uses a non-densifying largest-target shift policy and certifies `10/10` with native shifted-tridiagonal preconditioner provenance, but remains performance-red on speed (`0.49x` speed, `2.75x` memory versus dense base, `163` iterations). J remains open because sparse-largest and broader generalized LOBPCG production gates are still red. |
| K B-orthogonal Lanczos refinement | `R/reference_generalized_lanczos.R`, `R/problem.R`, `R/solve_eigen.R`, `inst/benchmarks/bench-generalized-lobpcg.R`, `tests/testthat/test-generalized-lobpcg.R` | focused generalized Lanczos adversarial/agreement tests; installed focused benchmark rows; full suite | Partial reference. Explicit generalized-SPD `method = lanczos()` requests now use an honest `reference generalized SPD B-orthogonal Lanczos refinement` label when `B` has a dense, diagonal, or CSC SPD solve. Focused tests cover diagonal, dense, and sparse CSC SPD metrics, B-orthogonality, original generalized residual certificates, and agreement with LOBPCG certificates. Fresh installed focused quick strict benchmark rows on 2026-05-17 certify `diagonal_generalized_lanczos_ref_smallest` and `sparse_csc_generalized_lanczos_ref_smallest`, pass the generalized-Lanczos reference contract, and record diagonal / sparse-Cholesky metric-solve provenance. Native/block production promotion remains open. |
| L shift-invert | `R/transform_shift_invert.R`, `src/scalar_krylov.cpp`, `src/native_operators.cpp`, `tests/testthat/test-shift-invert.R`, `inst/benchmarks/bench-shift-invert.R` | shift-invert tests and benchmark evidence | Partial native. Dense standard, dense generalized SPD, diagonal standard, sparse symmetric-tridiagonal standard, and sparse/diagonal tridiagonal generalized shift-invert are native; general sparse standard and general sparse/diagonal generalized remain honest reference paths. Fresh installed non-quick strict evidence on 2026-05-17 certifies all five native benchmark rows in original coordinates, including `tridiagonal_thomas_generalized_native`. A fresh installed quick strict saved run also covers general sparse standard and sparse/diagonal generalized reference rows: both retain sparse-LU cache provenance, converge in original coordinates, keep nonnative labels, and correctly withhold `passed` because certificate scale is estimated. The `matrix_free_user_solve_reference` quick strict row now covers the user-supplied solve boundary: it certifies with exact scale metadata, records `factorization = user_solve` and `external_cache = TRUE`, and stays nonnative/reference-labelled. L remains partial because general sparse native shift-invert factorization ownership is not implemented. |
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

- H: production SVD is not promoted. Current tiny sparse Gram rows certify and
  pass memory but are speed-red on fresh installed evidence, while general
  sparse and matrix-free SVD still depend on staging paths and
  benchmark-negative retained restart work.
- I: randomized SVD remains a reference/prototype implementation outside the
  validated exact-low-rank dense regime.
- J: generalized SPD LOBPCG is not production-promoted. Dense auto now uses the
  native dense generalized LAPACK fallback, and sparse-smallest plus
  sparse-largest shifted-tridiagonal slices have focused certification
  evidence, but sparse-largest performance remains red and broader generalized
  production gates still fail.
- K: B-orthogonal Lanczos now has an explicit reference scalar refinement path
  for dense/diagonal/CSC SPD metric solves plus focused installed diagonal and
  sparse-CSC benchmark evidence, but native/block production promotion remains
  open.
- L: general sparse/native shift-invert factorization ownership is not complete
  beyond the native diagonal/tridiagonal standard and tridiagonal generalized
  slices. The general sparse standard, sparse/diagonal generalized, and
  user-supplied solve reference boundaries now have focused installed benchmark
  contracts with explicit cache/factorization provenance.
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
- G1 has fresh red benchmark evidence on `path_laplacian:1000`; sparse block
  `auto()` promotion has been demoted to a diagnostic opt-in rather than a
  production default. The append-only projected-matrix update removes the
  largest hotspot, but the post-fix scalar and block-candidate rows are still
  below release speed/parity gates.

## Stop Rule

Do not mark the thread goal complete until:

1. Every row above is either green with fresh evidence or explicitly scoped out
   by a revised PRD.
2. The strict benchmark gates for promoted solver families have been rerun from
   a clean installed package.
3. `R CMD check --no-manual` is green from a fresh tarball.
4. Mote issues for completed gates are closed or handed off with remaining
   scope recorded.
