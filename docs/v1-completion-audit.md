# eigencore V2 CRAN Completion Audit

Last refreshed: 2026-07-03

This is the final stop-rule checklist for the active V2 CRAN release goal. It
preserves the scoped V1 evidence base but treats completion as a CRAN release
candidate gate, not a broad solver-expansion milestone.
`V2 CRAN` is the internal release-boundary name used by the tracker and audit
files; the package version for this release candidate is `eigencore` 1.0.0.

## Objective

V2 CRAN release readiness is complete only when all of these deliverables are
true at the same time:

1. Package checks remain clean from a fresh source tarball.
2. Public solver paths either run native engine code or carry honest
   reference/oracle planner labels.
3. Native operator foundation covers the required built-in dense, CSC,
   diagonal, adjoint, centered, scaled, summed, composed, and crossprod paths
   without silent sparse densification.
4. Production SVD and randomized SVD are promoted only after benchmark gates
   pass.
5. Generalized SPD, shift-invert, and nonsymmetric surfaces have native
   coverage where claimed and honest labels where not claimed.
6. Release hardening includes tests, CRAN-like checks, installed-package
   benchmark evidence, sanitizer or valgrind-style evidence, docs, migration
   guidance, known limitations, and mote handoff.

## Prompt-to-Artifact Checklist

| Requirement | Primary artifact | Required evidence | Current state |
|---|---|---|---|
| Clean package checks | `eigencore_1.0.0.tar.gz` and fresh `eigencore.Rcheck/` output | `LC_ALL=C LANG=C R CMD build /Users/bbuchsbaum/code/eigencore`; `LC_ALL=C LANG=C R CMD check --as-cran --no-manual eigencore_1.0.0.tar.gz` | Green on 2026-07-03 via `rcmdcheck::rcmdcheck(path = ".", args = c("--as-cran", "--no-manual"), build_args = "--no-manual")` after the generalized-eigen release gate and exported-return-value docs: `0 ERROR`, `0 WARNING`, `1 NOTE` for CRAN new submission. |
| Full test suite | `tests/testthat/` | `Rscript -e 'pkgload::load_all("."); testthat::test_dir("tests/testthat", reporter="summary")'`; installed `tests/testthat.R` inside the source-tarball check | Green on 2026-07-03: the full source suite was green after the generalized-eigen gate, and the installed-package `testthat.R` pass in the fresh 1.0.0 CRAN-like check was OK. |
| Diff hygiene | whole repo | `git diff --check` | Green on 2026-07-03 after release metadata, generalized-eigen docs, benchmark, roxygen, and exported return-value documentation updates. |
| CRAN remote checks | GitHub Actions, R-hub, win-builder maintainer emails | GitHub Actions for the release candidate; `rhub::rhub_check(platforms = c("ubuntu-release", "windows"))`; win-builder R-release and R-devel submissions with `--no-manual` | GitHub Actions are green for commit 088d3bd: R-CMD-check, lint, coverage, and pkgdown all completed successfully. R-hub run 28686224852 completed successfully for Ubuntu release and Windows. win-builder R-release and R-devel submissions from 088d3bd were accepted on 2026-07-03; maintainer email results are still required before closing the CRAN cross-platform gate. |
| Planner honesty | `R/problem.R`, `R/solve.R`, result diagnostics, solver tests | Every nonnative public path has reference/oracle/fallback labels; no production claim points to an R prototype. | Mostly covered for current paths; re-audit after any solver promotion. |
| Native operator foundation | `R/operator_algebra.R`, `src/native_operators.cpp`, `tests/testthat/test-operator-algebra.R` | Built-in dense/CSC/diagonal transforms preserve native provenance without sparse densification. | Green for explicit built-ins; matrix-free centering remains an honest callback boundary. |
| G1 Hermitian native gate | `inst/benchmarks/bench-native-hermitian-gate.R` | Clean installed strict gate and saved artifacts. | Green for the promoted structured-tridiagonal default. Installed strict `path_laplacian:1000` on 2026-06-05 certified `20/20` with the native selected tridiagonal solver, median `0.00799s`, memory `579696` bytes, `6.15x` speed versus RSpectra, `1.07x` memory versus the best reference, and `22.58x` parity versus PRIMME. |
| H production SVD | `inst/benchmarks/bench-svd-surface.R` | Installed SVD surface gate and target contracts pass against certified references. | Green for the promoted non-quick tall/wide sparse H surface plus sparse smallest/interior target contracts. Fresh installed 2026-06-05 strict evidence certifies `20/20` triplets on `tall_sparse` (`100000 x 500`) and `wide_sparse` (`500 x 100000`), uses the bounded native Gram special case with `native_gram_eigensolver = "lapack_dsyevr"` and `materialized_gram = TRUE`, and passes speed/memory (`3.849574x` / `2.049985x` tall, `9.882995x` / `2.049985x` wide). Fresh 2026-06-06 source strict target evidence certifies `smallest_sparse:5000x500` through the native certified Gram special case and `interior_sparse:5000x500` through native full-subspace Golub-Kahan, with target-contract speed/memory gates green versus dense `base_*` rows. The quick `600x90` fixture is diagnostic only after a final strict rerun proved too noisy for release signoff. The `complex_dense` benchmark row covers the native dense complex SVD label outside the sparse H gate. |
| I randomized SVD | `inst/benchmarks/bench-randomized-rsvd.R` | Scoped randomized release gate passes only against certified `rsvd` baseline rows tagged `release_gate_required = TRUE`; native controller contracts also prove dense/sparse controller provenance. | Green for the scoped dense release row and sparse native-controller contract. Fresh 2026-06-06 strict source evidence certifies `exact_low_rank_dense:2000x500` against certified `rsvd` at about `2.55x`, with the dense QR native randomized controller, native fused sketch/projection, native projected-core SVD, native certificate diagnostics, and q=0 early stop active. The same strict run certifies `low_rank_sparse:2000x500` through `native_csc_randomized_controller` with `randomized_sparse_native_controller = TRUE` and a green controller contract; the sparse speed ratio is diagnostic, not a promoted release win. Quick small-size rows and uncertified-baseline slow/near-low-rank rows remain visible diagnostics; matrix-free and LU/none randomized control remain reference-labelled. |
| J generalized SPD LOBPCG | `inst/benchmarks/bench-generalized-lobpcg.R` | Sparse smallest/largest production gate passes with native generalized/B-orthogonal/shifted-tridiagonal provenance, while fallback/reference rows remain honestly labelled. | Green for the promoted sparse shifted-tridiagonal surface. Installed 2026-06-05 full strict saved evidence certifies `10/10` for both sparse-smallest and sparse-largest at `n = 1000`, passes speed/memory versus dense base, and keeps native, generalized-Lanczos-reference, constrained, matrix-free-B, and adversarial-B contract rows green. Dense generalized `auto()` remains a native dense LAPACK fallback boundary rather than a broader iterative promotion. Standard Hermitian LOBPCG/preconditioner promotion is closed as diagnostic/prototype-only because the native shifted-tridiagonal row certifies and beats external references but loses to the current eigencore default. |
| K generalized SPD Lanczos | `R/reference_generalized_lanczos.R`, generalized benchmark K rows | Native dense/diagonal transformed labels, explicit block transformed rows, native sparse tridiagonal metric-solve provenance, and general sparse-CSC reference labels remain explicit and pass B-orthogonality and certification contracts. | Green for the dense/diagonal native transformed slice, including the explicit block transformed row, plus sparse CSC metric-boundary diagnostics. The focused strict probe certifies `diagonal_generalized_block_lanczos_native_smallest` with `block_size = 2` and a green block-native contract. The tridiagonal sparse-B row uses an eigencore-owned native Thomas solve inside the reference-labelled Lanczos refinement; the non-tridiagonal sparse-CSC row remains `Matrix::Cholesky` reference provenance. Sparse-CSC block production Lanczos and arbitrary sparse-CSC native factorization are not claimed. |
| L shift-invert | `R/transform_shift_invert.R`, `inst/benchmarks/bench-shift-invert.R` | Native labels cover claimed factorization paths; general sparse/user-solve boundaries remain honest. | Green for the scoped V1 surface. `auto()` plus `target = nearest(sigma)` now plans and dispatches through an implicit `shift_invert(sigma)` transform for supported factorized Hermitian regimes and rejects matrix-free/unfactorized inputs without a solve. Installed 2026-06-05 strict saved evidence passes all eight native/reference contract rows, with dense/diagonal/tridiagonal native factorized Lanczos rows certifying in original coordinates and general sparse/user-solve rows retaining honest reference/cache provenance. Native general sparse LU and native ownership of user-supplied solve functions are explicit PRD non-goals unless a future PRD reopens them. |
| Nonsymmetric eigen | `R/reference_arnoldi.R`, `inst/benchmarks/bench-nonsymmetric.R` | Dense/native-CSC refined Arnoldi, native dense complex LAPACK, and native matrix-free callback projected-Arnoldi labels are distinguished. | Green for the scoped V1 compatibility surface, with a current V2 tranche for dense/sparse CSC refined Ritz extraction. Dense/sparse installed evidence from 2026-06-05 is green for the original V1 native rows, current source tests cover the refined-Ritz provenance, and the focused installed `matrix_free_nonnormal:30` row requires native callback Arnoldi, restart diagnostics, and exact right-residual certification. The current benchmark script also records a dense complex LAPACK row with exact right-residual certification. Adjoint-capable rows expose left vectors with left-residual and biorthogonality diagnostics; full Krylov-Schur or harmonic/interior extraction, matrix-free refined extraction, and native complex sparse/operator paths remain future scope. |
| Certificates and result contracts | `R/certification.R`, `tests/testthat/test-validation.R`, `tests/testthat/test-result-contracts.R` | Estimated-scale certificates cannot mark `passed`; diagnostics/accessors are stable across current families. | Covered for current paths; re-audit after new certificate constructors or solver promotions. |
| Benchmark manifest | `docs/v1-benchmark-manifest.md` | Every release surface maps to installed commands, saved artifact names, and current gate status. | Green for promoted solver gates. The 2026-06-07 all-surface quick-smoke sweep passed `post_v1_svd_hard_surface`, `post_v1_operator_sidecars`, `post_v1_randomized_svd_hard_surface`, `post_v1_generalized_preconditioned_surface`, `post_v1_shift_invert_boundaries`, and `post_v1_nonsymmetric_matrix_free_surface`; installed strict evidence remains the release signoff surface. |
| Documentation scope | `docs/v1-doc-scope-audit.md`, README, vignettes, migration docs, limitations docs | User-facing docs match the final solver surface and open limitations. | Green for scoped V1 after README/vignette refresh and release-boundary docs. Future solver promotions must re-open this row. |
| Sanitizer / valgrind-style evidence | `inst/validation/native-smoke.R`, sanitizer logs, release notes | At least one sanitizer or valgrind-equivalent native smoke path is green, with ASan/valgrind blockers explained if local environment cannot run them. | Scoped green locally. Final 2026-06-05 UBSan install plus installed native smoke passed; ASan is blocked locally by macOS/R `dlopen` interceptor ordering and valgrind is unavailable, both documented as environment boundaries. |
| Mote handoff | `.mote/` issue history | Completed gates are closed or handed off; remaining blockers have current evidence. | Final H and release-hardening close operations are recorded in mote history after validation. |

## Completion Decision

Current decision: **V2 CRAN release candidate pending final validation and mote closure**.

As of 2026-07-03, the remaining final-validation item is external win-builder
email evidence for the R-release and R-devel submissions from commit 088d3bd.
GitHub Actions and R-hub have already passed for that candidate.

The future solver families remain listed as V3 deferrals; they are not hidden
release blockers for the scoped CRAN surface.

## Stop Rule

Mark V2 CRAN complete only after:

1. Every checklist row is green with fresh evidence or explicitly scoped out by
   a revised PRD.
2. Every promoted solver family has a clean installed-package strict benchmark
   run and saved artifact names recorded in `benchmarks/RELEASES.md`.
3. `R CMD build`, `R CMD check --no-manual`, full testthat, and
   `git diff --check` are green after the final docs and benchmark updates.
4. win-builder R-release and R-devel maintainer emails are confirmed clean for
   the same release candidate.
5. `mote board` shows no active claims/reservations for completed work and all
   remaining blockers are handed off with current evidence.
