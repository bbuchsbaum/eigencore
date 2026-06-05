# eigencore V1 Completion Audit

Last refreshed: 2026-06-05

This is the final stop-rule checklist for the active V1 readiness goal. It is
not a release signoff. It exists to keep completion tied to evidence instead of
effort, intent, or a green proxy command.

## Objective

V1 is complete only when all of these deliverables are true at the same time:

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
| Clean package checks | `eigencore_0.0.0.9000.tar.gz` and `eigencore.Rcheck/` | `R CMD build .`; `LC_ALL=C LANG=C R CMD check --no-manual eigencore_0.0.0.9000.tar.gz` | Green in latest local runs, but must be rerun after final code, docs, and benchmark updates. |
| Full test suite | `tests/testthat/` | `Rscript -e 'pkgload::load_all("."); testthat::test_dir("tests/testthat", reporter="summary")'` | Green with four expected CRAN skips in latest local runs. |
| Diff hygiene | whole repo | `git diff --check` | Green in latest local runs. |
| Planner honesty | `R/problem.R`, `R/solve.R`, result diagnostics, solver tests | Every nonnative public path has reference/oracle/fallback labels; no production claim points to an R prototype. | Mostly covered for current paths; re-audit after any solver promotion. |
| Native operator foundation | `R/operator_algebra.R`, `src/native_operators.cpp`, `tests/testthat/test-operator-algebra.R` | Built-in dense/CSC/diagonal transforms preserve native provenance without sparse densification. | Green for explicit built-ins; matrix-free centering remains an honest callback boundary. |
| G1 Hermitian native gate | `inst/benchmarks/bench-native-hermitian-gate.R` | Clean installed strict gate and saved artifacts. | Green for the promoted structured-tridiagonal default. Installed strict `path_laplacian:1000` on 2026-06-05 certified `20/20` with the native selected tridiagonal solver, median `0.00799s`, memory `579696` bytes, `6.15x` speed versus RSpectra, `1.07x` memory versus the best reference, and `22.58x` parity versus PRIMME. |
| H production SVD | `inst/benchmarks/bench-svd-surface.R` | Installed SVD surface gate passes against certified references. | Red. Tiny-Gram rows certify and pass memory but miss speed. Retained IRLBA/LBD now has residual-augmented, cached-`A Q`, and BPRO native diagnostics that certify the H-shaped wide row directly and beat retained block-GK source-loaded, but they are still slower and higher-allocation than certified RSpectra and remain unpromoted. |
| I randomized SVD | `inst/benchmarks/bench-randomized-rsvd.R` | Randomized gate passes only against certified `rsvd` baseline rows. | Mixed, not promoted. Large exact-low-rank dense row is green; quick exact-low-rank and slow-decay/native sketch regimes remain open. |
| J generalized SPD LOBPCG | `inst/benchmarks/bench-generalized-lobpcg.R` | Sparse smallest/largest and broader generalized production gates pass with native labels. | Partial. Sparse-smallest shifted-tridiagonal row passes; sparse-largest certifies and passes memory but remains speed-red. |
| K generalized SPD Lanczos | `R/reference_generalized_lanczos.R`, generalized benchmark K rows | Native/block alternate has production benchmark evidence, or reference label remains explicit. | Partial reference only. Dense/diagonal/CSC metric solves are honestly labelled reference refinements. |
| L shift-invert | `R/transform_shift_invert.R`, `inst/benchmarks/bench-shift-invert.R` | Native labels cover claimed factorization paths; general sparse/user-solve boundaries remain honest. | Partial native. Dense, diagonal, tridiagonal, and tridiagonal generalized native slices are covered; general sparse native factorization remains open. |
| Nonsymmetric eigen | `R/reference_arnoldi.R`, `inst/benchmarks/bench-nonsymmetric.R` | Dense/native-CSC Arnoldi compatibility and matrix-free reference labels are distinguished. | Compatibility-grade green for current labels; production-grade fully native restarted Arnoldi remains open. |
| Certificates and result contracts | `R/certification.R`, `tests/testthat/test-validation.R`, `tests/testthat/test-result-contracts.R` | Estimated-scale certificates cannot mark `passed`; diagnostics/accessors are stable across current families. | Covered for current paths; re-audit after new certificate constructors or solver promotions. |
| Benchmark manifest | `docs/v1-benchmark-manifest.md` | Every release surface maps to installed commands, saved artifact names, and current gate status. | Present but not signoff because red rows remain. |
| Documentation scope | `docs/v1-doc-scope-audit.md`, README, vignettes, migration docs, limitations docs | User-facing docs match the final solver surface and open limitations. | Partial. Current docs are mapped; README/vignette refresh remains blocked by final gate decisions. |
| Sanitizer / valgrind-style evidence | `inst/validation/native-smoke.R`, sanitizer logs, release notes | At least one sanitizer or valgrind-equivalent native smoke path is green, with ASan/valgrind blockers explained if local environment cannot run them. | Partial. UBSan smoke is green; ASan is blocked locally by macOS/R `dlopen` interceptor ordering and valgrind is unavailable. |
| Mote handoff | `.mote/` issue history | Completed gates are closed or handed off; remaining blockers have current evidence. | In progress. Several p1 solver gates remain open or red. |

## Completion Decision

Current decision: **not V1 ready**.

The blocking requirements are H, I, J, K native production promotion, L
general sparse native shift-invert, production-grade nonsymmetric Arnoldi,
final sanitizer/valgrind-style coverage, final benchmark artifacts, and final
README/vignette refresh. Passing tests and `R CMD check` are necessary but not
sufficient because they do not cover the red benchmark gates or release
hardening blockers.

## Stop Rule

Do not mark V1 complete until:

1. Every checklist row is green with fresh evidence or explicitly scoped out by
   a revised PRD.
2. Every promoted solver family has a clean installed-package strict benchmark
   run and saved artifact names recorded in `benchmarks/RELEASES.md`.
3. `R CMD build`, `R CMD check --no-manual`, full testthat, and
   `git diff --check` are green after the final docs and benchmark updates.
4. `mote board` shows no active claims/reservations for completed work and all
   remaining blockers are handed off with current evidence.
