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
| H production SVD | `inst/benchmarks/bench-svd-surface.R` | Installed SVD surface gate passes against certified references. | Green for the promoted tall/wide sparse H surface. The executable `--h-candidate` gate targets the promoted `eigencore` SVD path, with retained BPRO/block-GK rows diagnostic. Warning-free installed 3-iteration tall/wide sparse probes certify all five requested triplets and pass speed/memory; tall uses native right-normal `implicit_normal_lanczos` without materializing the Gram, and wide uses the certified native left-Gram path (`1.194x` / `2.392x` tall, `1.223x` / `2.392x` wide). |
| I randomized SVD | `inst/benchmarks/bench-randomized-rsvd.R` | Scoped randomized gate passes only against certified `rsvd` baseline rows tagged `release_gate_required = TRUE`. | Green for the scoped V1 release row. Fresh strict installed evidence certifies `exact_low_rank_dense:2000x500` against certified `rsvd` at about `3.1x`, with native fused sketch/projection diagnostics active. Quick small-size rows and uncertified-baseline slow/near-low-rank rows remain visible diagnostics and are not promoted planner wins. |
| J generalized SPD LOBPCG | `inst/benchmarks/bench-generalized-lobpcg.R` | Sparse smallest/largest production gate passes with native generalized/B-orthogonal/shifted-tridiagonal provenance, while fallback/reference rows remain honestly labelled. | Green for the promoted sparse shifted-tridiagonal surface. Installed 2026-06-05 full strict saved evidence certifies `10/10` for both sparse-smallest and sparse-largest at `n = 1000`, passes speed/memory versus dense base, and keeps native, generalized-Lanczos-reference, constrained, matrix-free-B, and adversarial-B contract rows green. Dense generalized `auto()` remains a native dense LAPACK fallback boundary rather than a broader iterative promotion. |
| K generalized SPD Lanczos | `R/reference_generalized_lanczos.R`, generalized benchmark K rows | Reference label remains explicit and passes B-orthogonality, certification, and LOBPCG-agreement contracts. | Green for the scoped reference-refinement gate. Fresh 2026-06-05 saved and focused installed evidence certifies diagonal and sparse-CSC metric solves; native/block production Lanczos is future scope. |
| L shift-invert | `R/transform_shift_invert.R`, `inst/benchmarks/bench-shift-invert.R` | Native labels cover claimed factorization paths; general sparse/user-solve boundaries remain honest. | Green for the scoped V1 surface. Installed 2026-06-05 strict saved evidence passes all eight native/reference contract rows, with dense/diagonal/tridiagonal native factorized Lanczos rows certifying in original coordinates and general sparse/user-solve rows retaining honest reference/cache provenance. |
| Nonsymmetric eigen | `R/reference_arnoldi.R`, `inst/benchmarks/bench-nonsymmetric.R` | Dense/native-CSC Arnoldi compatibility and matrix-free reference labels are distinguished. | Green for the scoped V1 compatibility surface on fresh 2026-06-05 installed evidence; fully restarted matrix-free native Arnoldi is future scope. |
| Certificates and result contracts | `R/certification.R`, `tests/testthat/test-validation.R`, `tests/testthat/test-result-contracts.R` | Estimated-scale certificates cannot mark `passed`; diagnostics/accessors are stable across current families. | Covered for current paths; re-audit after new certificate constructors or solver promotions. |
| Benchmark manifest | `docs/v1-benchmark-manifest.md` | Every release surface maps to installed commands, saved artifact names, and current gate status. | Present for solver gates; final aggregate signoff still pending. |
| Documentation scope | `docs/v1-doc-scope-audit.md`, README, vignettes, migration docs, limitations docs | User-facing docs match the final solver surface and open limitations. | Partial. Current docs are mapped; README/vignette refresh remains blocked by final gate decisions. |
| Sanitizer / valgrind-style evidence | `inst/validation/native-smoke.R`, sanitizer logs, release notes | At least one sanitizer or valgrind-equivalent native smoke path is green, with ASan/valgrind blockers explained if local environment cannot run them. | Partial. UBSan smoke is green; ASan is blocked locally by macOS/R `dlopen` interceptor ordering and valgrind is unavailable. |
| Mote handoff | `.mote/` issue history | Completed gates are closed or handed off; remaining blockers have current evidence. | In progress. Several p1 solver gates remain open or red. |

## Completion Decision

Current decision: **not V1 ready**.

The blocking requirements are final sanitizer/valgrind-style coverage, final
aggregate benchmark artifacts, and final README/vignette refresh. Passing
tests and `R CMD check` are necessary but not sufficient because they do not
cover the remaining release-hardening blockers.

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
