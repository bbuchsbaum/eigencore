# V1 Documentation Scope Audit

Date: 2026-05-17

This audit maps the V1 documentation requirement to concrete files.

This audit is not a release signoff. It does not override the solver gates or
final benchmark artifacts named in `docs/v1-readiness-audit.md` and
`docs/v1-benchmark-manifest.md`.

## Source Inventory

| Area | Current files | Release role |
|---|---|---|
| README and quick start | `README.md`, `README.Rmd` | First-contact examples, planner labels, certificate shape, and current public surface. |
| Vignettes | `vignettes/eigencore.Rmd`, `vignettes/certificates.Rmd` | Longer usage path and certificate/diagnostic explanation exercised by package build/check. |
| Migration and limits | `docs/rspectra-migration.md`, `docs/known-limitations.md` | RSpectra compatibility contract, unsupported paths, and honest fallback labels. |
| Method selection | `docs/method-selection-and-workflows.md` | User-facing guidance for eigen, SVD, generalized, shift-invert, operator, and certificate workflows. |
| Release evidence | `docs/v1-readiness-audit.md`, `docs/v1-benchmark-manifest.md`, `benchmarks/RELEASES.md` | Gate inventory, installed benchmark commands, saved-artifact expectations, and current blocker record. |
| Native design notes | `docs/native-lobpcg.md`, `docs/native-generalized-spd-lobpcg.md`, `docs/native-block-lanczos.md`, `docs/hegelsvd_svd_acceleration.md` | Design contracts and implementation notes for native solver work. |

## Coverage Matrix

| Requirement | Documentation coverage | Current status |
|---|---|---|
| Quick start and first-contact API | `README.md`, `README.Rmd`, `vignettes/eigencore.Rmd` | Covered for current public paths. Refresh after any solver promotion or public label change. |
| Certificates and diagnostics | `vignettes/certificates.Rmd`, `docs/method-selection-and-workflows.md`, `docs/v1-readiness-audit.md` | Covered for current certificate fields and planner diagnostics. Re-audit after new certificate constructors or estimated-scale paths. |
| RSpectra migration | `docs/rspectra-migration.md`, `docs/known-limitations.md` | Covered for current shim behavior. Revisit after H, I, J, L, and native Arnoldi decisions. |
| Partial Hermitian eigen | `README.md`, `docs/method-selection-and-workflows.md`, `docs/native-block-lanczos.md`, `docs/v1-readiness-audit.md` | Covered with the green structured-tridiagonal G1 default and explicit block-Lanczos diagnostic caveat. |
| Partial SVD and randomized SVD | `README.md`, `docs/method-selection-and-workflows.md`, `docs/hegelsvd_svd_acceleration.md`, `docs/v1-benchmark-manifest.md` | Covered with the green H promoted surface and scoped I release gate. Final examples must avoid broad randomized performance claims beyond the exact-low-rank release row while public randomized control remains reference-labelled. |
| Generalized SPD and LOBPCG | `docs/method-selection-and-workflows.md`, `docs/native-generalized-spd-lobpcg.md`, `docs/native-lobpcg.md`, `docs/v1-readiness-audit.md` | Covered for current native slices; final production wording is blocked by J/K. |
| Shift-invert | `docs/method-selection-and-workflows.md`, `docs/rspectra-migration.md`, `docs/v1-readiness-audit.md`, `docs/v1-benchmark-manifest.md` | Covered for the scoped V1 surface. Dense standard/generalized, diagonal and symmetric-tridiagonal standard, and tridiagonal generalized-with-diagonal-B native labels are documented; general sparse standard and general sparse diagonal-generalized remain reference-labelled as the honest V1 boundary. |
| Nonsymmetric eigen | `docs/method-selection-and-workflows.md`, `docs/known-limitations.md`, `docs/v1-readiness-audit.md`, `docs/v1-benchmark-manifest.md` | Partial. Dense oracle, sparse native Arnoldi-cycle with restart controls, and matrix-free reference Arnoldi labels are documented; production-grade fully native restarted Arnoldi remains open. |
| Operator algebra and no densification | `docs/method-selection-and-workflows.md`, `docs/v1-readiness-audit.md` | Covered for explicit built-ins; matrix-free centering remains callback-boundary policy. |
| Benchmark and release reports | `docs/v1-benchmark-manifest.md`, `benchmarks/RELEASES.md`, `docs/v1-readiness-audit.md` | Partial. Final strict installed-package artifacts are still required. |
| Sanitizer and valgrind-style checks | `inst/validation/native-smoke.R`, `docs/v1-readiness-audit.md`, `benchmarks/RELEASES.md`, `docs/known-limitations.md` | Partial. UBSan smoke is recorded through a reusable installed-package smoke artifact; ASan/valgrind-equivalent release coverage is not green in this environment. |

## Blocking Documentation Gaps

- Final README and vignette examples must be refreshed after H, I, J, L, and
  native Arnoldi decisions, because those gates can change planner labels and
  recommended methods.
- Final benchmark reports must list saved artifacts from strict installed
  runs, not only source-loaded or quick diagnostic runs.
- The migration guide must be revisited after solver promotions or demotions
  so RSpectra-compatible behavior and warnings match the actual planner.
- The package currently relies on R package build/check vignette rendering for
  release documentation verification. There is no separate pkgdown/site signoff
  artifact in this audit.
- ASan or valgrind-equivalent native coverage is not green locally. The reusable
  native smoke script keeps the UBSan/ASan/valgrind command surface stable, but
  do not present sanitizer coverage as complete until that row is resolved in a
  suitable environment.

## Verification Commands

Use these checks whenever this audit is changed:

```sh
Rscript -e 'pkgload::load_all("."); testthat::test_file("tests/testthat/test-bench-smoke.R", reporter = "summary")'
Rscript -e 'pkgload::load_all("."); testthat::test_dir("tests/testthat", reporter = "summary")'
git diff --check
R CMD build .
LC_ALL=C LANG=C R CMD check --no-manual eigencore_0.0.0.9000.tar.gz
Rscript inst/validation/native-smoke.R --load-all
```

The final V1 benchmark commands remain the installed-package commands in
`docs/v1-benchmark-manifest.md`.

## Stop Rule

Do not treat this documentation audit as release signoff. V1 remains open until
every linked red solver and release-hardening row is either green with fresh
evidence or explicitly scoped out by a revised PRD, and the final strict
installed-package benchmark artifacts are recorded in `benchmarks/RELEASES.md`.
