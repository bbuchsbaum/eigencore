# V2 CRAN Documentation Scope Audit

Date: 2026-06-07

This audit maps the V2 CRAN documentation requirement to concrete files. The
file name is retained because existing benchmark smoke tests and release notes
refer to it as the documentation-scope artifact.

This audit is the documentation-scope companion to the solver and benchmark
evidence named in `docs/v1-readiness-audit.md` and
`docs/v1-benchmark-manifest.md`.

## Source Inventory

| Area | Current files | Release role |
|---|---|---|
| README and quick start | `README.md`, `README.Rmd` | First-contact examples, planner labels, certificate shape, V2 CRAN scope, and current public surface. |
| Vignettes | `vignettes/eigencore.Rmd`, `vignettes/certificates.Rmd` | Longer usage path and certificate/diagnostic explanation exercised by package build/check. |
| Migration and limits | `docs/rspectra-migration.md`, `docs/known-limitations.md` | RSpectra compatibility contract, unsupported paths, and honest fallback labels. |
| Method selection | `docs/method-selection-and-workflows.md` | User-facing guidance for eigen, SVD, generalized, shift-invert, operator, and certificate workflows. |
| Release evidence | `docs/v1-readiness-audit.md`, `docs/v1-benchmark-manifest.md`, `benchmarks/RELEASES.md` | Gate inventory, installed benchmark commands, saved-artifact expectations, and current blocker record. |
| Native design notes | `docs/native-lobpcg.md`, `docs/native-generalized-spd-lobpcg.md`, `docs/native-block-lanczos.md`, `docs/hegelsvd_svd_acceleration.md` | Design contracts and implementation notes for native solver work. |

## Coverage Matrix

| Requirement | Documentation coverage | Current status |
|---|---|---|
| Quick start and first-contact API | `README.md`, `README.Rmd`, `vignettes/eigencore.Rmd` | Covered for the scoped V2 CRAN public paths. Refresh after any future solver promotion or public label change. |
| Certificates and diagnostics | `vignettes/certificates.Rmd`, `docs/method-selection-and-workflows.md`, `docs/v1-readiness-audit.md` | Covered for current certificate fields and planner diagnostics. Re-audit after new certificate constructors or estimated-scale paths. |
| RSpectra migration | `docs/rspectra-migration.md`, `docs/known-limitations.md` | Covered for current shim behavior and V1 release boundaries. Revisit after future public label changes. |
| Partial Hermitian eigen | `README.md`, `docs/method-selection-and-workflows.md`, `docs/native-block-lanczos.md`, `docs/v1-readiness-audit.md` | Covered with the green structured-tridiagonal G1 default and explicit block-Lanczos diagnostic caveat. |
| Partial SVD and randomized SVD | `README.md`, `docs/method-selection-and-workflows.md`, `docs/hegelsvd_svd_acceleration.md`, `docs/v1-benchmark-manifest.md` | Covered with the green H promoted surface, scoped I dense release gate, and sparse randomized controller contract. Final examples must avoid broad randomized performance claims beyond the exact-low-rank release row. Dense QR randomized control is native for largest-target dense double requests, and sparse CSC QR randomized control is native with a certificate/provenance contract; sparse speed, matrix-free, LU/none, slow-decay, and broader adaptive randomized regimes remain diagnostic or reference-labelled. |
| Generalized SPD and LOBPCG | `docs/method-selection-and-workflows.md`, `docs/native-generalized-spd-lobpcg.md`, `docs/native-lobpcg.md`, `docs/v1-readiness-audit.md` | Covered for current native LOBPCG slices, native dense/diagonal transformed generalized Lanczos, the explicit block transformed generalized Lanczos contract, and the sparse-CSC generalized Lanczos reference boundary; sparse-CSC block promotion remains unclaimed. |
| Shift-invert | `docs/method-selection-and-workflows.md`, `docs/rspectra-migration.md`, `docs/v1-readiness-audit.md`, `docs/v1-benchmark-manifest.md` | Covered for the scoped V1 surface. Dense standard/generalized, diagonal and symmetric-tridiagonal standard, and tridiagonal generalized-with-diagonal-B native labels are documented; general sparse standard and general sparse diagonal-generalized remain reference-labelled as the honest V1 boundary. |
| Nonsymmetric eigen | `docs/method-selection-and-workflows.md`, `docs/known-limitations.md`, `docs/v1-readiness-audit.md`, `docs/v1-benchmark-manifest.md` | Covered for the scoped V1 compatibility surface plus the dense/sparse CSC refined-Ritz tranche. Dense and sparse CSC operators use native Arnoldi-cycle/native-refined-Ritz labels with restart controls; real matrix-free callback operators use the native callback Arnoldi/native-projected-Ritz label. Adjoint-capable rows expose left vectors with separate left-residual and biorthogonality diagnostics. Base complex dense inputs use native dense complex LAPACK labels; full Krylov-Schur or harmonic/interior extraction, matrix-free refined extraction, and native complex sparse/operator inputs remain future scope. |
| Operator algebra and no densification | `docs/method-selection-and-workflows.md`, `docs/v1-readiness-audit.md` | Covered for explicit built-ins; matrix-free centering remains callback-boundary policy. |
| Benchmark and release reports | `docs/v1-benchmark-manifest.md`, `benchmarks/RELEASES.md`, `docs/v1-readiness-audit.md` | Covered for promoted solver surfaces with final installed artifacts named. |
| Sanitizer and valgrind-style checks | `inst/validation/native-smoke.R`, `docs/v1-readiness-audit.md`, `benchmarks/RELEASES.md`, `docs/known-limitations.md` | Scoped locally. UBSan smoke is recorded through a reusable installed-package smoke artifact; ASan and valgrind are unavailable in this environment and documented as local boundaries. |

## Documentation Boundaries

- README and vignette examples match the current scoped V2 CRAN surface. Reopen this
  audit after any solver promotion or demotion that changes planner labels or
  recommended methods.
- Benchmark reports list saved artifacts from strict installed runs for the
  promoted solver surfaces; quick diagnostic runs are not release signoff.
- The migration guide is current for the scoped V1 planner. Revisit it after
  future promotions so RSpectra-compatible behavior and warnings stay current.
- The package currently relies on R package build/check vignette rendering for
  release documentation verification. There is no separate pkgdown/site signoff
  artifact in this audit.
- ASan or valgrind-equivalent native coverage is not available locally. The
  reusable native smoke script keeps the UBSan/ASan/valgrind command surface
  stable; current release docs present UBSan as the local sanitizer evidence
  and ASan/valgrind as environment boundaries.

## Verification Commands

Use these checks whenever this audit is changed:

```sh
Rscript -e 'pkgload::load_all("."); testthat::test_file("tests/testthat/test-bench-smoke.R", reporter = "summary")'
Rscript -e 'pkgload::load_all("."); testthat::test_dir("tests/testthat", reporter = "summary")'
git diff --check
R CMD build .
LC_ALL=C LANG=C R CMD check --no-manual eigencore_0.2.0.tar.gz
Rscript inst/validation/native-smoke.R --load-all
```

The final V1 benchmark commands remain the installed-package commands in
`docs/v1-benchmark-manifest.md`.

## Stop Rule

Treat this documentation audit as current only for the scoped V2 CRAN surface. Future
solver promotions must update README/vignettes, migration docs, the benchmark
manifest, release notes, and this audit in the same slice.
