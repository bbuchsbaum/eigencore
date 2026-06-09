# Post-V2 / V3 Benchmark Gates

Date: 2026-06-07

This is the benchmark truth surface for promoting eigencore beyond the V2 CRAN
release boundaries. The structured source is
`inst/benchmarks/post-v1-gate-manifest.R`; this document explains how to use it.
The contribution-facing interpretation is
`docs/contribution-methods-artifact.md`, which summarizes the current wins,
losses, and non-promoted solver boundaries.

The rule is simple: a solver surface is not promoted just because a prototype
exists. It needs installed-package benchmark evidence, exact or explicitly
qualified certificates, memory evidence, planner-label provenance, saved
artifacts, and a current owner in the structured gate manifest. When a
promotion program is closed as a no-promotion or non-goal decision, its
regression gate becomes umbrella-owned until a revised V3 PRD creates a new
owner.

## Installed-Package Pattern

Run post-release/V3 gates from an installed package:

```sh
R CMD INSTALL --library=/tmp/eigencore-bench-lib .
R_LIBS=/tmp/eigencore-bench-lib Rscript inst/benchmarks/<script>.R --strict --save
```

Saved artifacts belong under `inst/benchmarks/results/` with names of the form
`YYYYMMDD-<gate>-<artifact>.rds`. Curated release/contribution evidence should
also be summarized in `benchmarks/RELEASES.md`.

## Gates

| Gate | Script | Required regimes | Baselines | Promotion threshold |
|---|---|---|---|---|
| Hard SVD surface | `inst/benchmarks/bench-svd-surface.R` | tall/wide sparse, rank-deficient sparse, clustered dense, slow-decay dense, low-rank sparse | RSpectra, PRIMME, irlba, base LAPACK for small truth | Future-promotion gate only: certified requested rank, exact two-sided SVD residuals, no normal-equation default unless explicitly labelled, at least `1.10x` time-to-certified-answer, no worse memory than the best certified reference. Current PRD broad thick-restart promotion is closed as diagnostic-only. |
| Operator sidecars | `inst/benchmarks/bench-post-v1-operator-sidecars.R` | matrix-free SVD, matrix-free nonsymmetric, matrix-free generalized B | Current planner labels and exact certificates | Certificate passes, planner label exactly matches expected native/reference boundary, native provenance matches the documented surface. |
| Randomized SVD hard surface | `inst/benchmarks/bench-randomized-rsvd.R` | exact low rank, nearly low rank, slow decay, sparse low rank | rsvd, irlba | Certified baseline required, no worse accuracy/subspace error, at least `2.00x` time-to-certified-answer where randomized is planner-promoted. |
| Generalized/preconditioned SPD | `inst/benchmarks/bench-generalized-lobpcg.R` | sparse largest/smallest, `adversarial_explicit_spd_matrix_free_b_smallest`, ill-conditioned diagonal B, sparse CSC B | base dense generalized solve, generalized Lanczos reference | Certificate and B-orthogonality pass, no sparse densification, at least `1.25x` speed and `4.00x` memory versus dense baseline for promoted sparse rows. |
| Shift-invert boundaries | `inst/benchmarks/bench-shift-invert.R` | `sparse_general_reference`, `sparse_general_diagonal_b_reference`, `matrix_free_user_solve_reference` | Matrix::lu, user solve, dense base small truth | Original-coordinate certificate passes, `shift_invert_factorization_contract_v1` metadata is present, and native labels are allowed only where the contract provider is `eigencore_native_factorization`. |
| Nonsymmetric Arnoldi | `inst/benchmarks/bench-nonsymmetric.R` plus operator sidecar | matrix-free nonnormal, dense native Arnoldi, sparse native Arnoldi | RSpectra, base LAPACK small truth, current native callback Arnoldi | Boundary-regression gate for the current Arnoldi surface: right-residual certificate passes, restart diagnostics are present, dense/sparse CSC rows carry native refined-Ritz provenance, and matrix-free real-operator rows carry the native callback projected-Ritz label. This is not a full Krylov-Schur, harmonic/interior, or matrix-free refined-extraction promotion gate. |

## Closed No-Promotion Decisions

- General sparse block Hermitian Lanczos `auto()` promotion is closed under
  `bd-01KTEH48HH4X8G9Q69HHAJ983B`. The explicit block candidate remains
  runnable through `bench-native-hermitian-gate.R --block-candidate`, but the
  current non-quick `path_laplacian:1000`, `k = 20` gate certifies while failing
  speed and PRIMME parity. Keep sparse block auto-promotion disabled unless a
  future redesign supplies fresh installed evidence that beats certified
  references and the current promoted eigencore path.
- Native standard Hermitian LOBPCG/preconditioner promotion is closed under
  `bd-01KTEH4G1QPR4RT14B4G78PF1M`. The shifted-tridiagonal native LOBPCG row
  certifies and beats certified external references on path-Laplacian probes,
  but it fails the scalar-speed gate against the current eigencore default at
  `n = 200`, `1000`, and `2000`. Keep standard Hermitian LOBPCG labelled as a
  diagnostic/prototype path unless a future redesign beats the current default.
- Broader randomized SVD promotion is closed under
  `bd-01KTE8JFKPA90ZJTXK496SBMK4` for the current PRD. Dense exact-low-rank QR
  requests and sparse CSC QR requests now have native controller boundaries, but
  sparse speed, slow-decay, matrix-free, and LU/none regimes remain diagnostic
  or reference-control unless a future PRD reopens them.
- Sparse-CSC block generalized B-orthogonal Lanczos and arbitrary sparse-CSC
  native metric factorization are closed under
  `bd-01KTE8JVEPGA1EEQYERZS1V7S1` as current-PRD no-promotion boundaries. The
  covered contribution is native dense/diagonal transformed generalized
  Lanczos, including the focused block row, plus honest sparse metric
  provenance.
- Broad sparse/matrix-free thick-restarted SVD promotion is closed under
  `bd-01KTE8J9SF16Y1832D8HQQ9KEC` as diagnostic-only for the current PRD. The
  retained workspace and callback-boundary slices remain useful, but hard SVD
  candidate rows are still speed or memory red against certified references.
- Scalable refined/shift-invert sparse or matrix-free smallest/interior SVD
  promotion is closed under `bd-01KTEH6862GB19JJWX2M3FQP6T` for the current PRD.
  The current contribution is exact certified smallest/interior boundaries, not
  a scalable non-full-subspace production claim.

## V3 Deferrals And Closed Non-Goal Decisions

- Native general sparse shift-invert LU/factorized apply and native ownership
  of user-supplied solve functions are closed under
  `bd-01KTEH4ZPPBESDDPR0Z5MRSQCG` as explicit PRD non-goals. Keep the
  shift-invert boundary gate active as regression evidence for
  `Matrix::lu_reference_factorization`, `user_supplied_solve`, cache
  provenance, no dense sparse-rcond fallback, original-coordinate residuals,
  and the rule that native labels require `eigencore_native_factorization`.
- Full Krylov-Schur-style, harmonic/interior, matrix-free refined, and native
  complex sparse/operator nonsymmetric extraction remain future scope under
  `bd-01KTEH5JM64A4CBZG7ECBWT9WB`. The first refined-extraction tranche,
  `bd-01KTF6H41S9XDN286TR3V184P4`, promotes native refined Ritz extraction for
  dense and sparse CSC Arnoldi paths only. Keep the nonsymmetric benchmark as a
  boundary-regression gate for that native Arnoldi compatibility surface.
- Complex sparse Matrix, native complex block-operator kernels, and promoted
  complex matrix-free solver operators are closed under
  `bd-01KTEH60X91VZRSW7NGV65FBDR` as current-PRD non-goals. Dense complex
  LAPACK paths and exact certificates remain the current promoted complex
  contribution.
- Jacobi-Davidson, Davidson, full nonsymmetric Krylov-Schur or harmonic/interior
  workflows, scalable sparse/matrix-free interior SVD, native general sparse LU
  ownership, GraphBLAS/GPU/distributed/SLEPc/PRIMME plugins, and broad matrix
  ecosystem adapters are V3 scope. They are not V2 CRAN blockers without a new
  PRD and strict installed evidence.

## Smoke Profile

The profile runner validates the structured manifest and selects gates by tier:

```sh
Rscript inst/benchmarks/run-post-v1-gates.R --tier=smoke --dry-run
Rscript inst/benchmarks/run-post-v1-gates.R --tier=smoke
Rscript inst/benchmarks/run-post-v1-gates.R --tier=smoke --load-all
```

By default the smoke tier runs `post_v1_operator_sidecars`, which is fast enough
for local and CI use. It is not a performance signoff. It proves the current
boundary truth remains executable: matrix-free SVD now carries its native
callback-boundary label, the matrix-free nonsymmetric native callback boundary
and the matrix-free generalized-B native contract remain certified, and the
strict gate fails if provenance drifts.

To smoke every structured release surface from the source tree, run:

```sh
Rscript inst/benchmarks/run-post-v1-gates.R --tier=smoke --gates=post_v1_svd_hard_surface,post_v1_operator_sidecars,post_v1_randomized_svd_hard_surface,post_v1_generalized_preconditioned_surface,post_v1_shift_invert_boundaries,post_v1_nonsymmetric_matrix_free_surface --load-all
```

This all-surface smoke profile is contract/provenance strict, not a speed
promotion signoff. For example, SVD quick-smoke rows keep speed and memory
ratios in the output but treat them as diagnostics because the tiny fixtures are
fixed-overhead sensitive. Use the installed-package strict tier for performance
promotion evidence.

The structural smoke check remains:

```sh
Rscript -e 'pkgload::load_all("."); testthat::test_file("tests/testthat/test-bench-smoke.R", reporter="summary")'
```

## CI And Nightly Profiles

The GitHub Actions workflow is `.github/workflows/post-v1-benchmarks.yaml`.
It runs the smoke tier on a weekly schedule and supports manual dispatch with:

- `tier=smoke`: local/CI boundary-truth gate.
- `tier=strict`: installed-package strict profile using each gate's strict
  command.
- `tier=long`: nightly/weekly long profile using each gate's `long_command`,
  generally the strict command with expanded iteration counts.
- `gates=<id1,id2>`: optional comma-separated gate filter for targeted reruns.

Examples:

```sh
Rscript inst/benchmarks/run-post-v1-gates.R --tier=strict --gates=post_v1_operator_sidecars
Rscript inst/benchmarks/run-post-v1-gates.R --tier=long --dry-run
```

Use `--load-all` only for local source-tree smoke/debug runs. Installed-package
strict and long evidence should not use it.

Failures report the gate id, owner mote, tier, case set, baseline set, subject
where applicable, and command. That is the intended handoff unit for deciding
whether a solver can be promoted or needs more implementation work.
