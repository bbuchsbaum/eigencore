# V1 Benchmark Manifest

Date: 2026-06-05

This manifest is the runnable benchmark inventory for the V1 readiness audit.
It is not a release signoff. A row is release-grade only when it was run from
a clean installed package, saved its result artifacts, and passed the strict
gate for the promoted solver family.

The companion documentation scope audit is
`docs/v1-doc-scope-audit.md`. That audit maps user-facing docs, migration
docs, design notes, and release evidence to V1 requirements; it does not relax
the benchmark stop rule below.

The companion completion audit is `docs/v1-completion-audit.md`. It is the
single checklist to consult before any V1 completion claim.

Use an installed package path for release evidence, for example:

```sh
R CMD INSTALL --library=/tmp/eigencore-bench-lib .
R_LIBS=/tmp/eigencore-bench-lib Rscript inst/benchmarks/<script>.R ...
```

## Release Surfaces

| Surface | Script | Primary command | Saved artifacts | Current status |
|---|---|---|---|---|
| G1 Hermitian native gate | `inst/benchmarks/bench-native-hermitian-gate.R` | `R_LIBS=/tmp/eigencore-bench-lib Rscript inst/benchmarks/bench-native-hermitian-gate.R --strict --save --cases=path_laplacian:1000` | `*-native-hermitian-gate-rows.rds`, `*-native-hermitian-gate-summary.rds` | Green for the promoted structured-tridiagonal default. Installed strict evidence from 2026-06-05 saved `20260605-native-hermitian-gate-rows.rds` and `20260605-native-hermitian-gate-summary.rds`; `path_laplacian:1000` certified `20/20`, passed speed (`6.15x` versus RSpectra), memory (`1.07x` versus best reference), and PRIMME parity (`22.58x`). |
| G1 sparse Hermitian coverage | `inst/benchmarks/bench-hermitian-sparse.R` | `R_LIBS=/tmp/eigencore-bench-lib Rscript inst/benchmarks/bench-hermitian-sparse.R --strict --save --cases=<case:id>` | `*-hermitian-sparse-rows.rds`, `*-hermitian-sparse-summary.rds` | Diagnostic. Used with case filtering for sparse Hermitian reruns; current sparse block auto-promotion is disabled except diagnostic opt-in. |
| H production SVD | `inst/benchmarks/bench-svd-surface.R` | Executable production H candidate: `R_LIBS_USER=/tmp/eigencore-bench-lib Rscript inst/benchmarks/bench-svd-surface.R --quick --iterations=3 --h-candidate --methods=eigencore,RSpectra,PRIMME,irlba,rsvd --cases=tall_sparse:600x90,wide_sparse:90x600`; default `--h-candidate` gates the promoted `eigencore` subject while keeping retained BPRO/block-GK rows available as diagnostics when methods are not overridden. Retained diagnostics may still use `--subject=eigencore_irlba_lbd_retained_bpro`; normal-scout diagnostics may use `--methods=eigencore_irlba_lbd_normal_scout,eigencore_golub_kahan_one_sided,RSpectra --subject=eigencore_irlba_lbd_normal_scout`. | `*-svd-surface-rows.rds`, `*-svd-surface-gates.rds`, `*-svd-surface-memory.rds` | Green for the promoted tall/wide sparse H surface. Warning-free installed 2026-06-05 `eigencore` H probes certify all five requested triplets on `tall_sparse:600x90` and `wide_sparse:90x600`; the tall row uses native right-normal `implicit_normal_lanczos` without materializing the Gram, and the wide row uses the certified native left-Gram special case. The gate passes speed/memory (`1.194x` / `2.392x` tall, `1.223x` / `2.392x` wide). Retained IRLBA/LBD/BPRO diagnostics remain benchmark-visible but unpromoted. |
| I randomized SVD | `inst/benchmarks/bench-randomized-rsvd.R` | Strict scoped release gate: `R_LIBS_USER=/tmp/eigencore-bench-lib Rscript inst/benchmarks/bench-randomized-rsvd.R --iterations=1 --methods=eigencore_randomized,rsvd --strict`; quick diagnostics: `--quick --iterations=5 --methods=eigencore_randomized,rsvd` | `*-randomized-rsvd-rows.rds`, `*-randomized-rsvd-gates.rds` | Green for the scoped V1 release row. Fresh installed 2026-06-05 strict non-quick evidence certifies `exact_low_rank_dense:2000x500` for both eigencore and `rsvd`, reports `randomized_native_sketch = TRUE`, `randomized_sketch_kind = "native_fused_a_omega"`, and `randomized_projection_kind = "native_direct_qt_a"`, and passes speed/accuracy at about `3.1x` versus certified `rsvd`. Nearly-low-rank and slow-decay rows remain printed diagnostics with `release_gate_required = FALSE` because `rsvd` fails eigencore certification. Quick exact-low-rank dense and low-rank sparse rows are also non-blocking diagnostics because fixed public-result overhead makes the `2x` ratio unstable at those sizes. |
| J generalized SPD LOBPCG | `inst/benchmarks/bench-generalized-lobpcg.R` | Full installed scoped gate: `R_LIBS_USER=/tmp/eigencore-bench-lib Rscript inst/benchmarks/bench-generalized-lobpcg.R --iterations=1 --strict --save`; focused sparse slices may use `--cases=sparse_generalized_path_smallest:1000` or `--cases=sparse_generalized_path_largest:1000` with `--methods=eigencore_shifted_tridiagonal,base --subject=eigencore_shifted_tridiagonal --strict`; dense auto fallback boundary remains inspectable with `--cases=dense_generalized_partial_smallest:180,dense_generalized_partial_largest:180 --methods=eigencore_auto,base --strict` | `20260605-generalized-lobpcg-rows.rds`, `20260605-generalized-lobpcg-gates.rds`, `20260605-generalized-lobpcg-native-contracts.rds`, `20260605-generalized-lanczos-reference-contracts.rds`, `20260605-generalized-lobpcg-adversarial-b-contracts.rds` | Green for the promoted sparse shifted-tridiagonal generalized LOBPCG surface. Fresh installed saved evidence on 2026-06-05 certifies `10/10` for both sparse-smallest and sparse-largest at `n = 1000`, keeps native generalized/B-orthogonal/shifted-tridiagonal provenance visible, and passes speed/memory versus dense base at about `1.56x` / `8.8x` for smallest and `1.67x` / `25.0x` for largest. Dense generalized `auto()` stays on the native dense LAPACK fallback, and generalized Lanczos plus adversarial B rows remain contract/reference diagnostics rather than broader production claims. |
| K generalized SPD B-orthogonal Lanczos | `inst/benchmarks/bench-generalized-lobpcg.R` focused K rows plus focused tests | Installed diagonal probe: `R_LIBS=/tmp/eigencore-bench-lib Rscript inst/benchmarks/bench-generalized-lobpcg.R --quick --strict --iterations=1 --cases=diagonal_generalized_lanczos_ref_smallest --methods=eigencore_lanczos_reference,eigencore,base --subject=eigencore_lanczos_reference`; installed sparse-CSC probe: same command with `--cases=sparse_csc_generalized_lanczos_ref_smallest`; source test: `Rscript -e 'pkgload::load_all("."); testthat::test_file("tests/testthat/test-generalized-lobpcg.R", reporter="summary")'` | `*-generalized-lanczos-reference-contracts.rds` when saved with the broader generalized benchmark script | Partial reference with installed benchmark evidence. Explicit generalized-SPD `lanczos()` requests now have an honest reference scalar B-orthogonal refinement for dense, diagonal, and CSC SPD metric solves; focused tests and installed K rows cover B-orthogonality, original-coordinate certification, LOBPCG agreement, planner label, diagonal and sparse-Cholesky metric-solve provenance, and no dense fallback for the diagonal/CSC cases. No native/block production benchmark gate exists yet. |
| M standard preconditioned LOBPCG | `inst/benchmarks/bench-lobpcg-preconditioned.R` | `R_LIBS=/tmp/eigencore-bench-lib Rscript inst/benchmarks/bench-lobpcg-preconditioned.R --strict --save` | `*-lobpcg-preconditioned-rows.rds`, `*-lobpcg-preconditioned-gates.rds` | Green for the path-Laplacian release surface. Native shifted-tridiagonal rows beat certified references on `n = 200, 1000, 2000`, `k = 5`. Broader non-Laplacian policy is future scope. |
| L shift-invert contract | `inst/benchmarks/bench-shift-invert.R` | Native label surface: `R_LIBS=/tmp/eigencore-bench-lib Rscript inst/benchmarks/bench-shift-invert.R --iterations=1 --save`; quick full label surface: `R_LIBS=/tmp/eigencore-bench-lib Rscript inst/benchmarks/bench-shift-invert.R --quick --strict --iterations=1 --save`; user-solve boundary: `R_LIBS=/tmp/eigencore-bench-lib Rscript inst/benchmarks/bench-shift-invert.R --quick --strict --iterations=1 --cases=matrix_free_user_solve_reference` | `20260517-shift-invert-rows.rds`, `20260517-shift-invert-contracts.rds` | Green for current labels. Fresh installed non-quick strict evidence on 2026-05-17 certifies dense standard, dense generalized, diagonal standard, sparse symmetric-tridiagonal standard, and sparse/diagonal tridiagonal-generalized native rows in original coordinates; the generalized tridiagonal row uses `tridiagonal_thomas_generalized_native`. A fresh installed quick strict saved run also covers general sparse standard and sparse/diagonal generalized reference rows: both converge with estimated-scale certificates that correctly do not mark `passed`, retain sparse-LU cache provenance, and remain nonnative. The installed quick strict `matrix_free_user_solve_reference` row certifies in original coordinates with exact scale metadata, records `factorization = user_solve`, `external_cache = TRUE`, and `label_kind = user_solve`, and remains nonnative/reference-labelled. L remains partial because general sparse native shift-invert is not implemented. |
| Nonsymmetric compatibility | `inst/benchmarks/bench-nonsymmetric.R` | `R_LIBS=/tmp/eigencore-bench-lib Rscript inst/benchmarks/bench-nonsymmetric.R --iterations=1 --save --strict`; quick strict smoke: `--quick --strict --iterations=1` | `20260517-nonsymmetric-rows.rds`, `20260517-nonsymmetric-contracts.rds` | Green for current compatibility labels. Dense nonsymmetric and dense `eigs(..., which = "LI")` rows now certify right residuals through the same native Arnoldi-cycle plus native-Ritz compatibility label as sparse real- and imaginary-target rows. Installed non-quick strict evidence on 2026-05-17 certifies `dense_native_arnoldi_lm`, `dense_native_arnoldi_li`, `dense_eigs_native_arnoldi_li`, `sparse_native_arnoldi_lr`, and `sparse_native_arnoldi_li`, all with `native_arnoldi_label = TRUE`, `ritz_extraction_native = TRUE`, and exact right-residual certificates. Dense rows use full dense subspace compatibility; sparse CSC rows keep the bounded subspace/restart policy. Matrix-free Arnoldi remains reference-labelled and production-grade fully native restarted Arnoldi remains open. |

## Diagnostic And Baseline Surfaces

| Surface | Script | Purpose |
|---|---|---|
| Block Hermitian prototype | `inst/benchmarks/bench-block-hermitian-prototype.R` | Development diagnostics for block Hermitian variants before promotion. |
| G1 candidate baseline | `inst/benchmarks/bench-g1-candidate-baseline.R` | Rebuilds `inst/benchmarks/baselines/g1_candidate_pre.csv` for before/after G1 candidate comparisons. |
| Performance baseline | `inst/benchmarks/bench-performance-baseline.R` | Broad comparative smoke across eigen, SVD, and generalized cases; not a strict release gate. |
| SVD tall-skinny | `inst/benchmarks/bench-svd-tallskinny.R` | Focused SVD development surface retained for historical tall-skinny diagnostics. |
| Tiny Gram eigensolvers | `inst/benchmarks/bench-tiny-gram-eigensolvers.R` | Tracks tiny dense Gram solve choices used by bounded SVD special cases. |
| Reference suite | `benchmarks/reference_suite.R` | Legacy reference comparisons; use the explicit `inst/benchmarks` release surfaces for V1 gate evidence. |

## Stop Rule

Do not use this manifest to close V1 by itself. V1 benchmark signoff still
requires:

1. Every promoted solver family to pass its strict installed-package gate.
2. Red rows above to be fixed, demoted, or explicitly scoped out by a PRD
   revision.
3. Saved artifacts from the final runs to be listed in `benchmarks/RELEASES.md`
   with the command that produced them.
4. A fresh `R CMD build` and `R CMD check --no-manual` after the final
   benchmark/doc updates.
