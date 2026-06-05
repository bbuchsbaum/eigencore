# eigencore Benchmark Releases

This file records curated release benchmark summaries. Raw benchmark runs are
written under `inst/benchmarks/results/` and are intentionally ignored because
they are machine-dependent.

## Methodology

- Metric: time to certified answer.
- Timer: `bench::mark`, with certificate recomputation after each method.
- Eigen reference category: sparse path graph Laplacian, smallest eigenpairs.
- SVD reference category: tall-skinny sparse matrix, largest singular triplets.
- Reported columns include median time, allocated memory, maximum residual,
  maximum backward error, orthogonality loss, converged count, seed, and
  eigencore package version.

## Unreleased

- Benchmark plumbing added.
- Native operator provenance now survives `adjoint()` for built-in dense,
  CSC, and diagonal operators. The adjoint wrapper continues to use the native
  parent apply closures, records `fused = "adjoint"`, preserves
  `metadata$native = TRUE`, and carries an explicit transposed source/matrix
  when available. This is a small operator-fusion step: planner and diagnostic
  surfaces can distinguish native adjoints from R callback operators while
  centered/scaled/composed native fusion remains open.
- Built-in scalar, row, and column scaling now fuse to native-backed operators
  when the parent has an explicit dense, CSC, or diagonal backing. Dense inputs
  materialize a scaled dense operator, CSC inputs preserve sparse `dgCMatrix`
  storage instead of densifying, and diagonal scalar scaling remains
  diagonal-backed. Matrix-free operators continue to use the existing R-level
  transform path.
- Built-in explicit operator algebra now fuses more transforms into native-backed
  operators when it can do so without densifying sparse inputs: dense and CSC
  sums/compositions, explicit dense/CSC `crossprod_operator()`, and dense
  centering return native-backed operators with `metadata$fused` provenance.
  Sparse CSC centering now uses a native centered-block apply wrapper that
  combines CSC apply with row/column low-rank corrections without materializing
  the centered dense matrix; it is marked as a fused centered operator but not
  as a plain `dgCMatrix` solver kernel. Matrix-free centering remains on the
  honest callback path.
- The former native source monolith has been split into focused translation
  units. `src/small_dense.cpp` provides standalone LAPACK/scalar shells: full and selected
  symmetric eigensolves, dense symmetry checking, generalized-SPD dense
  eigensolve, dense SVD, tridiagonal eigensolve, bidiagonal SVD, and the small
  tridiagonal solve helper. A second low-risk split moved shared native
  stage/history structs, timing helpers, LP64 index checks, apply-status
  errors, and the compact basis-combination helper into
  `src/eigencore_common.h`. The native Arnoldi cycle wrappers and
  projected-Ritz extraction wrapper now live in `src/arnoldi.cpp`, the dense
  symmetric Rayleigh-Ritz projection wrapper plus compiled scalar and block
  Golub-Kahan Ritz implementations now live in `src/projection.cpp`, and the
  registered native orthogonalization wrappers
  plus their private basis-workspace helper now live in
  `src/orthogonalization.cpp`. The native operator structs and built-in
  dense, CSC, diagonal, R-callback, and factorized shift-invert apply functions
  now live in `src/native_operators.h` and `src/native_operators.cpp`. The
  R-facing native operator block-apply/check wrappers now live alongside those
  kernels in `src/native_operators.cpp`. Dense/R-facing certificate wrappers
  and built-in native-operator certificate wrappers now live in
  `src/certificates.cpp`; `src/certificates.h` exposes the internal cached-
  `A v` helper needed by retained SVD. The native LOBPCG run helper, standard,
  generalized, matrix-free-B, and shifted-tridiagonal LOBPCG wrappers now live
  in `src/lobpcg.cpp`. Scalar native Lanczos, shift-invert Lanczos, and
  scalar Golub-Kahan staging wrappers now live in `src/scalar_krylov.cpp`,
  with `src/scalar_krylov.h` exposing the internal Golub-Kahan run helper used
  by retained SVD callers. The block Golub-Kahan basis builders and their
  dense/CSC registered wrappers now live in `src/block_golub_kahan_basis.cpp`,
  with `src/block_golub_kahan_basis.h` exposing the internal basis builder used
  by retained SVD callers. The block Golub-Kahan fit/retained-cycle wrappers
  and one-sided retained IRLBA/LBD wrappers now live in `src/retained_svd.cpp`.
  The block/thick-restart Lanczos helpers and registered dense/CSC wrappers now
  live in `src/block_lanczos.cpp`.
  The CSC left/right Gram SVD special-case kernels and their fast-result
  wrappers now live in `src/gram_svd.cpp`. The `.Call` registration table and
  `R_init_eigencore()` now live in
  `src/init.cpp`, with exact-arity extern declarations for the existing
  registered symbols. The obsolete `src/dense_block_apply.cpp` compatibility
  stub has been removed.
- Native LOBPCG wrappers now share one internal allocate/run/pack helper after
  each registered entry point validates its storage-specific inputs. This
  removes duplicated work-buffer/result packing code across dense, CSC,
  diagonal, generalized, matrix-free-`B`, and shifted-tridiagonal entry points
  without changing the `.Call` registration surface or claiming new solver
  capability.
- The block reference solvers now share `reference_block_solver_skeleton.R`
  helpers for restart-control validation, basis acceptance, initial basis
  construction, lock counting, restart keep counts, and convergence-history
  rows. Block Hermitian subspace iteration and block Golub-Kahan still own
  their solver-specific expansion and Ritz extraction, but common restart
  bookkeeping now has one implementation. Scalar reference Lanczos and scalar
  Golub-Kahan now also share scalar maxit/subspace validation through the same
  internal skeleton helper, while their recurrence loops and Ritz extraction
  remain solver-local.
- Added `docs/v1-readiness-audit.md` as the release-hardening checklist that
  maps each PRD/plan gate to current evidence, verification commands, and
  blockers. It records that V1 is not release-ready while H, I, J, K, sparse
  shift-invert, nonsymmetric Arnoldi, and release-hardening checks remain open.
- Added `docs/v1-benchmark-manifest.md` as the benchmark inventory for the V1
  audit. It maps each release surface to the installed-package command, saved
  artifacts, gate meaning, and current red/green status so benchmark evidence
  is not scattered only through release-note prose. The manifest is not a
  signoff artifact by itself; red gates still need fixes, demotions, or PRD
  scope changes before V1.
- Added `docs/v1-completion-audit.md` as the final stop-rule checklist for the
  active V1 readiness goal. It restates the deliverables, maps each prompt
  requirement to concrete artifacts and evidence, and records the current
  decision as not V1 ready while G1, H, I, J, K, sparse-native L,
  production-grade nonsymmetric Arnoldi, sanitizer/valgrind-style coverage,
  final benchmark artifacts, and final README/vignette refresh remain open.
- Added migration-facing release-hardening docs:
  `docs/rspectra-migration.md` for the RSpectra shim contract and
  `docs/known-limitations.md` for the current non-V1 solver and validation
  gaps. Added `docs/method-selection-and-workflows.md` to cover the broader
  partial eigen, SVD, generalized SPD, shift-invert, operator, certificate,
  diagnostics, and method-selection workflow. These docs intentionally keep
  reference, prototype, and oracle-labelled paths out of production performance
  claims.
- Added `docs/v1-doc-scope-audit.md` to map README, vignette, migration,
  method-selection, native-design, benchmark, and release-evidence docs to V1
  documentation requirements. It records the remaining doc blockers explicitly:
  README/vignette refresh after solver-gate decisions, final strict benchmark
  artifacts, and unresolved ASan/valgrind-equivalent release coverage.
- Fresh installed-package Hermitian gate evidence was collected on 2026-05-17
  using
  `R_LIBS=/tmp/eigencore-bench-lib Rscript inst/benchmarks/bench-native-hermitian-gate.R --quick --include-dense --strict --save`.
  The quick strict run certified both rows and saved
  `inst/benchmarks/results/20260517-native-hermitian-gate-rows.rds` plus
  `inst/benchmarks/results/20260517-native-hermitian-gate-summary.rds`.
  This is certification smoke only: the printed quick speed/memory gates were
  false (`0.098x` and `0.627x` speed ratios versus the best reference;
  `0.703x` and `0.461x` memory ratios). A non-quick
  `bench-native-hermitian-gate.R --include-dense --strict --save` attempt was
  stopped after roughly 18 minutes with no gate output in this local session,
  so final G1 release signoff still requires fresh non-quick strict gate
  evidence from an installed package.
- Partial sanitizer-style evidence was collected on 2026-05-17. A temporary
  UBSan-only install succeeded with
  `PKG_CXXFLAGS='-fsanitize=undefined -fno-sanitize-recover=undefined'`
  and `PKG_LIBS='-fsanitize=undefined'`, followed by a native smoke covering
  dense Hermitian eigen, sparse CSC Hermitian Lanczos, dense SVD, and dense
  shift-invert; the smoke printed `UBSan native smoke passed`. ASan builds
  compiled but failed during package load because the ASan interceptors are
  loaded too late through R's `dlopen` path on this macOS setup, even when
  `DYLD_INSERT_LIBRARIES` was set. `valgrind` is not installed locally. Treat
  the sanitizer row as partial until ASan or valgrind-equivalent coverage is
  green in a suitable environment.
- The native sanitizer smoke is now a reusable installed-package artifact at
  `inst/validation/native-smoke.R`. It covers dense Hermitian eigen, sparse CSC
  Hermitian Lanczos, dense generalized SPD LOBPCG, dense SVD, sparse CSC SVD,
  dense shift-invert, and native tridiagonal shift-invert, and prints
  `eigencore native smoke passed` on success. Use
  `Rscript -e 'source(system.file("validation/native-smoke.R", package = "eigencore"))'`
  after a sanitizer/valgrind-style install; `Rscript inst/validation/native-smoke.R --load-all`
  is the source-checkout smoke.
- Hermitian benchmark reruns are now case-filterable and progress-visible.
  `bench-native-hermitian-gate.R` and `bench-hermitian-sparse.R` accept
  `--cases=<case>` or `--cases=<case:n>` using stable ids such as
  `path_laplacian:10000`, and print the case id before entering each expensive
  benchmark. This does not change gate math; it exists so stalled non-quick
  release reruns can be bisected and diagnosed. A fresh installed quick check
  with `--cases=path_laplacian:200` printed progress and ran only that row.
- The first filtered non-quick installed Hermitian rerun is red:
  `R_LIBS=/tmp/eigencore-bench-lib Rscript inst/benchmarks/bench-native-hermitian-gate.R --strict --save --cases=path_laplacian:1000`
  certified all 20 requested eigenpairs and passed the memory gate
  (`1.78x` memory ratio versus the best certified reference), but failed the
  release speed gate (`0.0064x` versus RSpectra) and PRIMME parity
  (`0.025x`). The row shows the hotspot clearly:
  `stage_projected_solve_seconds` about `6.77s` and
  `stage_projection_update_seconds` about `6.75s` in a `7.04s` total median.
  Treat G1 as red until the projected-subspace/update cost is fixed or the
  promotion policy is corrected.
- Sparse block Hermitian auto-promotion has been corrected after that red
  evidence: default `auto()` now keeps sparse `dgCMatrix` Hermitian requests on
  scalar Lanczos, while explicit `lanczos(block > 1)`, dense full-subspace
  block selection, benchmark candidates, and the diagnostic
  `options(eigencore.promote_sparse_block_lanczos = TRUE)` path remain
  available for development and gate reruns.
- A post-demotion installed rerun of the same row confirms the new default
  sparse route is scalar (`block_size = 1`) and still red:
  `R_LIBS=/tmp/eigencore-bench-lib Rscript inst/benchmarks/bench-native-hermitian-gate.R --strict --save --cases=path_laplacian:1000`
  certified all 20 requested eigenpairs, but failed speed (`0.0086x` versus
  RSpectra), memory (`0.636x` versus the best certified reference), and PRIMME
  parity (`0.0336x`). The hotspot remains projected-subspace work:
  `stage_projected_solve_seconds` about `4.93s` and
  `stage_projection_update_seconds` about `4.91s` in a `5.41s` total median.
  This keeps G1 red even after the sparse block-promotion correction.
- The projected-matrix update in native thick-restart Lanczos now updates only
  the appended basis block during expansion, instead of recomputing the whole
  active projection after every accepted column. On an installed
  `path_laplacian:1000` scalar rerun this reduced total median time from about
  `5.41s` to `0.327s`, `stage_projected_solve_seconds` from about `4.93s` to
  `0.105s`, and `stage_projection_update_seconds` from about `4.91s` to
  `0.095s`. The row still certifies all 20 pairs but remains red against the
  release gate: speed `0.138x` versus RSpectra, memory `0.636x`, and PRIMME
  parity `0.513x`. Larger scalar subspaces (`120`, `160`, `240`, `320`) were
  slower in single-run probes. The post-fix block candidate also remains red on
  this row: `0.488s` median, memory gate green, speed `0.093x` versus RSpectra,
  and PRIMME parity `0.356x`.
- Certificate scale-estimate hardening is now explicitly tested across current
  certificate entry points. `new_certificate(scale_is_estimate = TRUE)`,
  matrix-free eigen certificates, matrix-free SVD certificates, and
  residual-backed generalized eigen certificates all withhold `passed` even
  when residual/backward-error convergence is true. This keeps stochastic
  Hutchinson norm estimates from being counted as release-grade certificates.
- Result diagnostics now have a cross-path contract test. Dense Hermitian,
  native sparse Lanczos, generalized LOBPCG, dense shift-invert, nonsymmetric
  dense oracle, randomized SVD, and RSpectra-compatible shims all expose stable
  method/plan/certificate/restart/stage diagnostics through `diagnostics()`.
  The test found and fixed a `diagnostics()` edge case where results without a
  list-valued `restart` could error while looking for restart preconditioner
  metadata.
- Generalized SPD B-orthogonal Lanczos now has an honest reference refinement
  slice for explicit `method = lanczos()` requests when `B` has a dense,
  diagonal, or CSC SPD solve. The planner label is
  `reference generalized SPD B-orthogonal Lanczos refinement`; it performs
  shift-invert-free B-orthogonal Rayleigh-Ritz extraction and certifies in the
  original generalized coordinates. Focused tests cover diagonal, dense, and
  sparse CSC SPD metrics, B-orthogonality, and agreement with generalized LOBPCG
  certificates. This advances K from absent to partial reference coverage; it
  is not a native/block production promotion.
- The generalized LOBPCG benchmark harness now includes a focused
  `eigencore_lanczos_reference` method row and a
  `generalized_lanczos_reference_contract` table. A fresh installed quick
  strict probe:
  `R_LIBS=/tmp/eigencore-bench-lib Rscript inst/benchmarks/bench-generalized-lobpcg.R --quick --strict --iterations=1 --cases=diagonal_generalized_lanczos_ref_smallest --methods=eigencore_lanczos_reference,eigencore,base --subject=eigencore_lanczos_reference`
  certified the reference Lanczos row (`nconv = 2/2`), passed the K contract,
  recorded `diagonal solve for B`, and kept the row honestly nonnative.
  A matching installed quick strict sparse-CSC metric probe with
  `--cases=sparse_csc_generalized_lanczos_ref_smallest` also certifies
  `2/2`, passes the K reference contract, records
  `sparse Cholesky solve for B`, and remains nonnative while the native LOBPCG
  comparison row certifies.
- The shift-invert benchmark now covers the general sparse reference boundary
  explicitly. A fresh installed quick strict saved run covers the native rows
  plus `sparse_general_reference` and
  `sparse_general_diagonal_b_reference`; both reference rows retain
  `Matrix::lu` cache provenance, certify convergence in original coordinates,
  stay nonnative/reference-labelled, and correctly withhold
  `certificate_passed` because their norm scale is estimated. Saved artifacts:
  `inst/benchmarks/results/20260517-shift-invert-rows.rds` and
  `inst/benchmarks/results/20260517-shift-invert-contracts.rds`. This
  documents the L boundary without promoting general sparse native
  factorization.
- The shift-invert benchmark also has an installed quick strict
  `matrix_free_user_solve_reference` row for user-supplied shifted solves. It
  certifies in original coordinates with exact operator-scale metadata, records
  `factorization = user_solve`, `external_cache = TRUE`, and
  `label_kind = user_solve`, and remains nonnative/reference-labelled rather
  than claiming native factorization ownership.
- The SVD surface benchmark now uses the shared case-filter/progress helpers.
  Cases have stable ids such as `tall_sparse:600x90` and
  `tall_sparse:100000x500`, `--cases=` accepts either the stable id or the case
  name, and each selected row prints progress before entering the expensive
  method loop. A fresh installed quick H probe:
  `R_LIBS=/tmp/eigencore-bench-lib Rscript inst/benchmarks/bench-svd-surface.R --quick --h-candidate --iterations=1 --save --cases=tall_sparse:600x90`
  keeps H red. The gate subject `eigencore_block_golub_kahan_retained`
  certified all five requested singular triplets but failed speed
  (`0.0188x` versus the best certified reference) and memory (`0.0619x`).
  The fastest certified eigencore diagnostic rows on this tiny tall sparse case
  were Gram/implicit-normal paths around `0.00039s`, but these are bounded
  diagnostic paths, not the general H production answer.
  Fresh default tiny-Gram probes further correct the previous Track A promotion
  note: `--subject=eigencore --methods=eigencore,RSpectra,irlba,rsvd` on
  `tall_sparse:600x90` and `wide_sparse:90x600` certifies all five eigencore
  triplets and passes the memory gate, but fails speed (`0.54x` tall,
  `0.58x` wide). Treat the older repeated 90-by-600 green result as historical
  narrow evidence, not current V1 release signoff.
- The native CSC Gram fast-result path now covers the tall/right-Gram special
  case as well as the previous wide/left-Gram case. It builds the
  certificate, plan, restart diagnostics, and classed SVD result in native code
  only when the Gram certificate already passes; failed Gram certificates still
  return to the existing Golub-Kahan fallback. A direct installed comparison on
  `tall_sparse:600x90` shows lower overhead for `svd_partial()` fast-result
  construction (about `604us`) than the R-assembled `solve(svd_problem())`
  route (about `692us`) with the same backward error. This does not close H:
  fresh 3-iteration installed quick probes still report `0.40x` tall and
  `0.63x` wide speed versus the best certified reference, while memory remains
  green.
- The randomized-rsvd benchmark now uses the shared case-filter/progress
  helpers. Cases have stable ids such as `exact_low_rank_dense:120x80`,
  `slow_decay_dense:140x90`, and `exact_low_rank_dense:2000x500`, and
  `--cases=` accepts either the stable id or the case name. A fresh installed
  quick probe:
  `R_LIBS=/tmp/eigencore-bench-lib Rscript inst/benchmarks/bench-randomized-rsvd.R --quick --iterations=1 --save --cases=exact_low_rank_dense:120x80,slow_decay_dense:140x90`
  keeps I red on quick evidence. On `exact_low_rank_dense:120x80`,
  `eigencore_randomized` and `rsvd` both certified and eigencore was faster
  (`1.85x` in the latest saved quick probe), but still below the randomized
  `2x` release speed gate. On `slow_decay_dense:140x90`, eigencore certified
  all eight requested singular triplets while `rsvd` failed eigencore
  certification. The gate now requires `baseline_certified = TRUE` before
  speed/parity can pass, so this row records `baseline_certified = FALSE`,
  `speed_gate = FALSE`, and `passed = FALSE` instead of treating an
  uncertified reference as a valid parity target.
- A fresh installed non-quick exact-low-rank randomized row remains green:
  `R_LIBS=/tmp/eigencore-bench-lib Rscript inst/benchmarks/bench-randomized-rsvd.R --iterations=1 --save --cases=exact_low_rank_dense:2000x500`
  certified all 50 requested singular triplets for both eigencore and `rsvd`,
  passed accuracy, and reported `2.97x` time-to-certified-answer speed versus
  certified `rsvd` with the `left_gram_eigen` randomized core solver. Treat
  this as a validated exact-low-rank regime, not a full I promotion: quick
  small rows and slow-decay/native sketch-projection gates remain open.
- A source-loaded quick experiment replacing the tiny projected-core
  `left_gram_eigen` helper with base `svd()` for the
  `exact_low_rank_dense:120x80` core was rejected. Isolated micro-timing on the
  16-by-80 projected core favored `svd()`, but the full benchmark row worsened
  from about `1.84x` to about `1.63x` versus certified `rsvd` and allocated
  more. Keep `left_gram_eigen`; the remaining I work needs native
  sketch/projection or planning changes, not this small-core swap.
- A fresh installed quick sparse randomized probe:
  `R_LIBS=/tmp/eigencore-bench-lib Rscript inst/benchmarks/bench-randomized-rsvd.R --quick --iterations=1 --cases=low_rank_sparse:140x90 --methods=eigencore_randomized,rsvd`
  certifies both eigencore and `rsvd`, with eigencore faster at `1.87x` but
  still below the `2x` randomized gate. A direct CSC projected-core experiment
  (`Q' A` instead of `t(A' Q)`) was rejected after the full row worsened to
  `0.71x` versus certified `rsvd` and inflated solver allocation, so the CSC
  randomized path keeps the existing `A' Q` projection until a real native
  sparse sketch/projection kernel is available.
- G1 block Hermitian candidate baseline captured in
  `inst/benchmarks/baselines/g1_candidate_pre.csv`; regenerate it with
  `inst/benchmarks/bench-g1-candidate-baseline.R --save`. The current baseline
  certifies the sparse Laplacian, dense Hermitian, clustered, and
  ill-conditioned diagonal block-candidate rows. Strict mode now fails if any
  block-candidate row is missing, errors, or returns without a passing
  certificate.
- G1 quick tuning grid captured in
  `inst/benchmarks/baselines/g1_tuning_quick.csv`.
- H SVD surface benchmark now supports `--subject=<method>` and the
  `--h-candidate` preset, so projected Golub-Kahan candidate rows can be gated
  directly against external references without changing the default `eigencore`
  SVD release gate. The preset selects `eigencore_golub_kahan_projected` as the
  subject, keeps the plain Golub-Kahan row for projected-stop comparison, and
  errors if the requested gate subject is absent. The surface now also emits
  `svd-surface-memory` diagnostics that split total, solver, and certificate
  allocation gaps against the best certified reference. A quick projected
  subject check certifies sampled
  wide-sparse and clustered-dense rows, but remains well below the `1.5x` SVD
  speed gate and memory parity; H is still a performance milestone. The native
  Golub-Kahan kernel now reports `native_workspace_bytes` and returns compact
  R-visible bases sized to realized iterations, which reduced projected solver
  allocation in the quick H check to `464576` bytes for wide sparse and
  `147664` bytes for clustered dense while keeping certificates green. The
  default non-diagnostic path now goes further and returns compact native Ritz
  fits instead of Krylov bases; the same quick H rows allocate `187600` bytes
  for wide sparse and `112544` bytes for clustered dense with
  `basis_returned = FALSE`. Dense Golub-Kahan reorthogonalization now uses
  BLAS `dgemv` projections; CSC keeps scalar projection after the quick sparse
  row showed BLAS overhead without iteration savings. Native Golub-Kahan rows
  now populate `stage_apply_seconds`, `stage_recurrence_seconds`,
  `stage_reorthogonalization_seconds`, and `stage_projected_solve_seconds`;
  current H samples show reorthogonalization, not operator application or
  certification, is the primary sparse SVD hotspot. The scalar Golub-Kahan
  path now uses an adaptive DGKS-style second reorthogonalization pass and
  reports `reorthogonalization_passes`; sampled projected rows remain
  certified, with wide sparse around `0.0028s` and clustered dense around
  `0.0012s`. The SVD surface rows now also report normalized H hotspot
  diagnostics, including accounted native stage time, reorthogonalization time
  fraction, seconds per pass, passes per iteration, native seconds per matvec,
  projected-check cost, and projected-stop savings fractions.
- Native Golub-Kahan now completes exact zero singular triplets after
  rank-deficient breakdown before certification. In the quick H candidate
  surface, `rank_deficient_sparse` now reports the projected subject as
  certified with `nconv = requested = 6`; the H gate still fails on speed, not
  missing requested triplets.
- An unexported reference block Golub-Kahan thick-restart SVD oracle now
  exercises the production H contract before native C++ promotion: forced
  restart, clustered singular subspaces, locking diagnostics, and exact zero
  singular triplet completion. This does not change the release gate subject or
  H status. A new internal native block-GK Ritz kernel computes selected
  singular Ritz slices from a right basis and cached `A V`, with tests covering
  active-column windows, target selection, and original-coordinate
  certification. Dense and CSC internal native basis-cycle staging paths now
  feed that kernel from native block right bases and cached `A V`; these are
  full-subspace scaffolding tests, not release gate subjects. The SVD surface
  benchmark can now expose that native block basis-cycle as
  `eigencore_block_golub_kahan_cycle` under `--h-candidate`, while excluding it
  from external release-reference gates. A quick installed-package probe
  certified `rank_deficient_sparse` and `clustered_dense`, with the block-cycle
  row faster than scalar/projected Golub-Kahan on both cases, but still behind
  the best certified external reference on time and memory. On `wide_sparse`
  the same row was faster than scalar/projected Golub-Kahan but failed the
  requested certificate, so H promotion still requires a real thick-restart and
  convergence policy rather than this full-subspace staging cycle. The staging
  row now has a conservative adaptive subspace retry policy that records
  attempted subspaces and total native work. On the quick `wide_sparse` row it
  turned the previous uncertified block-cycle result into a certified result
  after three attempts, but the cost rose to about `0.0084s` and `1.15MB`,
  much worse than the scalar/projected candidates and certified references.
  This is useful convergence evidence, not a promotion candidate: H still needs
  restart/reuse of retained subspaces instead of rebuilding larger bases from
  scratch. The native block basis scratch buffer was also widened to the
  operator dimension used by fallback orthogonalization work vectors. The
  adaptive cycle now seeds later attempts with previous Ritz vectors plus a
  fresh random block. On the same `wide_sparse` fixture this cuts total native
  matvecs from `167` to `73`, but the installed quick row still reports about
  `0.0090s` and `1.31MB`; the remaining loss is attempt materialization and
  lack of a retained native restart workspace, not only operator-call count.
  A lean internal comparator,
  `eigencore_block_golub_kahan_cycle_lean`, now restarts from Ritz vectors
  without the extra random block. On the same fixture it still certifies and
  remains faster than cold adaptive restarts, but it needs more matvecs than
  the default Ritz-plus-random start; keep the default speed path until a real
  retained native restart workspace lands.
- H SVD surface rows now expose adaptive restart-efficiency diagnostics:
  attempted subspaces, maximum attempted subspace, maximum restart start width,
  warm-started attempt count, certified attempt, final-attempt matvecs,
  final-attempt orthogonalization passes, and total orthogonalization passes.
  These fields are the benchmark contract for retained-restart work: a future
  native thick-restart candidate should reduce attempt materialization and
  total work without hiding failed intermediate certificates.
- Native block-GK basis generation now has an internal cached-start path for
  retained Ritz vectors. The `eigencore_block_golub_kahan_cycle_cached`
  benchmark row restarts from Ritz vectors and their exact cached `A V`,
  avoiding the initial operator apply on warm attempts. This is still a staging
  comparator, not H promotion, but it is the first direct retained-state step
  instead of another from-scratch attempt variant.
- The H staging surface also includes
  `eigencore_block_golub_kahan_cycle_cached_random`, which caches the retained
  Ritz-vector `A V` prefix while keeping the random exploration tail used by
  the current default restart. On the source-loaded wide sparse probe
  (`90 x 600`, rank `5`, seed `701`), both the default and cached-random rows
  certify in three attempts with `73` native apply calls; cached-random is
  therefore diagnostic only, not a promotion. Restart-block construction now
  avoids one `cbind()` copy in the default Ritz-plus-random path, reducing the
  same probe's allocation slightly (`~1.31MB` to `~1.30MB`).
- Native block-GK staging now has a compact fit path: the adaptive cycle keeps
  full right bases and cached `A V` inside C for Ritz extraction and returns
  only selected singular triplets plus diagnostics to R. On the same
  source-loaded wide sparse probe, the default block-GK candidate still
  certifies in three attempts with `73` native apply calls, but R-visible
  allocation drops from about `1.30MB` to `0.37MB`; the lean comparator drops
  to about `0.28MB`. This closes a materialization gap, not the remaining
  speed/thick-restart gap.
- Compact block-GK rows now report native iteration and Ritz extraction stage
  timings. On the same source-loaded wide sparse probe, the default row reports
  about `0.0066s` in native basis iteration and about `0.0026s` in Ritz
  extraction, confirming that the next H speed work should primarily attack the
  retained native restart/workspace path rather than R result materialization.
- Compact block-GK fit work arrays now use uninitialized native allocation and
  rely on the basis runner's single explicit zeroing pass. The same wide sparse
  probe remains certified with `73` native apply calls on the default row, and
  sampled median time dropped to about `0.0106s`; this is a cleanup win, not an
  H gate closure.
- The H restart surface now includes
  `eigencore_block_golub_kahan_cycle_residual`, which uses right residual
  vectors as the warm-restart exploration tail. On the source-loaded wide
  sparse probe it certifies, but needs `87` native apply calls versus `73` for
  the default Ritz-plus-random row, so it is evidence against residual-tail
  restarts as the next promotion path. Existing block-GK candidates continue to
  use the native built-in certificate path; cached `A v` is only used when the
  residual-tail diagnostic explicitly needs residual vectors.
- Native SVD certificates now have a cached-`A v` path for dense, CSC, and
  diagonal built-in operators. Compact block-GK rows use the selected Ritz
  `A v` block already formed during extraction, avoiding a redundant left-side
  operator apply during solver-internal certification while preserving native
  certificate formulas. On the same wide sparse probe, the default block-GK
  row remains certified with `73` apply calls and R-visible allocation drops
  modestly (`~0.37MB` to `~0.36MB`); this is an incremental certification
  cleanup, not the retained-restart speed breakthrough.
- The H benchmark surface now includes
  `eigencore_block_golub_kahan_retained`, the first native retained restart
  candidate. It constructs Ritz-plus-random restart blocks inside C, returns
  compact selected triplets, and keeps certification on the existing native
  cached-`A v` path. It now also performs native per-attempt cached-`A v`
  certification and can stop on the first certified retained subspace. On the
  source-loaded wide sparse probe it certifies on attempt 3 with the same `73`
  apply calls as the R adaptive block-GK cycle, while R-visible allocation is
  lower (`~0.23MB` versus `~0.36MB` in the sampled row); this is a staging
  candidate, not H closure, because cached `A V_keep` still needs benchmark
  proof under its certificate gate, and locking plus production thick-restart
  policy are still open.
- Cached-`A V_keep` retention inside the retained block-GK cycle is now
  attempted behind a certificate gate. The native path applies the same
  Cholesky-QR normalization to retained `V_keep` and its cached operator images;
  if the cached attempt does not certify, the wrapper reruns the same
  deterministic retained restart without cached `A V_keep` and records
  `fallback_method = "retained_uncached_after_cached_av_failure"`. On one wide
  sparse probe the cached path certifies; on the benchmark seed it falls back
  and still returns a certified row, so this is a controlled optimization rather
  than H closure. The primary `eigencore_block_golub_kahan_retained` H row now
  keeps this cache disabled so it remains the stable uncached retained baseline;
  `eigencore_block_golub_kahan_retained_cached` tracks the guarded cached
  experiment separately.
- SVD benchmark rows now expose retained-restart diagnostics directly:
  `retained_restart`, `retained_restart_native`, `retained_av_cache`,
  `native_attempt_certification`, and `native_early_stop`.
- Explicit benchmark method requests are now contractual. If `--methods=...`
  names an SVD or eigen method that is not available in the loaded eigencore
  namespace, the helper errors instead of silently intersecting the request
  with available methods. This keeps H probes from dropping
  `eigencore_irlba_lbd_retained_native` when a stale installed package is
  accidentally used.
- The retained IRLBA/LBD H benchmark candidate now uses a single retained
  native fixed-work attempt before the certified fallback. The current native
  scaffold does not yet update a coupled augmented recurrence across attempts,
  so repeated seeded attempts added diagnostic overhead without satisfying the
  release gate. Rows now expose this boundary directly through
  `irlba_lbd_restart_state_kind = "ritz_subspace_only"`,
  `irlba_lbd_recurrence_available = FALSE`,
  `irlba_lbd_augmented_recurrence = FALSE`,
  `irlba_lbd_retained_seed_strategy`, `irlba_lbd_retained_from_scout`,
  `irlba_lbd_retained_padding`, and
  `irlba_lbd_retained_fixed_work_attempts`.
- Added `eigencore_irlba_lbd_normal_scout` as an H diagnostic benchmark
  method. It runs bounded matrix-free normal scouts at 8/12/16/20 steps,
  uses the selected scout only as a warm start for certified one-sided
  Golub-Kahan/LBD polish, and records scout steps, side, matvecs, operator
  matvecs, polish cost, and the fact that scout certificates are not trusted.
  The H-shaped wide sparse probe certifies only after the final SVD
  certificate and is slower/higher-allocation than direct one-sided GK and
  RSpectra, so this is rejection evidence rather than an H promotion.
- Shift-invert factorization-cache keys now include the problem structure and
  standard/generalized `B` boundary in addition to the operator fingerprint and
  sigma. This closes the cache invalidation contract for structure changes
  while the shift-invert path remains honestly labelled as reference until the
  native hot loop lands.
- Hermitian shift-invert now has native factorized Lanczos hot loops for
  dense standard, dense generalized-SPD, diagonal standard, sparse
  symmetric-tridiagonal standard, and sparse/diagonal tridiagonal generalized
  problems. Standard dense
  `shift_invert(sigma)` plans factor `A - sigma I` with LAPACK `dgetrf`; dense
  generalized plans factor dense `B` with LAPACK `dpotrf`, factor
  `A - sigma B` with `dgetrf`, and apply
  `R solve(A - sigma B, R' y)` inside the native scalar Lanczos recurrence.
  Diagonal and sparse tridiagonal standard plans factor `A - sigma I` with a
  native Thomas recurrence and apply shifted solves inside native scalar
  Lanczos. Tridiagonal generalized plans with diagonal `B` factor
  `A - sigma B` with the same native Thomas recurrence and apply
  `sqrt(B) solve(A - sigma B, sqrt(B) y)`. These paths
  recover `lambda = sigma + 1 / mu` and certify residuals on the original
  problem. General sparse, general sparse diagonal-metric generalized, and
  user-solve shift-invert paths remain honestly reference-labelled until
  broader native non-densifying factorization ownership lands; the user-solve
  benchmark boundary records its external cache provenance explicitly.
- Shift-invert now has a release-hardening benchmark surface at
  `inst/benchmarks/bench-shift-invert.R`. It reports method labels,
  native/reference factorization-cache provenance, original-coordinate
  certificate status, and contract rows for dense standard, dense generalized
  SPD, diagonal standard, sparse tridiagonal standard, and sparse/diagonal
  tridiagonal-generalized cases. A fresh
  installed non-quick `--iterations=1 --save` pass certified all requested
  rows (`6/6`) and passed the contract gate after the generalized tridiagonal
  promotion: dense standard used `dense_lu_native`, dense generalized used
  `dense_lu_generalized_native`, diagonal and sparse tridiagonal standard used
  `tridiagonal_thomas_native`, and sparse/diagonal tridiagonal generalized used
  `tridiagonal_thomas_generalized_native`. L remains partial because general
  sparse native shift-invert is still not implemented.
- Dense nonsymmetric eigen oracle results received right-residual
  certificates even when LAPACK returned complex eigenpairs. This was the
  compatibility bridge before dense explicit nonsymmetric paths moved onto the
  native Arnoldi-cycle plus native-Ritz label.
- Nonsymmetric eigen now has a release-hardening benchmark surface at
  `inst/benchmarks/bench-nonsymmetric.R`. It checks real non-normal direct
  `eig_partial()`, complex-pair direct `eig_partial()`, and
  RSpectra-compatible `eigs(..., which = "LI")` rows, and the contract records
  right-residual certification, non-orthogonality semantics, planner label
  honesty, restart diagnostics, and native Arnoldi/Ritz provenance.
- Nonsymmetric matrix-free and sparse real-spectrum auto paths initially gained
  an honest `reference Arnoldi (prototype/oracle fallback)` route. This avoided
  pretending a dense oracle could service matrix-free operators or silently
  densifying sparse general matrices when a Krylov path was available. The
  path carried nonnative restart diagnostics and right-residual certification;
  it was compatibility evidence, not the native restarted Arnoldi
  implementation required for V1 promotion.
- Dense and sparse CSC nonsymmetric real- and imaginary-target auto paths now run a native Arnoldi
  cycle and native projected Ritz extraction instead of the R-level reference
  cycle. The compatibility label is
  `native Arnoldi cycle + native Ritz extraction (compatibility)`: basis
  expansion, two-pass orthogonalization, the projected Hessenberg `dgeev`
  solve, and Ritz-vector formation run in native code, while exact
  right-residual certification remains in the existing R certificate layer.
  The native compatibility default uses a bounded larger subspace
  (`min(n, max(k + 8, 9k))`) after the installed non-quick
  `sparse_native_arnoldi_lr:80` row exposed that the old `2k + 4` default did
  not certify. The planner now wires a restart budget into native sparse
  Arnoldi auto paths, and restart loops retain the best attempt by certificate
  pass state, convergence count, and backward error instead of blindly
  returning the last attempt. Dense explicit nonsymmetric `eigs(..., which =
  "LI")` now shares this native compatibility route. Matrix-free nonsymmetric operators remain
  reference-labelled, and production-grade fully native restarted Arnoldi
  remains open. Installed non-quick strict evidence from
  `R_LIBS=/tmp/eigencore-bench-lib Rscript inst/benchmarks/bench-nonsymmetric.R --iterations=1 --save --strict`
  certifies `dense_native_arnoldi_lm`, `dense_native_arnoldi_li`,
  `dense_eigs_native_arnoldi_li`, `sparse_native_arnoldi_lr`, and
  `sparse_native_arnoldi_li`; all rows report `native_arnoldi_label = TRUE`,
  `arnoldi_native = TRUE`, and `ritz_extraction_native = TRUE`. Dense explicit
  native Arnoldi uses the full dense subspace by default (`80` matvecs for
  `dense_native_arnoldi_lm`, `40` for the dense `LI` rows), while sparse CSC
  rows keep the bounded restart policy (`72` matvecs for
  `sparse_native_arnoldi_lr`, `56` for `sparse_native_arnoldi_li`). Saved
  installed artifacts:
  `inst/benchmarks/results/20260517-nonsymmetric-rows.rds` and
  `inst/benchmarks/results/20260517-nonsymmetric-contracts.rds`.
- A follow-up quick installed strict smoke after restart-control wiring also
  passed on 2026-05-17 with the same restart diagnostics for the sparse native
  Arnoldi row.
- Sparse CSC `eigs(..., which = "LI")` now uses the same native Arnoldi-cycle
  plus native-Ritz compatibility path. Complex Ritz-vector certification uses
  the explicit source matrix for right-residual checks so the solver does not
  fall back to the real-only native CSC apply wrapper. A fresh installed
  non-quick strict saved run added and certified the `sparse_native_arnoldi_li`
  row with `nconv = 2/2`, `matvecs = 56`, `restart_count = 3`,
  `certified_attempt = 4`, and `ritz_extraction_native = TRUE`.
- The nonsymmetric benchmark row now also exposes
  `stage_arnoldi_cycle_seconds`, `stage_ritz_extraction_seconds`, and
  `ritz_extraction_native`. These make the current compatibility boundary
  auditable: sparse CSC real- and imaginary-target rows now report
  `ritz_extraction_native = TRUE`, so production-grade native restarted Arnoldi
  remains open on restart policy, matrix-free/native callback support, and final performance
  gates rather than on the projected Ritz extraction itself.
- The randomized SVD milestone now has an explicit `rsvd` parity benchmark
  surface at `inst/benchmarks/bench-randomized-rsvd.R`. It compares
  `eigencore_randomized` against `rsvd` using oracle singular-value error,
  left/right subspace error, true SVD certificate fields, and
  time-to-certified-answer gates.
- Randomized SVD benchmark dispatch now uses `oversample = 10`, matching the
  PRD default and `rsvd`'s default `p = 10`, instead of scaling oversampling
  with the requested rank. The dense explicit randomized path also projects
  with direct `Q' A` rather than `A' Q` plus transpose. On the installed
  package large exact-low-rank dense parity row (`2000 x 500`, rank `50`,
  three iterations), eigencore certifies with identical singular-value and
  subspace accuracy and reports about `1.13x` speed versus certified `rsvd`.
  Randomized SVD now also reuses the cached projected matrix for the
  `A' u - sigma v` side of its certificate, so only the `A v - sigma u` side
  needs a full operator apply after the small SVD. The same installed-package
  parity row now reports about `1.18x` speed versus certified `rsvd`; the
  randomized `2x` release speed gate remains open.
- Randomized SVD now has a conservative certified early-stop for
  QR-normalized range finding: after the initial `q = 0` sketch it computes the
  same residual/backward-error certificate and returns immediately only if that
  certificate passes. On the installed-package large exact-low-rank dense
  parity row (`2000 x 500`, rank `50`, three iterations), this reaches about
  `3.08x` speed versus certified `rsvd`, with the accuracy gate still green.
  The quick slow-decay dense row remains a performance gap, so the broader
  randomized release gate is not closed.
- Historical pre-demotion G1 note: native block Hermitian Lanczos was promoted
  for previously benchmark-proven regimes. Fresh 2026-05-17 installed evidence
  invalidated sparse auto-promotion, so the default `eigencore` sparse path no
  longer claims these rows as current release evidence. Representative older
  installed-package gate rows:
  path Laplacian `n = 1000, k = 20` certified with block size `2`, about
  `1.42x` to `1.49x` faster than certified RSpectra, memory ratio about
  `1.79x`, and PRIMME parity above `5x`;
  path Laplacian `n = 10000, k = 20` certifies with block size `4` and
  `max_subspace = 320`, about `2.77x` to `2.99x` faster than certified
  RSpectra, memory ratio about `1.65x`, and PRIMME parity above `3.8x`;
  dense Hermitian `n = 200, k = 20` certifies through the small-dense
  full-subspace native path, about `1.42x` faster than certified RSpectra,
  memory ratio about `1.08x`, and PRIMME parity above `9x`.
  The native block path includes structured projected recurrence, capped Ritz
  padding, sparse CSC Ritz-vector residual application, in-solver certificate
  reuse in the benchmark harness, selected-Ritz workspace shrink, selected
  dense LAPACK extraction for small-dense full-subspace fallback, tiny-block
  native projection updates, one-pass reorthogonalization for `n >= 64`,
  upper-triangle projected-matrix copies, read-only restart Ritz residual
  norming that preserves `A * V` blocks for retained Ritz vectors, guarded
  native final subspace polishing with residual reuse, residual-refresh and
  Rayleigh-Ritz fallback paths, best-so-far restart snapshot recovery, and
  native malloc-based workspaces where zero-fill is unnecessary.
- Native standard LOBPCG with the shifted tridiagonal preconditioner passes
  `bench-lobpcg-preconditioned.R --strict` on path Laplacians
  `n = 200, 1000, 2000`, `k = 5`, with failed references recorded as
  uncertified rows. The selected `eigencore_lobpcg_tridiagonal` row certified
  all requested pairs, used the native shifted-tridiagonal preconditioner, and
  beat scalar eigencore by about `24x`, `154x`, and `76x` on the three sampled
  path-Laplacian rows while beating the best certified external reference by
  about `2.5x`, `9.1x`, and `7.9x`.
- Native generalized SPD LOBPCG now has built-in dense, diagonal, and CSC A/B
  prototype slices for explicit `lobpcg()` plans. Native shifted-diagonal and
  shifted-tridiagonal typed preconditioners are accepted by the generalized
  native path; broader generalized preconditioning and promotion gates remain
  open.
- The generalized LOBPCG adversarial B bank now checks largest and smallest
  native paths for ill-conditioned diagonal B, sparse CSC B, and explicitly SPD
  matrix-free B, with `allow_dense_fallback = "never"` and B-orthonormality
  checked against dense oracle values.
- Generalized SPD benchmark rows now expose the promotion-critical restart
  diagnostics: `native`, `native_kernels`, `generalized`,
  `orthogonalization_native`, `orthogonalization_methods`, `q_rank_final`,
  `constrained`, and `constraints_rank`. The quick sparse generalized smoke
  check now verifies that eigencore rows certify through the native
  B-orthogonal path rather than silently using dense/reference fallbacks.
  The strict generalized LOBPCG gate now also emits native-contract rows for
  bare, shifted-diagonal preconditioned, shifted-tridiagonal preconditioned
  sparse-smallest, and constrained generalized paths, and requires their
  native/preconditioner/constraint diagnostics to pass. It now also emits
  adversarial B contract rows for ill-conditioned diagonal, sparse CSC, and
  explicitly SPD matrix-free B cases across largest and smallest targets,
  including exact expected native B-orthogonalization methods.
- A source-loaded non-quick strict generalized LOBPCG run showed that dense
  generalized `auto()` should not promote to native iterative LOBPCG yet. Dense
  generalized `auto()` now routes to the native dense generalized SPD LAPACK
  fallback, while explicit `method = lobpcg()` still exercises the native
  generalized LOBPCG slice for diagnostics and requested iterative runs. Treat
  the quick strict run as contract smoke only; J remains open until the
  non-quick sparse and broader production gates certify and pass speed/memory
  thresholds.
- The generalized LOBPCG benchmark now supports stable `--cases=` filtering,
  `--methods=`, `--subject`, and per-case progress messages. It also separates
  dense `eigencore_auto` fallback boundary rows from native LOBPCG contract
  rows. Fresh installed focused non-quick evidence shows dense auto fallback
  rows certify as non-performance-gated boundary checks, and
  `sparse_generalized_path_smallest:500` passes when gated on the typed
  shifted-tridiagonal subject (`10/10`, `1.04x` speed, `2.79x` memory versus
  dense base). The sparse-largest benchmark now also includes a target-aware
  shifted-tridiagonal preconditioner row whose largest-target shift is
  non-densifying and scale-aware. Fresh installed non-quick evidence certifies
  `sparse_generalized_path_largest:500` at `10/10` with native generalized
  kernels and shifted-tridiagonal provenance, and passes memory (`2.75x`
  versus dense base), but the row remains performance-red on speed (`0.49x`
  speed, `163` iterations). J remains partial because sparse-largest and
  broader generalized production gates are still red.
