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
- G1 native block Hermitian Lanczos is promoted for benchmark-proven regimes.
  The default `eigencore` path now passes both strict release scripts with dense
  diagnostics included:
  `Rscript inst/benchmarks/bench-native-hermitian-gate.R --include-dense --strict`
  and
  `Rscript inst/benchmarks/bench-hermitian-sparse.R --include-dense --strict`.
  Representative installed-package gate rows:
  path Laplacian `n = 1000, k = 20` certifies with block size `2`, about
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
  uncertified rows.
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
  bare, shifted-diagonal preconditioned, and constrained generalized paths, and
  requires their native/preconditioner/constraint diagnostics to pass.
