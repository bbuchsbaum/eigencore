# eigencore V1 Plan — Refinement of `prd.json`

This document refines `prd.json` after two review passes. It does **not**
replace the PRD. The PRD defines scope and acceptance; this plan defines the
*engineering discipline* required to reach that scope honestly, and sequences
the work so that claims in `vision.md` and `mission.md` hold at every
checkpoint.

Governing principle:

> Treat the current R code as the **mathematical specification layer**,
> not the production solver layer. Freeze the R prototypes as reference
> oracles. Move production work into a native engine driven by a frozen
> operator ABI.

## 1. Three-layer architecture

The codebase is split into three strictly separated layers. Nothing in Layer 1
may call Layer 2 in production paths. Nothing in Layer 2 may be exported.

- **Layer 1 — R public API** (`R/problem.R`, `R/classes.R`, `R/compatibility.R`,
  `R/results.R`, `R/certification.R` public surface).
  Problem objects, target/method descriptors, plan objects, certificate
  objects, result accessors, RSpectra shims. No iteration, no numerics
  beyond small dense solves.

- **Layer 2 — R reference solvers** (`R/reference_*.R`, not exported).
  Rename `lanczos_hermitian` → `reference_lanczos_hermitian`, `golub_kahan_svd`
  → `reference_golub_kahan_svd`. These exist only as executable math and as
  oracles in adversarial tests. Production code paths never dispatch here.

- **Layer 3 — native engine** (`src/`).
  Operators, workspaces, orthogonalization, projection, solvers, residual
  and certificate kernels. Consumes the frozen operator ABI, owns all hot
  loops, performs no R allocations inside iteration.

**Rule:** any code path reachable from `eig_partial()` / `svd_partial()` /
`solve()` must either run Layer 3 or emit a planner label that identifies the
path as a reference or oracle fallback. This preserves planner honesty as the
engine grows.

## 2. Frozen operator ABI

Before adding solvers, freeze the native operator contract:

```cpp
enum class Transpose  { None, Adjoint };
enum class Structure  { General, Hermitian, Diagonal, Symmetric, Triangular };
enum class ScalarType { F64, C128 };           // V1 ships F64 only

struct Workspace;                              // opaque, caller-owned

struct BlockOperator {
    int64_t    rows;
    int64_t    cols;
    ScalarType scalar;
    Structure  structure;
    bool       has_adjoint;
    double     frobenius_upper;                // exact if known, else NaN
    double     two_norm_upper;                 // exact if known, else NaN
    void*      impl;

    // Y := alpha * op(A) X + beta * Y   (op = None | Adjoint)
    // X is (inner x block_cols) with leading dim ldx
    // Y is (outer x block_cols) with leading dim ldy
    // No allocation permitted inside this call.
    int (*apply)(void* impl, Transpose op, int64_t block_cols,
                 const double* X, int64_t ldx,
                 double alpha, double beta,
                 double* Y, int64_t ldy,
                 Workspace* ws);
};
```

V1 built-in implementations required:

- `DenseColumnMajorOperator`          (BLAS `dgemm` / `dsymm`)
- `CSCOperator`                       (block SpMM; already a stub)
- `CSROperator`                       (transpose path of CSC; optional V1)
- `DiagonalOperator`                  (scalar axpy)
- `ScaledOperator` / `CenteredOperator` (fused; no R closures)
- `SumOperator` / `ComposedOperator`  (native dispatch, no temporaries)
- `RCallbackOperator`                 (slow path; planner warns when used
                                       inside a block-native solver)

Contract obligations:

- `apply` never allocates. All workspace comes from `Workspace*`.
- `frobenius_upper` / `two_norm_upper` are exact or NaN; estimates go into a
  separate field with a `norm_bound_type` tag (see §5).
- `has_adjoint = false` forbids SVD solvers at plan time (hard error, not
  silent fallback).

## 3. Solver delivery sequence (native-first)

Deliver the native engine in order. Each step ships with adversarial tests
(see §7) and a benchmark row.

1. **Native MGS2** + **CholQR2** with workspace reuse, B-inner-product
   variant for generalized problems.
2. **Native residual / certificate kernels** for eigen, generalized eigen,
   SVD. These run in the native layer on native buffers; R sees only a
   finished certificate struct.
3. **Native Hermitian block Lanczos with thick restart and locking.**
   Block size ≥ 1 configurable; default chosen by planner. Ritz extraction
   is the canonical one: solve the projected problem in the basis, then
   `V_k = Q · S_k` *without* independent post-orthogonalization. If
   orthogonality loss exceeds tolerance, re-extract — never patch with
   independent MGS. A scalar `block = 1` implementation is an acceptable
   staging step only if planner labels say so; it does not satisfy the block
   milestone or the performance gate.
4. **Native Golub–Kahan bidiagonalization with thick restart.** True SVD
   residuals `(||A v − σ u||, ||A* u − σ v||)` computed natively. No
   `A^T A` materialization. No independent MGS on `U` and `V` after
   extraction.
5. **Native randomized SVD** (range finder + subspace iteration), with
   honest certificate (`certificate_type = "estimated"` unless refinement
   pass runs and lowers residual below tolerance). Use
   `docs/hegelsvd_svd_acceleration.md` as the design note for mining
   HegelSVD's adaptive QB/PCS strategy, planner regimes, and benchmark cases
   without inheriting its dense-only assumptions. Treat `vendor/rsvd` as the
   minimum public-package bar: eigencore's native randomized path must match or
   improve `rsvd::rsvd()` singular-value/subspace accuracy on the randomized
   benchmark bank and beat it on time-to-certified-answer in the regimes where
   the planner is allowed to choose randomized SVD. The enforcement surface is
   `inst/benchmarks/bench-randomized-rsvd.R`, which reports oracle
   singular-value error, left/right subspace error, true SVD certificates, and
   speed ratios against `rsvd`. As of 2026-04-28, the R reference randomized
   path uses the PRD/rsvd-style `oversample = 10` benchmark default and a
   direct dense `Q' A` projection, moving the large exact-low-rank dense parity
   row from slower than `rsvd` to roughly `1.13x` faster on
   time-to-certified-answer. The randomized certificate now also reuses the
   cached projected matrix for the right residual (`A' u - sigma v`), reducing
   the initial certificate to one full post-SVD operator apply plus cheap
   projected algebra; the same parity row moved to roughly `1.18x` faster than
   certified `rsvd`. A conservative certified early-stop now checks the `q = 0`
   range-finder result for QR-normalized randomized SVD and exits only when the
   residual/backward-error certificate already passes. On the installed
   large exact-low-rank dense parity row this reaches roughly `3.08x` versus
   certified `rsvd`, closing the `2x` gate for that randomized-planner regime.
   The broader randomized release gate remains open for slow-decay and other
   non-exact cases and is expected to require native sketch/projection kernels
   and/or stronger adaptive randomized planning.
6. **Generalized SPD LOBPCG** for dense, sparse, and matrix-free operators.
   This is the primary scalable V1 path for `A x = lambda B x`, especially
   smallest eigenpairs and preconditioned problems. It uses block iteration,
   B-orthogonal Rayleigh-Ritz extraction, optional constraints/deflation, and
   native residual/certificate kernels. The planner maps largest generalized
   eigenpairs through the sign-flipped problem when appropriate; it never
   densifies sparse `A` or `B` silently. The design contract is tracked in
   `docs/native-generalized-spd-lobpcg.md`.
7. **Generalized SPD B-orthogonal block Lanczos** as an alternate/refinement
   path after LOBPCG. Cholesky paths treat `B` as an operator or factorization
   cache; never invert `R` explicitly.
8. **Shift-invert** with user-supplied solve / factorization-cache
   operator. Fails loud if no solve is provided; never silently dense.
9. **Standard symmetric LOBPCG refinements** for graph Laplacians and
   preconditioned standard Hermitian problems, sharing the generalized LOBPCG
   core.

Until a step lands natively, its planner label is
`"reference <method> (prototype — not production)"`. The existing honest-label
machinery is the template.

### Current status checkpoint

As of the current native Hermitian and certificate work, eigencore has:

- native dense, CSC, and diagonal block-apply operators;
- native MGS2 / CholQR2 / dense B-CholQR2 orthogonalization kernels;
- native dense eigen/SVD certificate diagnostics;
- native standard built-in-operator certificate diagnostics for dense, CSC,
  and diagonal operators;
- a certified `native scalar thick-restart Hermitian Lanczos` path for dense
  double matrices and `dgCMatrix` operators, with native locking metadata and
  sparse non-densification;
- a promoted `native block Hermitian Lanczos (thick restart, locking)` path for
  dense double matrices and `dgCMatrix` operators in benchmark-proven regimes,
  with native block operator application, thick restart, native locking
  metadata, adaptive block/subspace controls, and scalar fallback for small `k`
  or unproven regimes;
- a certified native Golub-Kahan staging path for dense double and `dgCMatrix`
  SVD, with adaptive subspace growth metadata.
- an internal reference LOBPCG spike showing that a shifted sparse solve
  preconditioner can reduce Laplacian smallest-eigenpair iteration count by an
  order of magnitude; this is evidence for the native preconditioner interface,
  not a production solver path.
- an explicit `lobpcg()` method descriptor and
  `shifted_cholesky_preconditioner()` helper. The current solver route is
  labeled `reference LOBPCG prototype` and runs in R, but it returns the normal
  eigencore result/certificate object.
- built-in preconditioners now carry typed metadata (`kind`, `native`,
  `factorization`, `shift`) and LOBPCG results report preconditioner call
  counts. This makes the planner and benchmark gates inspectable before the
  native LOBPCG loop replaces the R prototype.
- a `native standard Hermitian LOBPCG prototype` now supports dense double and
  `dgCMatrix` standard Hermitian problems, with optional native shifted
  tridiagonal preconditioning. The planner uses this label only for supported
  built-in cases; generalized SPD and opaque preconditioners remain on
  reference labels.

The Hermitian path supports native targets `largest`, `smallest`,
`largest_magnitude`, and `smallest_magnitude`; unsupported targets such as
`nearest()` route away from the native path with an honest planner label.

Milestone G is now split cleanly: G0 was the scalar-native staging path, and G1
is the promoted block-native Hermitian path for benchmark-proven regimes. The
quick `k = 5` smoke gates remain diagnostic rather than release gates; in
`--quick --strict` mode the Hermitian benchmark scripts enforce certification
only, while non-quick `--strict` enforces the speed/memory/parity release gate:

- sparse Hermitian path Laplacian `n = 200, k = 5`: eigencore certifies all
  requested pairs with the scalar-native path and the tuned default
  `max_subspace = 3*k + 20`, but is still slower than RSpectra and
  allocates more memory on the quick gate. In the current source-loaded quick
  gate, scalar eigencore certifies but reaches only about `0.35x` of RSpectra's
  median time and about `0.66x` of PRIMME's median time;
- the block Hermitian path now supports explicit restart and locking for
  dense double and `dgCMatrix` inputs, and explicit
  `lanczos(block > 1, max_subspace < n)` can certify the quick path-Laplacian
  case. The benchmark harness now reports a subject gate for
  `eigencore_block_candidate`. The promoted auto policy keeps scalar Lanczos for
  small `k`, uses `block = 2` for `k >= 16` medium sparse rows and small dense
  full-subspace regression rows, and uses `block = 4` with at least `16*k`
  restart space for large sparse `k >= 16` rows, with a four-vector capped Ritz
  pad. The block path now
  forms block-recurrence residuals before reorthogonalization, maintains a
  structured projected problem instead of recomputing the full projection at
  each restart, uses allocation-free CholQR2 block acceptance when rank permits,
  reports native stage timings, and uses fair solver/certificate/total
  benchmark accounting. Additional tuning now caps the thick-restart Ritz pad,
  uses the sparse CSC block apply for Ritz-vector residuals, reuses in-solver
  certificates in the benchmark harness, uses `dsyevd` only for larger
  projected dense solves, and avoids zero-filling native workspaces that are
  immediately overwritten. It also uses native loops for tiny block projection
  updates and one reorthogonalization pass for `n >= 64` when certificates stay
  tight. It also copies only the upper triangle into the
  compact projected matrix passed to LAPACK and runs a native final subspace
  polish in the returned subspace to satisfy the orthogonality certificate
  without R-tracked allocations. The polish first reuses the existing locked
  residuals when the returned basis already satisfies the certificate
  orthogonality tolerance; otherwise it tries a cheap Cholesky-orthonormalization
  plus residual refresh and falls back to final Rayleigh-Ritz when needed.
  Final polish now preserves any genuine locked prefix when the restart budget
  is exhausted, so a noncertified tail cannot rotate away already certified
  vectors. The main restart residual check also leaves `A * Ritz` blocks
  intact, so retained Ritz vectors can carry `AV_active` forward without
  reconstructing it from residuals. Native restart history now records, per
  Rayleigh-Ritz cycle, the active subspace size, selected count, lock count,
  wanted-pair convergence count, and max backward error; this is part of the
  G1 diagnostic contract.
  Certificate `passed` now includes
  an orthogonality gate, not just per-pair residual convergence. The native
  lock path rejects converged Ritz vectors that are numerically dependent on
  already locked vectors, so duplicate locked vectors no longer earn residual
  convergence while failing orthogonality. Default small-dense block runs with
  `n <= getOption("eigencore.block_dense_full_subspace_max_n", 256)` take an
  honestly labeled native full-subspace LAPACK/Rayleigh-Ritz path, closing the
  small dense, clustered, and ill-conditioned diagonal rows in the G1 baseline;
  explicit bounded `max_subspace` still exercises restart. The
  installed-package quick gates now default to the sparse path-Laplacian
  release case. Dense Hermitian rows remain available as opt-in diagnostics
  with `--include-dense`, but no longer block the sparse G1 gate. The quick
  path row still fails the RSpectra speed bar while passing memory and PRIMME
  parity. A k=20, n=1000 installed staging row with the smaller block-2
  controls certifies, passes memory/PRIMME parity, and reaches about `1.45x` to
  `1.50x` of RSpectra's median time on local installed runs. A larger
  `n = 10000, k = 20` path-Laplacian probe now certifies with the adaptive
  block-4 candidate (`max_subspace = 320`), reaches about `3.0x` of a certified
  robust RSpectra reference row, and passes PRIMME parity in the benchmark
  harness. The full dense Hermitian regression row now certifies and passes the
  RSpectra speed, PRIMME parity, and memory gates. The default `eigencore`
  strict gates pass for `bench-native-hermitian-gate.R --include-dense --strict`
  and `bench-hermitian-sparse.R --include-dense --strict`;
- preconditioned reference LOBPCG on the same quick Laplacian case certifies
  in about 10 iterations. The current gate shows it is about 1.8x faster than
  scalar eigencore, slightly faster than PRIMME on some runs, still slower than
  RSpectra, and allocates materially more memory because the path is R-level
  and factorization-backed;
- a native shifted tridiagonal preconditioner now covers path-Laplacian-style
  tridiagonal gates. It preserves the 10-iteration convergence behavior and
  reduces memory versus the Cholesky-backed reference setup, but the current
  native standard LOBPCG prototype now keeps the hot iteration in C++ for this
  case. The staging gate now covers path Laplacians at
  `n = 200, 1000, 2000` with `k = 5`; it certifies, beats scalar eigencore,
  beats the best certified external reference on the sampled run, and passes
  the restored `0.25` memory threshold. A direct `n = 10000, k = 20`
  path-Laplacian probe with the native shifted-tridiagonal preconditioner
  certifies 20/20 pairs in about `4.8s` locally, so that solver is a strong M
  milestone candidate but does not satisfy G1's block-Lanczos promotion
  contract. Benchmark rows and the LOBPCG gate record which typed
  preconditioner candidate was selected, and failed references are captured as
  uncertified rows instead of aborting the gate;
- native generalized SPD LOBPCG remains deliberately unpromoted for production
  use. A design contract now exists, residual-backed generalized certificate
  helpers preserve the original-coordinate formula, and native dense, diagonal,
  and CSC A/B loop slices are wired for explicit `lobpcg()` plans. Matrix-free
  `B`, generalized preconditioners, and production promotion remain open;
- sparse tall-skinny SVD `500 x 80, rank = 5`: eigencore certifies but is
  still materially slower than RSpectra and `irlba`.
- the SVD surface benchmark is now explicit at
  `inst/benchmarks/bench-svd-surface.R`, covering tall/wide sparse,
  rank-deficient sparse, clustered dense, slow-decay dense, and low-rank sparse
  cases against RSpectra, PRIMME, `irlba`, `rsvd`, and base where applicable.
  The dense double path can now use the same bounded certified Gram special
  case as sparse CSC, and the Gram solver uses native selected symmetric
  eigensolves instead of a full eigensystem for largest singular values. This
  also fixed a native `dsyevr` buffer bug in the selected-eigen wrapper. The
  saved 2026-04-26 quick snapshots show eigencore certifies every subject row,
  with current speed ratios versus the best certified reference roughly in the
  `0.27x` to `0.70x` range on the latest quick run and one sparse memory gate
  passing. H remains a
  performance gap, not a correctness gap; the release answer still requires a
  production thick-restarted Golub-Kahan path rather than relying on Gram
  special cases. The benchmark now also reports an internal
  `eigencore_golub_kahan` candidate row. That row is excluded from external
  release-reference gates, but it makes promotion readiness visible: it is
  faster than the Gram subject on the current low-rank sparse row, and after
  fixing planner controls so default Golub-Kahan remains adaptive rather than
  accidentally fixed-budget, it certifies tall/wide random sparse rows too.
  Its adaptive work accounting now counts every from-scratch retry, not only
  the final successful subspace, and benchmark rows expose retry attempts,
  final-vs-total iterations/matvecs, and native iteration/Ritz timing. A shared
  planner/native initial-subspace heuristic now starts wide rectangular
  Golub-Kahan runs at a larger budget; on the quick wide sparse row this removed
  the previous second from-scratch attempt, reducing the candidate from `200`
  total matvecs to `120`. The native candidate can also record opt-in sampled
  prefix diagnostics from the final basis
  (`options(eigencore.golub_kahan_prefix_diagnostics = TRUE)`), including the
  first sampled prefix that certifies and the final-subspace iteration/matvec
  overshoot. This tells us whether the next C++ loop should stop earlier inside
  one attempt or whether the current subspace size is genuinely needed, without
  charging normal benchmark timing for extra diagnostic certificates. A staged
  native projected-stop path is also available behind
  `options(eigencore.golub_kahan_projected_stop = TRUE)`: the C++ loop checks a
  selected bidiagonal SVD residual estimate at coarse intervals and can stop
  before filling the planned subspace, while the normal result certificate still
  decides whether adaptive growth is needed. A quick probe on the wide sparse
  fixture reduced the candidate from `60` iterations / `120` matvecs to `50` /
  `100` while preserving the final certificate; tall sparse did not stop early
  under the same criterion. This remains opt-in until quick and full SVD gates
  show it improves time-to-certified-answer. The SVD surface benchmark can now
  include an explicit `eigencore_golub_kahan_projected` row with
  `--projected-stop`, so normal release gates stay unchanged while on/off
  candidate comparisons are saved in the same surface table. It also accepts
  `--subject=<method>`, so H candidate rows such as
  `eigencore_golub_kahan_projected` can be gated directly against external
  references without changing the default `eigencore` release gate. The
  `--h-candidate` preset now selects the projected Golub-Kahan row as the gate
  subject, includes the plain Golub-Kahan row for on/off comparison, includes
  external references, and fails early if a requested `--subject` is missing
  from the selected methods. The same surface now prints and saves
  `svd-surface-memory` diagnostics, separating total, solver, and certificate
  allocation gaps against the best certified reference so H memory work can
  target the right phase instead of relying on one aggregate memory ratio. The
  first quick projected-stop surface run certified the projected candidate rows
  and showed useful but uneven movement: wide sparse dropped from about `60/120`
  iterations/matvecs to `50/100`, clustered dense from `44/88` to `24/48`,
  slow-decay dense from `44/88` to `36/72`, while tall, rank-deficient sparse,
  and low-rank sparse did not stop early. The stop policy now also records
  projected-check counts and time, skips checks that cannot save at least
  `max(5, rank)` iterations, and conservatively disables projected stopping for
  high-aspect tall sparse operators until evidence shows it helps there. The
  quick comparison table reports that disabled reason explicitly; on the current
  quick fixture, tall sparse has zero projected checks, clustered/slow-decay
  dense and wide sparse still stop early, and rank-deficient/low-rank sparse
  run with zero useful projected checks. This is an improvement over the first
  blind opt-in path, but the one-iteration quick table still shows timing noise
  on cases without iteration savings, so promotion requires repeated quick and
  full-surface evidence rather than this heuristic alone. The Ritz extraction
  path now accepts an explicit active iteration count in native code, so the R
  wrapper no longer copies `U[, used]` and `V[, used]` before extraction or
  prefix diagnostics. On the quick surface this cut the internal Golub-Kahan
  candidate's memory materially, for example tall sparse from roughly `594kB`
  to `370kB`, wide sparse from roughly `896kB` to `562kB`, clustered dense from
  roughly `298kB` to `233kB`, and low-rank sparse from roughly `203kB` to
  `174kB`, while preserving certificates. The native Golub-Kahan entry points
  now also keep the planned `maxit` basis in native C++ workspace and return
  only the realized `iterations` prefix to R, with `native_workspace_bytes`
  exposed in restart metadata and benchmark rows. A quick H candidate check
  after that change still certified the projected rows and reduced R-visible
  solver allocation further: wide sparse `519920 -> 464576` bytes and
  clustered dense `176768 -> 147664` bytes on the sampled cases. The remaining
  memory failure is now explicitly the gap between native workspace/result
  materialization and the external reference allocations, not certificate
  recomputation; the next H allocation work should reuse the per-iteration
  `v/z/u/u_prev` scratch and avoid copying or retaining basis columns that are
  no longer needed for the final Ritz result. The next pass added compact
  native Golub-Kahan "fit" entry points that run the small Ritz projection
  before crossing back into R, so the default non-diagnostic path returns only
  final `d/u/v` plus metadata instead of the full Krylov basis. Prefix
  diagnostics still force the basis-returning path. On the same quick H sampled
  cases, projected solver allocation dropped again to `187600` bytes for wide
  sparse and `112544` bytes for clustered dense, with `basis_returned = FALSE`
  and certificates still passing. H is still not closed: those rows remain below
  the speed gate and above PRIMME/RSpectra memory in the quick fixture, so the
  remaining work is algorithmic speed and eliminating smaller R-visible result
  allocations rather than certificate cost. Dense Golub-Kahan now uses BLAS
  `dgemv` block projections for two-pass reorthogonalization, while CSC keeps
  the scalar projection loop because the sampled wide sparse row regressed with
  BLAS call overhead. A tighter projected-stop cadence was tested and rejected:
  it increased projected checks without reducing iterations on the sampled
  cases. The 3-iteration quick H check still certifies the sampled rows, with
  projected medians around `0.0037s` for wide sparse and `0.0011s` for clustered
  dense, so this is only an incremental kernel cleanup, not an H gate closure.
  Native Golub-Kahan now reports `apply`, `recurrence`, `reorthogonalization`,
  and `projected_solve` stage timings into the benchmark rows. The same quick
  H sample shows the real hotspot: projected wide sparse spends roughly
  `0.0020s` in reorthogonalization versus `0.00012s` in operator application,
  while clustered dense spends roughly `0.00009s` in reorthogonalization and
  `0.00016s` in application. This makes the next H closure criterion explicit:
  the production SVD path must reduce or restructure orthogonalization through
  block/thick-restart Golub-Kahan rather than further shaving certificate or
  callback overhead. As an interim scalar-path improvement, native Golub-Kahan
  now uses a DGKS-style adaptive second reorthogonalization pass: every vector
  gets one projection pass, and the second pass runs only if the first pass
  shrinks the candidate vector enough to signal orthogonality loss. The quick H
  sampled rows remain certified and now report `reorthogonalization_passes`:
  projected wide sparse used `99` passes for `50` iterations and improved to
  about `0.0028s`, with reorthogonalization around `0.0011s`; clustered dense
  used `47` passes for `24` iterations and remained around `0.0012s`. This is
  useful headroom, but still not an H gate closure. The SVD benchmark surface
  now records normalized H hotspot metrics: accounted native stage time,
  reorthogonalization time fraction, seconds per reorthogonalization pass,
  passes per iteration, native seconds per matvec, projected-check cost, and
  projected-stop savings fractions. These
  fields make block/thick-restart SVD progress measurable by orthogonalization
  work avoided, not only by noisy wall-clock medians. A full non-quick projected-stop SVD
  surface run on 2026-04-26 saved
  `inst/benchmarks/results/20260426-svd-surface-rows.rds`,
  `inst/benchmarks/results/20260426-svd-surface-gates.rds`, and
  `inst/benchmarks/results/20260426-svd-projected-stop-comparison.rds`.
  That repeated run shows projected stopping is useful where it actually
  reduces native work: clustered dense and slow-decay dense improved by about
  `1.24x` (`100/200` to `80/160` iterations/matvecs), and wide sparse improved
  by about `1.23x` (`180/360` to `160/320`). The high-aspect tall sparse,
  low-rank sparse, and rank-deficient sparse cases record the disabled reason
  or no useful projected checks. Release gates still fail overall: wide sparse
  is the only current SVD surface row meeting the speed gate, and it still
  fails memory; the other rows remain speed/memory failures. The
  rank-deficient sparse certification gap was closed by completing exact zero
  singular triplets after Golub-Kahan breakdown, so H is now blocked by
  speed/memory rather than requested-rank accounting on that row. It is
  still much slower and more memory-hungry there, so it is not ready to replace
  the default planner path. An unexported
  `reference_block_golub_kahan_thick_restart_svd()` oracle now defines the next
  block/thick-restart SVD contract with restart, locking, clustered-subspace,
  and exact-zero-singular-triplet tests. This is deliberately a reference
  artifact only; it does not change planner labels or H completion status. The
  first native-facing block-GK kernel,
  `eigencore_block_golub_kahan_ritz`, now computes selected singular Ritz
  slices from a right basis and cached `A V`, including active-column windows
  and target taxonomy checks. A matching internal dense/CSC native basis-cycle
  staging path now builds block right bases and cached `A V` with native
  BLAS-3 reorthogonalization and certifies through the native Ritz kernel on
  full-subspace tests. The SVD surface can now expose that staging path as an
  internal `eigencore_block_golub_kahan_cycle` row under `--h-candidate`, while
  keeping it out of external-reference gates. A quick installed-package probe
  showed the correct shape but not release readiness: it certified
  `rank_deficient_sparse` and `clustered_dense` and was faster than the scalar
  Golub-Kahan rows, but it was still slower/more memory-heavy than the best
  certified external reference, and it failed the `wide_sparse` certificate
  despite reducing wall time. It is scaffolding for the future thick-restart
  loop, not a production solver. The staging path now has an adaptive
  subspace-growth mode that records each attempted subspace and total
  work. This closes the immediate "wide sparse is fast but uncertified"
  diagnostic hole: on the quick `wide_sparse` row it certifies after growing
  through three basis attempts. The result is intentionally still not
  promotable because the from-scratch retry cost is worse than the scalar
  projected candidate and the best certified references. The next H step must
  reuse retained vectors through a real thick-restart loop; adaptive rebuilding
  is only a convergence oracle and benchmark diagnostic. A Ritz-seeded variant
  now carries previous Ritz vectors into later adaptive attempts, cutting the
  same quick wide-sparse diagnostic from `167` to `73` native matvecs while
  preserving certification. It still does not improve the release surface
  because R-visible attempt materialization and missing retained native
  workspace dominate wall time and allocation. This narrows the production
  requirement: implement retained-subspace restart in native code, not more
  from-scratch adaptive cycles.
- shift-invert via `shift_invert(sigma)` is now wired through a reference
  Hermitian Lanczos path on the inverted operator, with three honest planner
  labels: `reference Hermitian Lanczos shift-invert (dense LU)` for dense
  double sources, `reference Hermitian Lanczos shift-invert (sparse LU)` for
  `dgCMatrix`/`dsCMatrix` sources, and `reference Hermitian Lanczos
  shift-invert (user solve)` for matrix-free `A` with a user-supplied solve
  function. Eigenvalues are returned in original coordinates with full
  operator-side certification; near-singular shifted operators and zero-
  magnitude inverted eigenvalues are rejected with actionable errors;
  generalized SPD shift-invert is rejected at plan time with a roadmap note.
  Milestone L's "smallest/interior works on at least one real Laplacian"
  exit criterion is satisfied by `tests/testthat/test-shift-invert.R` (1D
  Laplacian, smallest 4 eigenvalues). The native shift-invert hot loop and
  generalized SPD shift-invert remain open for V1 promotion.
- a property-based test grid now runs under `NOT_CRAN=true` covering Hermitian
  and SVD certificate residual contracts across 5 spectrum patterns
  (`uniform`, `clustered`, `exponential`, `geometric`, `two_cluster`),
  multiple seeds, and multiple sizes; total ~230 assertions. The grid is
  skipped on CRAN to keep CI under budget but runs on every local
  `testthat::test_dir()` and on the GH Actions check matrix.
- a SVD adversarial-and-honest-label bank now asserts that planner labels
  match the kernel that actually runs across dense `auto`, dense
  `golub_kahan()`, sparse CSC `auto`, and matrix-free reference paths, and
  exercises `vectors = "left"|"right"|"none"` modes (which previously leaked
  `result$values` into `result$v` via `$` partial-matching — fixed in
  `R/solve.R`). Sparse CSC `auto()` can now take a certified Gram SVD special
  case for small rectangular problems, with residuals recomputed in original
  coordinates. That path handles exact rank deficiency through deterministic
  zero-singular-vector completion and carries a certification-gated native
  Golub-Kahan fallback when the Gram certificate is weaker. The planner and
  benchmark rows expose this policy through `fallback_policy`,
  `runtime_fallback`, `fallback_attempted`, `fallback_used`,
  `gram_max_backward_error`, and `fallback_max_backward_error`.

Do not mark the Hermitian solver milestone complete until the block path and
speed/memory gate below pass on reproducible benchmark scripts.

### Current PRD alignment and attack surfaces

The project is now past the original skeleton/prototype phase. The architecture
and correctness machinery are credible: operator/problem/plan/result/certificate
boundaries exist, planner labels are mostly honest, sparse non-densification is
guarded, native dense/CSC/diagonal kernels exist, and the test bank catches
certificate and planner regressions. The codebase is still not V1-disruptive
because the PRD's performance promise is not yet defensible.

Working status against the sequenced milestones:

| Milestone | Status | Notes |
|---|---|---|
| A | mostly done | Public/reference/native layering is established; reference fallbacks must remain honestly labeled. |
| B | mostly done | Dense, CSC, diagonal, adjoint, and explicit built-in scaling preserve native block-apply provenance; composed/centered native fusion remains incomplete. |
| C | substantially done | Adversarial, property, oracle, and benchmark-smoke tests are in place and expanding. |
| D | mostly done | Typed certificates and target taxonomy are implemented enough for current paths. |
| E | mostly done | Dense fallback is memory-budgeted; sparse densification is rejected in production paths. |
| F | largely done for current paths | Native ortho and certificate kernels exist, but not every future solver path is fully native-certificate-backed. |
| G0 | done | Native scalar Hermitian staging path exists and certifies on dense/CSC cases. |
| G1 | done | Promoted native block Hermitian Lanczos runs by default in benchmark-proven regimes; strict Hermitian sparse and dense regression gates pass against certified RSpectra/PRIMME references. |
| H | staged, not complete | Native Golub-Kahan exists as a staging path; block-GK restart comparators now include cached Ritz-vector `A V` paths, compact native fit extraction, restart-efficiency diagnostics, and a first native retained-restart candidate that constructs Ritz-plus-random restarts inside C. Production thick-restart SVD and SVD performance gates remain open. |
| I | prototype | Randomized SVD has reference implementation, normalizers, and certified refinement; native approximate engine remains open. |
| J | partial | Native generalized SPD LOBPCG slices exist for built-in `B`, explicitly SPD matrix-free `B`, dense constraints, and typed shifted-diagonal / shifted-tridiagonal preconditioners; strict benchmark rows now gate bare, shifted-diagonal, shifted-tridiagonal sparse-smallest, constrained, and adversarial B native-contract diagnostics; the adversarial B bank covers largest/smallest ill-conditioned diagonal, sparse CSC, and explicitly SPD matrix-free B without dense fallback; broader generalized preconditioning and promotion remain open. |
| K | not complete | B-orthogonal block Lanczos is still a later generalized-SPD refinement path. |
| L | reference-complete, native-open | Shift-invert works through honest reference paths; native hot loop/factorization-cache production path remains open. |
| M | partial | Standard Hermitian LOBPCG and tridiagonal preconditioner staging are useful but not final release surfaces. |
| N | not started as release hardening | CRAN/sanitizer/valgrind, docs, migration guide, and release benchmark reports remain ahead. |

PRD truth check:

- **Trustworthy/certified results:** strong and improving. Certificates are now
  the central design surface, not a debug add-on.
- **Fastest path:** not yet true across the whole PRD. G1 Hermitian promotion
  is the first benchmark-backed win, but RSpectra and/or `irlba` still beat
  eigencore on important SVD and unpromoted regimes, even when eigencore
  certifies more tightly.
- **Block-native engine:** partially true. Built-in block operators and block
  Hermitian solver loops exist for promoted regimes, but V1 still needs the
  same native block discipline across SVD, generalized SPD, and operator-fusion
  paths.
- **Default SVD not normal equations:** directionally true, but the final V1
  answer must be native thick-restarted Golub-Kahan. The Gram path is an
  explicit, inspectable, bounded special case, not the default SVD doctrine.
- **Generalized SPD first-class:** partially true at the API/certificate level;
  not yet complete at native production-solver level.
- **Inspectable automation:** mostly true; new solver work must continue to
  update `plan_solver()`, dispatch, result diagnostics, and benchmarks together.

Primary attack surfaces, in order:

1. **H production SVD.** Replace staging/Gram-special-case dependence with a
   native thick-restarted Golub-Kahan path that wins sparse SVD time-to-certified
   answer against RSpectra/`irlba` on the PRD benchmark subset. Current block-GK
   comparators show retained native restart workspace is the next real
   algorithmic surface; thinner R-level restarts trade speed for memory but do
   not close the gate. The H benchmark rows now expose attempted subspaces,
   restart start width, warm-start count, certified attempt, final-attempt work,
   and total orthogonalization work; future retained-restart patches must move
   those fields in the right direction, not just lower wall-clock noise. The
   cached Ritz-start comparators reuse exact `A V` for retained vectors; the
   random-tail variant confirms that prefix caching alone is not enough while
   restarts are still rebuilt from R, so the bridge remains a proper native
   retained-restart workspace. Compact native block-GK fit extraction removes
   full-basis return materialization from the staging path and should remain the
   default benchmark surface for H diagnostics. Its stage timings show native
   basis iteration remains the dominant cost, with Ritz extraction secondary;
   optimize retained restart/workspace reuse before spending effort on small
   projected SVD tuning. Minor compact-fit workspace cleanup has removed
   avoidable double initialization, but the structural H gap is still retained
   restart/workspace reuse. A residual-tail restart diagnostic certifies but
   costs more native apply calls than Ritz-plus-random on the wide sparse H
   probe, so residual-tail restarts should not be promoted without new evidence.
   Cached-`A v` native certificates remove a redundant left residual apply from
   compact block-GK internal certification, but only modestly reduce allocation;
   the main H path still runs through retained native restart/workspace design.
   The first retained native restart candidate is exposed only as
   `eigencore_block_golub_kahan_retained` on the H benchmark surface. It builds
   Ritz-plus-random restart blocks inside C, returns compact selected triplets,
   and on the wide sparse probe certifies with the same `73` apply calls while
   reducing R-visible allocation versus the R adaptive cycle. It is a staging
   row, not H closure: it does not yet do per-attempt native certification,
   locking, or a full production thick-restart policy.
2. **J generalized SPD LOBPCG promotion.** Broaden generalized
   preconditioning beyond the typed shifted-diagonal and certified
   shifted-tridiagonal sparse-smallest case, keep the benchmark
   B-orthogonality/native-path diagnostics green, keep the adversarial B
   benchmark contract green, and promote the native path only after sparse
   no-densification gates and the remaining production checks pass.
3. **L native shift-invert.** Move from reference inverted-operator Lanczos to
   factorization-aware native transforms with cached solves and original
   problem residual certification.
4. **Operator fusion.** Native centered/scaled/composed operators are needed
   for matrix-free PCA/SVD and for keeping the mathematical API elegant without
   paying R callback overhead.
5. **Release hardening.** Keep benchmark claims tied to reproducible
   `bench::mark()` scripts; do not weaken the PRD's "faster and unambiguously
   better" bar to fit current results.

### G1 implementation clarifications

The G1 solver is the first production Hermitian block path, so promotion is not
just a label change. The implementation must satisfy these extra constraints:

- Planner promotion must update both `plan_solver()` and solver dispatch. The
  plan must record the chosen block size, restart budget, and subspace budget
  well enough that `solve.eigencore_eigen_problem()` runs the same path the
  planner reported. `auto()` must not silently fall back to scalar Lanczos after
  planning a block path. Explicit `lanczos(block = 1)` remains the scalar
  staging/debug path.
- The current R-facing native MGS2, CholQR2, Rayleigh-Ritz, and dense symmetric
  eigen entry points are `.Call` wrappers and may allocate R objects. G1 must
  first factor the allocation-free C++ kernels needed by the native block loop,
  with all solver temporaries supplied by a setup-time workspace.
- Native locking and convergence tests must use the same backward-error scale
  definition as certificates: the operator norm source plus the identity `B`
  scale used by `eigen_backward_scale()`. A shortcut such as
  `tol * max(abs(theta), 1)` is acceptable only as a development diagnostic, not
  as the production lock criterion.
- The former `native block Hermitian Lanczos thick-restart candidate` counted
  as promoted G1 only after the adversarial bank and RSpectra/PRIMME benchmark
  gates passed with the default `eigencore` path.
- Correctness comparisons for iterative paths use residuals scaled by the
  certificate denominator, e.g. `max(tol * scale, 100 * eps * scale)`, and use
  subspace-distance assertions for clustered or repeated eigenvalues rather than
  per-vector identity.
- The G1 benchmark gate is median time-to-certified-answer at least `1.25x`
  faster than RSpectra on the Hermitian suite, with the memory gate still
  applied against certified references. PRIMME remains a required V1 reference
  point: if PRIMME is installed, G1 must at least pass parity against PRIMME or
  leave an explicit G1 performance gap recorded in this plan and
  `benchmarks/RELEASES.md`.
- Hermitian benchmark helpers must expose the reference set being gated and
  enforce strict failure behavior consistently across
  `bench-hermitian-sparse.R` and `bench-native-hermitian-gate.R`.

### G1 execution record

G1 was a promotion milestone, not just an implementation milestone. The
candidate became the default production Hermitian path only after it won the
correctness and benchmark gates.

**G1.0 Freeze and measure the current state**

- Keep the explicit block candidate benchmark subject (`eigencore_block`) on
  the same bounded controls that would be used after promotion: block size,
  `max_subspace`, restart budget, target, and tolerance. The helper must never
  force `max_subspace = n` for the promotion gate.
- Save a source-loaded baseline under
  `inst/benchmarks/baselines/g1_candidate_pre.csv` covering at least:
  path Laplacian smallest eigenpairs, dense Hermitian regression matrices,
  clustered spectra, and ill-conditioned diagonal spectra.
- Regenerate that baseline with
  `inst/benchmarks/bench-g1-candidate-baseline.R --save`. Its `--strict` mode
  checks case coverage, method-error-free rows, and full certification for the
  current block-candidate baseline.
- Record per-case `method`, `nconv`, `certificate_passed`, `iterations`,
  `matvecs`, `restarts`, `ortho_passes`, `locking_events`, `block_size`,
  memory, and median time-to-certified-answer. A row that does not certify is
  a failing row regardless of speed.

Exit: baseline committed, quick and full benchmark scripts can select scalar,
explicit block, RSpectra, and PRIMME as separate methods without ambiguity.

**G1.1 Make the candidate mathematically contract-complete**

- Keep `R/reference_block_lanczos.R` as the permanent unexported oracle, and
  compare native results against it on dense and CSC fixtures using
  certificate-scaled residual tolerances and subspace-distance assertions for
  clustered or repeated eigenvalues.
- Add white-box restart/locking tests: forced restart with
  `max_subspace < n`, forced soft locking on a well-separated cluster,
  duplicate Ritz-value rejection, lock-vector orthogonality against the active
  block, and no sparse densification.
- Treat native/native-oracle mismatches as regressions if either path returns
  `passed = TRUE` but the subspace distance or backward-error contract fails.
- Verify target semantics for `largest`, `smallest`, `largest_magnitude`, and
  `smallest_magnitude`; keep `nearest()` routed away from native Lanczos until
  shift-invert lands.

Exit: `tests/testthat/test-block-lanczos-oracle.R`,
`tests/testthat/test-block-lanczos-thick-restart.R`, and the adversarial bank
are green with the native block path enabled.

**G1.2 Convert the native candidate into a real hot-loop kernel**

- Move all solver temporaries into a setup-time workspace struct. The inner
  restart/extension loop must not call `R_alloc`, `Rf_*`, or allocate R
  objects; `.Call` allocation is allowed only for final result packaging.
- Preserve BLAS-3 operator application (`A * V_block`) for dense and CSC
  built-ins. Sparse inputs must stay on CSC block apply and never materialize a
  dense copy.
- Replace the current full-reorthogonalized block Krylov extension with the
  Hermitian block recurrence/projection structure needed for performance:
  block three-term recurrence before restart, arrow-plus-block-tridiagonal
  after thick restart, and dense `dsyev`/`dsyevd` only on the small projected
  problem.
- Keep the certificate-scale lock criterion:
  `||r_i|| <= tol * eigen_backward_scale(norm_A, 1, theta_i, v_i)`.
- Surface native counters for iteration count, block applies, reorth passes,
  restarts, locking events, final active subspace, locked count, and restart
  convergence history.

Exit: native block path has allocation-free hot-loop discipline, stable
counters, and no correctness loss against G1.1 tests.

**G1.3 Tune only after correctness is fixed**

Tune in this order, one isolated commit per tuning knob:

1. Planner block size: compare `b = 2`, `ceil(k / 4)`, `ceil(k / 2)`, and
   cap `8`, separately for `k < 8`, `8 <= k <= 32`, and larger `k`.
2. Subspace budget: compare `3*k + 20`, `4*k + 20`, and `max(3*k + 20,
   6*b + 20)` on sparse Laplacian and dense clustered cases.
3. Reorthogonalization: keep two-pass MGS2 as the correctness baseline, then
   test selective second pass when measured orthogonality loss stays within
   certificate tolerance.
4. Projection/RR cadence: avoid recomputing more projected state than needed
   between restart decisions.
5. Dense apply microkernels: use symmetric structure where it measurably helps,
   but only after it preserves certificate results on dense Hermitian fixtures.

Exit: the block path is faster than scalar on the Hermitian suite and within
striking distance of the RSpectra/PRIMME gate before any planner promotion.

**G1.4 Promote planner and dispatch**

- Promote `auto()` to the block path only in regimes that pass the benchmark
  gate. Initial policy should be conservative: `k >= 4`, native dense/CSC
  source, supported target, and a size threshold proven by the gate. Keep
  scalar thick restart for small `k`, small `n`, unsupported block regimes, and
  explicit `lanczos(block = 1)`.
- Update the planner label to the production string only after G1.1-G1.3 pass:
  `native block Hermitian Lanczos (thick restart, locking)`.
- Ensure `solve.eigencore_eigen_problem()` consumes `plan$controls` exactly:
  block size, `max_subspace`, and restart budget must match what the planner
  printed.
- Keep an explicit development escape hatch for the scalar native path, but do
  not let `auto()` report block and execute scalar.

Exit: `result$method`, `result$plan$method`, and actual native entry point agree
for promoted block runs; scalar remains honestly labeled.

**G1.5 Enforce the release gate**

- `inst/benchmarks/bench-native-hermitian-gate.R --strict` and
  `inst/benchmarks/bench-hermitian-sparse.R --strict` must fail on any
  uncertified eigencore row, `nconv < k`, missing reference row, RSpectra speed
  ratio below `1.25`, PRIMME parity below `1.0` when PRIMME is installed, or
  memory ratio below `1.0`.
- Gate only certified reference rows. If RSpectra or PRIMME fails to certify a
  case, record that in the row and exclude it from speed comparison rather than
  silently treating it as a win.
- Save passing rows and gate summaries to `benchmarks/RELEASES.md`, including
  package version, seed, platform, BLAS/LAPACK note when available, and exact
  benchmark command.

Exit: Hermitian sparse and dense regression gates pass reproducibly; G1 can be
marked done.

**G1.6 Documentation cleanup after promotion**

- Update `prd.json` only after G1.5 passes: move the promoted block path into
  `implemented_path` and remove the promoted items from `not_yet_satisfied`.
- Update this plan's status checkpoint and milestone table. If any PRIMME
  parity gap remains, G1 is not complete; record it as an explicit blocker
  rather than weakening the gate.
- Add a short release note explaining when `auto()` chooses block vs scalar and
  how users can request `lanczos(block = 1)` for scalar debugging.

## 4. Dense-fallback policy

Dense fallback is a correctness tool, not a performance path. Policy:

- `allow_dense_fallback` defaults to `"auto"` with a **memory** budget
  (not a dimension budget):
  `projected_peak_bytes <= getOption("eigencore.dense_fallback_mb", 256) * 1e6`.
- `"never"` forbids it. `"always"` permits it with a warning.
- When dense fallback runs, the certificate records `method = "dense oracle"`
  and the plan records the rejected native path + reason.
- **Sparse operators never silently densify.** If a sparse input would
  require a dense copy above budget, the planner errors with an estimated
  byte size and suggests a method change.

## 5. Certificate typing

Every certificate carries explicit provenance:

```r
certificate <- list(
  passed              = logical(1),
  tolerance           = numeric(1),
  certificate_type    = c("residual_backward_error", "uncomputed", ...),
  norm_bound_type     = c("frobenius_exact", "frobenius_metadata",
                          "identity_exact", "frobenius_hutchinson_estimate",
                          ...),
  scale_is_estimate   = logical(1),
  max_residual        = numeric(1),
  max_backward_error  = numeric(1),
  max_orthogonality_loss = numeric(1),
  failed_indices      = integer(),
  scale               = numeric(1),
  notes               = character()
)
```

Rules:

- `certificate_type = "residual_backward_error"` means the certificate is based
  on explicit residual recomputation or residuals returned by a native solver
  and scaled by the shared backward-error denominator.
- Hutchinson-style stochastic norm estimates produce
  `norm_bound_type = "frobenius_hutchinson_estimate"` and
  `scale_is_estimate = TRUE`. In that case `passed` is withheld even if every
  returned residual is below tolerance.
- Backward-error denominators currently use exact Frobenius norms for explicit
  dense matrices, exact Frobenius metadata for built-in sparse/diagonal
  operators, identity metadata for standard eigenproblem B = I, and labeled
  Hutchinson estimates only for matrix-free operators without norm metadata.
- Dense and built-in operator paths share one scale definition; only the norm
  source differs.

## 6. Target taxonomy

Replace the current 4-member target set with an explicit taxonomy:

```r
largest_algebraic()      # = current largest()
smallest_algebraic()     # = current smallest()
largest_magnitude()      # = current largest_magnitude()
smallest_magnitude()     # NEW — correct mapping for ARPACK "SM"
largest_real()           # NEW
smallest_real()          # NEW
largest_imaginary()      # NEW
smallest_imaginary()     # NEW
both_ends(k_low, k_high) # NEW — for interval-ish workloads
nearest(sigma)           # unchanged
```

Current `largest()` / `smallest()` become thin aliases for backward
compatibility, with a deprecation note when called inside
`eigs_sym()`/`eigs()` — those shims must map ARPACK `which` codes exactly:

| ARPACK | maps to                |
|--------|------------------------|
| LA     | largest_algebraic()    |
| SA     | smallest_algebraic()   |
| LM     | largest_magnitude()    |
| SM     | smallest_magnitude()   |
| LR     | largest_real()         |
| SR     | smallest_real()        |
| LI     | largest_imaginary()    |
| SI     | smallest_imaginary()   |
| BE     | both_ends()            |

The current `SM → smallest()` is a compatibility bug and is fixed here.

## 7. Test discipline — adversarial bank before engine swap

Build the correctness backstop *before* the C++ solvers land, so each native
implementation slots behind a fixed contract. Required test families:

- **Oracle tests**: diagonal known spectra, orthogonally generated SVD with
  prescribed σ, generalized SPD synthetic pair with known spectrum.
- **Adversarial tests**:
  - clustered eigenvalues (gap ≤ 1e-8)
  - nearly repeated σ
  - rank-deficient rectangular
  - graph Laplacian with known nullspace
  - ill-conditioned B
  - non-normal nonsymmetric
  - poorly scaled rows / columns
  - shift-invert near singular σ
- **Property tests**:
  - `check_adjoint()` randomized across operator types
  - residual recomputation equals reported residual (operator ≡ dense)
  - certificate agreement between dense and operator paths
  - subspace-angle convergence for clustered cases (not per-vector identity)
- **Regression benchmark smoke** (fast): one case per solver family,
  enforces a wall-clock ceiling using `bench::mark()` with warmup, not
  `system.time()`.

Each of these runs against the reference R solvers today. When a native
solver lands, the same bank runs against it unchanged.

## 8. Public API surface contract

Every exported name is a compatibility commitment. Current NAMESPACE exports
~50 symbols, several of which back stubs or internal helpers.

Export discipline:

- **Keep exported**: `linear_operator`, `as_operator`, `adjoint`, `center`,
  `scale_rows`, `scale_cols`, `compose`, `crossprod_operator`,
  `symmetric_operator`, `check_adjoint`, `eigen_problem`, `svd_problem`,
  `plan_solver`, `solve` (S3), `eig_partial`, `svd_partial`, `eigs`,
  `eigs_sym`, `svds`, `certificate`, `diagnostics`, `values`, `vectors`,
  `left_vectors`, `right_vectors`, `residuals`, `backward_error`, all target
  constructors (per §6), `auto`, `lanczos`, `golub_kahan`, `randomized`,
  `shift_invert`, `hermitian`, `general`, `euclidean`.
- **Un-export until implementation lands** (`@keywords internal`):
  `mgs2`, `cholqr2`, `b_orthogonalize`, `rayleigh_ritz`,
  `orthogonality_loss`, `operator_scale`, `operator_sum`,
  `validate_eigen_accuracy`, `validate_svd_accuracy`,
  `benchmark_eigen_methods`, `benchmark_svd_methods`.
- **Add when implemented**: `plot_convergence`, `reorthogonalize`,
  `preconditioner`, `factorization_cache`.

Benchmark and validation helpers become non-exported utilities under a
`tools/` or `inst/internal/` harness.

## 9. V1 acceptance gates (sharpened)

In addition to the PRD's criteria, V1 must pass:

- **No R-level iteration loops** in any non-reference solver path.
- **No silent dense fallback** above the memory budget (§4).
- **No default SVD through `A^T A`** for any operator class.
- **No `passed = TRUE`** from a stochastic norm estimate. Such cases must carry
  `scale_is_estimate = TRUE` and an explanatory note until a refinement path
  can produce a non-stochastic scale.
- **Sorted singular values** by target, always.
- **Planner honesty** — every result's `method` string equals the path
  that actually ran.
- **`eigs_sym()` `SM` maps to smallest-magnitude.**
- **Result object** always has `nconv`, `residuals`, `method`, `target`,
  `iterations`, `matvecs`, `certificate`, `warnings`.

Performance gates remain as in `prd.json` but measured with `bench::mark()`
over a fixed seed, not elapsed wall clock, and measured as **time-to-certified-answer**.
For G1 specifically, RSpectra carries the PRD's `1.25x` release threshold, while
PRIMME is a required parity reference when available unless the remaining gap is
explicitly documented as not yet G1-complete.

## 10. Sequenced milestones (replacing PRD Milestones 1–9 in order of work)

| # | Deliverable                                                   | Exit                                                       |
|---|---------------------------------------------------------------|------------------------------------------------------------|
| A | Layer split + rename R solvers to `reference_*`               | No public path reaches reference code without a label      |
| B | Frozen `BlockOperator` ABI + dense / CSC / diagonal natives   | ABI documented; ops pass `check_adjoint()` property tests  |
| C | Adversarial test bank (against reference solvers)             | Bank green on reference; named cases cover §7 list         |
| D | Certificate typing + target taxonomy + `SM` shim fix          | Dual-path cert agreement; ARPACK-compatible `which` codes  |
| E | Dense fallback policy + memory budget                         | Sparse never silently densifies; planner rejects w/ reason |
| F | Native MGS2 / CholQR2 + native certificate kernels            | Reference bank green against native ortho + cert           |
| G0 | Native scalar Hermitian thick-restart staging path          | Honest planner label; certified dense/CSC path; sparse non-densification; unsupported targets routed honestly |
| G1 | Native Hermitian block Lanczos (thick restart, locking)     | Reference bank green; beats RSpectra by PRD threshold; PRIMME parity when available; Hermitian memory gate green |
| H | Native Golub–Kahan (thick restart, true SVD residuals)        | Reference bank green; beats RSpectra on SVD subset         |
| I | Native randomized SVD + refinement + honest estimate cert     | HegelSVD-derived regimes ported; matches or improves `rsvd::rsvd()` accuracy; beats `rsvd` on time-to-certified-answer in randomized-planner regimes; 2× deterministic SVD on selected approximate cases; honest degraded certificates |
| J | Generalized SPD LOBPCG (dense/sparse/matrix-free)             | Largest/smallest generalized SPD paths certify without sparse densification; benchmark rows prove native B-orthogonal execution; adversarial B bank green |
| K | Generalized SPD B-orthogonal Lanczos alternate/refinement     | Shift-invert-free generalized path passes adversarial B and agrees with LOBPCG certificates |
| L | Shift-invert with user solve / factorization cache            | Smallest/interior works on at least one real Laplacian     |
| M | Standard symmetric LOBPCG refinements                         | Graph Laplacian test converges with fewer iterations      |
| N | Release hardening: vignettes, migration guide, CRAN/sanitizer | All PRD acceptance criteria + §9 gates green               |

A–E are *prerequisites* to any native solver work. F–M are the engine; G0 is
only the scalar Hermitian staging checkpoint, and G1 is the full Hermitian
block/performance milestone. N is the release.

## 11. Things to preserve exactly

- The problem / target / method / plan / result / certificate decomposition.
- `linear_operator()` with block apply + adjoint + structure.
- Certificate-first philosophy (already the strongest part of the codebase).
- RSpectra shims as the migration surface.
- The existing honest planner labels — they are the template for §9 gate.

## 12. Things explicitly out of scope for V1

Already listed in `prd.json` `v1_product_scope.non_goals`. Unchanged. Reiterated:
no GPU, no distributed memory, no DelayedArray / HDF5 integration, no PRIMME /
SLEPc wrapping, no full nonsymmetric generalized eig, no contour integral
methods, no nonlinear eigenproblems.

## 13. How this plan is enforced

- Every PR that adds a solver path also updates `plan_solver()` labels, solver
  dispatch, the certificate type fields, and the adversarial bank. No
  exceptions.
- Any export added without a corresponding native-or-reference
  implementation is a review-blocker.
- Benchmark claims in README / vignettes cite a reproducible `bench::mark()`
  script under `inst/benchmarks/` with pinned seed and package versions.

---

*This plan refines `prd.json` and `vision.md` but does not override them. If
a conflict arises, `vision.md` wins on goals, `prd.json` wins on scope, and
this plan wins on engineering method.*
