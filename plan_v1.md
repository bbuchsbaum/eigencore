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
   without inheriting its dense-only assumptions.
6. **Generalized SPD LOBPCG** for dense, sparse, and matrix-free operators.
   This is the primary scalable V1 path for `A x = lambda B x`, especially
   smallest eigenpairs and preconditioned problems. It uses block iteration,
   B-orthogonal Rayleigh-Ritz extraction, optional constraints/deflation, and
   native residual/certificate kernels. The planner maps largest generalized
   eigenpairs through the sign-flipped problem when appropriate; it never
   densifies sparse `A` or `B` silently.
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
- a `native block Hermitian Lanczos prototype` path for dense double matrices
  and `dgCMatrix` operators, with block native operator application and honest
  no-restart/no-locking metadata;
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

The Hermitian path supports native targets `largest`, `smallest`,
`largest_magnitude`, and `smallest_magnitude`; unsupported targets such as
`nearest()` route away from the native path with an honest planner label.

This is not Milestone G completion. The implementation is scalar (`block = 1`),
not block-native, and the quick benchmark gates currently show certification
but not performance leadership:

- sparse Hermitian path Laplacian `n = 200, k = 5`: eigencore certifies all
  requested pairs with the scalar-native path and the tuned default
  `max_subspace = 3*k + 20`, but is still slower than RSpectra and
  allocates more memory on the quick gate;
- the block Hermitian prototype certifies the same quick Laplacian case only
  when allowed to use the full `n`-dimensional subspace; it gives excellent
  residuals but is slower than the scalar path and therefore remains a
  prototype, not the G1 solution;
- preconditioned reference LOBPCG on the same quick Laplacian case certifies
  in about 10 iterations. The current gate shows it is about 1.8x faster than
  scalar eigencore, slightly faster than PRIMME on some runs, still slower than
  RSpectra, and allocates materially more memory because the path is R-level
  and factorization-backed;
- a native shifted tridiagonal preconditioner now covers path-Laplacian-style
  tridiagonal gates. It preserves the 10-iteration convergence behavior and
  reduces memory versus the Cholesky-backed reference setup, but the current
  R-level LOBPCG loop still dominates enough allocation that the memory gate
  remains open. Benchmark rows and the LOBPCG gate now record which typed
  preconditioner candidate was selected;
- sparse tall-skinny SVD `500 x 80, rank = 5`: eigencore certifies but is
  still materially slower than RSpectra and `irlba`.

Do not mark the Hermitian solver milestone complete until the block path and
speed/memory gate below pass on reproducible benchmark scripts.

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
| G1 | Native Hermitian block Lanczos (thick restart, locking)     | Reference bank green; beats RSpectra/PRIMME on Hermitian subset by time-to-certified-answer and memory gate |
| H | Native Golub–Kahan (thick restart, true SVD residuals)        | Reference bank green; beats RSpectra on SVD subset         |
| I | Native randomized SVD + refinement + honest estimate cert     | HegelSVD-derived regimes ported; 2× deterministic SVD on selected approximate cases; honest degraded certificates |
| J | Generalized SPD LOBPCG (dense/sparse/matrix-free)             | Largest/smallest generalized SPD paths certify without sparse densification; adversarial B bank green |
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

- Every PR that adds a solver path also updates `plan_solver()` labels, the
  certificate type fields, and the adversarial bank. No exceptions.
- Any export added without a corresponding native-or-reference
  implementation is a review-blocker.
- Benchmark claims in README / vignettes cite a reproducible `bench::mark()`
  script under `inst/benchmarks/` with pinned seed and package versions.

---

*This plan refines `prd.json` and `vision.md` but does not override them. If
a conflict arises, `vision.md` wins on goals, `prd.json` wins on scope, and
this plan wins on engineering method.*
