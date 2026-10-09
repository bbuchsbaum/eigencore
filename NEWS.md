# eigencore (development version)

## Behaviour change: two-norm backward errors

* Certificate backward errors now use the standard normwise 2-norm
  definition, `||A x - lambda B x|| / ((||A||_2 + |lambda| ||B||_2) ||x||)`
  for eigenpairs and `sqrt(||A v - s u||^2 + ||A^H u - s v||^2) / ||A||_2`
  for singular triplets, instead of dividing by Frobenius norms (up to
  `sqrt(rank)` times larger, so certificates were correspondingly more
  lenient). Reported backward errors are larger than before, and results
  whose residuals only met the Frobenius-scaled tolerance now fail.
* The denominator is exact where cheap (diagonal matrices, a full computed
  spectrum, `metadata$two_norm`) and otherwise a lower bound on `||A||_2`:
  the largest column norm, `||A x|| / ||x||` of the certified vectors, or,
  only when it could flip a failing pair, a short deterministic Lanczos
  estimate memoised per operator (it does not use the random-number stream).
  A lower bound over-states the backward error, so `passed` stays sound.
* `norm_bound_type` values are now `"two_norm_exact"`,
  `"two_norm_lower_bound"` and `"identity_exact"` (replacing
  `frobenius_exact`, `frobenius_metadata` and
  `frobenius_hutchinson_estimate`); new fields `norm_source` and
  `norm_values` record where the scale came from. The stochastic Hutchinson
  norm estimate is gone, so matrix-free operators without norm metadata now
  certify (`scale_is_estimate` is always `FALSE` for built-in
  certificates). The planner accordingly routes matrix-free smallest and
  nearest (interior) SVD targets to the native matrix-free Golub-Kahan
  routes without requiring `metadata$frobenius_norm` (the plan control
  `requires_nonestimated_norm_scale` is now `FALSE`); before, a callback
  operator without metadata got the largest-target route or an
  unsupported-interior error. Native solvers use the same lower bounds for their internal
  convergence scales; the implicit Gram SVD kernel runs at `tol / 2` and,
  if the exact certificate still fails, retries once at a tighter tolerance.
* Certificates for complex eigenvectors of real sparse or matrix-free
  operators apply the operator to real and imaginary parts instead of
  densifying the source.

## Correctness fixes

* `center(rows = TRUE, columns = TRUE)` now double centers correctly. Row
  means were taken from the uncentered matrix, so the grand mean was
  subtracted twice on the dense, callback, and native CSC paths.
* Column-centered `dgCMatrix` operators carry their exact Frobenius norm, so
  their SVD certificates can pass instead of always reporting an estimate.
* Complex Frobenius norms no longer drop imaginary parts, and sparse sources
  are no longer densified to compute a certificate norm.
* Complex Hermitian eigenproblems and Hermitian-definite pencils use the
  `zheev`-based kernels, so repeated eigenvalues keep an orthonormal
  (B-orthonormal) eigenbasis.
* `eig_full(A, B)` with a symmetric but indefinite or singular `B` falls back
  to QZ instead of failing in `dpotrf` (unless `structure = hermitian()` is
  requested explicitly).
* The Gram SVD zero threshold is relative to the matrix scale instead of
  `max(1, d)`, so tiny-norm matrices no longer return zero singular values.
* NA, NaN, and Inf matrix inputs are rejected when the operator is built.
* `k`/`rank` are validated once (whole number in `1..n`), eigenproblems
  require a square operator, and `both_ends(k_low, k_high)` must match `k`.
* `seed =` in `eig_partial()`/`svd_partial()` restores the global random
  stream on exit.
* `shifted_tridiagonal_preconditioner()` accepts symmetric storage and reads
  the bands without an R-level loop.
* Native kernels use 64-bit offsets for basis and certificate indexing, so
  problems with more than 2^31 basis or vector entries no longer overflow.
* Block Lanczos and block Golub-Kahan always run two Cholesky-QR passes
  (previously only for n < 64), keeping new blocks orthonormal when the
  residual block is ill-conditioned.
* Scalar Lanczos and Golub-Kahan stop with a clear error on non-finite
  values; NaN no longer passes the native symmetry and positive-diagonal
  checks, and is no longer hidden in maximum backward-error summaries.
* The native Gram SVD start vector is no longer an exact alternating sign
  pattern (orthogonal to the constant vector), and its zero threshold is
  scale invariant.

## Iteration limits, subspace sizes and nonsymmetric targets

* **Behaviour change:** `maxit` in `eig_partial()`, `solve()` and
  `plan_solver()` is now what its documentation always said, an iteration
  limit. It used to set the Krylov subspace size. It now caps thick-restart
  cycles (scalar, block, generalized and shift-invert callback Lanczos),
  Krylov-Schur restarts (nonsymmetric Arnoldi), restart cycles (reference
  Arnoldi), LOBPCG iterations, and Lanczos steps on unrestarted routes. The
  resolved limit is recorded in `plan$controls$iteration_limit` and
  `plan$controls$iteration_limit_kind`. To set the subspace size, use
  `max_subspace` on the method descriptor: `lanczos()`, `golub_kahan()`, the
  new `auto(max_subspace =)` (honoured by whichever Krylov route the planner
  picks, including shift-invert and the SVD Golub-Kahan / implicit-Gram
  routes), and the new `shift_invert(max_subspace =)`. Code that passed
  `maxit = m` to get an `m`-dimensional subspace should pass
  `method = auto(max_subspace = m)` or `lanczos(max_subspace = m)` instead.
  If `lanczos(max_restarts =)` or `lobpcg(maxit =)` disagrees with `maxit`,
  that is now an error. `lobpcg()` now defaults to `maxit = NULL`, which
  means the solve's `maxit` or the `eigencore.lobpcg_maxit` option (200).
* RSpectra shims: `opts$ncv` now maps to `auto(max_subspace = ncv)` (or
  `lanczos(max_subspace = ncv)` with `initvec`), so the route is the same as
  without it. In `eigs()`/`eigs_sym()`, `opts$maxitr` maps to `maxit` instead
  of being ignored. `svds()` still ignores `maxitr` and warns about it.
* `eigs()` computes left eigenvectors only when asked, as `RSpectra::eigs()`
  does. The new `left = FALSE` argument skips the adjoint Arnoldi solve and
  its certificate, which roughly halves the cost. `eig_partial()` and
  `solve()` take a new `left_vectors = c("auto", "none", "compute")`
  argument. The default `"auto"` keeps the two-sided certificate.
* The default scalar thick-restart Lanczos subspace is now ARPACK-like
  (`max(2k + 1, 20)`, capped at `n`) instead of `3k + 20`.
* Thick-restart Lanczos and implicit-Gram SVD results are sorted by target
  order. Previously values came back in lock order (for example 9, 7, 9).
* Nonsymmetric `smallest_magnitude()` is supported. Dense and sparse inputs
  run through shift-invert Arnoldi at `sigma = 0`. If `A` is singular, the
  shift is perturbed and the change is recorded. Matrix-free operators use
  the Krylov-Schur smallest-magnitude ranking directly.
* Nonsymmetric `nearest(sigma)` and `shift_invert(sigma)` run Krylov-Schur
  Arnoldi on a factorised `A - sigma I`: dense LAPACK QR, sparse LU with
  AMD ordering, or a user `solve`. Eigenvalues are recovered as
  `lambda = sigma + 1/theta`. Right and left residuals are certified on the
  original `A`, and transposed solves reuse the forward factorisation. As a
  result, `eigs(A, k, sigma =)` now works for nonsymmetric `A` with a real
  `sigma`.
* Factorization labels now match the code. Tridiagonal shift-invert and
  metric solves use pivoted LU (`dgttrf`/`dgttrs`), so the labels changed
  from `tridiagonal_thomas_*` to `tridiagonal_lu_*`, and
  `native_sparse_tridiagonal_thomas` became `native_sparse_tridiagonal_lu`.
  Shift-invert routes that run the native thick-restart Lanczos callback are
  now labelled `native thick-restart ... Lanczos shift-invert (... solve
  callback)` instead of `reference ...`.

## Multithreaded sparse kernels

* Native sparse (`dgCMatrix`) products now run on OpenMP threads: `A^T X`
  as a per-column gather, `A X` as a gather over a CSR copy (or per-thread
  row slabs for tall matrices) built once per solve, the centred and
  centred-scaled operators, and the implicit `A^T A` / `A A^T` operators of
  sparse `svds()` / `svd_partial()`. Sparse products are bitwise identical
  for every thread count. The thread count is `getOption("eigencore.threads")`;
  the default is 1 under `R CMD check` or `_R_CHECK_LIMIT_CORES_`, otherwise
  `OMP_NUM_THREADS` or the processor count capped at 8. Builds without
  OpenMP (Apple clang) stay serial. See `?"eigencore-threads"`.
* During a multithreaded sparse solve, a spinning-thread BLAS (OpenBLAS
  pthreads, FlexiBLAS) is switched to one thread and restored afterwards,
  and the Lanczos reorthogonalisation runs on eigencore's own OpenMP
  kernels, so BLAS and OpenMP threads do not compete. Results with one and
  several threads agree to rounding.
* Multi-column sparse products use row-major panels, and the centred-scaled
  sparse operator no longer applies one column at a time. R-level block
  applies no longer duplicate `Y` when `beta = 0` and accept `Y = NULL`.

## RSpectra compatibility

* `eigs()`, `eigs_sym()` and `svds()` follow the RSpectra signatures:
  `sigma` (nearest eigenvalues), function inputs (`n`/`args`,
  `Atrans`/`dim`), `eigs_sym(lower =)` reading one triangle, default
  `which = "LM"` for `eigs_sym()`, decreasing value order from `eigs_sym()`,
  and `opts$tol`, `ncv`, `retvec`, `initvec`, `center`, `scale`. Unused or
  unknown `opts` entries and unknown `which` codes are reported instead of
  being ignored, `nu`/`nv` are honoured, and non-convergence warns.

# eigencore 1.3.0 (2026-08-25)

## Certified positive-semidefinite geometry

* Add a frozen, additive real-double PSD surface: `psd_tolerance()`,
  `psd_policy()`, `psd_identity()`, `psd_factor()`, `psd_gram_factor()`,
  `psd_laplacian()`, `psd_capabilities()`, `psd_spectrum()`, `psd_rank()`,
  `psd_nullity()`, `psd_apply()`, `psd_operator()`, `psd_reduce()`,
  `psd_lift()`, `psd_solve()`, `psd_gram()`, `psd_orthonormalize()`, and
  `psd_reduced_operator()`. Identity, diagonal, dense spectral, and dense Gram
  factors expose complete numerical actions with typed validation,
  provenance, certificates, work, retained-memory, and RDS integrity records.
* Define the singular geometry explicitly: the source is a seminorm in its
  original coordinates, while canonical reduction represents the metric on
  `image(K)` or the quotient by `null(K)`. Roots, inverse roots,
  pseudoinverses, complementary projectors, reduction/lift, block
  orthonormalization, and reduced operators satisfy independently certified
  postconditions. `psd_solve()` is a strict equation solver and rejects RHS
  null components; pseudoinverse application remains total.
* Freeze separate Frobenius-relative tolerances for symmetry, positivity,
  numerical rank, and RHS compatibility. Classification retains original and
  repaired spectra, exact threshold categories, and repair defects. No
  `max(1, scale)` floor is inserted, so zero-absolute-tolerance classification
  remains invariant under finite positive rescaling.

## Structural sparse paths

* Add non-densifying `dgCMatrix` Gram factors that certify form and Gram
  actions from a supplied factor while withholding unsupported roots, rank,
  projectors, reduction, and solves.
* Add `dgCMatrix`/`dsCMatrix` graph-Laplacian validation for symmetry,
  non-positive off-diagonals, zero row sums, optional canonical diagonal
  repair, and connected-component algebraic nullity. The structural route
  exposes form and Gram actions without claiming a complete spectrum.
* Generic sparse symmetric matrices and opaque callbacks remain unsupported
  PSD-factor sources. Capability failures are typed and occur before callback
  probing, dense conversion, or undocumented approximation. The existing
  generalized-eigen `B`/`metric=` API remains SPD-only; singular PSD problems
  use explicit image-space reduction.

## Certification and downstream conformance

* Add adversarial numerical coverage for exact threshold neighbors, scales
  from `1e-12` through `1e12`, repeated eigenspaces, null contamination,
  permutation and block metamorphisms, mutation adequacy, strict-solve
  compatibility, serialization integrity, and fail-closed capability paths.
* Certify 50,000-row sparse Gram and path-Laplacian fixtures without dense
  `n` by `n` state, recording retained bytes, cumulative R allocation,
  isolated-process peak RSS, and reusable-action timing separately.
* Add a provenance-required installed-consumer gate using exported APIs only.
  gprocrustes O/SO and quotient-space checks and DKGE root, projector,
  orthonormalization, and K-Procrustes differential checks pass. A
  pre-existing DKGE tiny-scale threshold divergence remains explicit rather
  than weakening eigencore's relative policy. rfugw adoption is deferred
  because its candidate spectral auxiliaries do not require the PSD factor
  algebra.

# eigencore 1.2.0 (2026-08-25)

## Reusable solver workflows

* `plan_solver()` now returns a schema-versioned executable plan containing
  the problem, original method descriptor, selected route, canonical controls,
  execution arguments, all 20 route-affecting package options, operator
  identity/revision records, serialization capability, and retained-memory
  metadata. `solve(plan)` validates and executes that frozen route without
  invoking the planner. `replan = TRUE` is the explicit opt-in to a fresh
  decision under current policy; execution arguments cannot be overridden.
* `linear_operator()` gains `operator_id`, `revision`, and `portable` provenance.
  Matrix-backed and built-in composite operators derive deterministic portable
  identities from dimensions, structure, transformation parameters, and
  values. Callback operators without explicit provenance receive opaque
  session-local identities and restored cross-session plans fail before the
  callback is invoked. `operator_identity()` exposes the typed records.
* Eigen and SVD results distinguish `planned_method` from `actual_method`
  (`method` remains the compatibility alias for the latter), report
  `fallback_used` and a typed `fallback_reason`, and retain the exact immutable
  plan that governed execution. Legacy result fields remain available.
* `work()` returns a schema-versioned `eigencore_work` record with separate
  logical block-call and column counts for the operator, adjoint, metric,
  preconditioner, and current-certificate phases, plus iterations, restarts,
  and phase timings. Callback counts are captured at the operator boundary,
  including native callback cycles; `matvecs` remains unchanged as
  `legacy_matvecs`. Routes whose native controller cannot yet prove a split
  use `NA` and `complete = FALSE` rather than inventing equivalent work.
* `restart_state()` constructs immutable, schema-versioned basis states only
  from freshly certified eigen or SVD results, while `retained_bytes()` reports
  their exact version-1 R retention. `solve(plan, restart_state = ...)` admits
  basis reuse only for standard real Hermitian Lanczos on dense double, CSC,
  and callback operators. Exact-revision reuse may retain a deterministic,
  target-safe fitted start block, but never recurrence, locks, cached operator
  actions, convergence, or certificates. Changed revisions and lineages use
  only the public basis; incompatible and unsupported routes fail before the
  current operator is applied.

## Sparse operator performance

* `scale_cols(center(A, columns = TRUE), weights)` now recognizes a real
  `dgCMatrix` and builds one fused native operator for
  `(A - 1 mu^T) D`. Forward and adjoint block application, including
  `alpha`/`beta` accumulation, stay in C++; native Golub-Kahan consumes the
  sparse matrix, means, and weights directly, without materializing the
  centered matrix or crossing an R callback boundary inside the iteration.
  One CSC moments pass supplies exact centered-and-scaled Frobenius metadata,
  so the resulting two-sided SVD certificate reports
  `norm_bound_type = "frobenius_metadata"` and
  `scale_is_estimate = FALSE`. Zero columns, zero/negative/extreme finite
  weights, transpose algebra, planner provenance, and strict installed-package
  execution are covered by the v1.2 sparse-PCA gate.

# eigencore 1.1.0 (2026-08-24)

## New features

* `eig_partial()` and `solve()` gain an `initial_subspace` argument: a public,
  certified warm-start seam for standard real Hermitian Lanczos — the native
  paths on explicit dense double and `dgCMatrix` operators, the native
  matrix-free callback path selected by `lanczos(block > 1)`, and the scalar
  matrix-free reference path selected by `lanczos(block = 1)`. The supplied
  subspace is orthonormalized at the solver boundary and fitted to the method's
  start block: accepted directions are augmented deterministically when short,
  and when the accepted rank exceeds the block width the block is a seeded
  random rotation of the full accepted basis, so every supplied direction
  contributes generic weight (a k-column continuation subspace handed to a
  scalar method warm-starts all k targets, not just the first). Rank detection
  is invariant to column scaling, and the fitted start block is orthonormal.
  The start is treated only as a hint — every solve recomputes projected
  quantities, residuals, orthogonality, convergence, and a fresh
  current-operator certificate. A fully supplied subspace already invariant at
  the requested tolerance is discarded to a cold start: residual certification
  proves eigenpair accuracy, not that an invariant block contains the requested
  extremal pairs. `initial_subspace = NULL` preserves cold behavior exactly.
  Unsupported generalized, shift-invert, and dense-fallback plans error.
  Results, plans, and `diagnostics()` report start provenance plus exact
  operator block calls, operator columns, and certification columns. Reusable
  restart-state objects and generalized/transformed promotion remain future
  work.
* `lanczos()` gains `check_stride`. The default `0L` preserves full-sweep
  convergence checks exactly; a positive stride lets native block
  thick-restart paths, including real Hermitian matrix-free callbacks, check
  and stop within a sweep without additional operator applications.

## Performance

* Operator identities and workflow tokens use a native 128-bit structural
  hash over the data buffers instead of hashing `serialize()` output, about
  10x faster (dense 1500 x 1500 source: 0.047 s -> 0.005 s; 4000 x 4000:
  0.54 s -> 0.04 s). Equal values hash equal (`-0`/`0`, NaN payloads),
  attributes count regardless of order, and digests are the same across
  sessions and platforms. **Identity format change:** built-in operator
  identities, plan tokens and restart-state tokens all change value. Plans
  and restart states now record `serialization$hash_format`; ones saved by an
  earlier version are rejected with code `identity_format_changed` and a
  message asking to re-plan, rather than a generic identity mismatch.
  PSD factors (`psd_factor()`, `psd_gram_factor()`, ...) record the same
  `serialization$hash_format`; a factor persisted by an earlier version
  fails with an `eigencore_psd_corrupt_state` condition of code
  `identity_format_changed` asking to re-factor, instead of a generic
  integrity-token error. The unused legacy FNV entry point
  `eigencore_stable_raw_hash` is removed.
* Dense matrix operator construction screens finiteness and symmetry in one
  native pass instead of an R-level `all(is.finite())` plus a separate
  symmetry scan (dense `Matrix` classes are no longer screened twice).
  `as_operator()` at n = 4000: 0.17-0.21 s -> 0.058 s; `plan_solver()`:
  0.20 s -> 0.10 s (the rest is the identity hash). Errors and the
  relative symmetry tolerance are unchanged.
* `compose()`, `crossprod_operator()` and operator sums no longer form a
  dense product twice, and they materialise a product only when it is cheap
  and memory-safe: always for a product with at most 65,536 entries; a
  larger dense product only when it holds no more entries than its factors
  (for `crossprod_operator()`, no more columns than rows) and costs at most
  2^30 multiply-adds; a larger sparse product only when a structural
  nonzero bound stays within 4x the factors' nonzeros and below a quarter
  of its entries, so a sparse product is never densified. Otherwise the
  result is a lazy composition (e.g. a 3000 x 5 times 5 x 3000 dense
  composition: 0.12 s and 69 MB -> 0.001 s and no stored product).
* Sparse tridiagonal shift-invert now parses and validates the three matrix
  bands once per solve and reuses that immutable representation for planning,
  shift perturbation, factorization, and certification. The native kernel
  forms and returns only the requested Ritz vectors instead of exposing its
  full Krylov basis to R. On the installed `path_laplacian:1000`, `k = 20`
  release row this reduced cumulative R allocation from 3.82 MB to 1.51 MB
  while retaining a 20/20 original-coordinate certificate. The G1 gate now
  uses retained result size for its bounded memory envelope and reports
  `bench::mem_alloc` separately as diagnostic evidence, because R allocation
  counters do not observe native C/C++ working heaps in any compared engine.
* New production `auto` route for largest-target partial SVD: a native
  implicit normal-equations (Gram) thick-restart Lanczos that runs on
  \eqn{A^T A} or \eqn{A A^T} as an operator, without materializing the Gram
  matrix. It covers dense matrices and sparse operators whose smaller side
  exceeds the explicit-Gram caps — regimes that previously fell back to a
  full LAPACK SVD (dense) or a single unrestarted Golub-Kahan sweep
  (sparse). Representative same-machine speedups at `tol = 1e-8` with
  certificates still passing: dense 4000x1000 `k = 10` ~12x, dense
  2000x2000 `k = 10` ~80x, sparse 20000x5000 `k = 10` and `k = 50` ~2x.
  Results keep the exact two-sided residual certificate in original
  coordinates, and uncertified results still fall back to the native
  Golub-Kahan path.
* The scalar Golub-Kahan kernel now uses BLAS (dgemv) classical
  Gram-Schmidt reorthogonalization on the sparse CSC, matrix-free, and
  retained-restart paths, matching the dense path; previously these used
  scalar loops.
* The projected Golub-Kahan convergence check now runs `dbdsqr` directly on
  the projected bidiagonal, tracking only the last row of the left singular
  vectors, instead of a dense `dgesvd` with full vectors; scratch buffers
  are reused across checks. The check drops from O(iter^3) plus per-check
  allocations to O(iter^2).
* The scalar Lanczos convergence estimate (used by the shift-invert paths)
  computes eigenvalues with `dsterf` and recovers only the selected
  eigenvectors with `dstevr` instead of a full `dstev` decomposition each
  iteration, with scratch reuse across iterations.
* The CSC sparse matrix-multiply kernel now processes wide blocks in
  cache-friendly chunks of 10 columns instead of a strided generic loop.
* `as_operator()` no longer forces a full copy of an already-double dense
  matrix (`storage.mode<-` is now conditional), removing an O(m*n) copy
  plus GC churn from every problem construction — operator construction on
  a 2000x2000 input drops from ~52ms to ~6ms.
* The scalar thick-restart Hermitian Lanczos default subspace grows from
  `3k+20` to `max(3k+20, 5k)`: unchanged for small `k`, and at `k = 30`
  on a general sparse operator it cuts operator applications by ~25%.

## Portability

* The package again compiles on R < 4.4: guarded compatibility typedefs for
  `La_INT`/`La_LGL` and declarations for the complex QZ drivers
  (`zggev`, `zgges`) were added for headers that predate them, and
  `crossprod`/`tcrossprod` are now imported from Matrix so sparse-matrix
  dispatch does not rely on the base generics added in R 4.4.

## Bug fixes

* Generalized result contracts are now consistent across dense pencils, QZ,
  transformed sparse pencils, and GSVD. Every classified result exposes a
  `classification_policy`; GSVD records exact structural-zero semantics and
  documents its length-`n` `Inf`/`NA` value layout. Dense general-pencil left
  vectors now carry an original-coordinate adjoint-residual certificate and
  `W^H B V` diagnostics, and `diagnostics()` consistently exposes the same
  `left_vectors` matrix returned by `left_vectors()`.
* Benchmark/validation timing helpers no longer error on R builds without
  memory profiling (e.g. r-devel Linux fedora, configured without
  `--enable-memory-profiling`). `bench::mark()` memory measurement is now
  requested only when `capabilities("profmem")` is `TRUE`; timing still runs
  otherwise and `mem_alloc` is reported as `NA`. Benchmark smoke tests also use
  `testthat::skip_on_cran()`. Backported in 1.0.1 and completed in 1.0.2.

# eigencore 1.0.2

* Complete the no-memory-profiling fix from 1.0.1. The CRAN benchmark
  vignette and the installed README benchmark now also request
  `bench::mark()` memory measurements only when `capabilities("profmem")` is
  true. The vignette's summary tables also report unavailable results instead
  of failing when a benchmark regime has no successful methods.
* Add a regression test that audits every shipped benchmark entry point for
  unconditional `memory = TRUE` and exercises the package timing helper on the
  current R build.

# eigencore 1.0.1

* Fix an R CMD check ERROR on R builds without memory profiling (for example
  the r-devel Linux fedora flavors, which are configured without
  `--enable-memory-profiling`). The internal benchmark/validation timing
  helpers passed `memory = TRUE` to `bench::mark()`, which calls
  `utils::Rprofmem()` and aborts with "memory profiling is not available on
  this system" on those platforms. Memory measurement is now requested only
  when `capabilities("profmem")` is `TRUE`; elsewhere timing still runs and
  `mem_alloc` is reported as `NA`.
* Benchmark smoke tests now use `testthat::skip_on_cran()` so they are skipped
  on CRAN as intended (the previous `Sys.getenv("CRAN")` guard never fired on
  the CRAN check farm).

# eigencore 1.0.0

First CRAN release.

* `svd_partial()` and `eig_partial()` compute the top-*k* singular triplets or
  eigenpairs of large dense, sparse (CSC), diagonal, banded/tridiagonal, and
  matrix-free operators through native C++ kernels.
* Every result carries a numerical certificate: residuals for both singular
  relations, a backward-error bound, orthogonality loss, a labeled norm bound,
  and a single `passed` flag. Bounds that can only be estimated (for example
  stochastic norm estimates on centered sparse operators) are reported as
  estimates and never produce an unqualified `passed`.
* Operator algebra — `center()`, `scale_cols()`, `compose()`,
  `crossprod_operator()`, `linear_operator()` — solves centered, scaled, and
  composed problems without forming dense matrices.
* Transparent method selection: `plan_solver()` reports the chosen kernel
  before a solve, and `fit$method` names the path that actually ran. Problem
  classes without a production kernel carry explicit `reference` labels.
* RSpectra-compatible wrappers `eigs()`, `eigs_sym()`, and `svds()` accept the
  same `which` codes and additionally return certificates.
* Benchmarked against 'RSpectra', 'irlba', and 'PRIMME'; reproduce with
  `Rscript inst/benchmarks/bench-readme.R`.
* Generalized eigen support: `eig_full()` for dense SPD and general pencils,
  `generalized_schur()` and `generalized_svd()` for dense QZ/GSVD, partial
  sparse general pencils with nonsingular diagonal `B` via transformed native
  Arnoldi, left eigenvectors and conditioning diagnostics on supported dense
  paths, and `pencil_norm_scaled` alpha/beta classification. Sparse SPD partial
  paths remain under `eig_partial()` / LOBPCG / B-orthogonal Lanczos; general
  sparse QZ and non-diagonal sparse `B` are explicit unsupported boundaries.
  The real dense GSVD path currently requires a linked LAPACK that provides the
  deprecated `dggsvd` routine.
* The exported API is stable as of 1.0.0; breaking changes from here follow
  semantic versioning.
