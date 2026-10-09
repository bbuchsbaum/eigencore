# eigencore test assurance

## Native code assurance

Date: 2026-10-09. Tree: branch `assure-native` (from `test/assurance` at
`294eee5`). Scope: the C++ under `src/` after the 2026-10 review refactors
(exception-based unwind wrapper, OpenMP CSC kernels with CSR / row-slab
caches, native composite operators, Krylov-Schur Arnoldi, the LOBPCG
rewrite, the completeness probe, the native identity hash). Machine: Linux,
4 shared cores, R 4.3.3 (Debian build, not instrumented), OpenBLAS, gcc 13,
valgrind 3.22.

### How to reproduce

ASan + UBSan (R itself stays uninstrumented; only the package is compiled
with the sanitizers through the `src/Makevars` hooks, which `make` reads from
the environment, and the ASan runtime is preloaded):

```sh
R CMD build --no-build-vignettes --no-manual .
EIGENCORE_SANITIZER_CXXFLAGS='-fsanitize=address,undefined -fno-omit-frame-pointer -fno-sanitize-recover=undefined -g' \
EIGENCORE_SANITIZER_LIBS='-fsanitize=address,undefined' \
  R CMD INSTALL --no-test-load -l "$SANLIB" eigencore_1.3.0.tar.gz

export R_LIBS="$SANLIB" NOT_CRAN=true
export LD_PRELOAD="$(gcc -print-file-name=libasan.so)"
export ASAN_OPTIONS=detect_leaks=0:abort_on_error=1:allocator_may_return_null=1
export UBSAN_OPTIONS=print_stacktrace=1:halt_on_error=1
EIGENCORE_TEST_THREADS=4 Rscript -e 'library(eigencore); testthat::test_dir("tests/testthat", package = "eigencore", load_package = "installed")'
Rscript inst/validation/native-smoke.R
Rscript inst/validation/native-entry-points.R
```

- `allocator_may_return_null=1` is required: the C10 unwind tests request
  absurd sizes on purpose and expect `std::bad_alloc` / an R error, which ASan
  would otherwise turn into an abort (`allocation-size-too-big`).
- `detect_leaks=1` is usable: R does not leak at exit, so LeakSanitizer
  reports are actionable (verified with a deliberately leaking `.so`). The
  only report in a full run is one byte leaked by `sed`, which the R
  front-end shell script runs under the same `LD_PRELOAD`.
- `EIGENCORE_TEST_THREADS=<n>` (read by
  `tests/testthat/helper-native-assurance.R`) runs the whole suite with
  `options(eigencore.threads = n)`.
- The RSS-based C10 leak test is skipped under ASan (its quarantine keeps
  freed blocks resident by design; LeakSanitizer covers the same question).
- CI: `.github/workflows/sanitizers.yaml` (manual + weekly) does exactly this
  on `ubuntu-latest` with 1 and 4 threads and a LeakSanitizer job.

valgrind memcheck (plain `-O2 -g` build):

```sh
R -d "valgrind --error-exitcode=1 --track-origins=yes" --vanilla \
  -f inst/validation/native-entry-points.R
```

`inst/validation/native-entry-points.R` reaches 109 of the 133 registered
`.Call` entry points at small sizes (public API plus the package's internal
wrappers for routes the planner does not pick at these sizes), including the
error paths: malformed CSC (`i` out of range, inconsistent `p`), non-finite
input, singular shift-invert (`perturb sigma`), the unwind selftest
(error / `bad_alloc` / `length_error` / huge vector / R allocation failure /
interrupt), singular tridiagonal solve, and R callback errors inside native
Lanczos, Krylov-Schur, Golub-Kahan and implicit-Gram solves. The 24 entries it
does not reach have no R caller at these sizes or none at all (e.g.
`eigencore_lanczos_csc`, `eigencore_block_lanczos_dense/csc`,
`eigencore_dense_symmetric_eigen_dsyevx_selected`,
`eigencore_dense_complex_generalized_svd`, the `_cached` block Golub-Kahan
variants); the test suite calls several of them directly.

Thread determinism: `Rscript inst/validation/thread-determinism.R 1 2 4`.

### Results

| Check | Configuration | Result |
|---|---|---|
| ASan + UBSan, full suite (first run, before this branch's changes) | threads 1 and 4 | 0 sanitizer reports; 8107 expectations, 1 failure = stale `api-surface` snapshot (F2) |
| ASan + UBSan, full suite (final, with F1 fix and new tests) | threads 1 | 0 sanitizer reports; 8213 expectations, 0 failures, 12 skips |
| ASan + UBSan, full suite (final) | threads 4 | 0 sanitizer reports; 8213 expectations, 0 failures, 12 skips |
| ASan + LeakSanitizer, full suite | threads 2, `detect_leaks=1` | 0 leaks attributable to eigencore (only `sed`'s 1 byte) |
| ASan + UBSan + `float-cast-overflow` + `float-divide-by-zero` | threads 2, full suite + entry-point script | 0 reports (these two checks are not part of `-fsanitize=undefined` in gcc; `float-divide-by-zero` ran in recover mode, so it would have printed even benign IEEE divisions: none occurred) |
| ASan + UBSan, `native-smoke.R`, `native-entry-points.R` | threads 2 and 4 | pass, 0 reports |
| valgrind memcheck, `native-entry-points.R` | threads 2 | `ERROR SUMMARY: 0 errors` |
| valgrind memcheck, test suite minus the four heaviest files (`t4a-omp`, `bench-smoke`, `t4b-svd`, `t3-api`) | threads 2 | `ERROR SUMMARY: 0 errors`, `definitely lost: 0 bytes`; 6934 expectations, 1 failure that is a valgrind artifact (F7) |
| valgrind memcheck, `test-native-assurance.R` + `test-t3-unwind.R` + `test-native-hardening.R` + entry-point script, after the F1 fix | threads 2 | `ERROR SUMMARY: 0 errors`, `definitely lost: 0 bytes`; 456 expectations, 0 failures |
| Full suite, plain build (final) | threads 1 / 4 | 8415 expectations, 0 failures, 11 skips, identical outcome at both thread counts (ASan counts are lower because the 200-expectation RSS test is skipped there) |

No memory error, undefined behaviour, data race symptom or leak was found in
the native code. One error-path defect was found and fixed (F1 below).

### Findings

**F1 (fixed): R callback errors lost their message in every matrix-free
native solver.** `eigencore_r_operator_apply()` evaluated the user's apply
closure with `R_tryEval()` and returned status `-8` on error, so thick-restart
Lanczos (Hermitian `auto()`), Krylov-Schur Arnoldi, Golub-Kahan, the
implicit-Gram SVD, LOBPCG with an operator `B`, and the completeness probe all
failed with `native ... failed with status=-8`, while the actual message was
only printed to stderr. Reproducer:

```r
calls <- 0L
op <- linear_operator(c(50, 50), apply = function(X, ...) {
  calls <<- calls + 1L; if (calls > 3L) stop("boom"); X
}, apply_adjoint = function(X, ...) X, structure = hermitian())
eig_partial(op, 3L)   # before: "native block thick-restart Lanczos failed with status=-8"
                      # after:  "R operator callback failed: boom"
```

Fix: evaluate with `R_tryEvalSilent()` and raise
`R operator callback failed: <message>` through the exception-based unwind
(all callers run below an `eigencore_call` entry point, so the solver's
buffers are released before the R error). Regression test:
`test-native-assurance.R` ("errors raised inside R operator callbacks keep
their message ...").

**F2 (fixed, test-only): stale `api-surface` snapshot.** The frozen
certificate field list predated `norm_source`/`norm_values` (C12),
`residual_passed`/`target_passed` and `target_completeness` (C50), so the
snapshot test failed whenever `NOT_CRAN=true` (CI, coverage, sanitizer runs).

**F3 (not a bug, documented): sanitizer configuration.** Running the suite
under ASan needs `allocator_may_return_null=1` (see above) and skipping the
RSS-growth assertion; both are now encoded in the CI workflow and the test.

**F4 (open, route quality, not memory): `eigs_sym(sparse, which = "BE")`.**
Sparse `both_ends` goes to the R reference Hermitian Lanczos
("target unsupported by native path") and, on a 60 x 60 path Laplacian plus
diagonal, `k = 3`, returns 0 converged pairs with a warning and a failed
certificate (backward errors up to 3e-4); the dense input certifies. Honest
(uncertified) but a convergence gap.

**F5 (open, route quality): explicit `lanczos()` on a matrix-free Hermitian
operator** takes the R-level prototype Lanczos ("native hot loop not yet
implemented"), stops after 29 applies on n = 60 and fails its certificate;
`maxit = 500` and `max_subspace = 30` do not change the number of applies.
`auto()` (C53) routes the same operator to the native thick-restart kernel and
certifies. The `lanczos()` method descriptor should route there too, or honour
`maxit`.

**F6 (open, message quality):** a singular
`shifted_tridiagonal_preconditioner()` (e.g. a path Laplacian with shift 0)
fails LOBPCG with the bare `native CSC LOBPCG failed with status=-5`.

**F7 (valgrind artifact; low-severity portability note).** Under valgrind,
`test-identity-hash.R` ("C45: identity is identical across separate R
sessions") fails: the sparse operator's identity differs from the one computed
in a native child session. Cause: `as_operator(<dgCMatrix>)` stores column
means / centred sums of squares computed with `long double` accumulators in
its metadata, and the identity hashes that metadata; valgrind models x87
`long double` with 64-bit precision, so e.g. the mean `1/3` comes out one ulp
different (`0x1.5555555555556p-2` vs `0x1.5555555555555p-2`). This is a
valgrind limitation, not a memory error, but it shows that the persisted
identity of a sparse operator depends on `long double` arithmetic of the
platform (80-bit x86, 128-bit aarch64 Linux, 64-bit on macOS arm64 / MSVC),
so in rare double-rounding cases a plan or restart state persisted on one
platform will not match the same matrix on another. Hashing only the source
slots (or computing the moments in double with compensated summation) would
make identities platform-independent.

### OpenMP race review (static)

All OpenMP regions (`EIGENCORE_OMP`) were read: `reorth_gemm_tn` /
`reorth_gemm_nn_minus` (`block_lanczos.cpp`), the CSC forward gather / slab
scatter / column-chunk paths, the adjoint gather, and the randomized
projection kernel (`native_operators.cpp`).

- Every output entry has exactly one writer: `C[j, col]` per loop index in
  `reorth_gemm_tn`; disjoint 512-row blocks in `reorth_gemm_nn_minus`;
  disjoint nonzero-balanced row ranges in the CSR gather; per-slab row ranges
  (including the shared row-major `panel`, whose rows are partitioned the same
  way) in the slab scatter; one output column per iteration in the adjoint and
  projection gathers. No reductions, atomics or critical sections are needed
  and none are used.
- The forward gather's shared column panel / `skip` flags are filled in an
  `omp for` whose implicit barrier precedes the gather that reads them.
- The lazily built caches (CSR copy, row slabs, panel, `skip`) are built and
  resized only on the main thread before a parallel region starts, and are
  per operator; no parallel region applies an operator or calls the R API,
  and kernels inside regions neither allocate nor throw.
- Thread-count changes between calls rebuild the slabs for the new part
  count; the slab loop strides over parts, so it is also correct when the
  runtime grants fewer threads than requested.
- Global state (`g_eigencore_threads`, BLAS quiesce/restore, call depth) is
  touched only on the main thread. A body that `longjmp`s past
  `eigencore_call_leave()` (an unprotected small R allocation failing) would
  leave the call depth non-zero and BLAS at one thread for the rest of the
  session; this is the documented C10 residual risk.

### Thread determinism

`inst/validation/thread-determinism.R` (threads 1, 2, 4): all raw CSC /
centred-scaled applies (forward, adjoint, cached-CSR wide matrix, row-slab
tall matrix, repeated applies) are bitwise identical across thread counts.
The implicit-Gram and Golub-Kahan SVDs, `svds`, Krylov-Schur Arnoldi, LOBPCG
and randomized SVD are bitwise identical as well. Lanczos eigensolves (CSC,
block) differ between 1 and >1 threads by rounding only (values <= 2.4e-14
relative, vector subspaces <= 4e-8 at `tol = 1e-8`), because with more than
one thread the reorthogonalisation switches from BLAS `dgemm` to the OpenMP
kernels, as documented in `?eigencore-threads`; 2 and 4 threads agree
bitwise. The full suite passes identically with 1 and 4 threads (plain and
ASan builds).

Note: `svds()` has no `seed` argument and consumes the global RNG stream, so
successive unseeded calls differ; seed the RNG when comparing (C26).

### Coverage

`covr::package_coverage(type = "tests")` with `NOT_CRAN=true` (no line
exclusions are configured), R and C++ (gcov). "Changed code" = lines added or
modified in `git diff 83ca827..HEAD -- R src` that gcov/covr count.

| | 83ca827 (pre-review) | 294eee5+fixes (before new tests) | after new tests |
|---|---|---|---|
| overall | 85.38% (26393/30912) | 86.33% (31326/36285) | 87.27% (31654/36272) |
| R | 86.80% | 87.25% | 87.42% |
| C++ | 83.12% | 85.03% | 87.05% |
| changed lines, all | - | 89.50% (7243/8093) | 90.97% (7362/8093) |
| changed lines, R | - | 89.82% | 90.46% |
| changed lines, C++ | - | 89.27% (4286/4801) | 91.31% (4384/4801) |

Per C++ file (line %):

| file | 83ca827 | before | after |
|---|---|---|---|
| arnoldi.cpp | 91.4 | 88.6 | 90.5 |
| block_golub_kahan_basis.cpp | 85.5 | 85.7 | 85.7 |
| block_lanczos.cpp | 77.5 | 78.6 | 90.7 |
| certificates.cpp | 92.3 | 92.5 | 92.5 |
| completeness_probe.cpp | - | 95.6 | 95.6 |
| eigencore_common.h | 83.3 | 87.1 | 94.1 |
| eigencore_unwind.cpp | - | 78.4 | 81.6 |
| gram_svd.cpp | 95.7 | 95.1 | 95.1 |
| identity_hash.cpp | - | 92.7 | 100.0 |
| lobpcg.cpp | 79.6 | 83.8 | 83.8 |
| native_operators.cpp | 89.6 | 86.8 | 87.6 |
| orthogonalization.cpp | 85.4 | 86.1 | 86.1 |
| projection.cpp | 86.5 | 86.9 | 86.9 |
| retained_svd.cpp | 78.1 | 78.1 | 78.4 |
| scalar_krylov.cpp | 75.2 | 83.2 | 83.2 |
| small_dense.cpp | 79.3 | 84.3 | 84.3 |

New tests (`tests/testthat/test-native-assurance.R`): non-restarted block
Lanczos entries (dense/CSC, both targets, argument errors, malformed CSC);
composite spec validation (missing/unknown node type, bad dense leaf, empty
product, non-conformable product, non-finite sum weights, nesting depth),
cleared (deserialised) composite pointers, rank-one leaves with `beta != 0`
and adjoint; identity hash of calls, pairlists, symbols, expressions and S4
objects; Krylov-Schur breakdown on an invariant start subspace (random
restart path); overflow/underflow-safe residual norms at 1e+-150 scale (dense
and CSC, block and single-vector Lanczos); repeated eigenvalues at the
selection boundary; LOBPCG when the block exhausts a tiny search space
(standard and generalized); retained IRLBA under all five
reorthogonalisation policies and both targets; nonsymmetric shift-invert with
sparse transposed (left-vector) solves, `certify = FALSE`, and the
`factorization =` / generalized guard rails; R callback errors through every
matrix-free native solver (F1).

Largest uncovered changed regions that remain, and why:

- LAPACK non-convergence fallbacks: `tridiagonal_selected_eigenpairs`
  (`dstevr` -> `dstev`, scalar_krylov.cpp:1185), `eigencore_bidiagonal_svd`
  (`dbdsdc` -> `dbdsqr`, small_dense.cpp:1995), `eigencore_dense_svd`
  (`dgesdd` -> `dgesvd`, small_dense.cpp:1465), `eigencore_tridiagonal_eigen`
  (small_dense.cpp:1778). They need a failing LAPACK (fault injection) to
  reach.
- `eigencore_blas_probe` FlexiBLAS branch (native_operators.cpp:507): needs
  FlexiBLAS.
- `projected_eigen_selected` `dsyevr` tie fallback (block_lanczos.cpp:275):
  exact repeated eigenvalues at the boundary did not make bisection return
  extra values in practice.
- Retained IRLBA augmented-projection fallbacks (retained_svd.cpp:2180-2730):
  defensive branches for breakdown/rank loss; the C42 sweeps never take them.
- R side: shift-invert retry/fallback branches
  (`R/transform_shift_invert.R`, 89 lines), Arnoldi certificate operator
  split for complex vectors (`R/reference_arnoldi.R:289`).

### Not done

- ThreadSanitizer was not run (R is not TSan-instrumented and libgomp needs a
  TSan-built runtime to avoid false positives); the race review above is
  static.
- valgrind was not run on the four heaviest test files (`t4a-omp`,
  `bench-smoke`, `t4b-svd`, `t3-api`); they ran clean under ASan/UBSan.
