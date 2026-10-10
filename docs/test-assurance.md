# eigencore test assurance

Two complementary checks of the package after the 2026-10 review refactors:

1. [Oracle-based differential sweep](#oracle-based-differential-sweep): "certified implies correct" against independent references across the input space.
2. [Native code assurance](#native-code-assurance): sanitizers, valgrind, thread determinism and coverage of the C++ code.

## Oracle-based differential sweep

Date: 2026-10-09. Branch `assure-oracle` (from `test/assurance` at `294eee5`).

The per-bug regression tests in `tests/testthat` pin the defects found by the
review (`docs/review-2026-10.md`). They do not say much about inputs nobody
thought of. The oracle sweep fills that gap. It generates problems across the
input space, solves each one with eigencore, and checks the result against a
dense oracle (base `eigen()`/`svd()`, or a Cholesky reduction for SPD
pencils). The property it checks hardest is **"certified ⇒ actually
correct"**.

### Files

| File | Role |
|---|---|
| `tests/testthat/helper-oracle.R` | Case generator, oracle, invariant checks, subprocess runner, summary tables |
| `tests/testthat/test-oracle-sweep.R` | Runs the level's cases. Fails on any hard violation and prints a one-line reproducer per case. Also has harness self-tests (a doctored result must trip each class of check) |
| `tests/testthat/test-oracle-determinism.R` | Same seed gives identical results. `eigencore.threads` 1 vs 4 agree to 1e-12. `seed =` leaves the global RNG alone |
| `tests/testthat/test-oracle-shims.R` | `eigs_sym`/`eigs`/`svds` against RSpectra on well-separated spectra: same values in RSpectra's order |
| `tests/testthat/test-oracle-regressions.R` | Minimised reproducers for every violation found (O1–O18). Known, unfixed items are `skip("known: O<n>")` |
| `tools/oracle_sweep.R` | Stand-alone parallel runner. Writes `records.rds`/`.csv` and `summary.md` |
| `.github/workflows/oracle-sweep.yaml` | Runs the `extended` level weekly and on dispatch, then uploads the results |

### Levels

| Level | Case ids | Size cap | Default when | Runtime here (4 shared cores) |
|---|---|---|---|---|
| `cran` | 1–40 | n ≤ 40 | `NOT_CRAN` unset (R CMD check on CRAN). Runs in-process | ~3 s |
| `ci` | 1–400 | n ≤ 160 | `NOT_CRAN=true` (devtools, CI). 2 callr workers | ~35 s |
| `extended` | 1–15000 | n ≤ 300 | `EIGENCORE_ORACLE_LEVEL=extended` | ~15 min with 3 workers |

```sh
EIGENCORE_ORACLE_LEVEL=ci  Rscript -e 'testthat::test_file("tests/testthat/test-oracle-sweep.R")'
Rscript tools/oracle_sweep.R --level extended --workers 3 --out oracle-out
Rscript tools/oracle_sweep.R --ids 1:200,4711 --workers 2
```

The levels are nested prefixes of a single id stream. A case is a pure
function of its integer id: `oracle_case(id)` draws every axis from a private
RNG seeded by the id, and the global RNG is saved and restored. To re-run a
failing case on its own:

```r
source("tests/testthat/helper-oracle.R"); library(eigencore)
str(oracle_run_case(4711))   # record: status, certified, hard/soft findings, metrics
case <- oracle_case(4711); prob <- oracle_build(case)   # the input object and dense truth
```

At the `ci` and `extended` levels, cases run in batches of 20 inside `callr`
subprocesses, with one eigencore thread and one BLAS thread per worker. If a
batch crashes or hangs (900 s), every case in it is re-run alone with a
180 s limit. A case that still crashes or hangs is reported with its id.

### Generator coverage

| Axis | Values |
|---|---|
| Family (20-slot cycle) | real symmetric 30%, SVD 25%, nonsymmetric 15%, RSpectra shims 15%, SPD pencils 10%, complex 5% |
| Storage | base dense, `dgCMatrix` (built exactly from triplets), `dsCMatrix` upper and lower, `dgTMatrix`, `dgRMatrix`, `ddiMatrix`, `dgeMatrix`, `dsyMatrix` (both built from slots, because `Matrix::Matrix()` silently symmetrises tiny nonsymmetric input), matrix-free `linear_operator` (no norm metadata), and composites: `crossprod_operator(X)`, `crossprod_operator(center(X))`, `crossprod_operator(scale_cols(X, w))`, sparse `crossprod_operator`, `compose(A, D)` (nonsymmetric), `center(X)` dense and sparse, `scale_cols(X, w)` and `compose(X, Y)` (SVD). Shims also take function input (`n=`, `Atrans=`/`dim=`) |
| Structure | real symmetric, real nonsymmetric (prescribed `V J V^-1` with known `V`, including 2×2 complex-pair blocks and Jordan blocks), random real matrices, complex Hermitian, complex general, generalized SPD pencils `A x = λ B x` (dense, sparse, diagonal `B`), rectangular tall/wide/square |
| Spectrum | random, clustered (relative spread 1e-4), exactly repeated (multiplicities 2–8), near-repeated (1e-10 splits), power-law decay, tiny (1e-12), very tiny (1e-20, ‖A‖ < eps), huge (1e12), indefinite (symmetric about 0), singular/rank-deficient, defective (Jordan blocks), complex pairs |
| Size | n in {5, 6, 8, 12, 20, 30, 50, 80, 120, 200, 300}, weighted toward small |
| k | 1; small (2–6); near n (n, n−1, n−2) |
| Target | `largest`, `smallest`, `largest_magnitude`, `smallest_magnitude`, `both_ends`, `nearest(σ)`, `largest_real`/`smallest_real`, `largest_imaginary`/`smallest_imaginary` (nonsymmetric, complex). SVD: `largest`, `smallest`, `nearest(σ)` |
| Method | `auto()`, `lanczos(block = 1/2/3)`, `lobpcg()`, `shift_invert(σ)` (85% paired with `nearest(σ)`), `golub_kahan()`, `randomized()` |
| tol | 1e-8 (60%), 1e-10, 1e-6, 1e-12 |
| Shims | `eigs_sym` (`LM/SM/LA/SA/BE`, `sigma`, `lower = TRUE/FALSE` with the other triangle filled with junk, `dsCMatrix`), `eigs` (`LM/SM/LR/SR/LI/SI`, `sigma`), `svds` (`center`, `scale`, `nu`/`nv` in {k, 0, 1}); `opts$tol`, `opts$ncv`, `opts$maxitr`, `opts$retvec` |

### Invariants

**Hard** (a test failure):

* **(a) No crash.** Each case returns a result or an error. Error messages are
  matched against a list of *expected* "unsupported input" messages, such as
  `shift_invert() currently requires a Hermitian eigenproblem`. Internal-looking
  errors (`subscript out of bounds`, `status=-`, LAPACK argument errors,
  `max_subspace must be at least` when no subspace was requested, and so on)
  and any unmatched message count as *unexpected*. A dead or hung subprocess
  is a failure.
* **(b) Soundness**, checked whenever `certificate$passed` is TRUE:
  * The harness computes its own backward error, `‖Ax − λBx‖ / ((‖A‖₂ + |λ|‖B‖₂)‖x‖)`.
    For the SVD it is `‖[Av − σu; Aᴴu − σv]‖ / ‖A‖₂`. Here ‖A‖₂ comes from a
    dense SVD or eigendecomposition, and 0/0 is taken as 0. This must be at
    most `tol·1.05 + 64·n·eps`.
  * Regardless of `passed`, the reported per-pair backward error must not be
    below the true one. C12 promises an over-estimate.
  * Vectors must be orthonormal to within `2·max(tol, √eps)`: B-orthonormal
    for pencils, and both U and V for the SVD.
  * Every certified value must lie within its rigorous residual bound of a
    true eigenvalue. Hermitian: `‖r‖/‖x‖`. Pencils: `‖r‖_{B⁻¹}/‖x‖_B`. SVD:
    the Jordan–Wielandt bound. Nonsymmetric: Bauer–Fike with κ(V) taken from
    the generator's exact `V`, or Chatelin's bound `tᵐ ≤ κε(1+t)^{m−1}` for
    Jordan blocks of size m.
  * If `target_completeness` is `exact`, `probed`, `repaired` or
    `inertia_verified`, the returned multiset must equal the oracle's target
    set, including multiplicity. Scores are compared in sorted order, within
    `2√k·max‖r‖ + 256·n·eps·‖A‖`, so ties at the target boundary are allowed.
* **(e) API.** A certificate may not pass with fewer than k values. Also
  checked: `nconv ≤ k`, the dimensions of vectors, U and V (and `nu`/`nv` for
  `svds`), and that values are ordered per target (`both_ends`: low end
  ascending, then high end descending; `eigs_sym`/`svds`: decreasing). When
  RSpectra is installed and its answer matches the oracle, a certified shim
  result must match RSpectra's values (compared by score, so ties may differ).

**Soft** (counted and reported in the summary, never a failure):

* **(c)** Uncertified results, with their value error against the oracle.
* A certified result whose set differs from the oracle while the certificate
  makes **no** completeness claim (`not_checked`, or SVD, which has no
  completeness field). The per-pair certificate is still true here, but the
  set is wrong. The known classes are listed below as O6–O10 and O16–O18.
* Expected "unsupported" errors, RSpectra disagreements in either direction,
  and differences in order alone.

**(d) Determinism** (in `test-oracle-determinism.R`):

* Every non-shim sweep case, solved twice with the same seed, gives identical
  values, vectors and `passed`.
* Results agree to 1e-12 (relative) between 1 and 4 threads, on 2500×2500
  sparse symmetric problems (largest, smallest, block Lanczos), on a sparse
  nonsymmetric problem, and on sparse, centred and crossprod/scale_cols SVD and
  eigen problems.
* `seed =` leaves `.Random.seed` unchanged across nine routes, including
  LOBPCG, randomized, Golub–Kahan, Arnoldi and sparse Lanczos.

A self-test feeds `oracle_check()` doctored results (a wrong set with a
completeness claim, a bad pair, a wrong order, a short result, an
under-reported backward error). It checks that each one trips its invariant
and that a wrong set *without* a claim is reported as soft.

### What the sweep found

Results of running the `extended` level on this machine (R 4.3.3, OpenBLAS).

* Before the fixes: the first 400-case runs found O1–O5. The first extended
  run, made after those were fixed, then flagged 67 cases. Of these, 32 were
  harness classification of clear "unsupported" messages and 2 were the
  harness's own `Matrix()` conversion. The other 33 were real: O12, O13 and
  O14 (24 cases).
* After the fixes: **0 hard violations**. Status counts: 13944 solved, 1056
  expected errors, 0 unexpected errors, 0 crashes, 0 hangs. 12977 of the
  13944 solved cases were certified.

Each violation was minimised into `test-oracle-regressions.R`.

| Id | Finding (sweep case) | Class | Status |
|---|---|---|---|
| O1 | Certificate scale floored at `eps`: the backward error of operators with ‖A‖₂ < eps was reported up to 1e8 times too small. Certified pairs had a true backward error of 1e-5 to 0.67 (cases 17, 88, 244, 284, 342: `crossprod_operator` of 1e-12 data, 1e-20 spectra) | certified but wrong | **fixed**: floor only at `xmin` (`R/certificate_norms.R`, `R/certification.R`) |
| O2 | `passed = TRUE` with fewer than k pairs: shift-invert on {1,1,7,7,7,9,9,9}, k = 8, returned 1, 7, 9 (case 188) | certified but wrong | **fixed**: `withhold_short_certificate()` after every plan dispatch (`R/solve.R`) |
| O3 | An explicit `shift_invert(σ)` with target `smallest()`/`largest()`/`smallest_magnitude()` returned the eigenvalues nearest σ, labelled and certified as the requested target (cases 191, 264, 327, 376, 107, 347) | certified but wrong | **fixed**: such targets are an error, and `eig_partial(method = shift_invert(σ))` defaults to `nearest(σ)` (`R/problem.R`, `R/solve.R`) |
| O4 | A complex general matrix at 1e-12 scale was classified Hermitian. `isSymmetric.matrix()` falls back to an *absolute* `all.equal` when mean\|x\| < tol, so the matrix was solved with zheev and its backward error under-reported (cases 129, 309) | wrong route, under-reported | **fixed**: relative test (`R/operator.R`) |
| O5 | The reference Golub–Kahan absolute breakdown threshold returned no triplets on a 1e-12 diagonal. The empty U/V then reached `dsyrk` with `ldc = 0`, giving `BLAS/LAPACK routine 'DSYRK ' gave error code -10` (a longjmp out of C++) (case 282) | crash-class error | **fixed**: breakdown relative to a running ‖A‖ bound, and `k = 0` guard in `gram_upper_dsyrk_cert` (`R/reference_golub_kahan.R`, `src/certificates.cpp`) |
| O11 | `svds(nu = 0 / nv = 0 / nu = 1)` was never certified and warned "only 0 eigenvalue(s) converged" | API / drop-in | **fixed**: always solve with both sides, then trim (`R/compatibility.R`) |
| O12 | The randomized SVD wide-core Gram path ranked `nearest(σ)` by \|d² − σ\|, returning the wrong set in the wrong order. Its zero-singular-value floor was absolute (cases 502, 5250, …) | wrong set / API order | **fixed** (`R/reference_golub_kahan.R`) |
| O13 | Implicit `smallest_magnitude` on a singular sparse nonsymmetric matrix failed with `native Krylov-Schur Arnoldi operator apply failed with status=-8`: `Matrix::solve()` refused the LU at apply time (case 6743) | unexpected error | **fixed**: singular pivot ratio detected at factor time, so the implicit route perturbs σ (`R/transform_shift_invert.R`) |
| O14 | `eigs_sym()` re-sorted values and vectors but left `certificate$residuals`/`backward_error`/`converged` in solver order | API | **fixed** (`R/compatibility.R`) |
| O15 | Unclear errors: `svd_partial(method = lanczos())` gave "Invalid eigencore plan (dispatch_unavailable): planned_method." (159 cases). `k = n` on Lanczos routes gave "max_subspace must be at least k + 1" when no subspace had been requested (235 cases) | unclear error | **fixed**: messages now name the cause and the remedy (`R/problem.R`, `R/reference_lanczos.R`) |
| O6 | `nearest(σ)` and shift-invert routes miss copies of repeated eigenvalues while certified (case 194) | wrong set, `not_checked` | **fixed**: `nearest`/`smallest_magnitude` are inertia-counted under `completeness = "auto"` when the cost gate passes (proof), else probed on the route's (A − σI)⁻¹ or on (A − σI)² (evidence), with repair (`R/completeness_hermitian.R`) |
| O7 | `both_ends`/`nearest` on the reference Lanczos route miss copies (cases 341, 357, 1544) | wrong set, `not_checked` | **fixed**: every non-exact Hermitian route is checked; `both_ends` is counted (or probed) per end (`R/completeness_hermitian.R`) |
| O8 | SVD routes (prototype/native GK, implicit Gram) miss copies of repeated singular values (case 215) | wrong set, no completeness field | **fixed**: SVD target completeness (`R/completeness_svd.R`): augmented-matrix inertia count or Gram complement probe with repair; matrix-free `nearest()` SVD stays `not_checked` (`passed = FALSE`) |
| O9 | Matrix-free `smallest_magnitude` (`eigs_sym(f, which = "SM", ncv = 20)`) misses a multiple zero eigenvalue (case 133) | wrong set, `not_checked` | **fixed**: small matrix-free operators are materialised and counted; larger ones are probed on A² (`R/completeness_hermitian.R`) |
| O10 | Sparse `both_ends` takes the unrestarted reference Lanczos and does not converge even at n = 40 (honest, but RSpectra solves it) | quality | **fixed**: native both-ends route, two thick-restart solves merged by Rayleigh–Ritz (`R/both_ends_lanczos.R`) |
| O16 | Nonsymmetric Arnoldi certifies pairs outside the `largest/smallest_imaginary` and `largest_magnitude` target set on random matrices (cases 2692, 2783, 4543, 7072, 13198) | wrong set, `not_checked` | **fixed**: deflated Krylov-Schur completeness probe on the complement of the returned Schur basis (`R/completeness_nonsym.R`); all five cases are `repaired` to the correct set |
| O17 | LOBPCG with magnitude targets on indefinite problems misses the other end of the spectrum (cases 8847, 10647, 13154) | wrong set, `not_checked` | **fixed**: the count proves the other end missing and the deflated-complement repair (in the B-transformed space for generalized problems) finds it |
| O18 | Matrix-free nonsymmetric `smallest_magnitude` (no shift-invert available) certifies pairs that are not the smallest (cases 518, 10320, 14913) | wrong set, `not_checked` | **fixed**: the nonsymmetric probe finds the smaller value and repairs the set (`R/completeness_nonsym.R`) |

The sweep also exposed several harness-side pitfalls, which the harness now
avoids:

* `Matrix::Matrix()` and dense→sparse coercions with a tiny scale silently
  symmetrise or drop entries. This is a Matrix-package behaviour that users
  can hit too.
* RSpectra's `SR`/`SM` order is not monotone in the target score.
* `eigs_sym` re-sorting (fixed as O14).

#### Certified results whose set is wrong without a completeness claim

In the final extended run there were 165 such cases out of 12977 certified
(1.3%). 151 of them have repeated, near-repeated, clustered, singular or
indefinite (±λ) spectra, which is the multiplicity class O6–O9. The other 14
are O16–O18. None of them claims completeness: every one reports
`target_completeness = "not_checked"` or is an SVD. The deterministic fix is
the inertia/completeness work planned for tranche 5 (gap 4), and the probe
should be extended to nearest, both_ends, smallest_magnitude, SVD and
nonsymmetric targets.

### Uncertified rate per configuration

These are quality signals from the final extended run (15000 cases). The last
column is the median of max|score error| / ‖A‖ over the uncertified runs.

| Family / method | Runs | Uncertified % | Certified, wrong set (no claim) | Median rel. value error when uncertified |
|---|---|---|---|---|
| gen / lanczos1 | 304 | 16.1 | 7 | 8.7e-05 |
| herm / lanczos1 | 753 | 15.1 | 21 | 3.8e-05 |
| svd / golub_kahan | 802 | 14.5 | 16 | 2.5e-04 |
| gen / lobpcg | 312 | 14.1 | 3 | 5.5e-06 |
| svd / randomized | 611 | 13.9 | 1 | 3.1e-02 |
| herm / lanczos2 | 502 | 12.5 | 16 | 1.6e-07 |
| herm / lobpcg | 561 | 12.5 | 3 | 5.5e-04 |
| herm / lanczos3 | 289 | 10.0 | 7 | 2.5e-05 |
| gen / shift_invert | 69 | 7.2 | 3 | 4.2e-08 |
| gen / auto | 658 | 6.5 | 2 | 3.5e-06 |
| nonsym / shift_invert | 130 | 6.2 | 0 | 8.7e-02 |
| herm / auto | 1729 | 5.2 | 19 | 6.4e-08 |
| nonsym / auto | 1867 | 5.1 | 11 | 1.6e-01 |
| svd / auto | 2079 | 4.6 | 18 | 7.8e-04 |
| nonsym / lanczos1 | 113 | 4.4 | 1 | 2.2e-01 |
| herm / shift_invert | 334 | 3.0 | 15 | 7.0e-09 |
| shim / (eigs, eigs_sym, svds) | 2197 | 2.0 | 22 | 1.5e-08 |
| complex / auto | 634 | 0.2 | 0 | – |

By storage, the highest uncertified rates are on **easy** inputs, which
points to route-quality problems rather than hard problems:

| Family / storage | Runs | Uncertified % |
|---|---|---|
| herm / `ddiMatrix` | 194 | 28.4 |
| svd / `ddiMatrix` | 155 | 28.4 |
| svd / `center(sparse)` | 194 | 17.5 |
| gen / `dgCMatrix` | 373 | 16.6 |
| svd / matrix-free | 420 | 16.4 |
| nonsym / `dgT`/`dgR`/`dgC` | 744 | 14–15 |
| herm / matrix-free | 351 | 14.2 |
| herm / dense | 851 | 4.7 |
| nonsym / dense | 784 | 0.0 |

Observations to follow up (not violations):

* Diagonal input is uncertified about 28% of the time. This is C58: the
  tridiagonal shift-invert route runs out of steps on clustered targets.
* `vectors = FALSE` / `"none"` results are never certified (the certificate
  needs vectors). The shims now work around this; the core API does not.
* Sparse generalized pencils with `auto()` refuse to densify even when `B` is
  a `ddiMatrix` with dense `A`, or the target is `largest` (82 expected
  errors).
* `eigs(which = "SM")` with a function input is often certified-wrong (O18).
  RSpectra usually fails to converge on the same problem instead.
* RSpectra disagreed with the oracle in 160 shim cases where eigencore was
  certified and right. The reverse (RSpectra right, eigencore uncertified)
  happened in 29 cases.

### After the strict `passed` change (extended level, 2026-10-09)

`certificate$passed` now requires the returned set to be verified
(`exact`, `inertia_verified`, `probed` or `repaired`), and every Hermitian,
nonsymmetric and SVD route carries a completeness check. Extended level,
15000 cases, 3 workers, 661 s: **0 hard violations** and **0 certified
results with a wrong set** in every configuration (before: 165 of 12977
certified results had a wrong set without a completeness claim). 12873 of
14022 solved cases certify (91.8%); every remaining uncertified case either
fails its residual certificate or reports an honest unverified status.

| config | runs | certified | uncertified_pct | wrong_set_certified | median_value_err_uncert |
|---|---|---|---|---|---|
| gen / lanczos1 |  304 |  246 | 19.1 | 0 | 1.4e-05 |
| svd / golub_kahan |  802 |  651 | 18.8 | 0 | 2.6e-09 |
| svd / randomized |  611 |  501 | 18.0 | 0 | 1.9e-03 |
| gen / lobpcg |  312 |  268 | 14.1 | 0 | 5.5e-06 |
| herm / lobpcg |  561 |  491 | 12.5 | 0 | 5.5e-04 |
| herm / lanczos1 |  753 |  670 | 11.0 | 0 | 2.4e-07 |
| nonsym / auto | 1867 | 1689 |  9.5 | 0 | 3.2e-05 |
| herm / lanczos2 |  502 |  455 |  9.4 | 0 | 3.5e-09 |
| svd / auto | 2079 | 1895 |  8.9 | 0 | 4.0e-16 |
| herm / lanczos3 |  289 |  266 |  8.0 | 0 | 3.2e-05 |
| nonsym / shift_invert |  130 |  121 |  6.9 | 0 | 3.5e-02 |
| gen / auto |  683 |  640 |  6.3 | 0 | 3.5e-06 |
| nonsym / lanczos1 |  113 |  106 |  6.2 | 0 | 1.6e-01 |
| gen / shift_invert |  122 |  117 |  4.1 | 0 | 4.2e-08 |
| herm / auto | 1729 | 1669 |  3.5 | 0 | 7.2e-12 |
| shim / shim | 2197 | 2131 |  3.0 | 0 | 5.6e-16 |
| herm / shift_invert |  334 |  324 |  3.0 | 0 | 7.0e-09 |
| complex / auto |  634 |  633 |  0.2 | 0 | NA |

### Workflow

`.github/workflows/oracle-sweep.yaml` runs every Monday at 04:17 UTC and on
`workflow_dispatch`, with inputs `level` and `ids`. It uses `ubuntu-latest`
and R release. It installs the package with RSpectra and callr, runs
`tools/oracle_sweep.R` with 2 workers, appends `summary.md` to the job
summary, and uploads `oracle-out/` (records and summary) as an artifact. The
job fails when any hard invariant is violated.

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
| valgrind memcheck, test suite minus the four heaviest files (`t4a-omp`, `bench-smoke`, `t4b-svd`, `t3-api`) | threads 2 | `ERROR SUMMARY: 0 errors`, `definitely lost: 0 bytes`; 6934 expectations, 1 failure that is a valgrind artifact (F7, since fixed) |
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

**F4 (fixed, route quality): `eigs_sym(sparse, which = "BE")`.**
Sparse `both_ends` went to the R reference Hermitian Lanczos ("target
unsupported by native path") and, on a 60 x 60 path Laplacian plus
diagonal, `k = 3`, returned 0 converged pairs with a warning and a failed
certificate (backward errors up to 3e-4); the dense input certified. Root
cause: the reference Lanczos does not restart, so its 3k + 20 step budget
cannot resolve the poorly separated ends. The default shim call was already
fixed on `main` by the native both-ends route ("native Hermitian Lanczos both
ends (two thick-restart solves)", `R/both_ends_lanczos.R`): `eigs_sym()` maps
`"BE"` to `both_ends(k %/% 2, k - k %/% 2)` (alternating ends, the extra
pair from the high end when `k` is odd, values returned in decreasing
order), identical to RSpectra 0.16 for k = 1..6 on the reproducer. One path
still fell back: a warm start (`opts = list(initvec =)`, which the shim
passes as `initial_subspace` with `lanczos()`) excluded the both-ends route
and again converged 0 of 3 pairs. Fix: the both-ends route consumes the
warm start (both thick-restart solves start from it; it is admitted by the
warm-start/restart-state seam via `plan_dispatches_native_warm_lanczos()`),
so every `"BE"` call on sparse, matrix-free or function input takes the
native route and reaches `inertia_verified` / `probed` completeness.
Regression tests: `test-assurance-findings.R` ("F4: ...", compares with
RSpectra when installed).

**F5 (fixed, route quality): explicit `lanczos()` on a matrix-free
Hermitian operator** took the R-level prototype Lanczos ("native hot loop not
yet implemented"), stopped after 29 applies on n = 60 and failed its
certificate; `maxit = 500` and `max_subspace = 30` did not help. Root cause:
the planner sent scalar (`block = 1`) matrix-free `lanczos()` to the
unrestarted reference Lanczos on purpose (only `block > 1` was opted into the
native callback kernel), and that prototype's whole iteration is its
subspace: 3k + 20 = 29 steps by default; the solve-level `maxit` was resolved
as "Lanczos steps" but could only cap that default (`min(29, 500)`), and
`max_subspace = 30` bought exactly one more step. Fix: explicit `lanczos()`
on a matrix-free Hermitian operator with a native target now routes to the
native thick-restart callback kernel as `auto()` does (C53), including warm
starts; the reference route remains only for targets without a native kernel
(e.g. `nearest()` without a factorization), and there `maxit` now sets the
step budget (`min(n, maxit)`) when `max_subspace` is not given (it still only
caps an explicit `max_subspace`). Reproducer after the fix: certified,
`inertia_verified`, 54-80 applies. Regression tests: `test-assurance-findings.R`
("F5: ...").

**F6 (fixed, message quality):** a singular
`shifted_tridiagonal_preconditioner()` (e.g. a path Laplacian with shift 0)
failed LOBPCG with the bare `native CSC LOBPCG failed with status=-5`.
Fix: the preconditioner factors `A + shift * I` once at construction (O(n),
`dgttrf` + `dgtcon`) and a singular system is an error naming
`shifted_tridiagonal_preconditioner()` and asking for a positive shift; the
native LOBPCG maps its status `-5` (only returned by the shifted
diagonal/tridiagonal preconditioner factor and solve) to the same kind of
message. Regression test: `test-assurance-findings.R` ("F6: ...").

**F7 (fixed; found as a valgrind artifact, a portability defect).** Under
valgrind, `test-identity-hash.R` ("C45: identity is identical across separate
R sessions") failed: the sparse operator's identity differed from the one
computed in a native child session. Cause: `as_operator(<dgCMatrix>)` stores
column sums, sums of squares, means and centred sums of squares computed with
`long double` accumulators (`eigencore_csc_column_moments`) in its metadata,
and the built-in identity hashed that metadata; valgrind models x87
`long double` with 64-bit precision, so e.g. the mean `1/3` came out one ulp
different (`0x1.5555555555556p-2` vs `0x1.5555555555555p-2`). The same
happens between platforms (80-bit x86, 128-bit aarch64 Linux, 64-bit on
macOS arm64 / MSVC), so a plan, restart state or PSD factor persisted on one
platform could fail to match the same matrix on another. Fix
(`R/workflow_contract.R`): the identity payload drops the derived metadata
(`derived_identity_metadata_keys()`: the four column moments and
`frobenius_norm`) and hashes the exact content of a Matrix source instead:
its class and every slot (`i`, `p`, `x`, `Dim`, `Dimnames`, `uplo`, ...)
except the `factors` cache, plus the operator's structure flags
(`storage`, `input_storage`, `symmetric_storage`). The native hash streams
the slots without copying, so the cost is unchanged (one pass over the
slots). `identity_hash_format()` is now `"eigencore-identity-hash-v3"`
(the hash algorithm itself is unchanged): plans, restart states and PSD
factors persisted under v2 fail with the typed `identity_format_changed`
error (`eigencore_plan_error`, `eigencore_restart_state_error`) asking for a
re-plan, never with a crash or a silent mismatch. Residual: centred
operators (`center()`) still hash their centring vectors, which are operator
parameters (not derived metadata) computed by `colMeans()`. Regression tests:
`test-assurance-findings.R` ("F7: ...": one-ulp changes of every derived
moment leave the identity unchanged, the `factors` cache is ignored, content
slots and `Dimnames` count, v2 plans and restart states are rejected with the
typed error).

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
