# Test assurance: oracle-based differential sweep

Date: 2026-10-09. Branch `assure-oracle` (from `test/assurance` at `294eee5`).

The per-bug regression tests in `tests/testthat` pin the defects found by the
review (`docs/review-2026-10.md`). They do not say much about inputs nobody
thought of. The oracle sweep fills that gap. It generates problems across the
input space, solves each one with eigencore, and checks the result against a
dense oracle (base `eigen()`/`svd()`, or a Cholesky reduction for SPD
pencils). The property it checks hardest is **"certified ⇒ actually
correct"**.

## Files

| File | Role |
|---|---|
| `tests/testthat/helper-oracle.R` | Case generator, oracle, invariant checks, subprocess runner, summary tables |
| `tests/testthat/test-oracle-sweep.R` | Runs the level's cases. Fails on any hard violation and prints a one-line reproducer per case. Also has harness self-tests (a doctored result must trip each class of check) |
| `tests/testthat/test-oracle-determinism.R` | Same seed gives identical results. `eigencore.threads` 1 vs 4 agree to 1e-12. `seed =` leaves the global RNG alone |
| `tests/testthat/test-oracle-shims.R` | `eigs_sym`/`eigs`/`svds` against RSpectra on well-separated spectra: same values in RSpectra's order |
| `tests/testthat/test-oracle-regressions.R` | Minimised reproducers for every violation found (O1–O18). Known, unfixed items are `skip("known: O<n>")` |
| `tools/oracle_sweep.R` | Stand-alone parallel runner. Writes `records.rds`/`.csv` and `summary.md` |
| `.github/workflows/oracle-sweep.yaml` | Runs the `extended` level weekly and on dispatch, then uploads the results |

## Levels

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

## Generator coverage

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

## Invariants

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

## What the sweep found

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

### Certified results whose set is wrong without a completeness claim

In the final extended run there were 165 such cases out of 12977 certified
(1.3%). 151 of them have repeated, near-repeated, clustered, singular or
indefinite (±λ) spectra, which is the multiplicity class O6–O9. The other 14
are O16–O18. None of them claims completeness: every one reports
`target_completeness = "not_checked"` or is an SVD. The deterministic fix is
the inertia/completeness work planned for tranche 5 (gap 4), and the probe
should be extended to nearest, both_ends, smallest_magnitude, SVD and
nonsymmetric targets.

Since `certificate$passed` requires a verified set (`require_verified_completeness()`
in `R/solve.R`), every Hermitian route is checked by
`R/completeness_hermitian.R`: inertia counts (a proof; ties at the target edge
are verified as ties, see `hermitian_completeness_tie()`), materialised counts
for small matrix-free operators, and otherwise the deflated-complement probe
(evidence) in the B-transformed standard space, with repair. A Hermitian
result that cannot be checked reports `passed = FALSE` with
`residual_passed = TRUE`.

## Uncertified rate per configuration

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

## After the strict `passed` change (extended level, 2026-10-09)

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

## Workflow

`.github/workflows/oracle-sweep.yaml` runs every Monday at 04:17 UTC and on
`workflow_dispatch`, with inputs `level` and `ids`. It uses `ubuntu-latest`
and R release. It installs the package with RSpectra and callr, runs
`tools/oracle_sweep.R` with 2 workers, appends `summary.md` to the job
summary, and uploads `oracle-out/` (records and summary) as an artifact. The
job fails when any hard invariant is violated.
