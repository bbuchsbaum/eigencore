# Benchmark smoke tests

This vignette runs five small eigenvalue and SVD cases with eigencore,
`RSpectra`, `irlba`, and base R. The matrices are intentionally small
and each method is timed three times so the article remains practical to
rebuild. Results describe one render on one machine; they are not
evidence of package rankings or scaling behavior.

``` r

library(eigencore)
```

## What is included?

Dense and sparse inputs exercise different code paths, and SVD has
different comparison methods than Hermitian eigenproblems. The smoke
suite covers these five cases:

``` r

benchmark_regimes <- data.frame(
  regime = c(
    "dense Hermitian",
    "sparse path Laplacian",
    "dense low-rank SVD",
    "tall sparse SVD",
    "wide sparse SVD"
  ),
  input = c("120 x 120 dense", "300 x 300 dgCMatrix", "180 x 70 dense",
            "320 x 60 dgCMatrix", "60 x 320 dgCMatrix"),
  compared_methods = c("eigencore, RSpectra, base",
                       "eigencore, RSpectra, base",
                       "eigencore, RSpectra, irlba, base",
                       "eigencore, RSpectra, irlba, base",
                       "eigencore, RSpectra, irlba, base")
)
knitr::kable(benchmark_regimes)
```

| regime                | input               | compared_methods                 |
|:----------------------|:--------------------|:---------------------------------|
| dense Hermitian       | 120 x 120 dense     | eigencore, RSpectra, base        |
| sparse path Laplacian | 300 x 300 dgCMatrix | eigencore, RSpectra, base        |
| dense low-rank SVD    | 180 x 70 dense      | eigencore, RSpectra, irlba, base |
| tall sparse SVD       | 320 x 60 dgCMatrix  | eigencore, RSpectra, irlba, base |
| wide sparse SVD       | 60 x 320 dgCMatrix  | eigencore, RSpectra, irlba, base |

Base R supplies the dense reference values. For sparse matrices it
converts the input with
[`as.matrix()`](https://rdrr.io/r/base/matrix.html), so its timing and
memory results apply only to these small cases. `irlba` appears only in
the SVD cases.

## What is recorded?

Each row records four things:

- median elapsed time from
  [`bench::mark()`](https://bench.r-lib.org/reference/mark.html);
- allocated memory reported by `bench`;
- relative error against a dense base truth for the requested values;
- an explicit residual/backward-error check on the returned vectors.

For eigencore, the planner label identifies the implementation used in
each case.

## Results from this render

### Elapsed time

The entries are median solver-call milliseconds; lower values are
better. An eigencore solver call includes its native certificate. The
external calls do not, and the shared residual check used elsewhere in
this vignette runs after the timer. These smoke timings therefore
describe call overhead, not a common time-to-certified-answer. The
installed benchmark scripts report both raw solver time and a common
eigencore-certified total. The `lowest_median` column is selected from
unrounded measurements and reports the smallest value in each row.

| case                  | lowest median | eigencore | RSpectra |   irlba | base R |
|:----------------------|:--------------|----------:|---------:|--------:|-------:|
| dense Hermitian       | RSpectra      |    1.6080 |   0.4697 | not run |  2.105 |
| sparse path Laplacian | eigencore     |    2.4110 |   7.9610 | not run | 13.810 |
| dense low-rank SVD    | RSpectra      |    1.5830 |   0.2256 |  0.3615 |  1.699 |
| tall sparse SVD       | RSpectra      |    0.5371 |   0.3700 |  0.7712 |  1.850 |
| wide sparse SVD       | RSpectra      |    0.5463 |   0.4213 |  1.0060 |  3.362 |

Median solver-call time in milliseconds from 3 iterations per method.
{.table style="width:100%;"}

In this render, eigencore recorded the lowest median in **1 of 5**
cases. The result is mixed, and the sub-millisecond rows are especially
sensitive to setup overhead and run-to-run noise.

### Why can larger cases look different?

The smoke cases above are small enough that fixed work matters:
dispatch, result construction, certification, and the dense eigensolve
for a bounded Gram matrix can be a substantial fraction of total time.
They are useful for checking behavior, but they are poor evidence about
scaling.

For a highly rectangular sparse SVD, eigencore can form the smaller Gram
problem and certify the resulting triplets in the original coordinates.
If the smaller dimension stays bounded while the long dimension grows,
the Gram eigensolve retains the same dimension while iterative methods
perform more expensive matrix-vector products. A crossover is therefore
plausible. The installed large endpoint cases show that eigencore can
win in this regime, but those endpoints do not locate a crossover
because their dimensions, density, rank, and random fixture all change
together. The controlled mode of `bench-svd-gram-cutoff.R` fixes the
small side, rank, density, and a nested sparse fixture while varying
only the long side, and reports raw solver and common certified-total
ratios separately.

In the 2026-07-15 installed five-iteration sweep with small side 90, the
certified-total ratio crossed above one by long side 5,000 in both
orientations. The raw-call ratio crossed there for the wide case and by
25,000 for the tall case. These are machine-specific release
measurements, not a universal crossover constant.

This does not mean that increasing either matrix dimension helps. The
current planner limits the smaller dimension to 512 for tall matrices
and 1024 for wide matrices, and also checks aspect ratio, requested-rank
fraction, and a memory budget. The repository’s cutoff evidence is mixed
outside the core regime: increasing the smaller side can remove the
advantage, particularly for tall matrices. The current full endpoint
cases are `100000 x 500` and its transpose; see the [benchmark
manifest](https://github.com/bbuchsbaum/eigencore/blob/main/docs/v1-benchmark-manifest.md)
for the installed evidence and exact commands.

The larger benchmark scripts pass the requested tolerance to RSpectra
and irlba, then recompute eigencore certificates for all external
results. They retain both raw solver time and common time through
certification. A row that still fails the requested certificate is
excluded from time-to-certified-answer claims rather than being treated
as a slow certified result.

### Allocated memory

These values are megabytes allocated during the timed expression. They
do not measure peak resident memory.

| case                  | eigencore | RSpectra |   irlba | base R |
|:----------------------|----------:|---------:|--------:|-------:|
| dense Hermitian       |   1.07900 |  0.06661 | not run | 0.4240 |
| sparse path Laplacian |   1.84800 |  0.01991 | not run | 3.1880 |
| dense low-rank SVD    |   1.31200 |  0.09405 |  0.5606 | 0.5942 |
| tall sparse SVD       |   0.01906 |  0.02502 |  0.1760 | 0.8979 |
| wide sparse SVD       |   0.02361 |  0.02246 |  0.1542 | 0.9098 |

Allocated memory in megabytes. {.table}

### Numerical checks

The summary below reports the largest error observed across the methods
in each case. Timing should be ignored for any case that has a failed
method or residual check.

| case | methods completed | max relative error | max backward error | residual checks |
|:---|:--:|---:|---:|:---|
| dense Hermitian | 3/3 | 2.77e-15 | 1.57e-11 | all pass |
| sparse path Laplacian | 3/3 | 8.33e-15 | 8.45e-10 | all pass |
| dense low-rank SVD | 4/4 | 1.33e-15 | 2.52e-15 | all pass |
| tall sparse SVD | 4/4 | 3.36e-15 | 1.49e-10 | all pass |
| wide sparse SVD | 4/4 | 3.94e-15 | 6.71e-10 | all pass |

Numerical checks across all methods in each case. {.table
style="width:100%;"}

### Eigencore planner paths

| case | planner label | status |
|:---|:---|:---|
| dense Hermitian | native dense Hermitian LAPACK fallback | ok |
| sparse path Laplacian | native tridiagonal Hermitian shift-invert (factorized Lanczos) | ok |
| dense low-rank SVD | native certified Gram SVD special case | ok |
| tall sparse SVD | native certified Gram SVD special case | ok |
| wide sparse SVD | native certified Gram SVD special case | ok |

## Limits of this comparison

The rows above are deliberately small. They keep the vignette runnable,
and they catch obvious problems: wrong target, broken sparse handling,
missing vectors, certification drift, and large timing regressions. They
do not prove large-scale speed claims.

For larger sparse matrices, base R is a truth oracle only while
`as.matrix(A)` is affordable. Studying scaling and memory behavior
requires installed-package runs, larger fixtures, repeated measurements,
and saved artifacts.

## Run the larger benchmark scripts

The scripts under `inst/benchmarks` run larger cases with stricter
checks and saved outputs. Run them from an installed package rather than
a mutable source tree:

``` sh
R CMD INSTALL --library=/tmp/eigencore-bench-lib .
R_LIBS=/tmp/eigencore-bench-lib \
  Rscript inst/benchmarks/bench-native-hermitian-gate.R \
  --strict --save --cases=path_laplacian:1000
```

For the SVD surface, compare eigencore against `RSpectra`, `irlba`, and
base rows where those baselines are meaningful:

``` sh
R_LIBS=/tmp/eigencore-bench-lib \
  Rscript inst/benchmarks/bench-svd-surface.R \
  --iterations=1 --h-candidate \
  --methods=eigencore,RSpectra,irlba,base \
  --cases=tall_sparse,wide_sparse \
  --subject=eigencore --strict --save
```

Those scripts add strict pass/fail gates, larger sparse fixtures, saved
RDS artifacts, memory ratios, and planner-provenance checks. Use
`docs/v1-benchmark-manifest.md` for the current release-gate inventory
and `docs/post-v1-benchmark-gates.md` for future-promotion gates.
