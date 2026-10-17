# Multithreading in eigencore's sparse kernels

eigencore parallelises its native sparse (CSC, `dgCMatrix`) operator
applications with OpenMP: `A^T X` products (a gather per column), `A X`
products (a gather over a cached CSR copy, or a scatter over cached
per-thread row slabs for tall matrices with short rows), the centred and
centred-scaled sparse operators, and the implicit normal operators
`A^T A` / `A A^T` behind sparse
[`svds()`](https://bbuchsbaum.github.io/eigencore/reference/svds.md) /
[`svd_partial()`](https://bbuchsbaum.github.io/eigencore/reference/svd_partial.md).
Every output entry is accumulated by one thread in the same order and
with the same arithmetic as the serial loop, so sparse products are
bitwise identical for any thread count. While such a solve runs, the
Lanczos reorthogonalisation also uses these threads (see below):
single-threaded BLAS calls on fixed blocks of rows, shared among the
threads, whose results do not depend on the thread count either, so
results computed with one and with several threads differ only by
rounding, and results with 2, 4 or 8 threads agree exactly.

## Options

- `eigencore.threads`:

  Number of threads, read at the start of every native call. When unset,
  a default fixed at package load is used: `1` under `R CMD check`
  (`_R_CHECK_PACKAGE_NAME_` set) or when `_R_CHECK_LIMIT_CORES_` is set
  (CRAN's core limit); otherwise the first value of `OMP_NUM_THREADS`
  when set; otherwise the number of processors available to OpenMP
  (respecting the CPU affinity mask), further limited to the physical
  cores (`parallel::detectCores(logical = FALSE)`) and to the cgroup CPU
  quota of a container, and capped at 8. Builds without OpenMP (for
  example Apple clang without `libomp`) always use one thread.

- `eigencore.adaptive_threads`:

  `TRUE` (default) or `FALSE`. While enabled, eigencore watches the
  parallel efficiency of its own OpenMP regions and lowers the number of
  threads it actually uses when the machine is oversubscribed (see
  below). `FALSE` always uses `eigencore.threads` threads, as in
  eigencore 1.3.0.

- `eigencore.csr_cache_mb`:

  Memory cap in MB (default 4096) for the per-operator CSR copy or row
  slabs (about 12 bytes per nonzero) that the parallel `A X` builds once
  per solve. When the copy would be larger, `A X` stays serial for
  blocks of up to 10 columns.

- `eigencore.parallel_min_work`:

  Minimum work per thread of a parallel region (default 65536), in
  multiply-add units: a sparse product counts 4 units per nonzero and
  right-hand side, the reorthogonalisation 1 per basis entry and column.
  A region gets at most `floor(work / eigencore.parallel_min_work)`
  threads and runs serially below two, so small products (a 20000 x
  20000 matrix with 3 nonzeros per row, about 60 microseconds) use few
  threads or one: they cannot gain from more and, when any core is busy,
  every join waits for the slowest team member. Like the adaptive thread
  count it only changes how many threads run a region, never the result.
  `0` lets every region use all threads; `Inf` keeps eigencore's own
  kernels serial. `inst/benchmarks/bench-parallel-threshold.R` measures
  the crossover on a machine.

## Busy machines

OpenMP runtimes keep idle team threads spinning for a while (libgomp's
default `OMP_WAIT_POLICY`), which is fast on idle cores. When other
processes occupy the cores, every fork, join and barrier of a parallel
region instead waits a scheduler time slice for a team member that is
not running, and eigencore 1.3.0's short sparse kernels became 10-80
times slower with 4 threads than with 1. eigencore therefore times each
parallel region: when a region takes longer than the master thread's
share of the work times the team size (i.e. longer than running it
serially would have), twice in a row, the team is halved, and a larger
team is tried again after a pause that grows (0.25 s to 4 s) while the
machine stays busy. Only the number of threads in a region changes; the
algorithm and the work partition that determine rounding are chosen from
`eigencore.threads`, so results do not depend on the machine load. The
wait policy itself cannot be changed from a package (R loads the OpenMP
runtime at start-up); on machines that are routinely oversubscribed,
setting `OMP_WAIT_POLICY=passive` in the environment before starting R
also removes the spinning. A multithreaded OpenBLAS (pthreads build)
spins the same way in its own dense calls (including those made between
eigencore's native calls and by CHOLMOD); eigencore does not change it
on its own because OpenBLAS results depend on its thread count, so on a
busy machine also limit it (`OPENBLAS_NUM_THREADS=1`, or
`RhpcBLASctl::blas_set_num_threads(1)`).

## Interaction with a multithreaded BLAS

OpenBLAS (pthreads build) and FlexiBLAS keep idle worker threads
spinning after each call; OpenMP threads started meanwhile would compete
with them for cores. When `eigencore.threads` is above one, a native
sparse solve (for example
[`eigs_sym()`](https://bbuchsbaum.github.io/eigencore/reference/eigs_sym.md),
[`svds()`](https://bbuchsbaum.github.io/eigencore/reference/svds.md) or
[`svd_partial()`](https://bbuchsbaum.github.io/eigencore/reference/svd_partial.md)
on a `dgCMatrix`) therefore switches such a BLAS to one thread at its
first sparse product (also when that product is too small to run in
parallel) for the rest of that call and restores the previous setting
when the call returns (or fails); its reorthogonalisation then runs on
eigencore's threads. Sparse products applied one at a time from R-level
operators (for example inside a matrix-free
[`center()`](https://bbuchsbaum.github.io/eigencore/reference/center.md)
solve) stay serial while such a BLAS is multithreaded, so the
surrounding BLAS work keeps its threads. An OpenMP-built OpenBLAS shares
eigencore's thread pool; MKL, BLIS and other libraries are not changed.
With those, or when running several R workers in parallel, consider
`options(eigencore.threads = 1)` or limiting the BLAS threads (for
example with `RhpcBLASctl::blas_set_num_threads()`).

## Examples

``` r
old <- options(eigencore.threads = 1)
getOption("eigencore.threads")
#> [1] 1
options(old)
```
