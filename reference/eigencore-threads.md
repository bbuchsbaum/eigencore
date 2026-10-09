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
Lanczos reorthogonalisation also uses these threads (see below); its
OpenMP kernels do not depend on the thread count either, so results
computed with one and with several threads differ only by rounding.

## Options

- `eigencore.threads`:

  Number of threads, read at the start of every native call. When unset,
  a default fixed at package load is used: `1` under `R CMD check`
  (`_R_CHECK_PACKAGE_NAME_` set) or when `_R_CHECK_LIMIT_CORES_` is set
  (CRAN's core limit); otherwise the first value of `OMP_NUM_THREADS`
  when set; otherwise the number of processors reported by OpenMP,
  capped at 8. Builds without OpenMP (for example Apple clang without
  `libomp`) always use one thread.

- `eigencore.csr_cache_mb`:

  Memory cap in MB (default 4096) for the per-operator CSR copy or row
  slabs (about 12 bytes per nonzero) that the parallel `A X` builds once
  per solve. When the copy would be larger, `A X` stays serial for
  blocks of up to 10 columns.

## Interaction with a multithreaded BLAS

OpenBLAS (pthreads build) and FlexiBLAS keep idle worker threads
spinning after each call; OpenMP threads started meanwhile would compete
with them for cores. When a native sparse solve (for example
[`eigs_sym()`](https://bbuchsbaum.github.io/eigencore/reference/eigs_sym.md),
[`svds()`](https://bbuchsbaum.github.io/eigencore/reference/svds.md) or
[`svd_partial()`](https://bbuchsbaum.github.io/eigencore/reference/svd_partial.md)
on a `dgCMatrix`) runs a multithreaded sparse kernel, eigencore
therefore switches such a BLAS to one thread for the rest of that call
and restores the previous setting when the call returns (or fails); its
reorthogonalisation then runs on eigencore's OpenMP threads. Sparse
products applied one at a time from R-level operators (for example
inside a matrix-free
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
