#!/usr/bin/env Rscript
#
# bench-parallel-threshold.R -- calibrate getOption("eigencore.parallel_min_work")
# (P26) on this machine.
#
# Usage (against an installed eigencore; run on an otherwise idle machine and
# with the BLAS limited to one thread, as run-suite.R does):
#
#   OPENBLAS_NUM_THREADS=1 Rscript inst/benchmarks/bench-parallel-threshold.R \
#     --lib=.rlib --threads=4
#
# For random sparse square matrices with ~5 nonzeros per row and 1 or 10
# right-hand sides, it times one forward (A X) and one adjoint (A^T X) CSC
# apply with every region forced parallel (eigencore.parallel_min_work = 0)
# and forced serial (Inf), and prints the speedup against the region's work
# estimate (4 units per nonzero and right-hand side, eigencore_common.h).
# The work threshold should sit where the speedup reaches about 1.2-1.5: a
# region just above it must not be slower than serial even when a core is
# briefly busy. The reorthogonalisation kernels call single-threaded BLAS on
# fixed row slabs, so their per-thread speed equals the serial path and the
# same threshold applies to their n * k * b work.

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^--", name, "="), "", hit[[1L]]) else default
}
lib <- get_arg("lib", NULL)
if (!is.null(lib)) .libPaths(c(normalizePath(lib), .libPaths()))
threads <- as.integer(get_arg("threads", "4"))
sizes <- as.integer(strsplit(get_arg("n", "1000,2500,5000,10000,20000,40000,80000"), ",")[[1L]])
suppressPackageStartupMessages({
  library(Matrix)
  library(eigencore)
})
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) RhpcBLASctl::blas_set_num_threads(1)

per_apply_us <- function(A, X, transpose, min_work, nthreads, reps) {
  old <- options(eigencore.threads = nthreads, eigencore.parallel_min_work = min_work)
  on.exit(options(old))
  best <- Inf
  for (k in 1:7) {
    ms <- eigencore:::csc_apply_repeat(A, X, transpose = transpose, reps = reps)$ms
    best <- min(best, ms)
  }
  1000 * best
}

set.seed(20261015)
rows <- list()
for (n in sizes) {
  A <- methods::as(Matrix::rsparsematrix(as.double(n), as.double(n), nnz = 5 * n),
                   "CsparseMatrix")
  for (b in c(1L, 10L)) {
    X <- matrix(stats::rnorm(n * b), n, b)
    work <- 4 * length(A@x) * b
    reps <- as.integer(min(2000, max(10, 4e7 / work)))
    for (transpose in c(FALSE, TRUE)) {
      serial <- per_apply_us(A, X, transpose, Inf, threads, reps)
      parallel <- per_apply_us(A, X, transpose, 0, threads, reps)
      rows[[length(rows) + 1L]] <- data.frame(
        n = n, nnz = length(A@x), rhs = b, op = if (transpose) "A'X" else "AX",
        work = work, serial_us = round(serial, 1), parallel_us = round(parallel, 1),
        speedup = round(serial / parallel, 2))
    }
  }
}
out <- do.call(rbind, rows)
out <- out[order(out$work), ]
print(out, row.names = FALSE)
cat(sprintf("\nthreads = %d, load average %s\n", threads,
            paste(utils::head(strsplit(readLines("/proc/loadavg"), " ")[[1L]], 3), collapse = " ")))
ok <- out$work[out$speedup >= 1.2]
bad <- out$work[out$speedup < 1.0]
cat(sprintf("smallest work with speedup >= 1.2: %s; largest work with speedup < 1: %s\n",
            if (length(ok)) format(min(ok)) else "none",
            if (length(bad)) format(max(bad)) else "none"))
cat(sprintf("current eigencore.parallel_min_work: %s\n",
            format(getOption("eigencore.parallel_min_work", 262144))))
