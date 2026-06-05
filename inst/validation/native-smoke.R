#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
use_load_all <- "--load-all" %in% args
quiet <- "--quiet" %in% args

if (use_load_all) {
  if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("--load-all requires the pkgload package", call. = FALSE)
  }
  pkgload::load_all(".", quiet = TRUE)
} else {
  suppressPackageStartupMessages(library(eigencore))
}
suppressPackageStartupMessages(requireNamespace("Matrix"))

fail <- function(name, message) {
  stop(sprintf("%s failed: %s", name, message), call. = FALSE)
}

assert_certificate <- function(name, fit, passed = TRUE, max_backward_error = 1e-7) {
  cert <- eigencore::certificate(fit)
  if (!isTRUE(all(cert$converged))) {
    fail(name, "not all requested vectors converged")
  }
  if (isTRUE(passed) && !isTRUE(cert$passed)) {
    fail(name, "certificate did not pass")
  }
  if (!isTRUE(passed) && isTRUE(cert$passed)) {
    fail(name, "estimated-scale certificate unexpectedly passed")
  }
  if (!is.finite(cert$max_backward_error) ||
      cert$max_backward_error > max_backward_error) {
    fail(
      name,
      sprintf(
        "max backward error %.3g exceeds %.3g",
        cert$max_backward_error,
        max_backward_error
      )
    )
  }
  invisible(cert)
}

path_laplacian <- function(n) {
  Matrix::bandSparse(
    n,
    k = c(-1, 0, 1),
    diagonals = list(
      rep(-1, n - 1L),
      c(1, rep(2, n - 2L), 1),
      rep(-1, n - 1L)
    )
  )
}

dense_hermitian <- eigencore::eig_partial(
  diag(c(6, 3, 1)),
  k = 2L,
  target = eigencore::largest()
)
assert_certificate("dense Hermitian eigen", dense_hermitian)

sparse_hermitian <- eigencore::eig_partial(
  path_laplacian(35L),
  k = 3L,
  target = eigencore::smallest(),
  seed = 101
)
assert_certificate("sparse CSC Hermitian Lanczos", sparse_hermitian)

A_gen <- diag(c(2, 5, 9, 14))
B_gen <- diag(c(1, 2, 3, 4))
generalized <- eigencore::eig_partial(
  A_gen,
  B = B_gen,
  k = 2L,
  target = eigencore::smallest(),
  method = eigencore::lobpcg(maxit = 80L),
  seed = 102
)
assert_certificate("dense generalized SPD LOBPCG", generalized)

dense_svd <- eigencore::svd_partial(
  diag(c(8, 4, 2, 1)),
  rank = 2L,
  target = eigencore::largest()
)
assert_certificate("dense SVD", dense_svd)

sparse_svd_matrix <- Matrix::sparseMatrix(
  i = c(1L, 2L, 3L, 4L, 5L, 1L, 2L, 3L),
  j = c(1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L),
  x = c(9, 7, 5, 3, 1, 0.5, 0.4, 0.3),
  dims = c(20L, 8L)
)
sparse_svd <- eigencore::svd_partial(
  sparse_svd_matrix,
  rank = 2L,
  target = eigencore::largest(),
  seed = 103
)
assert_certificate("sparse CSC SVD", sparse_svd)

shifted <- eigencore::eig_partial(
  diag(c(1, 3, 7)),
  k = 1L,
  target = eigencore::nearest(2.8),
  method = eigencore::shift_invert(2.8)
)
assert_certificate("dense shift-invert", shifted)

tridiagonal_shifted <- eigencore::eig_partial(
  path_laplacian(35L),
  k = 2L,
  target = eigencore::nearest(0.02),
  method = eigencore::shift_invert(0.02),
  seed = 104,
  allow_dense_fallback = "never"
)
assert_certificate("native tridiagonal shift-invert", tridiagonal_shifted)

if (!quiet) {
  cat("eigencore native smoke passed\n")
}
