#!/usr/bin/env Rscript

library(eigencore)

set.seed(1)

cases <- list(
  eigen_dense = {
    A <- crossprod(matrix(rnorm(80 * 30), 80, 30))
    list(A = A, k = 5)
  },
  svd_dense = {
    A <- matrix(rnorm(120 * 40), 120, 40)
    list(A = A, rank = 5)
  }
)

cat("Eigen benchmark\n")
print(eigencore:::benchmark_eigen_methods(cases$eigen_dense$A, k = cases$eigen_dense$k, repeats = 3))
print(eigencore:::validate_eigen_accuracy(cases$eigen_dense$A, k = cases$eigen_dense$k))

cat("\nSVD benchmark\n")
print(eigencore:::benchmark_svd_methods(cases$svd_dense$A, rank = cases$svd_dense$rank, repeats = 3))
print(eigencore:::validate_svd_accuracy(cases$svd_dense$A, rank = cases$svd_dense$rank))

cat("\nSparse CSC block-apply smoke\n")
if (requireNamespace("Matrix", quietly = TRUE)) {
  A <- Matrix::rsparsematrix(2000, 1000, density = 0.002)
  X <- matrix(rnorm(1000 * 8), 1000, 8)
  op <- as_operator(A)
  native_time <- system.time(Y_native <- op$apply(X))[["elapsed"]]
  matrix_time <- system.time(Y_matrix <- as.matrix(A %*% X))[["elapsed"]]
  cat("  native CSC:", native_time, "s\n")
  cat("  Matrix %*%:", matrix_time, "s\n")
  cat("  max abs diff:", max(abs(Y_native - Y_matrix)), "\n")
}
