test_that("dense operator block apply matches matrix multiplication", {
  A <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 3)
  X <- matrix(c(1, 0, 2, 1), nrow = 2)
  op <- as_operator(A)

  expect_equal(op$apply(X), A %*% X)
  expect_equal(op$apply_adjoint(matrix(1, 3, 2)), t(A) %*% matrix(1, 3, 2))
})

test_that("dense operator block apply honors alpha beta Y contract", {
  A <- matrix(rnorm(12), nrow = 3)
  X <- matrix(rnorm(8), nrow = 4)
  Y <- matrix(rnorm(6), nrow = 3)
  op <- as_operator(A)

  expect_equal(op$apply(X, alpha = 2, beta = -0.5, Y = Y), 2 * A %*% X - 0.5 * Y)
})

test_that("adjoint operator swaps apply directions", {
  A <- matrix(rnorm(12), nrow = 3)
  op <- adjoint(as_operator(A))
  X <- matrix(rnorm(6), nrow = 3)

  expect_equal(op$apply(X), t(A) %*% X)
})

test_that("adjoint operator preserves native metadata for built-in kernels", {
  dense <- matrix(rnorm(12), nrow = 3)
  sparse <- Matrix::rsparsematrix(5, 3, density = 0.4)
  diagonal <- Matrix::Diagonal(x = c(3, 2, 1))

  dense_adj <- adjoint(as_operator(dense))
  sparse_adj <- adjoint(as_operator(sparse))
  diagonal_adj <- adjoint(as_operator(diagonal))

  expect_true(dense_adj$metadata$native)
  expect_true(sparse_adj$metadata$native)
  expect_true(diagonal_adj$metadata$native)
  expect_equal(dense_adj$metadata$fused, "adjoint")
  expect_equal(sparse_adj$metadata$fused, "adjoint")
  expect_equal(diagonal_adj$metadata$fused, "adjoint")
  expect_equal(dense_adj$metadata$source, t(dense))
  expect_s4_class(sparse_adj$metadata$matrix, "dgCMatrix")
  expect_equal(as.matrix(sparse_adj$metadata$matrix), as.matrix(Matrix::t(sparse)))
  expect_equal(diagonal_adj$metadata$storage, "adjoint:ddiMatrix")

  X_dense <- matrix(rnorm(9), nrow = 3)
  X_sparse <- matrix(rnorm(10), nrow = 5)
  X_diag <- matrix(rnorm(6), nrow = 3)
  expect_equal(dense_adj$apply(X_dense), t(dense) %*% X_dense)
  expect_equal(sparse_adj$apply(X_sparse), as.matrix(Matrix::t(sparse) %*% X_sparse))
  expect_equal(diagonal_adj$apply(X_diag), as.matrix(Matrix::t(diagonal) %*% X_diag))
})

test_that("dgCMatrix operator uses native CSC block apply", {
  A <- Matrix::sparseMatrix(
    i = c(1, 3, 2, 4),
    j = c(1, 1, 3, 4),
    x = c(2, -1, 3, 4),
    dims = c(4, 4)
  )
  X <- matrix(rnorm(12), nrow = 4)
  Y <- matrix(rnorm(12), nrow = 4)
  op <- as_operator(A)

  expect_equal(op$name, "sparse_csc_matrix")
  expect_true(op$metadata$native)
  expect_equal(op$metadata$storage, "dgCMatrix")
  expect_null(eigencore:::source_or_null(op))
  expect_equal(op$apply(X), as.matrix(A %*% X))
  expect_equal(op$apply_adjoint(X), as.matrix(Matrix::t(A) %*% X))
  expect_equal(op$apply(X, alpha = 2, beta = -0.25, Y = Y), 2 * as.matrix(A %*% X) - 0.25 * Y)
})

test_that("dgCMatrix norm metadata supports deterministic operator certificates", {
  A <- Matrix::sparseMatrix(i = 1:3, j = 1:3, x = c(5, 3, 1), dims = c(3, 3))
  op <- as_operator(A)
  v <- diag(3)[, 1:2]
  vals <- c(5, 3)
  cert <- eigencore:::certify_eigen_operator(op, vals, v)

  expect_equal(cert$scale, eigencore:::eigen_backward_scale(sqrt(35), 1, vals, v))
  expect_true(cert$passed)
})

test_that("ddiMatrix operator uses native diagonal block apply", {
  D <- Matrix::Diagonal(x = c(4, -2, 1))
  X <- matrix(rnorm(6), nrow = 3)
  Y <- matrix(rnorm(6), nrow = 3)
  op <- as_operator(D)

  expect_equal(op$name, "diagonal_matrix")
  expect_true(op$metadata$native)
  expect_equal(op$metadata$storage, "ddiMatrix")
  expect_null(eigencore:::source_or_null(op))
  expect_equal(op$apply(X), as.matrix(D %*% X))
  expect_equal(op$apply_adjoint(X), as.matrix(Matrix::t(D) %*% X))
  expect_equal(op$apply(X, alpha = -1.5, beta = 0.25, Y = Y), -1.5 * as.matrix(D %*% X) + 0.25 * Y)
})

test_that("unit ddiMatrix operator uses native diagonal block apply", {
  D <- Matrix::Diagonal(3)
  X <- matrix(rnorm(6), nrow = 3)
  op <- as_operator(D)

  expect_equal(op$apply(X), X)
  expect_equal(op$metadata$frobenius_norm, sqrt(3))
  expect_equal(op$metadata$two_norm_upper, 1)
})

test_that("native dense CSC and diagonal operators pass adjoint checks", {
  dense <- as_operator(matrix(rnorm(20), nrow = 5))
  sparse <- as_operator(Matrix::rsparsematrix(6, 4, density = 0.4))
  diagonal <- as_operator(Matrix::Diagonal(x = c(3, 2, 1, 0.5)))

  expect_true(check_adjoint(dense, seed = 10)$passed)
  expect_true(check_adjoint(sparse, seed = 11)$passed)
  expect_true(check_adjoint(diagonal, seed = 12)$passed)
})

test_that("native built-in apply functions report zero workspace allocations", {
  dense <- matrix(rnorm(20), nrow = 5)
  sparse <- Matrix::rsparsematrix(5, 4, density = 0.4)
  diagonal <- Matrix::Diagonal(x = c(3, 2, 1, 0.5))

  dense_counters <- eigencore:::native_apply_noalloc_check("dense", dense, matrix(rnorm(8), nrow = 4))
  sparse_counters <- eigencore:::native_apply_noalloc_check("csc", sparse, matrix(rnorm(8), nrow = 4))
  diagonal_counters <- eigencore:::native_apply_noalloc_check("diagonal", diagonal, matrix(rnorm(8), nrow = 4))

  expect_equal(unname(dense_counters), c(0L, 0L))
  expect_equal(unname(sparse_counters), c(0L, 0L))
  expect_equal(unname(diagonal_counters), c(0L, 0L))
})
