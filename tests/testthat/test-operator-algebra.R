test_that("operator_sum and operator_scale match dense algebra", {
  A <- matrix(c(1, 2, 3, 4), nrow = 2)
  B <- matrix(c(2, 0, 1, 3), nrow = 2)
  X <- matrix(c(1, -1, 2, 0), nrow = 2)

  op <- eigencore:::operator_sum(as_operator(A), eigencore:::operator_scale(B, 2))

  expect_equal(op$apply(X), (A + 2 * B) %*% X)
  expect_equal(op$apply_adjoint(X), t(A + 2 * B) %*% X)
  expect_true(check_adjoint(op, seed = 1)$passed)
})

test_that("compose and crossprod_operator match dense algebra", {
  A <- matrix(rnorm(15), nrow = 5)
  B <- matrix(rnorm(12), nrow = 3)
  X <- matrix(rnorm(8), nrow = 4)

  op <- compose(as_operator(A), as_operator(B))
  expect_equal(op$apply(X), A %*% B %*% X)

  Y <- matrix(rnorm(10), nrow = 5)
  expect_equal(op$apply_adjoint(Y), t(B) %*% t(A) %*% Y)
  expect_true(check_adjoint(op, seed = 2)$passed)

  cp <- crossprod_operator(A)
  Z <- matrix(rnorm(6), nrow = 3)
  expect_equal(cp$apply(Z), t(A) %*% A %*% Z)
  expect_equal(cp$metadata$fused, "crossprod")
  expect_equal(cp$structure$kind, "hermitian")
})

test_that("row and column scaling match dense algebra", {
  A <- matrix(1:12, nrow = 3)
  X <- matrix(rnorm(8), nrow = 4)
  row_weights <- c(2, -1, 0.5)
  col_weights <- c(1, 0.5, -2, 3)

  row_op <- scale_rows(A, row_weights)
  col_op <- scale_cols(A, col_weights)

  expect_equal(row_op$apply(X), (row_weights * A) %*% X)
  expect_equal(col_op$apply(X), sweep(A, 2, col_weights, `*`) %*% X)
  expect_true(check_adjoint(row_op, seed = 3)$passed)
  expect_true(check_adjoint(col_op, seed = 4)$passed)
})

test_that("built-in scaling fuses to native-backed operators without densifying sparse inputs", {
  dense <- matrix(rnorm(20), nrow = 5)
  sparse <- Matrix::rsparsematrix(5, 4, density = 0.35)
  diagonal <- Matrix::Diagonal(x = c(4, -2, 1, 0.5))
  X <- matrix(rnorm(12), nrow = 4)
  Y <- matrix(rnorm(10), nrow = 5)
  row_weights <- c(2, -1, 0.5, 3, -0.25)
  col_weights <- c(1.5, -2, 0.25, 4)

  dense_row <- scale_rows(dense, row_weights)
  dense_col <- scale_cols(dense, col_weights)
  sparse_row <- scale_rows(sparse, row_weights)
  sparse_col <- scale_cols(sparse, col_weights)
  diagonal_scaled <- eigencore:::operator_scale(diagonal, -2)

  expect_true(dense_row$metadata$native)
  expect_true(dense_col$metadata$native)
  expect_true(sparse_row$metadata$native)
  expect_true(sparse_col$metadata$native)
  expect_true(diagonal_scaled$metadata$native)
  expect_equal(dense_row$metadata$fused, "scale_rows")
  expect_equal(dense_col$metadata$fused, "scale_cols")
  expect_equal(sparse_row$metadata$fused, "scale_rows")
  expect_equal(sparse_col$metadata$fused, "scale_cols")
  expect_equal(diagonal_scaled$metadata$fused, "scalar_scale")
  expect_s4_class(sparse_row$metadata$matrix, "dgCMatrix")
  expect_s4_class(sparse_col$metadata$matrix, "dgCMatrix")
  expect_null(eigencore:::source_or_null(sparse_row))
  expect_null(eigencore:::source_or_null(sparse_col))

  expect_equal(dense_row$apply(X), (row_weights * dense) %*% X)
  expect_equal(dense_col$apply(X), sweep(dense, 2L, col_weights, `*`) %*% X)
  expect_equal(sparse_row$apply(X), as.matrix((Matrix::Diagonal(x = row_weights) %*% sparse) %*% X))
  expect_equal(sparse_col$apply(X), as.matrix((sparse %*% Matrix::Diagonal(x = col_weights)) %*% X))
  expect_equal(diagonal_scaled$apply(Y[seq_len(4), , drop = FALSE]),
               as.matrix((-2 * diagonal) %*% Y[seq_len(4), , drop = FALSE]))
  expect_true(check_adjoint(sparse_row, seed = 31)$passed)
  expect_true(check_adjoint(sparse_col, seed = 32)$passed)
})

test_that("column centering matches dense centered matrix", {
  A <- matrix(1:12, nrow = 3)
  X <- matrix(rnorm(8), nrow = 4)
  centered <- sweep(A, 2, colMeans(A), `-`)
  op <- center(A, columns = TRUE)

  expect_equal(op$apply(X), centered %*% X)

  Y <- matrix(rnorm(6), nrow = 3)
  expect_equal(op$apply_adjoint(Y), t(centered) %*% Y)
  expect_true(check_adjoint(op, seed = 5)$passed)
})

test_that("row centering matches dense centered matrix", {
  A <- matrix(1:12, nrow = 3)
  X <- matrix(rnorm(8), nrow = 4)
  centered <- sweep(A, 1, rowMeans(A), `-`)
  op <- center(A, rows = TRUE, columns = FALSE)

  expect_equal(op$apply(X), centered %*% X)
  expect_true(check_adjoint(op, seed = 6)$passed)
})
