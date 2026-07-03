matrix_free_test_operator <- function(A, apply_adjoint = TRUE, name = "matrix_free") {
  linear_operator(
    dim = dim(A),
    apply = function(X, alpha = 1, beta = 0, Y = NULL) {
      out <- alpha * (A %*% X)
      if (!is.null(Y) && beta != 0) {
        out <- out + beta * Y
      }
      out
    },
    apply_adjoint = if (isTRUE(apply_adjoint)) {
      function(X, alpha = 1, beta = 0, Y = NULL) {
        out <- alpha * (t(A) %*% X)
        if (!is.null(Y) && beta != 0) {
          out <- out + beta * Y
        }
        out
      }
    } else {
      NULL
    },
    name = name
  )
}

test_that("operator_sum and operator_scale match dense algebra", {
  A <- matrix(c(1, 2, 3, 4), nrow = 2)
  B <- matrix(c(2, 0, 1, 3), nrow = 2)
  X <- matrix(c(1, -1, 2, 0), nrow = 2)

  op <- eigencore:::operator_sum(as_operator(A), eigencore:::operator_scale(B, 2))

  expect_equal(op$apply(X), (A + 2 * B) %*% X)
  expect_equal(op$apply_adjoint(X), t(A + 2 * B) %*% X)
  expect_true(check_adjoint(op, seed = 1)$passed)
})

test_that("matrix-free composition sum and scalar scale use fallback algebra", {
  A <- matrix(c(1, 0, -2, 3, 4, 1, 2, -1, 0, 5, -3, 2), nrow = 4)
  B <- matrix(c(2, -1, 0, 3, 1, 4), nrow = 3)
  C <- matrix(c(0.5, 2, -1, 3, 0, 1, -2, 0, 4, -3, 2, 1), nrow = 4)
  X_comp <- matrix(c(1, -2, 0, 3), nrow = 2)
  X_sum <- matrix(seq(-1, 1, length.out = 6), nrow = 3)
  Y4 <- matrix(seq(-0.4, 1.1, length.out = 8), nrow = 4)
  Z4 <- matrix(seq(0.2, 1.7, length.out = 8), nrow = 4)
  W2 <- matrix(seq(-0.5, 0.6, length.out = 4), nrow = 2)

  op_a <- matrix_free_test_operator(A, name = "A")
  op_b <- matrix_free_test_operator(B, name = "B")
  op_c <- matrix_free_test_operator(C, name = "C")

  composed <- compose(op_a, op_b)
  summed <- eigencore:::operator_sum(op_a, op_c)
  scaled <- eigencore:::operator_scale(op_a, -1.5)

  expect_false(isTRUE(composed$metadata$native))
  expect_false(isTRUE(summed$metadata$native))
  expect_false(isTRUE(scaled$metadata$native))
  expect_null(eigencore:::source_or_null(composed))
  expect_equal(composed$apply(X_comp, alpha = 1.25, beta = -0.5, Y = Y4),
               1.25 * A %*% B %*% X_comp - 0.5 * Y4)
  expect_equal(composed$apply_adjoint(Z4, alpha = -2, beta = 0.25, Y = W2),
               -2 * t(B) %*% t(A) %*% Z4 + 0.25 * W2)
  expect_equal(summed$apply(X_sum), (A + C) %*% X_sum)
  expect_equal(summed$apply_adjoint(Z4), t(A + C) %*% Z4)
  expect_equal(scaled$apply(X_sum), (-1.5 * A) %*% X_sum)
  expect_equal(scaled$apply_adjoint(Z4), t(-1.5 * A) %*% Z4)
  expect_true(check_adjoint(composed, seed = 41)$passed)
  expect_true(check_adjoint(summed, seed = 42)$passed)
  expect_true(check_adjoint(scaled, seed = 43)$passed)
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

test_that("matrix-free row and column scaling preserve adjoint algebra", {
  A <- matrix(c(1, -2, 3, 0, 4, 2, -1, 5, 0, 3, 2, -4), nrow = 3)
  op <- matrix_free_test_operator(A)
  X <- matrix(seq(-1, 1, length.out = 8), nrow = 4)
  Y <- matrix(seq(2, -1, length.out = 6), nrow = 3)
  Z <- matrix(seq(-0.5, 0.7, length.out = 6), nrow = 3)
  W <- matrix(seq(1.2, -0.3, length.out = 8), nrow = 4)
  row_weights <- c(2, -1, 0.5)
  col_weights <- c(1, 0.5, -2, 3)

  row_op <- scale_rows(op, row_weights)
  col_op <- scale_cols(op, col_weights)

  expect_false(isTRUE(row_op$metadata$native))
  expect_false(isTRUE(col_op$metadata$native))
  expect_null(eigencore:::source_or_null(row_op))
  expect_null(eigencore:::source_or_null(col_op))
  expect_equal(row_op$apply(X, alpha = -0.75, beta = 0.5, Y = Y),
               -0.75 * (row_weights * A) %*% X + 0.5 * Y)
  expect_equal(row_op$apply_adjoint(Z), t(row_weights * A) %*% Z)
  expect_equal(col_op$apply(X), sweep(A, 2, col_weights, `*`) %*% X)
  expect_equal(col_op$apply_adjoint(Z, alpha = 1.5, beta = -0.25, Y = W),
               1.5 * t(sweep(A, 2, col_weights, `*`)) %*% Z - 0.25 * W)
  expect_true(check_adjoint(row_op, seed = 44)$passed)
  expect_true(check_adjoint(col_op, seed = 45)$passed)
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

test_that("built-in sums and compositions fuse to native explicit operators", {
  dense_a <- matrix(rnorm(12), nrow = 3)
  dense_b <- matrix(rnorm(12), nrow = 4)
  sparse_a <- Matrix::rsparsematrix(5, 4, density = 0.35)
  sparse_b <- Matrix::rsparsematrix(4, 3, density = 0.35)
  sparse_c <- Matrix::rsparsematrix(5, 4, density = 0.25)
  X_dense <- matrix(rnorm(6), nrow = 3)
  X_sparse <- matrix(rnorm(6), nrow = 3)
  X_sum <- matrix(rnorm(8), nrow = 4)

  dense_comp <- compose(as_operator(dense_a), as_operator(dense_b))
  sparse_comp <- compose(as_operator(sparse_a), as_operator(sparse_b))
  sparse_sum <- eigencore:::operator_sum(as_operator(sparse_a), as_operator(sparse_c))

  expect_true(dense_comp$metadata$native)
  expect_true(sparse_comp$metadata$native)
  expect_true(sparse_sum$metadata$native)
  expect_equal(dense_comp$metadata$fused, "compose")
  expect_equal(sparse_comp$metadata$fused, "compose")
  expect_equal(sparse_sum$metadata$fused, "sum")
  expect_s4_class(sparse_comp$metadata$matrix, "dgCMatrix")
  expect_s4_class(sparse_sum$metadata$matrix, "dgCMatrix")
  expect_null(eigencore:::source_or_null(sparse_comp))
  expect_null(eigencore:::source_or_null(sparse_sum))

  expect_equal(dense_comp$apply(X_dense), dense_a %*% dense_b %*% X_dense)
  expect_equal(sparse_comp$apply(X_sparse), as.matrix((sparse_a %*% sparse_b) %*% X_sparse))
  expect_equal(sparse_sum$apply(X_sum), as.matrix((sparse_a + sparse_c) %*% X_sum))
  expect_true(check_adjoint(sparse_comp, seed = 33)$passed)
  expect_true(check_adjoint(sparse_sum, seed = 34)$passed)
})

test_that("column centering matches dense centered matrix", {
  A <- matrix(1:12, nrow = 3)
  X <- matrix(rnorm(8), nrow = 4)
  centered <- sweep(A, 2, colMeans(A), `-`)
  op <- center(A, columns = TRUE)

  expect_equal(op$apply(X), centered %*% X)
  expect_true(op$metadata$native)
  expect_equal(op$metadata$fused, "center")

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
  expect_true(op$metadata$native)
  expect_equal(op$metadata$fused, "center")
  expect_true(check_adjoint(op, seed = 6)$passed)
})

test_that("matrix-free centering requires explicit means and matches centered algebra", {
  A <- matrix(c(1, 3, 5, 2, 4, 6, -1, 0, 2, 7, 8, 9), nrow = 3)
  op <- matrix_free_test_operator(A)
  row_means <- rowMeans(A)
  col_means <- colMeans(A)
  centered <- sweep(sweep(A, 2L, col_means, `-`), 1L, row_means, `-`)
  X <- matrix(seq(-1, 1, length.out = 8), nrow = 4)
  Y <- matrix(seq(0.5, -0.5, length.out = 6), nrow = 3)
  Z <- matrix(seq(-2, 2, length.out = 6), nrow = 3)
  W <- matrix(seq(1, -1, length.out = 8), nrow = 4)

  expect_error(center(op, columns = TRUE), "col_means must be supplied")
  expect_error(center(op, rows = TRUE, columns = FALSE), "row_means must be supplied")
  expect_error(center(op, columns = TRUE, col_means = col_means[-1]),
               "column dimension")
  expect_error(center(op, rows = TRUE, columns = FALSE, row_means = row_means[-1]),
               "row dimension")

  centered_op <- center(
    op,
    rows = TRUE,
    columns = TRUE,
    row_means = row_means,
    col_means = col_means
  )

  expect_false(isTRUE(centered_op$metadata$native))
  expect_null(eigencore:::source_or_null(centered_op))
  expect_equal(centered_op$apply(X, alpha = 1.5, beta = -0.25, Y = Y),
               1.5 * centered %*% X - 0.25 * Y)
  expect_equal(centered_op$apply_adjoint(Z, alpha = -0.5, beta = 0.75, Y = W),
               -0.5 * t(centered) %*% Z + 0.75 * W)
  expect_true(check_adjoint(centered_op, seed = 46)$passed)
})

test_that("sparse centering uses native low-rank correction without densifying", {
  A <- Matrix::rsparsematrix(7, 5, density = 0.35)
  X <- matrix(rnorm(15), nrow = 5)
  Y <- matrix(rnorm(21), nrow = 7)
  Z <- matrix(rnorm(21), nrow = 7)

  op <- center(A, rows = TRUE, columns = TRUE)
  expected <- sweep(sweep(as.matrix(A), 2L, Matrix::colMeans(A), `-`),
                    1L, Matrix::rowMeans(A), `-`)

  expect_true(op$metadata$native)
  expect_equal(op$metadata$fused, "center")
  expect_equal(op$metadata$storage, "centered_dgCMatrix")
  expect_true(op$metadata$low_rank_correction)
  expect_s4_class(op$metadata$base_matrix, "dgCMatrix")
  expect_null(op$metadata$matrix)
  expect_null(eigencore:::source_or_null(op))
  expect_false(eigencore:::has_native_kernel(op))

  expect_equal(op$apply(X), expected %*% X)
  expect_equal(op$apply(X, alpha = 1.5, beta = -0.25, Y = Y),
               1.5 * expected %*% X - 0.25 * Y)
  expect_equal(op$apply_adjoint(Z), t(expected) %*% Z)
  expect_true(check_adjoint(op, seed = 37)$passed)
})

test_that("crossprod fuses native explicit dense and sparse operators", {
  dense <- matrix(rnorm(20), nrow = 5)
  sparse <- Matrix::rsparsematrix(6, 4, density = 0.35)
  X <- matrix(rnorm(12), nrow = 4)

  dense_cp <- crossprod_operator(dense)
  sparse_cp <- crossprod_operator(sparse)

  expect_true(dense_cp$metadata$native)
  expect_true(sparse_cp$metadata$native)
  expect_equal(dense_cp$metadata$fused, "crossprod")
  expect_equal(sparse_cp$metadata$fused, "crossprod")
  expect_true(dense_cp$metadata$materialized_crossprod)
  expect_true(sparse_cp$metadata$materialized_crossprod)
  expect_s4_class(sparse_cp$metadata$matrix, "dgCMatrix")
  expect_equal(dense_cp$structure$kind, "hermitian")
  expect_equal(sparse_cp$structure$kind, "hermitian")

  expect_equal(dense_cp$apply(X), crossprod(dense) %*% X)
  expect_equal(sparse_cp$apply(X), as.matrix(Matrix::crossprod(sparse) %*% X))
  expect_true(check_adjoint(dense_cp, seed = 35)$passed)
  expect_true(check_adjoint(sparse_cp, seed = 36)$passed)
})

test_that("operator algebra rejects incompatible or unsupported operator contracts", {
  A <- matrix_free_test_operator(matrix(1:6, nrow = 3), name = "A")
  B <- matrix_free_test_operator(matrix(1:8, nrow = 4), name = "B")
  no_adjoint <- matrix_free_test_operator(diag(2), apply_adjoint = FALSE)
  bad_adjoint <- linear_operator(
    dim = c(2, 2),
    apply = function(X, alpha = 1, beta = 0, Y = NULL) alpha * X,
    apply_adjoint = function(X, alpha = 1, beta = 0, Y = NULL) alpha * (2 * X)
  )

  expect_error(compose(A, B), "Cannot compose")
  expect_error(eigencore:::operator_sum(), "At least one operator")
  expect_error(eigencore:::operator_sum(A, B), "same dimensions")
  expect_error(eigencore:::operator_scale(A, Inf), "finite numeric")
  expect_error(scale_rows(A, c(1, 2)), "row weights length")
  expect_error(scale_cols(A, c(1, 2, 3)), "column weights length")
  expect_error(crossprod_operator(no_adjoint), "apply_adjoint")
  expect_error(check_adjoint(no_adjoint), "apply_adjoint")
  expect_error(check_adjoint(bad_adjoint, trials = 2, seed = 47),
               "Adjoint check failed")
  expect_error(symmetric_operator(A), "must be square")
  expect_error(symmetric_operator(bad_adjoint), "Adjoint check failed")

  unchecked <- symmetric_operator(no_adjoint, validate = FALSE)
  expect_equal(unchecked$structure$kind, "hermitian")
  expect_match(unchecked$name, "^symmetric")
})
