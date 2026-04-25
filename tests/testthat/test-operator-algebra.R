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
