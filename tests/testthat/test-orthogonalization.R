test_that("mgs2 produces an orthonormal basis", {
  set.seed(10)
  X <- matrix(rnorm(30), nrow = 10)
  out <- eigencore:::mgs2(X)

  expect_equal(crossprod(out$Q), diag(ncol(out$Q)), tolerance = 1e-12)
  expect_lt(out$orthogonality, 1e-12)
  expect_equal(out$Q %*% out$R, X, tolerance = 1e-10)
})

test_that("mgs2 orthogonalizes against an existing basis", {
  set.seed(11)
  Q0 <- qr.Q(qr(matrix(rnorm(40), nrow = 10)))[, 1:2]
  X <- matrix(rnorm(30), nrow = 10)
  out <- eigencore:::mgs2(X, against = Q0)

  expect_lt(max(abs(crossprod(Q0, out$Q))), 1e-12)
  expect_lt(eigencore:::orthogonality_loss(out$Q), 1e-12)
})

test_that("native reorthogonalize_against matches BLAS algebra contract", {
  set.seed(110)
  Q <- qr.Q(qr(matrix(rnorm(40), nrow = 10)))[, 1:3]
  X <- matrix(rnorm(50), nrow = 10)

  once <- eigencore:::reorthogonalize_against(X, Q, passes = 1L)
  twice <- eigencore:::reorthogonalize_against(X, Q, passes = 2L)
  direct_once <- X - Q %*% crossprod(Q, X)
  direct_twice <- direct_once - Q %*% crossprod(Q, direct_once)

  expect_equal(once, direct_once, tolerance = 1e-12)
  expect_equal(twice, direct_twice, tolerance = 1e-12)
  expect_lt(max(abs(crossprod(Q, twice))), 1e-12)
})

test_that("native reorthogonalize_against handles empty and invalid inputs", {
  Q <- matrix(rnorm(20), nrow = 5)
  X <- matrix(rnorm(15), nrow = 5)

  expect_equal(eigencore:::reorthogonalize_against(X, matrix(0, 5, 0), passes = 2L), X)
  expect_equal(eigencore:::reorthogonalize_against(X, Q, passes = 0L), X)
  expect_error(
    eigencore:::reorthogonalize_against(X, matrix(0, 4, 1), passes = 1L),
    "nrow\\(Q\\) == nrow\\(X\\)"
  )
})

test_that("basis workspace reuses coefficient storage for reorthogonalization", {
  set.seed(111)
  Q <- qr.Q(qr(matrix(rnorm(40), nrow = 10)))[, 1:4]
  X <- matrix(rnorm(30), nrow = 10)
  workspace <- eigencore:::basis_workspace(rows = 10, basis_cols = 4, block_cols = 3)
  info <- eigencore:::basis_workspace_info(workspace)

  expect_equal(info$rows, 10)
  expect_equal(info$basis_cols, 4)
  expect_equal(info$block_cols, 3)
  expect_equal(info$allocation_count, 1)
  expect_equal(info$bytes_allocated, 4 * 3 * 8)

  plain <- eigencore:::reorthogonalize_against(X, Q, passes = 2L)
  with_workspace <- eigencore:::reorthogonalize_against(X, Q, passes = 2L, workspace = workspace)
  info_after <- eigencore:::basis_workspace_info(workspace)

  expect_equal(with_workspace, plain, tolerance = 1e-12)
  expect_equal(info_after$allocation_count, info$allocation_count)
  expect_equal(info_after$bytes_allocated, info$bytes_allocated)
})

test_that("basis workspace rejects oversized reorthogonalization requests", {
  workspace <- eigencore:::basis_workspace(rows = 5, basis_cols = 1, block_cols = 1)

  expect_error(
    eigencore:::reorthogonalize_against(
      matrix(rnorm(10), nrow = 5),
      matrix(rnorm(10), nrow = 5),
      workspace = workspace
    ),
    "too small"
  )
})

test_that("native mgs2 handles rank-deficient blocks", {
  X <- cbind(
    c(1, 0, 0, 0),
    c(0, 1, 0, 0),
    c(1, 1, 0, 0),
    c(0, 0, 0, 0)
  )
  out <- eigencore:::mgs2(X, tol = 1e-12)

  expect_equal(out$rank, 2L)
  expect_equal(dim(out$Q), c(4L, 2L))
  expect_equal(dim(out$R), c(2L, 4L))
  expect_equal(out$Q %*% out$R, X, tolerance = 1e-12)
  expect_lt(out$orthogonality, 1e-12)
})

test_that("native mgs2 matches QR reconstruction on full-rank blocks", {
  set.seed(15)
  X <- matrix(rnorm(48), nrow = 12)
  out <- eigencore:::native_mgs2(X)

  expect_equal(out$rank, ncol(X))
  expect_equal(crossprod(out$Q), diag(ncol(X)), tolerance = 1e-12)
  expect_equal(out$Q %*% out$R, X, tolerance = 1e-10)
})

test_that("cholqr2 produces an orthonormal basis", {
  set.seed(12)
  X <- matrix(rnorm(40), nrow = 10)
  out <- eigencore:::cholqr2(X)

  expect_equal(crossprod(out$Q), diag(ncol(out$Q)), tolerance = 1e-12)
  expect_equal(out$Q %*% out$R, X, tolerance = 1e-10)
})

test_that("native cholqr2 matches the QR contract and rejects rank deficiency", {
  set.seed(120)
  X <- matrix(rnorm(60), nrow = 12)
  out <- eigencore:::native_cholqr2(X)

  expect_equal(out$rank, ncol(X))
  expect_equal(crossprod(out$Q), diag(ncol(X)), tolerance = 1e-12)
  expect_equal(out$Q %*% out$R, X, tolerance = 1e-10)
  expect_equal(out$R[lower.tri(out$R)], rep(0, sum(lower.tri(out$R))), tolerance = 1e-14)

  bad <- cbind(diag(4)[, 1:2], diag(4)[, 1])
  expect_error(eigencore:::native_cholqr2(bad), "CholQR2 first pass")
})

test_that("b_orthogonalize normalizes in the B inner product", {
  set.seed(13)
  X <- matrix(rnorm(24), nrow = 8)
  D <- diag(seq(1, 3, length.out = 8))
  out <- eigencore:::b_orthogonalize(X, D)

  expect_equal(crossprod(out$Q, D %*% out$Q), diag(ncol(out$Q)), tolerance = 1e-12)
  expect_equal(out$Q %*% out$R, X, tolerance = 1e-10)
  expect_lt(out$orthogonality, 1e-12)
})

test_that("native B-CholQR2 gives a B-orthonormal basis", {
  set.seed(130)
  X <- matrix(rnorm(40), nrow = 10)
  B0 <- matrix(rnorm(100), nrow = 10)
  B <- crossprod(B0) + diag(10)
  out <- eigencore:::native_b_cholqr2(X, B)

  expect_equal(out$rank, ncol(X))
  expect_equal(crossprod(out$Q, B %*% out$Q), diag(ncol(X)), tolerance = 1e-12)
  expect_equal(out$Q %*% out$R, X, tolerance = 1e-10)
  expect_equal(out$R[lower.tri(out$R)], rep(0, sum(lower.tri(out$R))), tolerance = 1e-14)
})

test_that("native diagonal B-CholQR2 avoids dense B materialization", {
  set.seed(404)
  X <- matrix(rnorm(40), nrow = 10)
  weights <- seq(1, 2, length.out = 10)
  out <- eigencore:::native_diagonal_b_cholqr2(X, weights)

  expect_equal(out$rank, ncol(X))
  expect_equal(crossprod(out$Q, weights * out$Q), diag(ncol(X)), tolerance = 1e-12)
  expect_equal(out$BQ, weights * out$Q, tolerance = 1e-12)
  expect_equal(out$Q %*% out$R, X, tolerance = 1e-10)
})

test_that("B-orthogonalization respects an existing B-orthonormal basis", {
  set.seed(131)
  B0 <- matrix(rnorm(100), nrow = 10)
  B <- crossprod(B0) + diag(10)
  Q0 <- eigencore:::b_orthogonalize(matrix(rnorm(20), nrow = 10), B)$Q
  X <- matrix(rnorm(30), nrow = 10)
  out <- eigencore:::b_orthogonalize(X, B, against = Q0)
  projected <- X - Q0 %*% crossprod(Q0, B %*% X)

  expect_equal(crossprod(out$Q, B %*% out$Q), diag(ncol(out$Q)), tolerance = 1e-12)
  expect_lt(max(abs(crossprod(Q0, B %*% out$Q))), 1e-12)
  expect_equal(out$Q %*% out$R, projected, tolerance = 1e-10)
})

test_that("native orthogonality loss matches direct R formulas", {
  set.seed(14)
  Q <- matrix(rnorm(40), nrow = 10)
  B <- crossprod(matrix(rnorm(100), nrow = 10)) + diag(10)

  standard <- max(abs(crossprod(Q) - diag(ncol(Q))))
  generalized <- max(abs(crossprod(Q, B %*% Q) - diag(ncol(Q))))

  expect_equal(eigencore:::orthogonality_loss(Q), standard, tolerance = 1e-12)
  expect_equal(eigencore:::orthogonality_loss(Q, B = B), generalized, tolerance = 1e-12)
  expect_equal(eigencore:::orthogonality_loss(matrix(0, 10, 0)), 0)
})

test_that("rayleigh_ritz recovers invariant subspace Ritz pairs", {
  A <- diag(c(5, 4, 3, 2, 1))
  Q <- diag(5)[, 1:3]
  rr <- eigencore:::rayleigh_ritz(A, Q, target = largest())

  expect_equal(rr$values[1:3], c(5, 4, 3))
  expect_lt(max(rr$residuals[1:3]), 1e-12)
})

test_that("native symmetric rayleigh_ritz matches direct projected eigensolve", {
  set.seed(16)
  A0 <- matrix(rnorm(36), nrow = 6)
  A <- crossprod(A0)
  Q <- qr.Q(qr(matrix(rnorm(18), nrow = 6)))
  native <- eigencore:::native_rayleigh_ritz_symmetric(A, Q)
  projected <- crossprod(Q, A %*% Q)
  projected <- (projected + t(projected)) / 2
  oracle <- eigen(projected, symmetric = TRUE)
  idx <- order(native$values, decreasing = TRUE)

  expect_equal(native$projected, projected, tolerance = 1e-12)
  expect_equal(native$values[idx], oracle$values, tolerance = 1e-12)
  expect_equal(abs(crossprod(native$vectors[, idx], Q %*% oracle$vectors)), diag(ncol(Q)), tolerance = 1e-10)
})

test_that("rayleigh_ritz reports residuals for non-invariant trial subspaces", {
  set.seed(17)
  A0 <- matrix(rnorm(64), nrow = 8)
  A <- crossprod(A0)
  Q <- qr.Q(qr(matrix(rnorm(24), nrow = 8)))
  rr <- eigencore:::rayleigh_ritz(A, Q, target = largest())
  direct <- eigencore:::dense_eigen_residuals(A, rr$values, rr$vectors)

  expect_equal(rr$residuals, direct, tolerance = 1e-12)
  expect_true(any(rr$residuals > 1e-8))
})

test_that("generalized rayleigh_ritz works with B-orthonormal basis", {
  A <- diag(c(6, 4, 2))
  B <- diag(c(3, 2, 1))
  Q <- eigencore:::b_orthogonalize(diag(3), B)$Q
  rr <- eigencore:::rayleigh_ritz(A, Q, B = B, target = largest())

  expect_equal(unname(rr$values), c(2, 2, 2), tolerance = 1e-12)
  expect_lt(max(rr$residuals), 1e-12)
})
