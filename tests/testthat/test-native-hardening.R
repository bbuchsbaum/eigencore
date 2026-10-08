known_svd_matrix <- function(m, n, d, seed = 1L) {
  set.seed(seed)
  U <- qr.Q(qr(matrix(rnorm(m * length(d)), m, length(d))))
  V <- qr.Q(qr(matrix(rnorm(n * length(d)), n, length(d))))
  U %*% (d * t(V))
}

test_that("native scalar Lanczos and Golub-Kahan reject non-finite input", {
  set.seed(11)
  A <- crossprod(matrix(rnorm(400), 20, 20))
  A_bad <- A
  A_bad[3, 5] <- NaN
  A_bad[5, 3] <- NaN
  msg <- "non-finite value encountered \\(check input for NA/NaN/Inf\\)"

  expect_error(
    .Call("eigencore_lanczos_dense", A_bad, 10L, rnorm(20), 2L, 0L, 1e-8,
          PACKAGE = "eigencore"),
    msg
  )
  start_bad <- rnorm(20)
  start_bad[4] <- Inf
  expect_error(
    .Call("eigencore_lanczos_dense", A, 10L, start_bad, 2L, 0L, 1e-8,
          PACKAGE = "eigencore"),
    msg
  )

  B <- matrix(rnorm(30 * 12), 30, 12)
  B_bad <- B
  B_bad[7, 2] <- NA_real_
  expect_error(
    .Call("eigencore_golub_kahan_dense", B_bad, 8L, rnorm(12), 2L, 0L, 1e-8,
          FALSE, PACKAGE = "eigencore"),
    msg
  )
  start_bad <- rnorm(12)
  start_bad[1] <- NaN
  expect_error(
    .Call("eigencore_golub_kahan_dense", B, 8L, start_bad, 2L, 0L, 1e-8,
          FALSE, PACKAGE = "eigencore"),
    msg
  )
})

test_that("native Golub-Kahan alpha breakdown does not report a zero left vector", {
  set.seed(12)
  # Rank-2 matrix: a start vector with a null-space component makes the
  # third left vector vanish (alpha breakdown) after two full steps.
  A <- matrix(rnorm(30 * 2), 30, 2) %*% matrix(rnorm(2 * 20), 2, 20)
  out <- .Call("eigencore_golub_kahan_dense", A, 6L, rnorm(20), 1L, 0L, 1e-8,
               FALSE, PACKAGE = "eigencore")
  expect_equal(out$iterations, 2L)
  expect_equal(ncol(out$U), out$iterations)
  expect_equal(unname(sqrt(colSums(out$U^2))), rep(1, out$iterations),
               tolerance = 1e-10)
  expect_true(all(out$alpha > 0))
})

test_that("dense symmetry check treats non-finite entries as not symmetric", {
  A <- diag(4)
  expect_true(.Call("eigencore_dense_is_symmetric", A, 1e-10, PACKAGE = "eigencore"))
  A_nan <- A
  A_nan[1, 2] <- NaN
  A_nan[2, 1] <- NaN
  expect_false(.Call("eigencore_dense_is_symmetric", A_nan, 1e-10, PACKAGE = "eigencore"))
  A_inf <- A
  A_inf[3, 3] <- Inf
  expect_false(.Call("eigencore_dense_is_symmetric", A_inf, 1e-10, PACKAGE = "eigencore"))
})

test_that("diagonal B CholQR2 rejects NaN diagonal entries", {
  X <- matrix(rnorm(12), 4, 3)
  expect_error(
    .Call("eigencore_diagonal_b_cholqr2", X, c(1, NaN, 1, 1), FALSE,
          PACKAGE = "eigencore"),
    "B diagonal must be positive"
  )
})

test_that("sparse Gram SVD fast path is scale invariant", {
  set.seed(13)
  A <- Matrix::rsparsematrix(4000, 60, 0.2)
  expected <- svd(as.matrix(A), nu = 0, nv = 0)$d[1:3]
  for (s in c(1, 1e-10)) {
    fit <- svd_partial(A * s, rank = 3)
    expect_equal(fit$d / s, expected, tolerance = 1e-10)
  }
})

test_that("dense svd_partial recovers top singular values of a tiny-scale matrix", {
  # The C-side Gram zero tolerance is scale invariant; the dense path also
  # depends on the R-side gram_svd_zero_tolerance(), so only run this once
  # that helper no longer floors its scale at 1.
  skip_if(
    eigencore:::gram_svd_zero_tolerance(5e-10, 1e-8) >= 5e-10,
    "R-side gram_svd_zero_tolerance() still uses an absolute scale floor"
  )
  d <- 10^seq(0, -3, length.out = 60) * 5
  A <- known_svd_matrix(3000, 60, d)
  fit <- svd_partial(A * 1e-10, rank = 3)
  expect_equal(fit$d / 1e-10, d[1:3], tolerance = 1e-10)
})
