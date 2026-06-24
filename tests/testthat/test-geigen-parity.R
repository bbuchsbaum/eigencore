skip_if_no_geigen <- function() {
  skip_on_cran()
  skip_if_not(
    requireNamespace("geigen", quietly = TRUE),
    "geigen is not installed"
  )
}

geigen_oracle_values <- function(A, B, target, k) {
  oracle <- geigen::geigen(as.matrix(A), as.matrix(B), symmetric = TRUE)
  values <- Re(oracle$values)
  idx <- eigencore:::order_indices(values, target)
  unname(values[idx[seq_len(k)]])
}

expect_geigen_values <- function(fit, A, B, target, k, tolerance = 1e-8,
                                 unordered = FALSE) {
  expected <- geigen_oracle_values(A, B, target, k)
  actual <- unname(values(fit))
  if (unordered) {
    expected <- sort(expected)
    actual <- sort(actual)
  }
  expect_equal(actual, expected, tolerance = tolerance)
}

expect_b_orthonormal <- function(fit, B, tolerance = 1e-8) {
  V <- vectors(fit)
  B <- as.matrix(B)
  expect_equal(crossprod(V, B %*% V), diag(ncol(V)), tolerance = tolerance)
}

test_that("dense generalized SPD fallback matches geigen oracle", {
  skip_if_no_geigen()

  set.seed(2601)
  A0 <- matrix(rnorm(36), nrow = 6)
  B0 <- matrix(rnorm(36), nrow = 6)
  A <- crossprod(A0) + diag(6)
  B <- crossprod(B0) + diag(6)
  target <- smallest()

  fit <- eig_partial(A, B = B, k = 3L, target = target, tol = 1e-10)

  expect_equal(fit$plan$method, "native dense generalized SPD LAPACK fallback")
  expect_geigen_values(fit, A, B, target, k = 3L, tolerance = 1e-8)
  expect_true(certificate(fit)$passed)
  expect_b_orthonormal(fit, B, tolerance = 1e-10)
})

test_that("eig_full dense generalized SPD matches geigen oracle", {
  skip_if_no_geigen()

  A <- diag(c(1, 4, 9))
  B <- diag(c(1, 2, 3))
  target <- smallest()

  fit <- eig_full(A, B = B, tol = 1e-10)

  expect_equal(fit$method, eigencore:::native_dense_generalized_spd_full_label())
  expect_geigen_values(fit, A, B, target, k = 3L, tolerance = 1e-10)
  expect_true(certificate(fit)$passed)
  expect_b_orthonormal(fit, B, tolerance = 1e-10)
})

test_that("native generalized SPD LOBPCG matches geigen oracle", {
  skip_if_no_geigen()

  A <- diag(c(2, 5, 9, 14, 20, 27))
  B <- diag(c(1, 2, 3, 4, 5, 6))
  target <- largest()

  fit <- eig_partial(
    A,
    B = B,
    k = 2L,
    target = target,
    method = lobpcg(maxit = 80L),
    seed = 2602,
    tol = 1e-9
  )

  expect_equal(fit$method, eigencore:::native_generalized_lobpcg_label())
  expect_geigen_values(fit, A, B, target, k = 2L, tolerance = 1e-8)
  expect_true(certificate(fit)$passed)
  expect_b_orthonormal(fit, B, tolerance = 1e-8)
})

test_that("sparse diagonal-B generalized Lanczos matches geigen oracle", {
  skip_if_no_geigen()

  A <- Matrix::Diagonal(x = c(1, 4, 9, 16, 25, 36))
  B <- Matrix::Diagonal(x = c(1, 2, 3, 4, 5, 6))
  target <- smallest()

  fit <- eig_partial(
    A,
    B = B,
    k = 3L,
    target = target,
    method = lanczos(max_subspace = 6L),
    seed = 2603,
    tol = 1e-9,
    allow_dense_fallback = "never"
  )

  expect_equal(fit$method, eigencore:::native_generalized_lanczos_label())
  expect_geigen_values(fit, A, B, target, k = 3L, tolerance = 1e-8)
  expect_true(certificate(fit)$passed)
  expect_b_orthonormal(fit, B, tolerance = 1e-8)
  expect_true(fit$restart$native)
  expect_equal(
    fit$restart$metric_solve,
    "diagonal scaling similarity transform for B"
  )
})

test_that("dense generalized shift-invert matches geigen oracle", {
  skip_if_no_geigen()

  n <- 10L
  A <- symmetric_with_spectrum(seq_len(n), seed = 2604)
  B <- diag(seq(1, 2, length.out = n))
  sigma <- 3.2
  target <- nearest(sigma)

  fit <- eig_partial(
    A,
    B = B,
    k = 2L,
    target = target,
    method = shift_invert(sigma = sigma),
    tol = 1e-9
  )

  expect_equal(
    fit$method,
    eigencore:::native_dense_generalized_shift_invert_label()
  )
  expect_geigen_values(
    fit,
    A,
    B,
    target,
    k = 2L,
    tolerance = 1e-7,
    unordered = TRUE
  )
  expect_true(certificate(fit)$passed)
  expect_b_orthonormal(fit, B, tolerance = 1e-8)
})

test_that("sparse tridiagonal generalized shift-invert matches geigen oracle", {
  skip_if_no_geigen()

  n <- 30L
  A <- Matrix::bandSparse(
    n,
    k = c(-1, 0, 1),
    diagonals = list(rep(-1, n - 1L), rep(2.5, n), rep(-1, n - 1L))
  )
  B <- Matrix::Diagonal(x = seq(1, 2, length.out = n))
  oracle <- geigen::geigen(as.matrix(A), as.matrix(B), symmetric = TRUE)
  sigma <- sort(oracle$values)[4L] + 0.01
  target <- nearest(sigma)

  fit <- eig_partial(
    A,
    B = B,
    k = 3L,
    target = target,
    method = shift_invert(sigma = sigma),
    tol = 1e-8,
    allow_dense_fallback = "never"
  )

  expect_equal(
    fit$method,
    eigencore:::native_tridiagonal_generalized_shift_invert_label()
  )
  expect_geigen_values(
    fit,
    A,
    B,
    target,
    k = 3L,
    tolerance = 1e-6,
    unordered = TRUE
  )
  expect_true(certificate(fit)$passed)
  expect_b_orthonormal(fit, B, tolerance = 1e-8)
  expect_true(fit$restart$native)
  expect_true(fit$restart$generalized)
})
