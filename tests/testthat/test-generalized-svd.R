gsvd_adjoint <- function(x) {
  if (is.complex(x)) {
    Conj(t(x))
  } else {
    t(x)
  }
}

expect_gsvd_reconstruction <- function(fit, A, B, tolerance = 1e-9) {
  expect_equal(
    fit$U %*% fit$D1 %*% fit$zero_R %*% gsvd_adjoint(fit$Q),
    A,
    tolerance = tolerance
  )
  expect_equal(
    fit$V %*% fit$D2 %*% fit$zero_R %*% gsvd_adjoint(fit$Q),
    B,
    tolerance = tolerance
  )
  expect_equal(crossprod(fit$U), diag(nrow(A)), tolerance = tolerance)
  expect_equal(crossprod(fit$V), diag(nrow(B)), tolerance = tolerance)
  expect_equal(crossprod(fit$Q), diag(ncol(A)), tolerance = tolerance)
}

test_that("generalized_svd computes finite real dense GSVD values", {
  A <- diag(c(3, 4, 5))
  B <- diag(c(4, 3, 2))

  fit <- generalized_svd(A, B, tol = 1e-10)
  coords <- alpha_beta(fit)

  expect_s3_class(fit, "eigencore_gsvd_result")
  expect_identical(fit$method, eigencore:::native_dense_generalized_svd_label())
  expect_identical(fit$plan$method, fit$method)
  expect_match(fit$method, "dggsvd")
  expect_identical(fit$plan$controls$lapack_driver, "dggsvd")
  expect_false(fit$plan$controls$sparse_densified)
  expect_equal(coords$classification, rep("finite", 3L))
  expect_equal(coords$alpha^2 + coords$beta^2, rep(1, 3L), tolerance = 1e-12)
  expect_equal(sort(unname(values(fit))), sort(c(3 / 4, 4 / 3, 5 / 2)),
               tolerance = 1e-10)
  expect_true(certificate(fit)$passed)
  expect_gsvd_reconstruction(fit, A, B, tolerance = 1e-10)
})

test_that("generalized_svd accessors return exact GSVD factors", {
  A <- diag(c(3, 4, 5))
  B <- diag(c(4, 3, 2))

  fit <- generalized_svd(A, B, tol = 1e-10)

  expect_equal(left_vectors(fit), fit$U)
  expect_equal(right_vectors(fit), fit$V)
  expect_false(identical(left_vectors(fit), fit$undefined))
  expect_false(identical(right_vectors(fit), fit$values))
})

test_that("generalized_svd keeps small positive beta values finite", {
  scale <- 1e8
  A <- matrix(scale, 1, 1)
  B <- matrix(1, 1, 1)

  fit <- generalized_svd(A, B, tol = 1e-10)
  coords <- alpha_beta(fit)

  expect_gt(coords$beta, 0)
  expect_lt(coords$beta, sqrt(.Machine$double.eps))
  expect_identical(coords$classification, "finite")
  expect_true(is.finite(values(fit)))
  expect_equal(unname(values(fit) / scale), 1, tolerance = 1e-12)
  expect_true(certificate(fit)$passed)
  expect_gsvd_reconstruction(fit, A, B, tolerance = 1e-8)
})

test_that("generalized_svd covers geigen manual rectangular rank layout", {
  A <- matrix(c(1, 2, 3, 3, 2, 1, 4, 5, 6, 7, 8, 8),
              nrow = 2, byrow = TRUE)
  B <- matrix(1:18, byrow = TRUE, ncol = 6)

  fit <- generalized_svd(A, B, tol = 1e-7)
  coords <- alpha_beta(fit)

  expect_equal(fit$dimensions, c(m = 2, n = 6, p = 3))
  expect_equal(fit$rank, fit$k + fit$l)
  expect_equal(dim(fit$D1), c(nrow(A), fit$rank))
  expect_equal(dim(fit$D2), c(nrow(B), fit$rank))
  expect_equal(dim(fit$R), c(fit$rank, fit$rank))
  expect_equal(dim(fit$zero_R), c(fit$rank, ncol(A)))
  expect_true(any(coords$infinite))
  expect_true(any(coords$undefined))
  expect_true(certificate(fit)$passed)
  expect_gsvd_reconstruction(fit, A, B, tolerance = 1e-8)
})

test_that("generalized_svd optionally matches geigen gsvd alpha/beta output", {
  skip_if_not_installed("geigen")

  A <- diag(c(3, 4, 5))
  B <- diag(c(4, 3, 2))

  fit <- generalized_svd(A, B, tol = 1e-10)
  ref <- get("gsvd", envir = asNamespace("geigen"))(A, B)

  expect_equal(fit$k, ref$k)
  expect_equal(fit$l, ref$l)
  expect_equal(fit$alpha, ref$alpha, tolerance = 1e-10)
  expect_equal(fit$beta, ref$beta, tolerance = 1e-10)
})

test_that("generalized_svd rejects unsupported inputs without densifying", {
  expect_error(
    generalized_svd(Matrix::Diagonal(2), diag(2)),
    "base dense matrix"
  )
  expect_error(
    generalized_svd(diag(2), Matrix::Diagonal(2)),
    "base dense matrix"
  )
  expect_error(
    generalized_svd(matrix(1 + 1i, 1, 1), matrix(1 + 0i, 1, 1)),
    "native complex GSVD requires a complex LAPACK GSVD driver"
  )
  expect_error(
    generalized_svd(diag(2), matrix(1, 3, 3)),
    "same number of columns"
  )
  expect_error(
    generalized_svd(matrix(c(1, NA), 1, 2), matrix(1, 1, 2)),
    "finite"
  )
  expect_error(
    generalized_svd(diag(2), diag(2), tol = Inf),
    "single finite non-negative"
  )
  expect_error(
    generalized_svd(diag(2), diag(2), tol = NaN),
    "single finite non-negative"
  )
})
