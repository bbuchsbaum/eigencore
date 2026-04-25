test_that("eig_partial returns certified eigenpairs", {
  A <- diag(c(4, 1, 3, 2))
  fit <- eig_partial(A, k = 2, target = largest())

  expect_equal(values(fit), c(4, 3))
  expect_true(certificate(fit)$passed)
  expect_equal(fit$nconv, 2)
  expect_equal(fit$plan$method, "native dense Hermitian LAPACK fallback")
  expect_equal(fit$method, "native dense Hermitian LAPACK fallback")
})

test_that("native dense Hermitian fallback matches base eigen ordering", {
  set.seed(24)
  A0 <- matrix(rnorm(36), nrow = 6)
  A <- crossprod(A0)
  fit <- eig_partial(A, k = 3, target = smallest())
  oracle <- eigen(A, symmetric = TRUE)

  expect_equal(values(fit), sort(oracle$values)[1:3], tolerance = 1e-10)
  expect_true(certificate(fit)$passed)
  expect_match(fit$warnings, "native dense Hermitian LAPACK fallback")
})

test_that("generalized SPD eigenproblem is certified in original coordinates", {
  A <- diag(c(6, 4, 2))
  B <- diag(c(3, 2, 1))
  fit <- eig_partial(A, B = B, k = 2, target = largest())

  expect_equal(unname(values(fit)), c(2, 2), tolerance = 1e-10)
  expect_true(certificate(fit)$passed)
  expect_equal(fit$plan$method, "native dense generalized SPD LAPACK fallback")
  expect_match(fit$warnings, "native dense generalized SPD LAPACK fallback")
})

test_that("native generalized SPD fallback produces B-orthonormal vectors", {
  set.seed(26)
  A0 <- matrix(rnorm(25), nrow = 5)
  B0 <- matrix(rnorm(25), nrow = 5)
  A <- crossprod(A0) + diag(5)
  B <- crossprod(B0) + diag(5)
  fit <- eig_partial(A, B = B, k = 3, target = smallest())

  expect_true(certificate(fit)$passed)
  expect_equal(crossprod(vectors(fit), B %*% vectors(fit)), diag(3), tolerance = 1e-10)
  expect_lt(max(abs(A %*% vectors(fit) - B %*% vectors(fit) %*% diag(values(fit)))), 1e-8)
})

test_that("svd_partial returns sorted singular values and certificate", {
  A <- diag(c(5, 1, 3))
  fit <- svd_partial(A, rank = 2)

  expect_equal(values(fit), c(5, 3))
  expect_true(certificate(fit)$passed)
  expect_equal(fit$nconv, 2)
  expect_equal(fit$plan$method, "native dense LAPACK SVD fallback")
  expect_equal(fit$method, "native dense LAPACK SVD fallback")
})

test_that("native dense SVD fallback matches base SVD on rectangular inputs", {
  set.seed(25)
  A <- matrix(rnorm(35), nrow = 7)
  fit <- svd_partial(A, rank = 3, target = largest())
  oracle <- svd(A, nu = 3, nv = 3)

  expect_equal(values(fit), oracle$d[1:3], tolerance = 1e-10)
  expect_true(certificate(fit)$passed)
  expect_match(fit$warnings, "native dense LAPACK SVD fallback")
})

test_that("shift_invert is honest until implemented", {
  expect_error(
    eig_partial(diag(3), k = 1, method = shift_invert(0)),
    "not implemented"
  )
})

test_that("RSpectra-compatible shims expose core fields", {
  A <- diag(c(3, 2, 1))
  ef <- eigs_sym(A, k = 1)
  sf <- svds(A, k = 1)

  expect_equal(ef$values, 3)
  expect_equal(sf$d, 3)
  expect_s3_class(ef$certificate, "eigencore_certificate")
  expect_s3_class(sf$certificate, "eigencore_certificate")
})

test_that("RSpectra SM maps to smallest magnitude, not smallest algebraic", {
  A <- diag(c(-10, 2, 5))
  ef <- eigs_sym(A, k = 1, which = "SM")

  expect_equal(ef$values, 2)
})

test_that("prototype Lanczos solves matrix-free Hermitian operators", {
  vals <- c(6, 5, 4, 3, 2, 1)
  op <- linear_operator(
    dim = c(length(vals), length(vals)),
    apply = function(X, alpha = 1, beta = 0, Y = NULL) {
      out <- alpha * (vals * X)
      if (!is.null(Y) && beta != 0) {
        out <- out + beta * Y
      }
      out
    },
    apply_adjoint = function(X, alpha = 1, beta = 0, Y = NULL) {
      out <- alpha * (vals * X)
      if (!is.null(Y) && beta != 0) {
        out <- out + beta * Y
      }
      out
    },
    structure = hermitian(),
    name = "diagonal_matrix_free"
  )

  fit <- eig_partial(op, k = 2, method = lanczos(max_subspace = 6), seed = 123)

  expect_equal(values(fit), c(6, 5), tolerance = 1e-10)
  expect_true(certificate(fit)$passed)
  expect_equal(fit$nconv, 2)
  expect_equal(fit$method, "prototype Hermitian Lanczos")
})

test_that("auto uses Lanczos for matrix-free Hermitian operators", {
  vals <- c(4, 3, 2, 1)
  op <- linear_operator(
    dim = c(4, 4),
    apply = function(X, alpha = 1, beta = 0, Y = NULL) {
      out <- alpha * (vals * X)
      if (!is.null(Y) && beta != 0) out <- out + beta * Y
      out
    },
    apply_adjoint = function(X, alpha = 1, beta = 0, Y = NULL) {
      out <- alpha * (vals * X)
      if (!is.null(Y) && beta != 0) out <- out + beta * Y
      out
    },
    structure = hermitian()
  )

  fit <- eig_partial(op, k = 1, seed = 321)
  expect_equal(values(fit), 4, tolerance = 1e-10)
  expect_equal(fit$plan$method, "prototype Hermitian Lanczos")
})

test_that("auto uses native CSC-backed Lanczos for sparse Hermitian matrices", {
  A <- Matrix::sparseMatrix(
    i = c(1, 2, 3, 4),
    j = c(1, 2, 3, 4),
    x = c(8, 5, 3, 1),
    dims = c(4, 4)
  )
  op <- as_operator(A)
  fit <- eig_partial(A, k = 2, seed = 101)

  expect_null(eigencore:::source_or_null(op))
  expect_equal(fit$plan$method, "native CSC-backed prototype Hermitian Lanczos")
  expect_equal(fit$method, "native CSC-backed prototype Hermitian Lanczos")
  expect_match(fit$warnings, "native CSC block apply")
  expect_equal(values(fit), c(8, 5), tolerance = 1e-10)
  expect_true(certificate(fit)$passed)
})

test_that("prototype Golub-Kahan solves matrix-free rectangular SVD", {
  sing <- c(7, 5, 3, 1)
  m <- 6
  n <- length(sing)
  A <- rbind(diag(sing), matrix(0, m - n, n))
  op <- linear_operator(
    dim = c(m, n),
    apply = function(X, alpha = 1, beta = 0, Y = NULL) {
      out <- alpha * (A %*% X)
      if (!is.null(Y) && beta != 0) out <- out + beta * Y
      out
    },
    apply_adjoint = function(X, alpha = 1, beta = 0, Y = NULL) {
      out <- alpha * (t(A) %*% X)
      if (!is.null(Y) && beta != 0) out <- out + beta * Y
      out
    },
    name = "rectangular_matrix_free"
  )

  fit <- svd_partial(op, rank = 2, method = golub_kahan(max_subspace = 4), seed = 123)

  expect_equal(values(fit), c(7, 5), tolerance = 1e-10)
  expect_true(certificate(fit)$passed)
  expect_equal(fit$nconv, 2)
  expect_equal(fit$method, "prototype Golub-Kahan")
})

test_that("auto uses Golub-Kahan for matrix-free SVD", {
  A <- rbind(diag(c(4, 2, 1)), 0)
  op <- linear_operator(
    dim = dim(A),
    apply = function(X, alpha = 1, beta = 0, Y = NULL) {
      out <- alpha * (A %*% X)
      if (!is.null(Y) && beta != 0) out <- out + beta * Y
      out
    },
    apply_adjoint = function(X, alpha = 1, beta = 0, Y = NULL) {
      out <- alpha * (t(A) %*% X)
      if (!is.null(Y) && beta != 0) out <- out + beta * Y
      out
    }
  )

  fit <- svd_partial(op, rank = 1, seed = 456)
  expect_equal(values(fit), 4, tolerance = 1e-10)
  expect_equal(fit$plan$method, "prototype Golub-Kahan")
})

test_that("auto uses native CSC-backed Golub-Kahan for sparse rectangular SVD", {
  A <- Matrix::sparseMatrix(
    i = c(1, 2, 3, 4),
    j = c(1, 2, 3, 4),
    x = c(9, 6, 2, 1),
    dims = c(6, 4)
  )
  op <- as_operator(A)
  fit <- svd_partial(A, rank = 2, seed = 202)

  expect_null(eigencore:::source_or_null(op))
  expect_equal(fit$plan$method, "native CSC-backed prototype Golub-Kahan")
  expect_equal(fit$method, "native CSC-backed prototype Golub-Kahan")
  expect_match(fit$warnings, "native CSC block apply")
  expect_equal(values(fit), c(9, 6), tolerance = 1e-10)
  expect_true(certificate(fit)$passed)
})

test_that("Golub-Kahan vector modes withhold full certificate when needed", {
  A <- rbind(diag(c(3, 2, 1)), 0)
  op <- linear_operator(
    dim = dim(A),
    apply = function(X, alpha = 1, beta = 0, Y = NULL) alpha * (A %*% X),
    apply_adjoint = function(X, alpha = 1, beta = 0, Y = NULL) alpha * (t(A) %*% X)
  )

  fit <- svd_partial(op, rank = 1, method = golub_kahan(max_subspace = 3), vectors = "right", seed = 789)
  expect_null(left_vectors(fit))
  expect_false(certificate(fit)$passed)
  expect_match(certificate(fit)$notes, "both left and right vectors")
})
