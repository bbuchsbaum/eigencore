test_that("validate_eigen_accuracy compares eigencore against base eigen", {
  A <- diag(c(7, 2, 5, 1))
  fit <- eig_partial(A, k = 2)
  val <- eigencore:::validate_eigen_accuracy(A, k = 2, fit = fit)

  expect_true(val$passed)
  expect_equal(val$oracle_values, c(7, 5))
  expect_lt(val$max_value_abs_error, 1e-12)
  expect_true(val$certificate_agrees)
})

test_that("validate_svd_accuracy compares eigencore against base svd", {
  A <- diag(c(8, 3, 6, 1))
  fit <- svd_partial(A, rank = 2)
  val <- eigencore:::validate_svd_accuracy(A, rank = 2, fit = fit)

  expect_true(val$passed)
  expect_equal(val$oracle_values, c(8, 6))
  expect_lt(val$max_value_abs_error, 1e-12)
  expect_true(val$certificate_agrees)
})

test_that("dense and operator certificates use the same eigen scale for explicit operators", {
  A <- diag(c(4, 2, 1))
  v <- diag(3)[, 1:2]
  vals <- c(4, 2)
  dense <- eigencore:::certify_eigen(A, vals, v)
  op <- eigencore:::certify_eigen_operator(as_operator(A), vals, v)

  expect_equal(op$backward_error, dense$backward_error, tolerance = 1e-14)
  expect_equal(op$scale, dense$scale, tolerance = 1e-14)
  expect_equal(op$norm_bound_type, "frobenius_exact+identity_exact")
  expect_false(op$scale_is_estimate)
})

test_that("dense and operator certificates use the same SVD scale for explicit operators", {
  A <- diag(c(4, 2, 1))
  s <- svd(A, nu = 2, nv = 2)
  dense <- eigencore:::certify_svd(A, s$d[1:2], s$u[, 1:2], s$v[, 1:2])
  op <- eigencore:::certify_svd_operator(as_operator(A), s$d[1:2], s$u[, 1:2], s$v[, 1:2])

  expect_equal(op$backward_error, dense$backward_error, tolerance = 1e-14)
  expect_equal(op$scale, dense$scale, tolerance = 1e-14)
  expect_equal(op$norm_bound_type, "frobenius_exact")
  expect_false(op$scale_is_estimate)
})

test_that("built-in CSC and diagonal operator certificates use native diagnostics", {
  A <- Matrix::Diagonal(x = c(5, 3, 1))
  vals <- c(5, 3)
  vecs <- diag(3)[, 1:2]
  diag_cert <- eigencore:::certify_eigen_operator(as_operator(A), vals, vecs)

  expect_true(diag_cert$passed)
  expect_equal(diag_cert$norm_bound_type, "frobenius_metadata+identity_exact")
  expect_equal(diag_cert$residuals, c(0, 0), tolerance = 1e-14)

  S <- Matrix::sparseMatrix(i = c(1, 2, 3), j = c(1, 2, 3), x = c(6, 4, 2), dims = c(3, 3))
  csc_cert <- eigencore:::certify_eigen_operator(as_operator(S), c(6, 4), vecs)
  dense_cert <- eigencore:::certify_eigen(as.matrix(S), c(6, 4), vecs)

  expect_equal(csc_cert$residuals, dense_cert$residuals, tolerance = 1e-14)
  expect_equal(csc_cert$backward_error, dense_cert$backward_error, tolerance = 1e-14)
  expect_equal(csc_cert$norm_bound_type, "frobenius_metadata+identity_exact")
})

test_that("built-in CSC and diagonal SVD certificates use native diagnostics", {
  D <- Matrix::Diagonal(x = c(7, 3, 1))
  s <- svd(as.matrix(D), nu = 2, nv = 2)
  diag_cert <- eigencore:::certify_svd_operator(as_operator(D), s$d[1:2], s$u[, 1:2], s$v[, 1:2])

  expect_true(diag_cert$passed)
  expect_equal(diag_cert$norm_bound_type, "frobenius_metadata")
  expect_equal(diag_cert$residuals$combined, c(0, 0), tolerance = 1e-14)

  S <- Matrix::sparseMatrix(i = c(1, 2, 3, 4), j = c(1, 2, 3, 4), x = c(8, 5, 2, 1), dims = c(4, 4))
  ss <- svd(as.matrix(S), nu = 2, nv = 2)
  csc_cert <- eigencore:::certify_svd_operator(as_operator(S), ss$d[1:2], ss$u[, 1:2], ss$v[, 1:2])
  dense_cert <- eigencore:::certify_svd(as.matrix(S), ss$d[1:2], ss$u[, 1:2], ss$v[, 1:2])

  expect_equal(csc_cert$residuals$combined, dense_cert$residuals$combined, tolerance = 1e-14)
  expect_equal(csc_cert$backward_error, dense_cert$backward_error, tolerance = 1e-14)
  expect_equal(csc_cert$norm_bound_type, "frobenius_metadata")
})

test_that("matrix-free stochastic norm estimates withhold passed certificates", {
  vals <- c(4, 2, 1)
  op <- linear_operator(
    dim = c(3, 3),
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

  cert <- eigencore:::certify_eigen_operator(op, vals[1:2], diag(3)[, 1:2])

  expect_equal(cert$norm_bound_type, "frobenius_hutchinson_estimate+identity_exact")
  expect_true(cert$scale_is_estimate)
  expect_true(all(cert$converged))
  expect_false(cert$passed)
  expect_match(cert$notes, "stochastic norm estimate")
})

test_that("native column norms preserve certificate residual formulas", {
  set.seed(20)
  X <- matrix(rnorm(35), nrow = 7)
  A <- matrix(rnorm(49), nrow = 7)
  A <- crossprod(A)
  eig <- eigen(A, symmetric = TRUE)
  values <- eig$values[1:3]
  vectors <- eig$vectors[, 1:3]
  perturbed <- vectors
  perturbed[1, 1] <- perturbed[1, 1] + 1e-5
  residual_matrix <- A %*% perturbed - sweep(perturbed, 2L, values, `*`)

  expect_equal(eigencore:::col_norms(X), sqrt(colSums(X^2)), tolerance = 1e-14)
  expect_lt(max(abs(
    eigencore:::certify_eigen(A, values, perturbed)$residuals -
      sqrt(colSums(residual_matrix^2))
  )), 1e-12)
})

test_that("native dense eigen residuals match direct standard and generalized formulas", {
  set.seed(21)
  A <- crossprod(matrix(rnorm(36), nrow = 6))
  B <- crossprod(matrix(rnorm(36), nrow = 6)) + diag(6)
  vectors <- qr.Q(qr(matrix(rnorm(18), nrow = 6)))
  values <- c(3, 2, 1)

  standard_matrix <- A %*% vectors - sweep(vectors, 2L, values, `*`)
  generalized_matrix <- A %*% vectors - sweep(B %*% vectors, 2L, values, `*`)

  expect_equal(
    eigencore:::dense_eigen_residuals(A, values, vectors),
    sqrt(colSums(standard_matrix^2)),
    tolerance = 1e-12
  )
  expect_equal(
    eigencore:::dense_eigen_residuals(A, values, vectors, B = B),
    sqrt(colSums(generalized_matrix^2)),
    tolerance = 1e-12
  )
})

test_that("native dense eigen certificate diagnostics match R formula contract", {
  set.seed(211)
  A <- crossprod(matrix(rnorm(36), nrow = 6))
  B <- crossprod(matrix(rnorm(36), nrow = 6)) + diag(6)
  vectors <- qr.Q(qr(matrix(rnorm(18), nrow = 6)))
  values <- c(3, 2, 1)
  tol <- 1e-8

  standard <- eigencore:::native_dense_eigen_certificate(A, values, vectors, tol = tol)
  residuals <- eigencore:::dense_eigen_residuals(A, values, vectors)
  scale <- eigencore:::eigen_backward_scale(norm(A, type = "F"), 1, values, vectors)

  expect_equal(standard$residuals, residuals, tolerance = 1e-12)
  expect_equal(standard$scale, scale, tolerance = 1e-12)
  expect_equal(standard$backward_error, residuals / scale, tolerance = 1e-12)
  expect_equal(standard$orthogonality, eigencore:::orthogonality_loss(vectors), tolerance = 1e-12)
  expect_equal(standard$converged, standard$backward_error <= tol)

  generalized <- eigencore:::native_dense_eigen_certificate(A, values, vectors, B = B, tol = tol)
  generalized_residuals <- eigencore:::dense_eigen_residuals(A, values, vectors, B = B)
  generalized_scale <- eigencore:::eigen_backward_scale(norm(A, type = "F"), norm(B, type = "F"), values, vectors)

  expect_equal(generalized$residuals, generalized_residuals, tolerance = 1e-12)
  expect_equal(generalized$scale, generalized_scale, tolerance = 1e-12)
  expect_equal(generalized$backward_error, generalized_residuals / generalized_scale, tolerance = 1e-12)
  expect_equal(generalized$orthogonality, eigencore:::orthogonality_loss(vectors, B = B), tolerance = 1e-12)
})

test_that("certificate passed requires residual convergence and orthogonality", {
  cert <- eigencore:::new_certificate(
    tol = 1e-8,
    residuals = c(1e-12, 1e-12),
    backward_error = c(1e-12, 1e-12),
    orthogonality = 0.5,
    converged = c(TRUE, TRUE),
    scale = c(1, 1)
  )

  expect_false(cert$passed)
  expect_false(cert$orthogonality_passed)
  expect_equal(cert$failed_indices, integer(0))
  expect_match(paste(cert$notes, collapse = " "), "orthogonality loss")

  nonnormal_cert <- eigencore:::new_certificate(
    tol = 1e-8,
    residuals = c(1e-12, 1e-12),
    backward_error = c(1e-12, 1e-12),
    orthogonality = 0.5,
    converged = c(TRUE, TRUE),
    scale = c(1, 1),
    require_orthogonality = FALSE
  )

  expect_true(nonnormal_cert$passed)
  expect_true(nonnormal_cert$orthogonality_passed)
  expect_false(nonnormal_cert$orthogonality_required)
})

test_that("operator generalized eigen certificates use dense native diagnostics when available", {
  set.seed(212)
  A <- crossprod(matrix(rnorm(36), nrow = 6)) + diag(6)
  B <- crossprod(matrix(rnorm(36), nrow = 6)) + diag(6)
  eig <- eigencore:::dense_generalized_spd_eigen(A, B)
  values <- eig$values[1:3]
  vectors <- eig$vectors[, 1:3, drop = FALSE]

  cert <- eigencore:::certify_eigen_operator(
    as_operator(A),
    values,
    vectors,
    Bop = as_operator(B),
    tol = 1e-8
  )

  expect_true(cert$passed)
  expect_equal(cert$norm_bound_type, "frobenius_exact+frobenius_exact")
  expect_equal(cert$residuals, eigencore:::dense_eigen_residuals(A, values, vectors, B = B),
               tolerance = 1e-12)
  expect_lt(cert$max_orthogonality_loss, 1e-10)
})

test_that("residual-backed generalized operator certificates preserve original-coordinate contract", {
  set.seed(213)
  A <- crossprod(matrix(rnorm(36), nrow = 6)) + diag(6)
  B <- crossprod(matrix(rnorm(36), nrow = 6)) + diag(6)
  eig <- eigencore:::dense_generalized_spd_eigen(A, B)
  values <- eig$values[1:3]
  vectors <- eig$vectors[, 1:3, drop = FALSE]
  residuals <- eigencore:::dense_eigen_residuals(A, values, vectors, B = B)

  from_residuals <- eigencore:::certify_eigen_operator_residuals(
    as_operator(A),
    values,
    vectors,
    residuals,
    Bop = as_operator(B),
    tol = 1e-8
  )
  direct <- eigencore:::certify_eigen_operator(
    as_operator(A),
    values,
    vectors,
    Bop = as_operator(B),
    tol = 1e-8
  )

  expect_true(from_residuals$passed)
  expect_equal(from_residuals$residuals, direct$residuals, tolerance = 1e-12)
  expect_equal(from_residuals$backward_error, direct$backward_error, tolerance = 1e-12)
  expect_equal(from_residuals$orthogonality, direct$orthogonality, tolerance = 1e-12)
  expect_equal(from_residuals$norm_bound_type, direct$norm_bound_type)
})

test_that("native dense SVD residuals match direct two-sided formulas", {
  set.seed(22)
  A <- matrix(rnorm(35), nrow = 7)
  u <- qr.Q(qr(matrix(rnorm(21), nrow = 7)))
  v <- qr.Q(qr(matrix(rnorm(15), nrow = 5)))
  d <- c(4, 2, 1)

  left_matrix <- A %*% v - sweep(u, 2L, d, `*`)
  right_matrix <- crossprod(A, u) - sweep(v, 2L, d, `*`)
  left <- sqrt(colSums(left_matrix^2))
  right <- sqrt(colSums(right_matrix^2))
  residuals <- eigencore:::dense_svd_residuals(A, d, u, v)

  expect_equal(residuals$left, left, tolerance = 1e-12)
  expect_equal(residuals$right, right, tolerance = 1e-12)
  expect_equal(residuals$combined, sqrt(left^2 + right^2), tolerance = 1e-12)
})

test_that("native dense SVD certificate diagnostics match R formula contract", {
  set.seed(221)
  A <- matrix(rnorm(35), nrow = 7)
  u <- qr.Q(qr(matrix(rnorm(21), nrow = 7)))
  v <- qr.Q(qr(matrix(rnorm(15), nrow = 5)))
  d <- c(4, 2, 1)
  tol <- 1e-8
  diag <- eigencore:::native_dense_svd_certificate(A, d, u, v, tol = tol)
  residuals <- eigencore:::dense_svd_residuals(A, d, u, v)
  scale <- eigencore:::svd_backward_scale(norm(A, type = "F"), d)

  expect_equal(diag$left, residuals$left, tolerance = 1e-12)
  expect_equal(diag$right, residuals$right, tolerance = 1e-12)
  expect_equal(diag$combined, residuals$combined, tolerance = 1e-12)
  expect_equal(diag$scale, scale, tolerance = 1e-12)
  expect_equal(diag$backward_error, residuals$combined / scale, tolerance = 1e-12)
  expect_equal(diag$orthogonality, c(U = eigencore:::orthogonality_loss(u), V = eigencore:::orthogonality_loss(v)), tolerance = 1e-12)
  expect_equal(diag$converged, diag$backward_error <= tol)
})

test_that("cached-Av native SVD certificates match native full certificates", {
  set.seed(222)
  A <- matrix(rnorm(35), nrow = 7)
  u <- qr.Q(qr(matrix(rnorm(21), nrow = 7)))
  v <- qr.Q(qr(matrix(rnorm(15), nrow = 5)))
  d <- c(4, 2, 1)
  Av <- A %*% v
  dense <- eigencore:::native_dense_svd_certificate(A, d, u, v)
  dense_cached <- eigencore:::native_dense_svd_certificate_cached_av(A, d, u, v, Av)
  expect_equal(dense_cached$left, dense$left, tolerance = 1e-12)
  expect_equal(dense_cached$right, dense$right, tolerance = 1e-12)
  expect_equal(dense_cached$combined, dense$combined, tolerance = 1e-12)
  expect_equal(dense_cached$backward_error, dense$backward_error, tolerance = 1e-12)

  S <- Matrix::Matrix(A, sparse = TRUE)
  op <- eigencore:::as_operator(S)
  csc <- eigencore:::certify_svd_operator(op, d, u, v)
  csc_cached <- eigencore:::certify_svd_operator_cached_av(op, d, u, v, Av)
  expect_equal(csc_cached$residuals$left, csc$residuals$left, tolerance = 1e-12)
  expect_equal(csc_cached$residuals$right, csc$residuals$right, tolerance = 1e-12)
  expect_equal(csc_cached$backward_error, csc$backward_error, tolerance = 1e-12)
})

test_that("benchmark helpers return timing rows for base and eigencore", {
  A <- diag(c(5, 4, 3, 2, 1))
  eb <- eigencore:::benchmark_eigen_methods(A, k = 2, repeats = 1, include = c("eigencore", "base"))
  sb <- eigencore:::benchmark_svd_methods(A, rank = 2, repeats = 1, include = c("eigencore", "base"))

  expect_equal(vapply(eb, `[[`, character(1), "method"), c("eigencore", "base"))
  expect_equal(vapply(sb, `[[`, character(1), "method"), c("eigencore", "base"))
  expect_true(all(vapply(eb, function(x) is.numeric(x$median_seconds), logical(1))))
  expect_true(all(vapply(sb, function(x) is.numeric(x$median_seconds), logical(1))))
  expect_true(all(vapply(eb, function(x) !is.null(x$norm_bound_type), logical(1))))
  expect_true(all(vapply(sb, function(x) !is.null(x$norm_bound_type), logical(1))))
  expect_true(all(vapply(sb, function(x) !is.null(x$fallback_attempted), logical(1))))
  expect_true(all(vapply(sb, function(x) !is.null(x$fallback_used), logical(1))))
})

test_that("base SVD benchmark adapter keeps singular values conformable", {
  set.seed(502)
  A <- matrix(rnorm(35), nrow = 7)
  fit <- eigencore:::run_svd_method("base", A, rank = 3, tol = 1e-8)
  cert <- eigencore:::certify_svd_operator(as_operator(A), fit$d, fit$u, fit$v, tol = 1e-8)

  expect_equal(length(fit$d), 3L)
  expect_equal(dim(fit$u), c(7L, 3L))
  expect_equal(dim(fit$v), c(5L, 3L))
  expect_true(cert$passed)
})

test_that("matrix-free Golub-Kahan values match base SVD oracle through explicit twin", {
  A <- rbind(diag(c(6, 4, 2)), matrix(0, 2, 3))
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
    },
    metadata = list(frobenius_norm = norm(A, type = "F"))
  )
  fit <- svd_partial(op, rank = 2, method = golub_kahan(max_subspace = 3), seed = 42)
  oracle <- svd(A, nu = 2, nv = 2)

  expect_equal(values(fit), oracle$d[1:2], tolerance = 1e-10)
  expect_true(certificate(fit)$passed)
})
