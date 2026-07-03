test_that("generalized pencil values classify finite infinite and undefined pairs", {
  pencil <- eigencore:::generalized_pencil_values(
    alpha = c(2, 3, 0, 1e-13),
    beta = c(1, 0, 0, 1e-14),
    tol = 1e-10
  )

  expect_equal(pencil$classification, c("finite", "infinite", "undefined", "undefined"))
  expect_equal(pencil$values[[1L]], 2)
  expect_true(is.infinite(pencil$values[[2L]]))
  expect_true(is.na(pencil$values[[3L]]))
  expect_equal(pencil$finite, c(TRUE, FALSE, FALSE, FALSE))
  expect_equal(pencil$infinite, c(FALSE, TRUE, FALSE, FALSE))
  expect_equal(pencil$undefined, c(FALSE, FALSE, TRUE, TRUE))
  expect_equal(pencil$tolerance_policy, "per_pair_magnitude")

  empty <- eigencore:::generalized_pencil_values(numeric(), numeric())
  expect_length(empty$values, 0L)
  expect_identical(empty$classification, character())
  expect_identical(empty$finite, logical())
  expect_identical(empty$infinite, logical())
  expect_identical(empty$undefined, logical())
})

test_that("norm-scaled pencil classification survives joint rescaling", {
  # Uniformly tiny pencil: per-pair max(1, |alpha|, |beta|) scaling would
  # classify every pair as undefined because both coordinates fall below
  # tol * 1. Norm-scaled classification compares against the pencil norms
  # instead and keeps the labels of the unscaled pencil.
  alpha <- c(2, 3, 0) * 1e-12
  beta <- c(1, 0, 0) * 1e-12

  per_pair <- eigencore:::generalized_pencil_values(alpha, beta, tol = 1e-10)
  expect_equal(per_pair$classification, rep("undefined", 3L))

  scaled <- eigencore:::generalized_pencil_values(
    alpha, beta,
    tol = 1e-10,
    norm_A = 3e-12, norm_B = 1e-12
  )
  expect_equal(scaled$classification, c("finite", "infinite", "undefined"))
  expect_equal(scaled$tolerance_policy, "pencil_norm_scaled")
  expect_equal(scaled$values[[1L]], 2)
  expect_equal(scaled$norm_A, 3e-12)
  expect_equal(scaled$norm_B, 1e-12)
  expect_true(all(scaled$alpha_threshold >= 0))
  expect_true(all(scaled$beta_threshold >= 0))
})

test_that("norm-scaled classification falls back to per-pair for invalid norms", {
  alpha <- c(2, 0)
  beta <- c(1, 0)

  for (bad in list(
    list(norm_A = NULL, norm_B = 1),
    list(norm_A = NA_real_, norm_B = 1),
    list(norm_A = Inf, norm_B = 1),
    list(norm_A = -1, norm_B = 1)
  )) {
    pencil <- eigencore:::generalized_pencil_values(
      alpha, beta,
      norm_A = bad$norm_A, norm_B = bad$norm_B
    )
    expect_equal(pencil$tolerance_policy, "per_pair_magnitude")
    expect_equal(pencil$classification, c("finite", "undefined"))
  }
})

test_that("generalized pencil dense certificate certifies finite pairs only", {
  A <- diag(c(2, 5, 7))
  B <- diag(3)
  V <- diag(3)

  cert <- eigencore:::certify_dense_generalized_pencil(
    A,
    B,
    alpha = c(2, 3, 0),
    beta = c(1, 0, 0),
    vectors = V,
    tol = 1e-10,
    beta_tol = 1e-12
  )

  expect_s3_class(cert, "eigencore_certificate")
  expect_false(cert$passed)
  expect_equal(cert$certificate_type, "generalized_pencil_right_residual_backward_error")
  expect_equal(cert$converged, c(TRUE, FALSE, FALSE))
  expect_equal(cert$failed_indices, c(2L, 3L))
  expect_equal(cert$residuals[[1L]], 0, tolerance = 1e-14)
  expect_true(all(is.infinite(cert$residuals[2:3])))
  expect_true(all(is.infinite(cert$backward_error[2:3])))
  expect_false(cert$orthogonality_required)
  expect_match(paste(cert$notes, collapse = " "), "infinite generalized eigenvalues")
  expect_match(paste(cert$notes, collapse = " "), "undefined generalized eigenvalues")
})

test_that("generalized pencil operator certificate withholds passed under an estimated scale", {
  diag_a <- c(2, 6, 12)
  diag_b <- c(1, 2, 3)
  # Matrix-free diagonal operators with no norm metadata force the certificate
  # onto the Hutchinson Frobenius estimate, so scale_is_estimate is TRUE even
  # though the eigenpairs are exact (residuals ~ 0). This exercises the
  # passed-withholding branch that the dense/exact pencil tests never hit.
  matrix_free_diag <- function(d) {
    linear_operator(
      dim = c(length(d), length(d)),
      apply = function(X, alpha = 1, beta = 0, Y = NULL) {
        out <- d * X
        if (is.null(Y)) alpha * out else alpha * out + beta * Y
      },
      structure = hermitian(),
      name = "matrix_free_diag"
    )
  }

  cert <- eigencore:::certify_generalized_pencil_operator(
    matrix_free_diag(diag_a),
    matrix_free_diag(diag_b),
    alpha = diag_a,
    beta = diag_b,
    vectors = diag(3),
    tol = 1e-10
  )

  expect_true(all(cert$converged))
  expect_lt(max(cert$backward_error), 1e-10)
  expect_true(cert$scale_is_estimate)
  expect_false(cert$passed)
  expect_match(paste(cert$notes, collapse = " "), "stochastic norm estimate")
})

test_that("generalized pencil dense and operator certificates share scale contract", {
  A <- diag(c(2, 6, 12))
  B <- diag(c(1, 2, 3))
  V <- diag(3)
  alpha <- c(2, 6, 12)
  beta <- c(1, 2, 3)

  dense <- eigencore:::certify_dense_generalized_pencil(
    A,
    B,
    alpha = alpha,
    beta = beta,
    vectors = V,
    tol = 1e-10
  )
  operator <- eigencore:::certify_generalized_pencil_operator(
    as_operator(A),
    as_operator(B),
    alpha = alpha,
    beta = beta,
    vectors = V,
    tol = 1e-10
  )

  expect_true(dense$passed)
  expect_true(operator$passed)
  expect_equal(dense$residuals, operator$residuals, tolerance = 1e-12)
  expect_equal(dense$backward_error, operator$backward_error, tolerance = 1e-12)
  expect_equal(dense$scale, operator$scale, tolerance = 1e-12)
  expect_equal(dense$norm_bound_type, operator$norm_bound_type)
  expect_equal(dense$scale_is_estimate, operator$scale_is_estimate)
})

test_that("current generalized SPD paths share original-coordinate contract", {
  A <- diag(c(1, 4, 9, 16, 25))
  B <- diag(c(1, 2, 3, 4, 5))
  Bop <- as_operator(B)

  expect_current_contract <- function(fit, expected_values, expected_method) {
    expect_s3_class(fit, "eigencore_eigen_result")
    expect_equal(values(fit), expected_values, tolerance = 1e-8)
    expect_identical(fit$method, fit$plan$method)
    expect_identical(fit$method, expected_method)
    expect_true(all(c(
      "values", "vectors", "method", "plan", "certificate", "warnings",
      "residuals", "backward_error"
    ) %in% names(fit)))
    expect_true(certificate(fit)$passed)
    expect_equal(
      crossprod(vectors(fit), B %*% vectors(fit)),
      diag(length(expected_values)),
      tolerance = 1e-8
    )

    direct <- eigencore:::certify_eigen_operator(
      as_operator(A),
      values(fit),
      vectors(fit),
      Bop = Bop,
      tol = 1e-8
    )
    expect_equal(certificate(fit)$residuals, direct$residuals, tolerance = 1e-10)
    expect_equal(certificate(fit)$backward_error, direct$backward_error,
                 tolerance = 1e-10)
  }

  dense <- eig_partial(A, B = B, k = 2L, target = smallest(), tol = 1e-10)
  expect_current_contract(
    dense,
    c(1, 2),
    "native dense generalized SPD LAPACK fallback"
  )

  lobpcg <- eig_partial(
    A,
    B = B,
    k = 2L,
    target = smallest(),
    method = lobpcg(maxit = 80L),
    seed = 44,
    tol = 1e-9
  )
  expect_current_contract(
    lobpcg,
    c(1, 2),
    eigencore:::native_generalized_lobpcg_label()
  )

  lanczos <- eig_partial(
    A,
    B = B,
    k = 2L,
    target = smallest(),
    method = lanczos(max_subspace = 5L),
    seed = 44,
    tol = 1e-9
  )
  expect_current_contract(
    lanczos,
    c(1, 2),
    eigencore:::native_generalized_lanczos_label()
  )
  expect_true(lanczos$generalized)
  expect_true(lanczos$restart$native)

  shifted <- eig_partial(
    A,
    B = B,
    k = 2L,
    target = nearest(2.1),
    method = shift_invert(2.1),
    tol = 1e-9
  )
  expect_current_contract(
    shifted,
    c(2, 3),
    eigencore:::native_dense_generalized_shift_invert_label()
  )
  expect_true(shifted$restart$generalized)
  expect_true(shifted$transform$factorization_cache$native)
})
