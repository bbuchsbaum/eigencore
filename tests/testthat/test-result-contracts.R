expect_eigen_result_contract <- function(fit, expected_method = NULL) {
  expect_s3_class(fit, "eigencore_eigen_result")
  expect_equal(length(values(fit)), fit$requested)
  if (!is.null(vectors(fit))) {
    expect_equal(ncol(vectors(fit)), fit$requested)
  }
  expect_s3_class(certificate(fit), "eigencore_certificate")
  expect_equal(fit$residuals, certificate(fit)$residuals)
  expect_equal(fit$backward_error, certificate(fit)$backward_error)
  expect_equal(fit$nconv, sum(certificate(fit)$converged))
  expect_identical(fit$method, fit$plan$method)
  if (!is.null(expected_method)) {
    expect_identical(fit$method, expected_method)
  }

  diag <- diagnostics(fit)
  expect_named(
    diag,
    c(
      "residuals", "backward_error", "orthogonality", "nconv",
      "iterations", "matvecs", "preconditioner_calls",
      "convergence_history", "restart", "stage_seconds",
      "preconditioner", "locked", "method", "plan", "warnings"
    )
  )
  expect_identical(diag$method, fit$method)
  expect_identical(diag$plan, fit$plan)
  expect_equal(diag$nconv, fit$nconv)
  expect_equal(diag$residuals, fit$residuals)
  expect_equal(diag$backward_error, fit$backward_error)
  expect_true(is.numeric(diag$stage_seconds))
}

expect_svd_result_contract <- function(fit, expected_method = NULL) {
  expect_s3_class(fit, "eigencore_svd_result")
  expect_equal(values(fit), fit$d)
  expect_equal(length(values(fit)), fit$requested)
  if (!is.null(left_vectors(fit))) {
    expect_equal(ncol(left_vectors(fit)), fit$requested)
  }
  if (!is.null(right_vectors(fit))) {
    expect_equal(ncol(right_vectors(fit)), fit$requested)
  }
  expect_s3_class(certificate(fit), "eigencore_certificate")
  expect_equal(fit$residuals, certificate(fit)$residuals)
  expect_equal(fit$backward_error, certificate(fit)$backward_error)
  expect_equal(fit$nconv, sum(certificate(fit)$converged))
  expect_identical(fit$method, fit$plan$method)
  if (!is.null(expected_method)) {
    expect_identical(fit$method, expected_method)
  }

  diag <- diagnostics(fit)
  expect_named(
    diag,
    c(
      "residuals", "backward_error", "orthogonality", "nconv",
      "iterations", "matvecs", "preconditioner_calls",
      "convergence_history", "restart", "stage_seconds",
      "preconditioner", "locked", "method", "plan", "warnings"
    )
  )
  expect_identical(diag$method, fit$method)
  expect_identical(diag$plan, fit$plan)
  expect_equal(diag$nconv, fit$nconv)
  expect_equal(diag$residuals, fit$residuals)
  expect_equal(diag$backward_error, fit$backward_error)
  expect_true(is.numeric(diag$stage_seconds))
}

test_that("eigen result diagnostics are stable across current public paths", {
  dense <- eig_partial(diag(c(5, 3, 1)), k = 2L)
  expect_eigen_result_contract(dense, "native dense Hermitian LAPACK fallback")
  expect_true(certificate(dense)$passed)

  sparse <- eig_partial(
    Matrix::sparseMatrix(
      i = 1:4,
      j = 1:4,
      x = c(9, 7, 4, 1),
      dims = c(4L, 4L)
    ),
    k = 2L
  )
  expect_eigen_result_contract(sparse, "native scalar thick-restart Hermitian Lanczos")
  expect_true(certificate(sparse)$passed)
  expect_equal(sparse$restart$kind, "thick_restart")
  expect_equal(diagnostics(sparse)$restart, sparse$restart)

  A <- diag(c(1, 4, 9, 16))
  B <- diag(c(1, 2, 3, 4))
  generalized <- eig_partial(
    A,
    B = B,
    k = 2L,
    target = smallest(),
    method = lobpcg(maxit = 50L),
    seed = 41
  )
  expect_eigen_result_contract(generalized, eigencore:::native_generalized_lobpcg_label())
  expect_true(certificate(generalized)$passed)
  expect_true(generalized$restart$generalized)
  expect_true(generalized$restart$native)
  expect_equal(diagnostics(generalized)$preconditioner_calls, generalized$preconditioner_calls)

  shifted <- eig_partial(
    diag(c(1, 3, 7)),
    k = 1L,
    target = nearest(2.8),
    method = shift_invert(2.8)
  )
  expect_eigen_result_contract(shifted, eigencore:::native_dense_shift_invert_label())
  expect_true(certificate(shifted)$passed)
  expect_equal(shifted$transform$certification$problem, "original")
  expect_equal(diagnostics(shifted)$restart, shifted$restart)

  nonsymmetric <- eig_partial(
    rbind(c(0, -1), c(1, 0)),
    k = 2L,
    target = largest_imaginary()
  )
  expect_eigen_result_contract(
    nonsymmetric,
    eigencore:::native_arnoldi_label()
  )
  expect_equal(nonsymmetric$certificate$certificate_type, "right_residual_backward_error")
  expect_false(nonsymmetric$certificate$orthogonality_required)
  expect_true(certificate(nonsymmetric)$passed)
  expect_true(nonsymmetric$restart$native)
  expect_true(nonsymmetric$restart$ritz_extraction_native)
})

test_that("SVD result diagnostics are stable across current public paths", {
  dense <- svd_partial(diag(c(6, 3, 1)), rank = 2L)
  expect_svd_result_contract(dense, "native dense LAPACK SVD fallback")
  expect_true(certificate(dense)$passed)

  randomized <- svd_partial(
    diag(c(6, 3, 1, 0)),
    rank = 2L,
    method = randomized(oversample = 2L, n_iter = 0L, refine = TRUE),
    seed = 42
  )
  expect_svd_result_contract(randomized, "reference randomized SVD prototype")
  expect_true(certificate(randomized)$passed)
  expect_equal(randomized$restart$kind, "randomized_range_finder")
  expect_equal(diagnostics(randomized)$restart, randomized$restart)
})

test_that("compatibility shims expose the same diagnostics contract", {
  sym <- eigs_sym(diag(c(5, 3, 1)), k = 2L)
  expect_named(sym, c("values", "vectors", "nconv", "niter", "nops", "certificate", "diagnostics"))
  expect_equal(sym$diagnostics$method, "native dense Hermitian LAPACK fallback")
  expect_s3_class(sym$certificate, "eigencore_certificate")
  expect_identical(sym$diagnostics$plan$method, sym$diagnostics$method)

  general <- eigs(rbind(c(0, -1), c(1, 0)), k = 2L, which = "LI")
  expect_equal(general$diagnostics$method, eigencore:::native_arnoldi_label())
  expect_equal(general$certificate$certificate_type, "right_residual_backward_error")
  expect_false(general$certificate$orthogonality_required)
  expect_true(general$diagnostics$restart$native)
  expect_true(general$diagnostics$restart$ritz_extraction_native)

  sv <- svds(diag(c(6, 3, 1)), k = 2L)
  expect_named(sv, c("d", "u", "v", "nconv", "niter", "nops", "certificate", "diagnostics"))
  expect_equal(sv$diagnostics$method, "native dense LAPACK SVD fallback")
  expect_s3_class(sv$certificate, "eigencore_certificate")
  expect_identical(sv$diagnostics$plan$method, sv$diagnostics$method)
})

test_that("matrix-free nonsymmetric Arnoldi result diagnostics are shaped", {
  A <- rbind(
    c(4, 1, 0),
    c(0, 2, 1),
    c(0, 0, -1)
  )
  op <- linear_operator(
    dim = dim(A),
    apply = function(X, alpha = 1, beta = 0, Y = NULL) {
      out <- alpha * (A %*% X)
      if (!is.null(Y) && beta != 0) out <- out + beta * Y
      out
    },
    structure = general(),
    metadata = list(frobenius_norm = sqrt(sum(A^2)))
  )
  fit <- eig_partial(op, k = 2L, target = largest_real(), seed = 13)
  diag <- diagnostics(fit)

  expect_equal(diag$method, eigencore:::reference_arnoldi_label())
  expect_true(is.list(diag$restart))
  expect_identical(diag$restart$kind, "reference_arnoldi")
  expect_true(fit$certificate$passed)
})
