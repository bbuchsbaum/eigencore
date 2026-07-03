test_that("ARPACK which codes map to exact targets", {
  expect_identical(eigencore:::target_from_which("LA")$kind, "largest")
  expect_identical(eigencore:::target_from_which("SA")$kind, "smallest")
  expect_identical(eigencore:::target_from_which("LM")$kind, "largest_magnitude")
  expect_identical(eigencore:::target_from_which("SM")$kind, "smallest_magnitude")
  expect_identical(eigencore:::target_from_which("LR")$kind, "largest_real")
  expect_identical(eigencore:::target_from_which("SR")$kind, "smallest_real")
  expect_identical(eigencore:::target_from_which("LI")$kind, "largest_imaginary")
  expect_identical(eigencore:::target_from_which("SI")$kind, "smallest_imaginary")
  be <- eigencore:::target_from_which("BE", k = 5)
  expect_identical(be$kind, "both_ends")
  expect_equal(be$value$k_low, 2L)
  expect_equal(be$value$k_high, 3L)
})

test_that("order_indices honours magnitude-based targets", {
  x <- c(-3, 1, -2, 4)
  expect_equal(x[eigencore:::order_indices(x, largest_magnitude())], c(4, -3, -2, 1))
  expect_equal(x[eigencore:::order_indices(x, smallest_magnitude())], c(1, -2, -3, 4))
})

test_that("order_indices distinguishes real and imaginary targets", {
  z <- complex(real = c(1, 3, 2), imaginary = c(4, 0, -2))
  expect_equal(Re(z[eigencore:::order_indices(z, largest_real())]), c(3, 2, 1))
  expect_equal(Re(z[eigencore:::order_indices(z, smallest_real())]), c(1, 2, 3))
  expect_equal(Im(z[eigencore:::order_indices(z, largest_imaginary())]), c(4, 0, -2))
  expect_equal(Im(z[eigencore:::order_indices(z, smallest_imaginary())]), c(-2, 0, 4))
})

test_that("order_indices selects both algebraic ends", {
  x <- c(-5, -2, 0, 3, 9)
  idx <- eigencore:::order_indices(x, both_ends(k_low = 2, k_high = 1))
  expect_equal(x[idx], c(-5, -2, 9))
  expect_equal(eigencore:::target_to_rspectra_which(both_ends(1, 1)), "BE")
})

test_that("eigs_sym BE returns both algebraic ends", {
  A <- diag(c(-7, -3, 0.5, 2, 11))

  fit <- eigs_sym(A, k = 3, which = "BE")

  expect_equal(fit$values, c(-7, 11, 2), tolerance = 1e-10)
  expect_match(paste(fit$diagnostics$plan$reasons, collapse = "\n"), "target: both_ends")
  expect_equal(fit$diagnostics$method, "native dense Hermitian LAPACK fallback")
  expect_true(fit$certificate$passed)
})

test_that("eigs_sym SM selects smallest-magnitude, not smallest-algebraic", {
  A <- diag(c(-5, 0.1, 3, -0.01, 2))
  fit <- eigs_sym(A, k = 2, which = "SM")
  expect_equal(sort(abs(fit$values)), c(0.01, 0.1), tolerance = 1e-10)
})

test_that("complex dense Hermitian inputs use native dense complex certification", {
  A <- matrix(c(1, 1i, -1i, 2), 2, 2)

  fit <- eigs_sym(A, k = 2, which = "LA", tol = 1e-10)

  expect_equal(fit$values, c(2.61803398875, 0.38196601125), tolerance = 1e-10)
  expect_true(is.complex(fit$vectors))
  expect_equal(fit$diagnostics$method, eigencore:::native_dense_complex_hermitian_label())
  expect_true(fit$certificate$passed)
  expect_equal(fit$certificate$norm_bound_type, "frobenius_exact")
  expect_lt(fit$certificate$max_orthogonality_loss, 1e-10)
})

test_that("complex dense nonsymmetric inputs use native dense complex certification", {
  A <- matrix(c(0, 1i, 2, 0), 2, 2)

  fit <- eigs(A, k = 2, which = "LM", tol = 1e-10)

  expect_equal(Mod(fit$values), c(sqrt(2), sqrt(2)), tolerance = 1e-10)
  expect_true(is.complex(fit$values))
  expect_true(is.complex(fit$vectors))
  expect_equal(fit$diagnostics$method, eigencore:::native_dense_complex_general_label())
  expect_equal(fit$certificate$certificate_type, "right_residual_backward_error")
  expect_true(fit$certificate$passed)
  expect_false(fit$certificate$orthogonality_required)
})

test_that("complex dense SVD inputs use native dense complex certification", {
  A <- matrix(c(1 + 1i, 0, 0, 2), 2, 2)

  fit <- svds(A, k = 2L, tol = 1e-10)

  expect_equal(fit$d, c(2, sqrt(2)), tolerance = 1e-10)
  expect_true(is.complex(fit$u))
  expect_true(is.complex(fit$v))
  expect_equal(fit$diagnostics$method, eigencore:::native_dense_complex_svd_label())
  expect_true(fit$certificate$passed)
  expect_equal(fit$certificate$norm_bound_type, "frobenius_exact")
  expect_lt(fit$certificate$max_orthogonality_loss, 1e-10)
})

test_that("complex dense operator adjoints preserve conjugate transpose metadata", {
  A <- matrix(c(1, 1i, -1i, 2), 2, 2)

  op <- as_operator(A)
  adj <- adjoint(op)

  expect_identical(op$dtype, "complex")
  expect_true(isTRUE(op$metadata$native))
  expect_identical(op$metadata$native_operator_kernel, "dense_complex_zgemm")
  expect_equal(adj$metadata$source, Conj(t(A)))
})

test_that("complex ABI contract matches dense and operator certificate semantics", {
  A <- matrix(c(1, 1i, -1i, 2), 2, 2)
  op <- as_operator(A)
  eig <- eigen(A, symmetric = FALSE)
  sv <- svd(A)
  X <- matrix(c(1 + 2i, 2 - 1i, -1i, 3), 2, 2)

  expect_identical(op$dtype, "complex")
  expect_identical(op$metadata$storage, "complex_dense_matrix")
  expect_true(isTRUE(op$metadata$native))
  expect_identical(op$metadata$native_operator_kernel, "dense_complex_zgemm")
  expect_equal(apply_operator(op, X), A %*% X)
  expect_equal(apply_adjoint_operator(op, X), Conj(t(A)) %*% X)
  expect_equal(eigencore:::certificate_gram(eig$vectors),
               Conj(t(eig$vectors)) %*% eig$vectors)

  norm_info <- eigencore:::operator_norm_for_certificate_info(op)
  expect_equal(norm_info$value, norm(A, type = "F"))
  expect_equal(norm_info$norm_bound_type, "frobenius_exact")
  expect_false(norm_info$scale_is_estimate)

  dense_cert <- eigencore:::certify_eigen(A, eig$values, eig$vectors, tol = 1e-10)
  operator_cert <- eigencore:::certify_eigen_operator(
    op, eig$values, eig$vectors, tol = 1e-10
  )
  expect_true(dense_cert$passed)
  expect_true(operator_cert$passed)
  expect_equal(operator_cert$norm_bound_type, "frobenius_exact+identity_exact")
  expect_false(operator_cert$scale_is_estimate)
  expect_equal(operator_cert$residuals, dense_cert$residuals, tolerance = 1e-12)
  expect_equal(operator_cert$backward_error, dense_cert$backward_error,
               tolerance = 1e-12)

  dense_svd_cert <- eigencore:::certify_svd(A, sv$d, sv$u, sv$v, tol = 1e-10)
  operator_svd_cert <- eigencore:::certify_svd_operator(
    op, sv$d, sv$u, sv$v, tol = 1e-10
  )
  expect_true(dense_svd_cert$passed)
  expect_true(operator_svd_cert$passed)
  expect_equal(operator_svd_cert$norm_bound_type, "frobenius_exact")
  expect_false(operator_svd_cert$scale_is_estimate)
  expect_equal(operator_svd_cert$residuals$combined,
               dense_svd_cert$residuals$combined,
               tolerance = 1e-12)
  expect_equal(operator_svd_cert$backward_error,
               dense_svd_cert$backward_error,
               tolerance = 1e-12)
})

test_that("complex matrix-free eigen operators fail with actionable boundary", {
  A <- matrix(c(0, 1i, 2, 0), 2, 2)
  op <- linear_operator(
    dim = c(2, 2),
    apply = function(X, alpha = 1, beta = 0, Y = NULL) {
      out <- alpha * (A %*% X)
      if (!is.null(Y) && beta != 0) out <- out + beta * Y
      out
    },
    apply_adjoint = function(X, alpha = 1, beta = 0, Y = NULL) {
      out <- alpha * (Conj(t(A)) %*% X)
      if (!is.null(Y) && beta != 0) out <- out + beta * Y
      out
    },
    dtype = "complex",
    structure = general()
  )

  expect_error(
    eig_partial(op, k = 2, target = largest_magnitude(), tol = 1e-10),
    "Complex matrix-free eigen operators are future scope"
  )
})

test_that("dense eigs LI routes nonsymmetric complex pairs through native Arnoldi", {
  A <- rbind(
    c(0, -1, 0),
    c(1, 0, 0),
    c(0, 0, 0.5)
  )

  fit <- eigs(A, k = 2, which = "LI", tol = 1e-10)

  expect_equal(fit$values, c(1i, 0.5), tolerance = 1e-10)
  expect_true(is.complex(fit$vectors))
  expect_equal(fit$certificate$certificate_type, "right_residual_backward_error")
  expect_false(fit$certificate$orthogonality_required)
  expect_true(fit$certificate$passed)
  expect_equal(
    fit$diagnostics$method,
    eigencore:::native_refined_arnoldi_label()
  )
  expect_true(fit$diagnostics$restart$native)
  expect_true(fit$diagnostics$restart$ritz_extraction_native)
  expect_true(fit$diagnostics$restart$refined_extraction_native)
  expect_match(fit$diagnostics$warnings, "right residuals certified")
})

test_that("sparse eigs LI routes through native Arnoldi compatibility path", {
  A <- Matrix::bdiag(
    matrix(c(0, -2, 2, 0), 2, 2, byrow = TRUE),
    matrix(c(0, -1, 1, 0), 2, 2, byrow = TRUE),
    Matrix::Diagonal(2, c(0.5, 0.25))
  )
  A <- methods::as(A, "dgCMatrix")

  fit <- eigs(A, k = 2L, which = "LI", tol = 1e-10)

  expect_equal(Im(fit$values), c(2, 1), tolerance = 1e-8)
  expect_true(is.complex(fit$vectors))
  expect_equal(fit$diagnostics$method, eigencore:::native_refined_arnoldi_label())
  expect_true(fit$diagnostics$restart$native)
  expect_true(fit$diagnostics$restart$ritz_extraction_native)
  expect_true(fit$diagnostics$restart$refined_extraction_native)
  expect_equal(fit$certificate$certificate_type, "right_residual_backward_error")
  expect_false(fit$certificate$orthogonality_required)
  expect_true(fit$certificate$passed)
})

test_that("public exports exclude reference and planner internals", {
  exports <- getNamespaceExports("eigencore")

  expect_false(any(startsWith(exports, "reference_")))
  expect_false(any(grepl("(_label|_supported|_or_null)$", exports)))
  expect_false("target_from_which" %in% exports)
  expect_true(all(c("eigs", "eigs_sym", "svds", "both_ends") %in% exports))
})
