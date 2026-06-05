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

test_that("eigs_sym SM selects smallest-magnitude, not smallest-algebraic", {
  A <- diag(c(-5, 0.1, 3, -0.01, 2))
  fit <- eigs_sym(A, k = 2, which = "SM")
  expect_equal(sort(abs(fit$values)), c(0.01, 0.1), tolerance = 1e-10)
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
    eigencore:::native_arnoldi_label()
  )
  expect_true(fit$diagnostics$restart$native)
  expect_true(fit$diagnostics$restart$ritz_extraction_native)
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
  expect_equal(fit$diagnostics$method, eigencore:::native_arnoldi_label())
  expect_true(fit$diagnostics$restart$native)
  expect_true(fit$diagnostics$restart$ritz_extraction_native)
  expect_equal(fit$certificate$certificate_type, "right_residual_backward_error")
  expect_false(fit$certificate$orthogonality_required)
  expect_true(fit$certificate$passed)
})
