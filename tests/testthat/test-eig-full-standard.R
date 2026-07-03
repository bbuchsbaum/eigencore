# Coverage for the eig_full() standard (B = NULL) surface and its dimension
# validation. Added under review item bd-01KVWRKQGTE4HQE57DDS1W8E1V.

# Order complex spectra by (Re, Im) so two decompositions can be compared
# position-by-position regardless of the order each solver returns.
key_sort <- function(z) z[order(round(Re(z), 6), round(Im(z), 6))]

test_that("eig_full standard real symmetric uses the native Hermitian path", {
  set.seed(1)
  M <- matrix(rnorm(25), 5, 5)
  A <- crossprod(M)

  fit <- eig_full(A)

  expect_identical(fit$method, "native dense Hermitian LAPACK fallback")
  expect_equal(sort(values(fit)), sort(eigen(A, symmetric = TRUE)$values),
               tolerance = 1e-10)
  expect_true(certificate(fit)$passed)
  expect_equal(crossprod(vectors(fit)), diag(5), tolerance = 1e-10)
})

test_that("eig_full standard real general with real spectrum is honestly labeled as a base fallback", {
  A <- matrix(c(2, 0, 0, 1, 3, 0, 4, 5, 6), 3, 3) # upper triangular: eigenvalues 2, 3, 6

  fit <- eig_full(A)

  expect_identical(fit$method, "dense LAPACK general eigen oracle (base fallback)")
  expect_match(paste(fit$warnings, collapse = " "), "base dense general eigen fallback")
  expect_equal(sort(Re(values(fit))), c(2, 3, 6), tolerance = 1e-10)
  expect_true(certificate(fit)$passed)
})

test_that("eig_full standard real general with complex spectrum certifies right residuals", {
  A <- matrix(c(0, -1, 1, 0), 2, 2) # eigenvalues +/- i

  fit <- eig_full(A)

  expect_equal(sort(Im(values(fit))), c(-1, 1), tolerance = 1e-10)
  expect_equal(Re(values(fit)), c(0, 0), tolerance = 1e-10)
  expect_true(certificate(fit)$passed)
  res <- A %*% vectors(fit) - vectors(fit) %*% diag(values(fit))
  expect_lt(max(Mod(res)), 1e-10)
})

test_that("eig_full standard complex Hermitian uses the native complex Hermitian path", {
  H <- matrix(c(2 + 0i, 1 - 1i, 1 + 1i, 3 + 0i), 2, 2) # Hermitian

  fit <- eig_full(H)

  expect_identical(fit$method, eigencore:::native_dense_complex_hermitian_label())
  expect_true(all(abs(Im(values(fit))) < 1e-12)) # Hermitian -> real eigenvalues
  oracle <- eigen(H, symmetric = FALSE)
  expect_equal(sort(Re(values(fit))), sort(Re(oracle$values)), tolerance = 1e-10)
  expect_true(certificate(fit)$passed)
})

test_that("eig_full standard complex general uses the native complex general path", {
  set.seed(3)
  A <- matrix(complex(real = rnorm(9), imaginary = rnorm(9)), 3, 3)

  fit <- eig_full(A)

  expect_identical(fit$method, eigencore:::native_dense_complex_general_label())
  expect_equal(key_sort(values(fit)),
               key_sort(eigen(A, only.values = TRUE)$values),
               tolerance = 1e-8)
  expect_true(certificate(fit)$passed)
  expect_true(is.complex(vectors(fit)))
})

test_that("eig_full rejects non-square and dimension-mismatched inputs", {
  expect_error(eig_full(matrix(1:6, 2, 3)), "square")
  expect_error(eig_full(diag(3), B = diag(2)), "same dimension")
  expect_error(eig_full(diag(3), B = matrix(1:6, 2, 3)), "square")
})
