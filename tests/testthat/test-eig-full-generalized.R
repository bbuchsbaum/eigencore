test_that("eig_full returns certified dense real generalized SPD results", {
  A <- diag(c(1, 4, 9))
  B <- diag(c(1, 2, 3))

  fit <- eig_full(A, B = B, tol = 1e-10)
  oracle <- eigencore:::dense_generalized_spd_eigen(A, B)

  expect_s3_class(fit, "eigencore_eigen_result")
  expect_identical(fit$method, eigencore:::native_dense_generalized_spd_full_label())
  expect_identical(fit$plan$method, fit$method)
  expect_equal(values(fit), oracle$values, tolerance = 1e-10)
  expect_equal(values(fit), c(1, 2, 3), tolerance = 1e-10)
  expect_true(certificate(fit)$passed)
  expect_true(fit$generalized)
  expect_equal(fit$classification, rep("finite", 3L))
  expect_equal(crossprod(vectors(fit), B %*% vectors(fit)), diag(3),
               tolerance = 1e-10)
  expect_equal(abs(crossprod(vectors(fit), B %*% oracle$vectors)), diag(3),
               tolerance = 1e-10)
})

test_that("eig_full exposes alpha beta and beta-zero classification for real pencils", {
  A <- diag(c(2, 3, 0))
  B <- diag(c(1, 0, 0))

  fit <- eig_full(A, B = B, structure = general(), tol = 1e-10)

  expect_identical(fit$method, eigencore:::native_dense_generalized_pencil_full_label())
  expect_true(is.complex(fit$alpha))
  expect_true(is.complex(fit$beta))
  expect_equal(Re(fit$alpha), c(2, 3, 0), tolerance = 1e-12)
  expect_equal(Re(fit$beta), c(1, 0, 0), tolerance = 1e-12)
  expect_equal(fit$classification, c("finite", "infinite", "undefined"))
  expect_equal(values(fit)[[1L]], 2 + 0i, tolerance = 1e-12)
  expect_true(is.infinite(values(fit)[[2L]]))
  expect_true(is.na(values(fit)[[3L]]))
  expect_false(certificate(fit)$passed)
  expect_equal(certificate(fit)$converged, c(TRUE, FALSE, FALSE))
  expect_match(paste(certificate(fit)$notes, collapse = " "), "beta equal to zero")
})

test_that("eig_full supports dense complex Hermitian SPD pencils", {
  A <- diag(as.complex(c(2, 6)))
  B <- diag(as.complex(c(1, 2)))

  fit <- eig_full(A, B = B, tol = 1e-10)

  expect_identical(fit$method, eigencore:::native_dense_generalized_spd_full_label())
  expect_equal(values(fit), c(2, 3), tolerance = 1e-10)
  expect_true(certificate(fit)$passed)
  expect_true(is.complex(vectors(fit)))
  expect_equal(
    Conj(t(vectors(fit))) %*% B %*% vectors(fit),
    matrix(as.complex(diag(2)), nrow = 2),
    tolerance = 1e-10
  )
})

test_that("eig_full supports dense complex general pencils", {
  A <- diag(c(1 + 1i, 4 + 2i))
  B <- diag(as.complex(c(1, 2)))

  fit <- eig_full(A, B = B, structure = general(), tol = 1e-10)

  expect_identical(fit$method, eigencore:::native_dense_generalized_pencil_full_label())
  expect_equal(values(fit), c(1 + 1i, 2 + 1i), tolerance = 1e-10)
  expect_equal(fit$classification, rep("finite", 2L))
  expect_true(certificate(fit)$passed)
  expect_true(is.complex(vectors(fit)))
})

test_that("eig_full rejects sparse generalized inputs instead of densifying", {
  A <- Matrix::Diagonal(x = c(1, 4, 9))
  B <- Matrix::Diagonal(x = c(1, 2, 3))

  expect_error(
    eig_full(A, B = B, allow_dense_fallback = "always"),
    "must be a base dense matrix"
  )
})

test_that("eig_full reconstructs complex-conjugate eigenpairs for a real general pencil", {
  # Rotation-like A: eigenvalues +/- i. With B = I these are the generalized
  # eigenvalues and the real DGGEV path must reconstruct the complex-conjugate
  # eigenvectors from the packed real Schur vectors. This exercises the
  # conjugate-pair branch that every diagonal-only test above misses.
  A <- matrix(c(0, -1, 1, 0), 2, 2) # [[0, 1], [-1, 0]] in column-major
  B <- diag(2)

  fit <- eig_full(A, B = B, structure = general(), tol = 1e-10)

  expect_identical(fit$method, eigencore:::native_dense_generalized_pencil_full_label())
  expect_equal(fit$classification, rep("finite", 2L))
  expect_equal(sort(Im(values(fit))), c(-1, 1), tolerance = 1e-10)
  expect_equal(Re(values(fit)), c(0, 0), tolerance = 1e-10)
  expect_true(certificate(fit)$passed)

  V <- vectors(fit)
  expect_true(is.complex(V))
  conjugate_pair <- isTRUE(all.equal(V[, 1], Conj(V[, 2]))) ||
    isTRUE(all.equal(V[, 2], Conj(V[, 1])))
  expect_true(conjugate_pair)
  res <- A %*% V - (B %*% V) %*% diag(values(fit))
  expect_lt(max(Mod(res)), 1e-10)
})

test_that("eig_full real general pencil with SPD B matches the B^{-1}A spectrum", {
  set.seed(11)
  M <- matrix(rnorm(16), 4, 4)
  B <- crossprod(matrix(rnorm(16), 4, 4)) + diag(4)

  fit <- eig_full(M, B = B, structure = general(), tol = 1e-8)
  ref <- eigen(solve(B) %*% M, only.values = TRUE)$values

  ks <- function(z) z[order(round(Re(z), 6), round(Im(z), 6))]
  expect_equal(ks(values(fit)), ks(as.complex(ref)), tolerance = 1e-6)
  expect_true(all(fit$classification == "finite"))
  expect_true(certificate(fit)$passed)
})

test_that("eig_full complex non-diagonal general pencil matches the B^{-1}A spectrum", {
  set.seed(21)
  A <- matrix(complex(real = rnorm(9), imaginary = rnorm(9)), 3, 3)
  B <- matrix(complex(real = rnorm(9), imaginary = rnorm(9)), 3, 3) + diag(3) * 3

  fit <- eig_full(A, B = B, structure = general(), tol = 1e-8)
  ref <- eigen(solve(B) %*% A, only.values = TRUE)$values

  expect_identical(fit$method, eigencore:::native_dense_generalized_pencil_full_label())
  ks <- function(z) z[order(round(Re(z), 6), round(Im(z), 6))]
  expect_equal(ks(values(fit)), ks(ref), tolerance = 1e-6)
  expect_true(all(fit$classification == "finite"))
  expect_true(certificate(fit)$passed)
  expect_true(is.complex(vectors(fit)))
})

test_that("eig_full rejects dimension-mismatched dense generalized inputs", {
  expect_error(eig_full(diag(3), B = diag(2)), "same dimension")
  expect_error(eig_full(matrix(1:6, 2, 3), B = diag(2)), "square")
})
