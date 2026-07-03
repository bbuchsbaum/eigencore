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

test_that("eig_full classifies singular B with multiple infinite eigenvalues", {
  A <- diag(c(2, 3, 5, 7))
  B <- diag(c(1, 0, 0, 1))

  fit <- eig_full(A, B = B, structure = general(), tol = 1e-10)

  expect_identical(fit$classification, c("finite", "infinite", "infinite", "finite"))
  expect_equal(sum(fit$infinite), 2L)
  expect_equal(Re(values(fit)[fit$finite]), c(2, 7), tolerance = 1e-12)
  expect_true(all(is.infinite(values(fit)[fit$infinite])))

  V <- vectors(fit)
  homogeneous_residual <- sweep(A %*% V, 2L, fit$beta, `*`) -
    sweep(B %*% V, 2L, fit$alpha, `*`)
  expect_lt(max(Mod(homogeneous_residual)), 1e-10)
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

test_that("eig_full returns certified left generalized eigenvectors for real pencils", {
  set.seed(31)
  A <- matrix(rnorm(25), 5, 5)
  B <- crossprod(matrix(rnorm(25), 5, 5)) + diag(5)

  fit <- eig_full(A, B = B, structure = general(), tol = 1e-8)

  W <- left_vectors(fit)
  expect_true(is.complex(W))
  expect_equal(dim(W), c(5L, 5L))
  # Left generalized eigenvectors satisfy w^H A = lambda w^H B.
  left_residual <- Conj(t(W)) %*% A - diag(values(fit)) %*% (Conj(t(W)) %*% B)
  expect_lt(max(Mod(left_residual)), 1e-10)
  expect_true(certificate(fit)$passed)
})

test_that("eig_full returns left generalized eigenvectors for complex pencils", {
  set.seed(32)
  A <- matrix(complex(real = rnorm(9), imaginary = rnorm(9)), 3, 3)
  B <- matrix(complex(real = rnorm(9), imaginary = rnorm(9)), 3, 3) + diag(3) * 3

  fit <- eig_full(A, B = B, structure = general(), tol = 1e-8)

  W <- left_vectors(fit)
  left_residual <- Conj(t(W)) %*% A - diag(values(fit)) %*% (Conj(t(W)) %*% B)
  expect_lt(max(Mod(left_residual)), 1e-10)
  # ZGGEV path: left vectors are available, conditioning diagnostics are not.
  expect_false(fit$conditioning$available)
  expect_match(fit$conditioning$note, "ZGGEVX|DGGEVX|expert")
})

test_that("eig_full omits left vectors and residual certificate when vectors = FALSE", {
  set.seed(33)
  A <- matrix(rnorm(16), 4, 4)
  B <- diag(4)

  fit <- eig_full(A, B = B, structure = general(), vectors = FALSE)

  expect_null(vectors(fit))
  expect_null(fit$left_vectors)
  expect_false(certificate(fit)$passed)
})

test_that("real pencil conditioning diagnostics come from DGGEVX", {
  set.seed(34)
  A <- matrix(rnorm(16), 4, 4)
  B <- crossprod(matrix(rnorm(16), 4, 4)) + diag(4)

  fit <- eig_full(A, B = B, structure = general(), tol = 1e-8)

  cond <- fit$conditioning
  expect_true(cond$available)
  expect_match(cond$source, "dggevx")
  expect_length(cond$rconde, 4L)
  expect_length(cond$rcondv, 4L)
  # DGGEVX reciprocal condition numbers are nonnegative but, unlike the
  # standard-eigenproblem convention, are not bounded by 1: they measure
  # sensitivity relative to the chordal metric on (alpha, beta).
  expect_true(all(cond$rconde >= 0))
  expect_true(all(cond$rcondv >= 0))
  expect_gt(cond$abnrm, 0)
  expect_gt(cond$bbnrm, 0)
})

test_that("conditioning diagnostics flag sensitivity that right residuals miss", {
  # Defective pencil: A is a Jordan block, B = I. The double eigenvalue 1 is
  # infinitely ill-conditioned, so rconde collapses to machine epsilon even
  # though the right-residual backward errors are at machine precision and
  # the certificate passes.
  A <- matrix(c(1, 0, 1, 1), 2, 2)
  B <- diag(2)

  fit <- eig_full(A, B = B, structure = general(), tol = 1e-8)
  well <- eig_full(diag(c(1, 2)), B = diag(2), structure = general(), tol = 1e-8)

  expect_true(certificate(fit)$passed)
  expect_lt(max(fit$certificate$backward_error), 1e-12)
  expect_lt(max(fit$conditioning$rconde), 1e-12)
  expect_gt(min(well$conditioning$rconde), 0.1)
})

test_that("alpha/beta classification is invariant under joint pencil rescaling", {
  set.seed(35)
  A <- matrix(rnorm(16), 4, 4)
  B <- crossprod(matrix(rnorm(16), 4, 4)) + diag(4)

  base <- eig_full(A, B = B, structure = general())
  expect_true(all(base$classification == "finite"))

  for (c_scale in c(1e-12, 1e-8, 1e8, 1e12)) {
    scaled <- eig_full(c_scale * A, B = c_scale * B, structure = general())
    expect_identical(scaled$classification, base$classification)
  }
})

test_that("scaled pencil with true infinite eigenvalues keeps beta-zero labels", {
  A <- diag(c(2, 3, 5))
  B <- diag(c(1, 1, 0))

  for (c_scale in c(1e-10, 1, 1e10)) {
    fit <- eig_full(c_scale * A, B = c_scale * B, structure = general())
    expect_identical(fit$classification, c("finite", "finite", "infinite"))
  }
})

test_that("near-singular B classification respects the norm-scaled boundary", {
  # B is nonsingular with one tiny diagonal entry. When beta stays above
  # tol * norm(B) the huge eigenvalue is still classified as finite; once
  # beta falls below the threshold it is reported as infinite rather than
  # returned as an untrustworthy enormous finite number.
  A <- diag(c(2, 3, 4))

  above <- eig_full(A, B = diag(c(1, 1, 1e-6)), structure = general())
  expect_identical(above$classification, rep("finite", 3L))
  expect_equal(sort(Mod(values(above)))[3L], 4e6, tolerance = 1e-6)

  below <- eig_full(A, B = diag(c(1, 1, 1e-12)), structure = general())
  expect_identical(below$classification, c("finite", "finite", "infinite"))
})

test_that("clustered finite eigenvalues near an infinite one stay separated", {
  A <- diag(c(1, 1 + 1e-9, 7))
  B <- diag(c(1, 1, 0))

  fit <- eig_full(A, B = B, structure = general())

  expect_identical(fit$classification, c("finite", "finite", "infinite"))
  finite_vals <- Re(values(fit)[fit$finite])
  expect_equal(sort(finite_vals), c(1, 1 + 1e-9), tolerance = 1e-12)
})

test_that("alpha_beta exposes the classification policy metadata", {
  fit <- eig_full(diag(c(2, 3, 0)), B = diag(c(1, 0, 0)), structure = general())

  ab <- alpha_beta(fit)
  expect_identical(ab$classification, c("finite", "infinite", "undefined"))
  policy <- ab$classification_policy
  expect_identical(policy$policy, "pencil_norm_scaled")
  expect_true(policy$tolerance > 0)
  expect_true(policy$norm_A > 0)
  expect_true(policy$norm_B > 0)
  expect_true(all(policy$alpha_threshold >= 0))
  expect_true(all(policy$beta_threshold >= 0))
})

test_that("generalized_schur classification is invariant under joint rescaling", {
  A <- diag(c(2, 3, 5))
  B <- diag(c(1, 1, 0))

  for (c_scale in c(1e-10, 1, 1e10)) {
    qz <- generalized_schur(c_scale * A, c_scale * B)
    expect_identical(qz$classification, c("finite", "finite", "infinite"))
    expect_identical(qz$classification_policy$policy, "pencil_norm_scaled")
  }
})
