test_that("clustered Hermitian eigenvalues are certified by residuals and subspaces", {
  values <- c(10, 10 - 1e-9, 10 - 2e-9, 2, 1, 0.5)
  A <- symmetric_with_spectrum(values, seed = 31)
  fit <- eig_partial(A, k = 3, target = largest(), tol = 1e-8)
  oracle <- eigen(A, symmetric = TRUE)

  expect_equal(values(fit), oracle$values[1:3], tolerance = 1e-8)
  expect_lt(subspace_distance(vectors(fit), oracle$vectors[, 1:3]), 1e-6)
  expect_certificate_clean(fit)
})

test_that("nearly repeated singular values are judged by subspace accuracy", {
  s <- c(7, 7 - 1e-9, 7 - 2e-9, 2, 0.5)
  A <- rectangular_with_singular_values(s, m = 9, n = 6, seed = 32)
  fit <- svd_partial(A, rank = 3, target = largest(), tol = 1e-8)
  oracle <- svd(A, nu = 3, nv = 3)

  expect_equal(values(fit), oracle$d[1:3], tolerance = 1e-8)
  expect_lt(subspace_distance(left_vectors(fit), oracle$u[, 1:3]), 1e-6)
  expect_lt(subspace_distance(right_vectors(fit), oracle$v[, 1:3]), 1e-6)
  expect_certificate_clean(fit)
})

test_that("rank-deficient rectangular SVD returns finite certified triplets", {
  A <- rectangular_with_singular_values(c(6, 3, 0, 0), m = 8, n = 5, seed = 33)
  fit <- svd_partial(A, rank = 4, target = largest(), tol = 1e-8)

  expect_equal(values(fit)[1:2], c(6, 3), tolerance = 1e-10)
  expect_true(all(is.finite(values(fit))))
  expect_false(any(is.nan(values(fit))))
  expect_certificate_clean(fit)
})

test_that("graph Laplacian nullspace is recovered without sparse densification", {
  n <- 8L
  A <- Matrix::bandSparse(n, k = c(-1, 0, 1), diagonals = list(rep(-1, n - 1), c(1, rep(2, n - 2), 1), rep(-1, n - 1)))
  fit <- eig_partial(A, k = 1, target = smallest(), method = lanczos(max_subspace = n), seed = 34, tol = 1e-8)

  expect_equal(fit$plan$method, "native scalar thick-restart Hermitian Lanczos")
  expect_lt(abs(values(fit)), 1e-8)
  expect_lt(subspace_distance(vectors(fit), matrix(rep(1 / sqrt(n), n), ncol = 1)), 1e-5)
  expect_certificate_clean(fit)
})

test_that("native LOBPCG certifies clustered smallest dense eigenpairs", {
  values <- c(0, 1e-8, 2e-8, 1, 2, 4, 8, 16)
  A <- symmetric_with_spectrum(values, seed = 71)
  fit <- eig_partial(
    A,
    k = 3,
    target = smallest(),
    method = lobpcg(maxit = 120L),
    seed = 71,
    tol = 1e-8
  )
  oracle <- eigen(A, symmetric = TRUE)
  idx <- order(oracle$values)[1:3]

  expect_equal(fit$method, "native standard Hermitian LOBPCG prototype")
  expect_lt(max(abs(values(fit) - sort(oracle$values)[1:3])), 1e-8)
  expect_lt(subspace_distance(vectors(fit), oracle$vectors[, idx]), 1e-6)
  expect_certificate_clean(fit)
})

test_that("native preconditioned LOBPCG handles Laplacian near-nullspace clusters", {
  n <- 50L
  A <- Matrix::bandSparse(
    n,
    k = c(-1, 0, 1),
    diagonals = list(rep(-1, n - 1), c(1, rep(2, n - 2), 1), rep(-1, n - 1))
  )
  preconditioner <- shifted_tridiagonal_preconditioner(A, shift = 1e-4)
  fit <- eig_partial(
    A,
    k = 4,
    target = smallest(),
    method = lobpcg(maxit = 80L, preconditioner = preconditioner),
    seed = 72,
    tol = 1e-8
  )
  oracle <- eigen(as.matrix(A), symmetric = TRUE)
  idx <- order(oracle$values)[1:4]

  expect_equal(fit$method, "native standard Hermitian LOBPCG prototype")
  expect_true(fit$restart$preconditioner_native)
  expect_lt(max(abs(values(fit) - sort(oracle$values)[1:4])), 1e-8)
  expect_lt(subspace_distance(vectors(fit), oracle$vectors[, idx]), 1e-5)
  expect_certificate_clean(fit)
})

test_that("ill-conditioned generalized SPD problems remain B-orthonormal", {
  A <- diag(c(9, 6, 4, 2, 1))
  B <- diag(c(1, 1e-2, 1e-4, 1e-6, 1e-8))
  fit <- eig_partial(A, B = B, k = 3, target = smallest(), tol = 1e-8)

  expect_true(certificate(fit)$passed)
  expect_equal(crossprod(vectors(fit), B %*% vectors(fit)), diag(3), tolerance = 1e-8)
  expect_lt(max(abs(A %*% vectors(fit) - B %*% vectors(fit) %*% diag(values(fit)))), 1e-7)
})

test_that("non-normal nonsymmetric matrices with real spectra are certified", {
  V <- matrix(c(
    1, 1, 1,
    0, 1e-3, 1,
    0, 0, 1
  ), nrow = 3, byrow = TRUE)
  A <- V %*% diag(c(5, 3, 1)) %*% solve(V)
  fit <- eig_partial(A, k = 2, target = largest_magnitude(), tol = 1e-8)

  expect_equal(values(fit), c(5, 3), tolerance = 1e-8)
  expect_certificate_clean(fit)
})

test_that("poorly scaled SVD inputs remain finite and certified", {
  A <- diag(c(1e6, 1e2, 1, 1e-2, 1e-6))
  A <- diag(c(1e-3, 1, 1e3, 1e-2, 1e2)) %*% A
  fit <- svd_partial(A, rank = 3, target = largest(), tol = 1e-8)
  oracle <- svd(A, nu = 3, nv = 3)

  expect_equal(values(fit), oracle$d[1:3], tolerance = 1e-8)
  expect_true(all(is.finite(certificate(fit)$backward_error)))
  expect_certificate_clean(fit)
})

test_that("near-singular shift-invert requests fail loudly", {
  A <- diag(c(-1e-10, 0, 1e-10))
  expect_error(
    eig_partial(A, k = 1, target = nearest(0), method = shift_invert(0)),
    "singular or near-singular"
  )
})
