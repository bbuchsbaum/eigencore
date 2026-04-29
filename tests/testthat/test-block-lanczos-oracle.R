test_that("reference block thick-restart Lanczos matches diagonal spectra", {
  A <- diag(c(9, 5, 2, 1, -1, -3))
  set.seed(101)
  fit <- eigencore:::reference_block_lanczos_thick_restart_hermitian(
    A,
    k = 3L,
    target = largest(),
    block = 2L,
    max_subspace = 5L,
    max_restarts = 20L,
    tol = 1e-10
  )

  expect_equal(fit$values, c(9, 5, 2), tolerance = 1e-8)
  expect_true(certificate(fit)$passed)
  expect_gt(fit$restarts, 0L)
  expect_equal(fit$restart$kind, "block_subspace_iteration_thick_restart_reference")
})

test_that("reference block thick-restart Lanczos handles clustered subspaces", {
  vals <- c(10, 10 - 1e-9, 10 - 2e-9, 4, 2, 1, 0.5)
  A <- symmetric_with_spectrum(vals, seed = 102)
  set.seed(102)
  fit <- eigencore:::reference_block_lanczos_thick_restart_hermitian(
    A,
    k = 3L,
    target = largest(),
    block = 2L,
    max_subspace = 6L,
    max_restarts = 20L,
    tol = 1e-8
  )
  oracle <- eigen(A, symmetric = TRUE)

  expect_equal(fit$values, oracle$values[1:3], tolerance = 1e-8)
  expect_lt(subspace_distance(vectors(fit), oracle$vectors[, 1:3]), 1e-5)
  expect_true(certificate(fit)$passed)
})

test_that("reference block thick-restart Lanczos recovers sparse Laplacian nullspace", {
  n <- 10L
  A <- Matrix::bandSparse(
    n,
    k = c(-1, 0, 1),
    diagonals = list(rep(-1, n - 1L), c(1, rep(2, n - 2L), 1), rep(-1, n - 1L))
  )
  set.seed(103)
  fit <- eigencore:::reference_block_lanczos_thick_restart_hermitian(
    A,
    k = 1L,
    target = smallest(),
    block = 2L,
    max_subspace = n,
    max_restarts = 20L,
    tol = 1e-8
  )

  expect_lt(abs(values(fit)), 1e-8)
  expect_lt(subspace_distance(vectors(fit), matrix(rep(1 / sqrt(n), n), ncol = 1)), 1e-5)
  expect_true(certificate(fit)$passed)
})
