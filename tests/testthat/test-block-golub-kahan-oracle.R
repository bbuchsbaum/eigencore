test_that("reference block Golub-Kahan thick restart matches known singular values", {
  s <- c(9, 5, 2, 1, 0.25, 0.1)
  A <- rectangular_with_singular_values(s, m = 14L, n = 9L, seed = 701)

  set.seed(701)
  fit <- eigencore:::reference_block_golub_kahan_thick_restart_svd(
    A,
    rank = 3L,
    target = largest(),
    block = 2L,
    max_subspace = 5L,
    max_restarts = 20L,
    tol = 1e-10
  )

  expect_equal(fit$d, s[1:3], tolerance = 1e-8)
  expect_true(fit$certificate$passed)
  expect_gt(fit$restarts, 0L)
  expect_identical(fit$restart$kind, "block_golub_kahan_thick_restart_reference")
})

test_that("reference block Golub-Kahan supports wide rectangular operators", {
  s <- c(7, 4, 3, 1, 0.5)
  A <- rectangular_with_singular_values(s, m = 9L, n = 14L, seed = 704)

  set.seed(704)
  fit <- eigencore:::reference_block_golub_kahan_thick_restart_svd(
    A,
    rank = 4L,
    target = largest(),
    block = 2L,
    max_subspace = 6L,
    max_restarts = 20L,
    tol = 1e-9
  )
  oracle <- svd(A, nu = 4L, nv = 4L)

  expect_equal(fit$d, oracle$d[1:4], tolerance = 1e-8)
  expect_lt(subspace_distance(fit$u, oracle$u[, 1:4, drop = FALSE]), 1e-5)
  expect_lt(subspace_distance(fit$v, oracle$v[, 1:4, drop = FALSE]), 1e-5)
  expect_true(fit$certificate$passed)
})

test_that("native block Golub-Kahan Ritz kernel matches full-basis SVD slices", {
  s <- c(9, 6, 2, 0.8, 0.2, 0.05, 0.01)
  A <- rectangular_with_singular_values(s, m = 11L, n = 7L, seed = 705)
  V <- diag(1, ncol(A))
  AV <- A %*% V
  V_aug <- cbind(V, matrix(rnorm(ncol(A) * 2L), ncol(A), 2L))
  AV_aug <- cbind(AV, matrix(rnorm(nrow(A) * 2L), nrow(A), 2L))

  top <- eigencore:::native_block_golub_kahan_ritz(
    V_aug,
    AV_aug,
    rank = 3L,
    target = largest(),
    active_cols = ncol(A)
  )
  small <- eigencore:::native_block_golub_kahan_ritz(
    V_aug,
    AV_aug,
    rank = 2L,
    target = smallest(),
    active_cols = ncol(A)
  )
  oracle <- svd(A, nu = 3L, nv = 3L)

  expect_equal(top$d, oracle$d[1:3], tolerance = 1e-8)
  expect_equal(small$d, sort(oracle$d, decreasing = FALSE)[1:2], tolerance = 1e-8)
  expect_equal(top$Avectors, A %*% top$v, tolerance = 1e-10)
  expect_equal(top$Avectors, sweep(top$u, 2L, top$d, `*`), tolerance = 1e-10)

  cert <- eigencore:::certify_svd_operator(
    eigencore:::as_operator(A),
    top$d,
    top$u,
    top$v,
    tol = 1e-8
  )
  expect_true(cert$passed)
})

test_that("native block Golub-Kahan basis cycle certifies dense and CSC full subspaces", {
  s <- c(8, 5, 3, 1, 0.4, 0.1)
  A <- rectangular_with_singular_values(s, m = 10L, n = 6L, seed = 706)
  A_csc <- methods::as(Matrix::Matrix(A, sparse = TRUE), "dgCMatrix")

  for (A_in in list(A, A_csc)) {
    set.seed(706)
    fit <- eigencore:::native_block_golub_kahan_cycle_svd(
      A_in,
      rank = 4L,
      target = largest(),
      block = 2L,
      max_subspace = ncol(A),
      tol = 1e-8
    )

    expect_identical(fit$restart$kind, "block_golub_kahan_native_basis_cycle")
    expect_true(fit$restart$native)
    expect_gte(fit$restart$active_cols, 4L)
    expect_gt(fit$matvecs, 0L)
    expect_gt(fit$restart$ortho_passes, 0L)
    expect_equal(fit$d, s[1:4], tolerance = 1e-8)
    expect_true(fit$certificate$passed)
  }
})

test_that("reference block Golub-Kahan handles clustered singular subspaces", {
  s <- c(10, 10 - 1e-9, 10 - 2e-9, 3, 1, 0.1)
  A <- rectangular_with_singular_values(s, m = 16L, n = 10L, seed = 702)

  set.seed(702)
  fit <- eigencore:::reference_block_golub_kahan_thick_restart_svd(
    A,
    rank = 3L,
    target = largest(),
    block = 2L,
    max_subspace = 6L,
    max_restarts = 30L,
    tol = 1e-8
  )
  oracle <- svd(A, nu = 3L, nv = 3L)

  expect_equal(fit$d, oracle$d[1:3], tolerance = 1e-8)
  expect_lt(subspace_distance(fit$u, oracle$u[, 1:3, drop = FALSE]), 1e-5)
  expect_lt(subspace_distance(fit$v, oracle$v[, 1:3, drop = FALSE]), 1e-5)
  expect_true(fit$certificate$passed)
})

test_that("reference block Golub-Kahan completes exact zero singular triplets", {
  s <- c(8, 4, 1, 0, 0)
  A <- rectangular_with_singular_values(s, m = 12L, n = 9L, seed = 703)

  set.seed(703)
  fit <- eigencore:::reference_block_golub_kahan_thick_restart_svd(
    A,
    rank = 5L,
    target = largest(),
    block = 2L,
    max_subspace = 7L,
    max_restarts = 30L,
    tol = 1e-8
  )

  expect_length(fit$d, 5L)
  expect_equal(fit$d[1:3], s[1:3], tolerance = 1e-8)
  expect_equal(fit$d[4:5], c(0, 0), tolerance = 1e-10)
  expect_true(fit$zero_singular_completion)
  expect_true(fit$certificate$passed)
})
