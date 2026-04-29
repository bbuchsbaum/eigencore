test_that("native block thick-restart production path restarts and locks separated pairs", {
  set.seed(7)
  A <- diag(c(50, 40, 30, 20, 10, 5, 4, 3, 2, 1))

  fit <- eig_partial(
    A,
    k = 2L,
    target = largest(),
    method = lanczos(block = 2L, max_subspace = 4L, max_restarts = 20L),
    seed = 7,
    tol = 1e-7
  )

  expect_equal(fit$method, "native block Hermitian Lanczos (thick restart, locking)")
  expect_equal(values(fit), c(50, 40), tolerance = 1e-6)
  expect_true(certificate(fit)$passed)
  expect_gt(fit$restart$restarts_used, 0L)
  expect_gt(fit$restart$locking_events, 0L)
  expect_equal(fit$restart$locked_count, 2L)
  expect_equal(fit$restart$max_subspace, 4L)
  expect_equal(fit$restart$block, 2L)
  expect_equal(fit$operator_allocations, 0)
  expect_equal(fit$restart$operator_allocations, 0)
  expect_s3_class(fit$restart$history, "data.frame")
  expect_equal(fit$convergence_history, fit$restart$history)
  expect_equal(fit$restart$history$restart, seq_len(nrow(fit$restart$history)) - 1L)
  expect_true(any(fit$restart$history$locked_after > fit$restart$history$locked_before))
  expect_lte(utils::tail(fit$restart$history$max_backward_error, 1L), 1e-7)
})

test_that("native block thick-restart production path locks on certificate scale and executes planned controls", {
  A <- diag(c(25, 18, 7, 3, 1, -2, -4, -9))
  method <- lanczos(block = 2L, max_subspace = 4L, max_restarts = 30L)

  plan <- plan_solver(eigen_problem(A, target = largest()), k = 2L, method = method)
  fit <- eig_partial(A, k = 2L, target = largest(), method = method, seed = 17, tol = 1e-8)
  cert <- certificate(fit)

  expect_equal(fit$method, "native block Hermitian Lanczos (thick restart, locking)")
  expect_equal(plan$method, fit$plan$method)
  expect_equal(fit$plan$controls$block, 2L)
  expect_equal(fit$plan$controls$max_subspace, 4L)
  expect_equal(fit$plan$controls$max_restarts, 30L)
  expect_equal(fit$restart$block, fit$plan$controls$block)
  expect_equal(fit$restart$max_subspace, fit$plan$controls$max_subspace)
  expect_equal(fit$restart$max_restarts, fit$plan$controls$max_restarts)
  expect_equal(fit$locked, seq_len(fit$restart$locked_count))
  for (idx in fit$locked) {
    expect_lte(cert$residuals[idx], cert$tolerance * cert$scale[idx])
  }
  expect_true(cert$passed)
})

test_that("default block controls use full native subspace for small dense cases", {
  A <- diag(as.numeric(seq(80, 1)))
  method <- lanczos(block = 2L, max_restarts = 30L)
  plan <- plan_solver(eigen_problem(A, target = largest()), k = 5L, method = method)
  fit <- eig_partial(A, k = 5L, target = largest(), method = method, seed = 18, tol = 1e-8)

  expect_equal(plan$controls$max_subspace, nrow(A))
  expect_equal(fit$restart$kind, "block_full_subspace_dense_lapack")
  expect_equal(fit$restart$max_subspace, nrow(A))
  expect_equal(values(fit), seq(80, 76), tolerance = 1e-12)
  expect_true(certificate(fit)$passed)

  selected_smallest <- eigencore:::native_dense_symmetric_eigen_selected(A, 4L, smallest())
  expect_equal(selected_smallest[[1L]], 1:4, tolerance = 1e-12)
  expect_lt(max(abs(crossprod(selected_smallest[[2L]]) - diag(4L))), 1e-12)
})

test_that("native block thick-restart production path supports sparse CSC without densifying", {
  A <- Matrix::sparseMatrix(
    i = 1:8,
    j = 1:8,
    x = c(9, 7, 5, 3, 2, 1, -1, -2)
  )

  fit <- eig_partial(
    A,
    k = 2L,
    target = largest(),
    method = lanczos(block = 2L, max_subspace = 4L, max_restarts = 20L),
    seed = 8,
    tol = 1e-7
  )

  expect_equal(fit$method, "native block Hermitian Lanczos (thick restart, locking)")
  expect_equal(values(fit), c(9, 7), tolerance = 1e-6)
  expect_true(certificate(fit)$passed)
  expect_equal(fit$restart$locking, "in_native_loop")
  expect_equal(fit$plan$controls$block, 2L)
  expect_equal(fit$operator_allocations, 0)
  expect_equal(fit$restart$operator_bytes_allocated, 0)
})

test_that("native block thick-restart production path residuals match direct operator certificate", {
  n <- 120L
  A <- Matrix::bandSparse(
    n,
    k = c(-1, 0, 1),
    diagonals = list(rep(-1, n - 1L), c(1, rep(2, n - 2L), 1), rep(-1, n - 1L))
  )

  fit <- eig_partial(
    A,
    k = 5L,
    target = smallest(),
    method = lanczos(block = 2L, max_subspace = 36L, max_restarts = 100L),
    seed = 77,
    tol = 1e-8
  )
  direct <- eigencore:::certify_eigen_operator(
    as_operator(A),
    values(fit),
    vectors(fit),
    tol = 1e-8
  )

  expect_true(certificate(fit)$passed)
  expect_true(direct$passed)
  expect_lt(max(abs(certificate(fit)$residuals - direct$residuals)), 1e-12)
  expect_lt(max(abs(certificate(fit)$backward_error - direct$backward_error)), 1e-12)
})

test_that("native block final polish preserves genuine locked prefix on exhaustion", {
  n <- 60L
  A <- Matrix::bandSparse(
    n,
    k = c(-1, 0, 1),
    diagonals = list(rep(-1, n - 1L), c(1, rep(2, n - 2L), 1), rep(-1, n - 1L))
  )

  fit <- eig_partial(
    A,
    k = 4L,
    target = smallest(),
    method = lanczos(block = 2L, max_subspace = 20L, max_restarts = 10L),
    seed = 1,
    tol = 1e-8
  )
  locked_by_restart <- max(fit$restart$history$locked_after)

  expect_equal(locked_by_restart, 1L)
  expect_equal(fit$restart$locked_count, locked_by_restart)
  expect_gte(fit$nconv, locked_by_restart)
  expect_true(certificate(fit)$converged[[1L]])
})

test_that("native block thick-restart production path agrees with R oracle on clustered subspaces", {
  vals <- c(12, 12 - 1e-9, 11.999999, 3, 1, 0, -1, -2)
  A <- symmetric_with_spectrum(vals, seed = 211)

  set.seed(211)
  ref <- eigencore:::reference_block_lanczos_thick_restart_hermitian(
    A,
    k = 3L,
    target = largest(),
    block = 2L,
    max_subspace = 6L,
    max_restarts = 50L,
    tol = 1e-8
  )
  set.seed(211)
  fit <- eig_partial(
    A,
    k = 3L,
    target = largest(),
    method = lanczos(block = 2L, max_subspace = 6L, max_restarts = 50L),
    tol = 1e-8
  )

  expect_true(certificate(ref)$passed)
  expect_true(certificate(fit)$passed)
  expect_equal(values(fit), values(ref), tolerance = 1e-8)
  expect_lt(subspace_distance(vectors(fit), vectors(ref)), 1e-6)
  expect_gt(fit$restart$restarts_used, 0L)
  expect_gt(fit$restart$locking_events, 0L)
})

test_that("native block thick-restart falls back rather than returning duplicate low-residual Ritz directions", {
  A <- random_symmetric_with_spectrum(200L, pattern = "clustered", seed = 1L)
  truth <- sort(spectrum_pattern("clustered", 200L), decreasing = TRUE)[seq_len(5L)]

  fit <- eig_partial(
    A,
    k = 5L,
    target = largest(),
    method = lanczos(block = 2L, max_subspace = 35L, max_restarts = 100L),
    seed = 1L,
    tol = 1e-8
  )

  expect_true(isTRUE(fit$restart$fallback_used))
  expect_match(fit$warnings, "failed certification", fixed = TRUE)
  expect_true(certificate(fit)$passed)
  expect_lt(certificate(fit)$max_orthogonality_loss, 1e-8)
  expect_equal(values(fit), truth, tolerance = 1e-6)
  expect_lt(max(abs(crossprod(vectors(fit)) - diag(5L))), 1e-8)
})

test_that("native block thick-restart production path preserves target taxonomy", {
  A <- diag(c(-9, 7, -2, 0.5, 0.1, 3))
  method <- lanczos(block = 2L, max_subspace = 6L, max_restarts = 20L)

  set.seed(12)
  lm <- eig_partial(A, k = 2L, target = largest_magnitude(), method = method, tol = 1e-8)
  set.seed(12)
  sm <- eig_partial(A, k = 2L, target = smallest_magnitude(), method = method, tol = 1e-8)
  set.seed(12)
  sa <- eig_partial(A, k = 2L, target = smallest(), method = method, tol = 1e-8)
  set.seed(12)
  la <- eig_partial(A, k = 2L, target = largest(), method = method, tol = 1e-8)

  expect_equal(values(lm), c(-9, 7), tolerance = 1e-8)
  expect_equal(values(sm), c(0.1, 0.5), tolerance = 1e-8)
  expect_equal(values(sa), c(-9, -2), tolerance = 1e-8)
  expect_equal(values(la), c(7, 3), tolerance = 1e-8)
  expect_true(certificate(lm)$passed)
  expect_true(certificate(sm)$passed)
  expect_true(certificate(sa)$passed)
  expect_true(certificate(la)$passed)
})

test_that("native block thick-restart production path handles repeated wanted clusters when block covers multiplicity", {
  A <- diag(c(5, 5, 5, 1, 0, -1, -2))

  fit <- eig_partial(
    A,
    k = 3L,
    target = largest(),
    method = lanczos(block = 3L, max_subspace = 6L, max_restarts = 30L),
    seed = 55,
    tol = 1e-8
  )

  expect_equal(values(fit), c(5, 5, 5), tolerance = 1e-10)
  expect_true(certificate(fit)$passed)
  expect_lt(certificate(fit)$max_orthogonality_loss, 1e-12)
  expect_equal(fit$restart$locked_count, 3L)
})

test_that("native block thick-restart production path rejects duplicate locked vectors", {
  vals <- 10^seq(3, -3, length.out = 30)

  fit <- eig_partial(
    diag(vals),
    k = 3L,
    target = largest(),
    method = lanczos(block = 2L, max_subspace = 30L, max_restarts = 100L),
    seed = 1500,
    tol = 1e-8
  )

  expect_lt(certificate(fit)$max_orthogonality_loss, 1e-12)
  expect_equal(diag(crossprod(vectors(fit))), rep(1, 3), tolerance = 1e-10)
  expect_lt(max(abs(crossprod(vectors(fit)) - diag(3))), 1e-10)
  expect_equal(length(unique(signif(values(fit), 12))), 3L)
})

test_that("native block production path uses dense full-subspace path when max_subspace spans n", {
  A <- symmetric_with_spectrum(c(9, 6, 4, 2, 1, 0.5, 0.25, 0.1), seed = 1208)

  fit <- eig_partial(
    A,
    k = 3L,
    target = largest(),
    method = lanczos(block = 2L, max_subspace = nrow(A), max_restarts = 100L),
    seed = 1208,
    tol = 1e-8
  )

  expect_true(certificate(fit)$passed)
  expect_equal(values(fit), c(9, 6, 4), tolerance = 1e-10)
  expect_equal(fit$restart$kind, "block_full_subspace_dense_lapack")
  expect_equal(fit$restart$locking, "not_required_full_subspace")
  expect_equal(fit$matvecs, 0L)
  expect_equal(fit$restarts, 0L)
})
