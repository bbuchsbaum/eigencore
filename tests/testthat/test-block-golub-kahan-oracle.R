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
    start <- matrix(stats::rnorm(ncol(A) * 2L), nrow = ncol(A), ncol = 2L)
    basis <- eigencore:::native_block_golub_kahan_basis(
      A_in,
      max_subspace = ncol(A),
      block = 2L,
      start = start
    )
    compact_fit <- eigencore:::native_block_golub_kahan_fit(
      A_in,
      max_subspace = ncol(A),
      rank = 4L,
      target = largest(),
      block = 2L,
      start = start
    )
    compact_ref <- eigencore:::native_block_golub_kahan_ritz(
      basis$V,
      basis$AV,
      rank = 4L,
      target = largest(),
      active_cols = basis$active_cols
    )

    expect_false("U" %in% names(basis))
    expect_true(all(c("V", "AV", "active_cols", "active_left_cols") %in% names(basis)))
    expect_equal(dim(basis$V), c(ncol(A), ncol(A)))
    expect_equal(dim(basis$AV), c(nrow(A), ncol(A)))
    expect_equal(compact_fit$d, compact_ref$d, tolerance = 1e-10)
    expect_equal(abs(crossprod(compact_fit$v, compact_ref$v)),
                 diag(4L), tolerance = 1e-8)
    expect_equal(compact_fit$active_cols, basis$active_cols)
    expect_true(all(c("native_iteration", "ritz") %in% names(compact_fit$stage_seconds)))
    expect_gt(compact_fit$stage_seconds[["native_iteration"]], 0)
    expect_gt(compact_fit$stage_seconds[["ritz"]], 0)
    cached_basis <- eigencore:::native_block_golub_kahan_basis(
      A_in,
      max_subspace = ncol(A),
      block = 2L,
      start = basis$V[, 1:2, drop = FALSE],
      start_av = basis$AV[, 1:2, drop = FALSE]
    )
    expect_true(cached_basis$cached_start_used)
    expect_lt(cached_basis$matvecs, basis$matvecs)
    expect_equal(cached_basis$AV[, 1:2, drop = FALSE],
                 basis$AV[, 1:2, drop = FALSE],
                 tolerance = 1e-10)
    cached_prefix_basis <- eigencore:::native_block_golub_kahan_basis(
      A_in,
      max_subspace = ncol(A),
      block = 3L,
      start = cbind(basis$V[, 1:2, drop = FALSE], stats::rnorm(ncol(A))),
      start_av = basis$AV[, 1:2, drop = FALSE]
    )
    expect_true(cached_prefix_basis$cached_start_used)
    expect_lt(cached_prefix_basis$matvecs, basis$matvecs)
    expect_equal(cached_prefix_basis$AV[, 1:2, drop = FALSE],
                 basis$AV[, 1:2, drop = FALSE],
                 tolerance = 1e-10)

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
    expect_false(fit$restart$basis_returned)
    expect_true(all(c("native_iteration", "ritz") %in% names(fit$stage_seconds)))
    expect_gte(fit$restart$active_cols, 4L)
    expect_gt(fit$matvecs, 0L)
    expect_gt(fit$restart$ortho_passes, 0L)
    expect_equal(fit$d, s[1:4], tolerance = 1e-8)
    expect_true(fit$certificate$passed)
  }
})

test_that("retained block Golub-Kahan restart ABI fixes the native implementation contract", {
  set.seed(701)
  A <- Matrix::t(Matrix::rsparsematrix(600, 90, density = 0.03))
  abi <- eigencore:::native_block_golub_kahan_retained_restart_abi(
    A,
    rank = 5L,
    block = 2L,
    max_attempts = 4L,
    target = largest()
  )

  expect_s3_class(abi, "eigencore_block_golub_kahan_retained_restart_abi")
  expect_equal(abi$version, 1L)
  expect_true(abi$implemented)
  expect_equal(abi$native_storage, "dgCMatrix")
  expect_equal(abi$input_schema$initial_start, c(ncol(A), 2L))
  expect_equal(abi$input_schema$tail_layout, "column-major n x (block * restart_count)")
  expect_true(all(diff(abi$max_subspace_sequence) > 0L))
  expect_lte(max(abi$max_subspace_sequence), ncol(A))
  expect_true("Avectors" %in% abi$output_schema)
  expect_true(any(grepl("no R-side restart block construction", abi$invariants, fixed = TRUE)))
  expect_true(any(grepl("reorthogonalized by the native basis runner", abi$invariants, fixed = TRUE)))
  expect_true(any(grepl("transformed together by native QR normalization", abi$invariants, fixed = TRUE)))
  expect_true(any(grepl("certify or fall back", abi$invariants, fixed = TRUE)))
  expect_equal(unname(abi$entry_points[["csc"]]),
               "eigencore_block_golub_kahan_csc_retained_cycle")

  dense_abi <- eigencore:::native_block_golub_kahan_retained_restart_abi(
    as.matrix(A),
    rank = 5L,
    block = 2L,
    max_attempts = 2L
  )
  expect_equal(dense_abi$native_storage, "double_matrix")
})

test_that("native retained block Golub-Kahan cycle builds restart state inside native code", {
  set.seed(703)
  A <- Matrix::t(Matrix::rsparsematrix(600, 90, density = 0.03))

  set.seed(701)
  retained <- eigencore:::native_block_golub_kahan_retained_cycle_svd(
    A,
    rank = 5L,
    target = largest(),
    block = 2L,
    tol = 1e-8,
    retained_av_cache = TRUE
  )

  expect_identical(retained$restart$kind, "block_golub_kahan_native_retained_cycle")
  expect_true(retained$restart$thick_restart)
  expect_true(retained$restart$retained_restart)
  expect_true(retained$restart$retained_restart_native)
  expect_true(retained$restart$retained_av_cache)
  expect_true(retained$restart$retained_av_cache_attempted)
  expect_false(retained$restart$retained_av_cache_fallback)
  expect_true(retained$restart$native_attempt_certification)
  expect_equal(retained$restart$retained_restart_abi_version, 1L)
  expect_true(retained$certificate$passed)
  expect_gte(sum(retained$certificate$converged), 5L)
  expect_gt(retained$matvecs, 0L)
  expect_gt(retained$restart$native_workspace_bytes, 0)
  expect_true(is.data.frame(retained$restart$attempt_history))
  expect_gt(nrow(retained$restart$attempt_history), 1L)
  expect_true(any(retained$restart$attempt_history$cached_start_used))
  expect_true(all(retained$restart$attempt_history$warm_started[-1L]))
  expect_true(any(retained$restart$attempt_history$certificate_passed))
  expect_true(all(c("converged_count", "leading_converged_count") %in%
                    names(retained$restart$attempt_history)))
  expect_gte(max(retained$restart$attempt_history$converged_count), 5L)
  expect_gte(max(retained$restart$attempt_history$leading_converged_count), 5L)
  expect_true(all(is.finite(retained$restart$attempt_history$max_backward_error)))
  expect_true(all(is.finite(retained$restart$attempt_history$max_residual)))
  expect_equal(retained$restart$attempted_subspaces,
               retained$restart$attempt_history$max_subspace)
  expect_false(retained$restart$basis_returned)
  expect_true(all(c("native_iteration", "ritz", "restart") %in% names(retained$stage_seconds)))
  expect_gt(retained$stage_seconds[["native_iteration"]], 0)
  expect_gt(retained$stage_seconds[["ritz"]], 0)
})

test_that("native retained block Golub-Kahan cached AV path certifies after MGS2 restart normalization", {
  set.seed(702)
  A <- Matrix::t(Matrix::rsparsematrix(600, 90, density = 0.03))

  set.seed(701)
  retained <- eigencore:::native_block_golub_kahan_retained_cycle_svd(
    A,
    rank = 5L,
    target = largest(),
    block = 2L,
    tol = 1e-8,
    retained_av_cache = TRUE
  )

  expect_true(retained$certificate$passed)
  expect_gte(sum(retained$certificate$converged), 5L)
  expect_true(retained$restart$retained_av_cache_attempted)
  expect_false(retained$restart$retained_av_cache_fallback)
  expect_true(retained$restart$retained_av_cache)
  expect_true(any(retained$restart$attempt_history$cached_start_used))
  expect_gte(max(retained$restart$attempt_history$converged_count), 5L)
  expect_gte(max(retained$restart$attempt_history$leading_converged_count), 5L)
  expect_false(retained$restart$fallback_attempted)
  expect_false(retained$restart$fallback_used)
  expect_true(is.na(retained$restart$fallback_method))
  expect_true(is.na(retained$restart$fallback_max_backward_error))
  expect_true(is.na(retained$restart$retained_av_cache_failed_backward_error))
})

test_that("native retained block Golub-Kahan deflation is opt-in and certificate guarded", {
  set.seed(702)
  A <- Matrix::t(Matrix::rsparsematrix(600, 90, density = 0.03))

  set.seed(701)
  retained <- eigencore:::native_block_golub_kahan_retained_cycle_svd(
    A,
    rank = 5L,
    target = largest(),
    block = 2L,
    tol = 1e-8,
    retained_av_cache = TRUE,
    retained_deflation = TRUE
  )

  expect_true(retained$certificate$passed)
  expect_gte(sum(retained$certificate$converged), 5L)
  expect_true(retained$restart$retained_deflation)
  expect_true(retained$restart$retained_av_cache_attempted)
  expect_false(retained$restart$retained_av_cache_fallback)
  expect_type(retained$restart$retained_locked_count, "integer")
  expect_gte(retained$restart$retained_locked_count, 0L)
  expect_true("eigencore_block_golub_kahan_retained_deflated" %in%
                eigencore:::available_svd_methods())
})

test_that("native block Golub-Kahan basis cycle records adaptive subspace attempts", {
  set.seed(702)
  A <- Matrix::t(Matrix::rsparsematrix(600, 90, density = 0.03))

  set.seed(701)
  fixed <- eigencore:::native_block_golub_kahan_cycle_svd(
    A,
    rank = 5L,
    target = largest(),
    block = 2L,
    tol = 1e-8,
    adaptive = FALSE
  )
  set.seed(701)
  adaptive <- eigencore:::native_block_golub_kahan_cycle_svd(
    A,
    rank = 5L,
    target = largest(),
    block = 2L,
    tol = 1e-8
  )
  set.seed(701)
  cold_adaptive <- eigencore:::native_block_golub_kahan_cycle_svd(
    A,
    rank = 5L,
    target = largest(),
    block = 2L,
    tol = 1e-8,
    adaptive_start = "initial"
  )
  set.seed(701)
  lean_adaptive <- eigencore:::native_block_golub_kahan_cycle_svd(
    A,
    rank = 5L,
    target = largest(),
    block = 2L,
    tol = 1e-8,
    adaptive_start = "ritz_lean"
  )
  set.seed(701)
  cached_adaptive <- eigencore:::native_block_golub_kahan_cycle_svd(
    A,
    rank = 5L,
    target = largest(),
    block = 2L,
    tol = 1e-8,
    adaptive_start = "ritz_cached"
  )
  set.seed(701)
  cached_random_adaptive <- eigencore:::native_block_golub_kahan_cycle_svd(
    A,
    rank = 5L,
    target = largest(),
    block = 2L,
    tol = 1e-8,
    adaptive_start = "ritz_cached_random"
  )
  set.seed(701)
  residual_adaptive <- eigencore:::native_block_golub_kahan_cycle_svd(
    A,
    rank = 5L,
    target = largest(),
    block = 2L,
    tol = 1e-8,
    adaptive_start = "ritz_residual"
  )

  expect_false(fixed$certificate$passed)
  expect_true(adaptive$certificate$passed)
  expect_true(lean_adaptive$certificate$passed)
  expect_true(cached_adaptive$certificate$passed)
  expect_true(cached_random_adaptive$certificate$passed)
  expect_true(residual_adaptive$certificate$passed)
  expect_true(adaptive$restart$adaptive)
  expect_equal(adaptive$restart$adaptive_start, "ritz")
  expect_equal(cached_adaptive$restart$adaptive_start, "ritz_cached")
  expect_equal(cached_random_adaptive$restart$adaptive_start, "ritz_cached_random")
  expect_equal(residual_adaptive$restart$adaptive_start, "ritz_residual")
  expect_equal(lean_adaptive$restart$adaptive_start, "ritz_lean")
  expect_gt(adaptive$restart$attempts, 1L)
  expect_equal(adaptive$restart$attempted_subspaces[[1L]],
               adaptive$restart$initial_max_subspace)
  expect_gt(max(adaptive$restart$attempt_history$start_cols), 2L)
  expect_lt(max(lean_adaptive$restart$attempt_history$start_cols),
            max(adaptive$restart$attempt_history$start_cols))
  expect_true(any(adaptive$restart$attempt_history$warm_started))
  expect_true(any(cached_adaptive$restart$attempt_history$cached_start_used))
  expect_true(any(cached_random_adaptive$restart$attempt_history$cached_start_used))
  expect_true(any(residual_adaptive$restart$attempt_history$warm_started))
  expect_true(any(lean_adaptive$restart$attempt_history$warm_started))
  expect_true(any(adaptive$restart$attempt_history$certificate_passed))
  expect_equal(adaptive$matvecs, sum(adaptive$restart$attempt_history$matvecs))
  expect_lt(adaptive$matvecs, cold_adaptive$matvecs)
  expect_lte(residual_adaptive$matvecs, cold_adaptive$matvecs)
  expect_lte(cached_random_adaptive$matvecs, adaptive$matvecs)
  expect_lt(cached_adaptive$matvecs, lean_adaptive$matvecs)
  expect_lt(lean_adaptive$matvecs, cold_adaptive$matvecs)
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
    max_subspace = 8L,
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
