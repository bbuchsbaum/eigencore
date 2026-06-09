## Unit tests for R/reference_arnoldi.R primitives.
## These focus on internals that are not reachable through the native Arnoldi
## path exercised in test-adversarial.R and test-solvers.R.

# ---------------------------------------------------------------------------
# reference_arnoldi_target_supported
# ---------------------------------------------------------------------------

test_that("reference_arnoldi_target_supported accepts all documented kinds", {
  supported <- list(
    largest(), smallest(), largest_magnitude(),
    largest_real(), smallest_real(),
    largest_imaginary(), smallest_imaginary()
  )
  for (t in supported) {
    expect_true(
      eigencore:::reference_arnoldi_target_supported(t),
      label = paste("should be supported:", t$kind)
    )
  }
})

test_that("reference_arnoldi_target_supported rejects unsupported targets", {
  expect_false(eigencore:::reference_arnoldi_target_supported(nearest(1)))
  expect_false(eigencore:::reference_arnoldi_target_supported(smallest_magnitude()))
})

# ---------------------------------------------------------------------------
# native_arnoldi_default_max_subspace
# ---------------------------------------------------------------------------

test_that("native_arnoldi_default_max_subspace respects n ceiling and k floor", {
  # m >= k + 1; m <= n
  expect_equal(eigencore:::native_arnoldi_default_max_subspace(100L, 1L), 9L)
  expect_equal(eigencore:::native_arnoldi_default_max_subspace(100L, 8L), 72L)
  expect_equal(eigencore:::native_arnoldi_default_max_subspace(100L, 16L), 100L)
  expect_equal(eigencore:::native_arnoldi_default_max_subspace(5L, 3L), 5L)
  expect_equal(eigencore:::native_arnoldi_default_max_subspace(10L, 200L), 10L)
})

# ---------------------------------------------------------------------------
# native_matrix_free_arnoldi_available
# ---------------------------------------------------------------------------

test_that("native_matrix_free_arnoldi_available detects eligible operators", {
  A <- diag(c(3, 2, 1))
  dense_op <- eigencore:::as_operator(A)
  expect_false(eigencore:::native_matrix_free_arnoldi_available(dense_op))

  mf_op <- linear_operator(
    dim = dim(A),
    apply = function(X) A %*% X,
    dtype = "double"
  )
  expect_true(eigencore:::native_matrix_free_arnoldi_available(mf_op))
})

test_that("native_matrix_free_arnoldi_available rejects non-double operators", {
  A <- matrix(as.integer(diag(3L)), 3, 3)
  int_op <- linear_operator(
    dim = dim(A),
    apply = function(X) A %*% X,
    dtype = "integer"
  )
  expect_false(eigencore:::native_matrix_free_arnoldi_available(int_op))
})

test_that("native_matrix_free_arnoldi_available rejects non-square operators", {
  rect_op <- linear_operator(
    dim = c(4L, 3L),
    apply = function(X) matrix(rnorm(4 * ncol(X)), 4),
    dtype = "double"
  )
  expect_false(eigencore:::native_matrix_free_arnoldi_available(rect_op))
})

# ---------------------------------------------------------------------------
# reference_arnoldi_cycle: structural correctness
# ---------------------------------------------------------------------------

test_that("reference_arnoldi_cycle produces orthonormal basis and upper Hessenberg H", {
  set.seed(31)
  n <- 6L
  A <- matrix(rnorm(n * n), n, n)
  op <- eigencore:::as_operator(A)
  start <- rnorm(n); start <- start / sqrt(sum(start^2))
  m_req <- 5L
  cyc <- eigencore:::reference_arnoldi_cycle(op, start, m_req)
  m <- cyc$iterations

  Vm <- cyc$V[, seq_len(m), drop = FALSE]
  Hm <- cyc$H[seq_len(m), seq_len(m), drop = FALSE]

  # V columns are orthonormal
  expect_equal(crossprod(Vm), diag(m), tolerance = 1e-12)

  # H_m is upper Hessenberg: H[i,j] == 0 for i > j+1
  below_subdiag <- 0
  for (j in seq_len(m)) {
    n_below <- max(0L, m - j - 1L)
    if (n_below > 0L) {
      rows_below <- seq_len(n_below) + j + 1L
      below_subdiag <- max(below_subdiag, max(abs(Hm[rows_below, j])))
    }
  }
  expect_equal(below_subdiag, 0, tolerance = 1e-14)

  # Arnoldi factorization: A * V[:,1:m-1] = V[:,1:m] * H_extended columns 1:m-1
  # Columns 1..m-1: A*v_j = V*H_col_j exactly (residual only in last column)
  Hext <- cyc$H[seq_len(m + 1L), seq_len(m), drop = FALSE]
  Vext <- cyc$V[, seq_len(m + 1L), drop = FALSE]
  AV <- A %*% Vm
  VH <- Vext %*% Hext
  col_resids <- sqrt(colSums((AV - VH)^2))
  expect_equal(col_resids[seq_len(m - 1L)], rep(0, m - 1L), tolerance = 1e-12)
  expect_equal(col_resids[m], abs(cyc$H[m + 1L, m]), tolerance = 1e-12)
})

test_that("reference_arnoldi_cycle matvecs counter equals m", {
  set.seed(32)
  A <- matrix(rnorm(16), 4, 4)
  op <- eigencore:::as_operator(A)
  start <- c(1, 0, 0, 0)
  cyc <- eigencore:::reference_arnoldi_cycle(op, start, 3L)
  expect_equal(cyc$matvecs, cyc$iterations)
})

# ---------------------------------------------------------------------------
# reference_arnoldi_score / reference_arnoldi_score_better
# ---------------------------------------------------------------------------

test_that("reference_arnoldi_score maps NA and non-finite errors to Inf", {
  cert_na  <- list(passed = FALSE, converged = FALSE, max_backward_error = NA_real_)
  cert_inf <- list(passed = FALSE, converged = FALSE, max_backward_error = Inf)
  cert_nan <- list(passed = FALSE, converged = FALSE, max_backward_error = NaN)
  expect_equal(eigencore:::reference_arnoldi_score(cert_na)$max_backward_error,  Inf)
  expect_equal(eigencore:::reference_arnoldi_score(cert_inf)$max_backward_error, Inf)
  expect_equal(eigencore:::reference_arnoldi_score(cert_nan)$max_backward_error, Inf)
})

test_that("reference_arnoldi_score_better respects passed > nconv > error priority", {
  s_pass  <- eigencore:::reference_arnoldi_score(list(passed = TRUE,  converged = c(TRUE, TRUE),  max_backward_error = 1e-4))
  s_fail2 <- eigencore:::reference_arnoldi_score(list(passed = FALSE, converged = c(TRUE, TRUE),  max_backward_error = 1e-3))
  s_fail1 <- eigencore:::reference_arnoldi_score(list(passed = FALSE, converged = c(TRUE, FALSE), max_backward_error = 1e-6))
  s_fail0 <- eigencore:::reference_arnoldi_score(list(passed = FALSE, converged = c(FALSE, FALSE), max_backward_error = 1e-6))

  # passed beats non-passed regardless of other fields
  expect_true(eigencore:::reference_arnoldi_score_better(s_pass, s_fail0))
  expect_false(eigencore:::reference_arnoldi_score_better(s_fail0, s_pass))

  # more converged beats fewer converged (both not passed)
  expect_true(eigencore:::reference_arnoldi_score_better(s_fail2, s_fail1))
  expect_false(eigencore:::reference_arnoldi_score_better(s_fail1, s_fail2))

  # same nconv: smaller error wins
  expect_true(eigencore:::reference_arnoldi_score_better(s_fail1, s_fail0))
})

# ---------------------------------------------------------------------------
# match_left_eigenvectors
# ---------------------------------------------------------------------------

test_that("match_left_eigenvectors matches by nearest-eigenvalue greedily", {
  lvals <- c(5.05, 3.02, 1.0)
  Lmat  <- diag(3)
  colnames(Lmat) <- NULL
  wanted <- c(5, 3)
  res <- eigencore:::match_left_eigenvectors(lvals, Lmat, wanted)
  expect_equal(res$values, c(5.05, 3.02))
  expect_equal(res$match_distance, c(0.05, 0.02), tolerance = 1e-12)
  expect_equal(res$vectors, Lmat[, c(1L, 2L), drop = FALSE])
})

test_that("match_left_eigenvectors returns NULL when wanted exceeds available", {
  lvals <- c(5, 3)
  Lmat  <- diag(2)
  expect_null(eigencore:::match_left_eigenvectors(lvals, Lmat, c(5, 3, 1)))
})

test_that("match_left_eigenvectors returns empty result for zero-length wanted", {
  lvals <- c(5, 3, 1)
  Lmat  <- diag(3)
  res <- eigencore:::match_left_eigenvectors(lvals, Lmat, numeric(0L))
  expect_equal(length(res$values), 0L)
  expect_equal(ncol(res$vectors), 0L)
})

test_that("match_left_eigenvectors does not reuse the same left eigenvalue", {
  # Both wanted values are closest to left[1]; greedy should pick [1] for
  # first wanted and [2] for second.
  lvals <- c(5.0, 4.8)
  Lmat  <- cbind(c(1, 0), c(0, 1))
  res <- eigencore:::match_left_eigenvectors(lvals, Lmat, c(5.0, 4.9))
  expect_equal(res$values, c(5.0, 4.8))
})

# ---------------------------------------------------------------------------
# normalize_left_eigenvectors
# ---------------------------------------------------------------------------

test_that("normalize_left_eigenvectors achieves unit diagonal biorthogonality", {
  set.seed(41)
  R <- matrix(rnorm(9), 3, 3)
  L <- matrix(rnorm(9), 3, 3)
  Ln <- eigencore:::normalize_left_eigenvectors(L, R)
  gram_diag <- diag(crossprod(Ln, R))
  expect_equal(Re(gram_diag), rep(1, 3), tolerance = 1e-12)
})

test_that("normalize_left_eigenvectors leaves near-zero-gram columns unchanged", {
  R <- diag(c(0, 1, 1))
  L <- diag(3)
  Ln <- eigencore:::normalize_left_eigenvectors(L, R)
  # Column 1: gram = <L[,1], R[,1]> = 0 => near-zero, unchanged
  expect_equal(Ln[, 1L], L[, 1L])
  # Columns 2 and 3: normalized
  expect_equal(diag(crossprod(Ln[, 2:3, drop=FALSE], R[, 2:3, drop=FALSE])),
               c(1, 1), tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# arnoldi_left_eigen_contract: failure modes
# ---------------------------------------------------------------------------

test_that("arnoldi_left_eigen_contract returns supported=FALSE when adjoint unavailable", {
  op_no_adj <- linear_operator(
    dim = c(3L, 3L),
    apply = function(X) matrix(rnorm(3L * ncol(X)), 3L),
    structure = general(),
    dtype = "double"
  )
  set.seed(7)
  right_vecs <- matrix(rnorm(6L), 3L, 2L)
  res <- eigencore:::arnoldi_left_eigen_contract(op_no_adj, c(1, 2), right_vecs, largest())
  expect_false(res$supported)
  expect_false(is.null(res$reason))
  expect_null(res$vectors)
  expect_null(res$certificate)
})

test_that("arnoldi_left_eigen_contract returns supported=FALSE when right_vectors is NULL", {
  A   <- diag(c(3, 2, 1))
  op  <- eigencore:::as_operator(A)
  res <- eigencore:::arnoldi_left_eigen_contract(op, c(3, 2), NULL, largest())
  expect_false(res$supported)
  expect_match(res$reason, "right eigenvectors", fixed = TRUE)
})

# ---------------------------------------------------------------------------
# reference_arnoldi_general: correctness and restart tracking
# ---------------------------------------------------------------------------

test_that("reference_arnoldi_general returns correct values with native=FALSE label", {
  A <- diag(c(7, 5, 3, 1))
  A[1, 2] <- 2; A[2, 3] <- 1
  op <- eigencore:::as_operator(A)

  res <- eigencore:::reference_arnoldi_general(
    op, k = 2L, target = largest_real(), tol = 1e-10
  )

  expect_equal(sort(Re(res$values), decreasing = TRUE), c(7, 5), tolerance = 1e-8)
  expect_true(res$certificate$passed)
  expect_equal(res$restart$kind, "reference_arnoldi")
  expect_false(res$restart$native)
  expect_equal(res$restart$restart_count, 0L)
})

test_that("reference_arnoldi_general works for smallest_real target", {
  A <- diag(c(7, 5, 3, 1))
  A[3, 4] <- 0.5
  op <- eigencore:::as_operator(A)

  res <- eigencore:::reference_arnoldi_general(
    op, k = 2L, target = smallest_real(), tol = 1e-10
  )

  expect_equal(sort(Re(res$values)), c(1, 3), tolerance = 1e-8)
  expect_true(res$certificate$passed)
})

test_that("reference_arnoldi_general works for smallest_imaginary target", {
  # [[0,1],[-1,0]] has eigenvalues +i and -i
  A  <- matrix(c(0, 1, -1, 0), 2, 2)
  op <- eigencore:::as_operator(A)

  res <- eigencore:::reference_arnoldi_general(
    op, k = 1L, target = smallest_imaginary(), tol = 1e-10
  )

  expect_true(is.complex(res$values))
  expect_equal(Im(res$values[[1L]]), -1, tolerance = 1e-8)
})

test_that("reference_arnoldi_general works for largest_imaginary target", {
  A  <- matrix(c(0, 1, -1, 0), 2, 2)
  op <- eigencore:::as_operator(A)

  res <- eigencore:::reference_arnoldi_general(
    op, k = 1L, target = largest_imaginary(), tol = 1e-10
  )

  expect_true(is.complex(res$values))
  expect_equal(Im(res$values[[1L]]), 1, tolerance = 1e-8)
})

test_that("reference_arnoldi_general rejects unsupported targets", {
  A  <- diag(3)
  op <- eigencore:::as_operator(A)

  expect_error(
    eigencore:::reference_arnoldi_general(op, k = 1L, target = nearest(1)),
    "largest/smallest real-part and largest-magnitude"
  )
  expect_error(
    eigencore:::reference_arnoldi_general(op, k = 1L, target = smallest_magnitude()),
    "largest/smallest real-part and largest-magnitude"
  )
})

test_that("reference_arnoldi_general rejects non-square operators", {
  rect_op <- linear_operator(
    dim = c(4L, 3L),
    apply = function(X) matrix(0, 4, ncol(X)),
    dtype = "double"
  )
  expect_error(
    eigencore:::reference_arnoldi_general(rect_op, k = 1L, target = largest()),
    "square"
  )
})

test_that("reference_arnoldi_general restart history has correct shape and selects best", {
  set.seed(42)
  n  <- 12L
  A  <- matrix(0, n, n)
  diag(A) <- seq(n, 1)
  A[cbind(seq_len(n - 1L), seq_len(n - 1L) + 1L)] <- 0.5
  op <- eigencore:::as_operator(A)

  res <- eigencore:::reference_arnoldi_general(
    op, k = 3L, target = largest_real(), tol = 1e-14,
    maxit = 4L, max_restarts = 2L
  )

  hist <- res$restart$attempt_history
  expect_equal(nrow(hist), 3L)
  expect_equal(res$restart$restart_count, 2L)
  expect_true(res$restart$selected_attempt %in% seq_len(3L))
  # selected attempt matches the row with minimum backward error
  finite_errors <- hist$max_backward_error[is.finite(hist$max_backward_error)]
  if (length(finite_errors)) {
    expect_equal(
      res$restart$selected_attempt,
      hist$attempt[which.min(hist$max_backward_error)]
    )
  }
  # stage_seconds are named and non-negative
  expect_setequal(names(res$restart$stage_seconds), c("cycle", "ritz_extraction"))
  expect_true(all(res$restart$stage_seconds >= 0))
})

test_that("reference_arnoldi_general vectors=FALSE omits vector matrix", {
  A  <- diag(c(5, 3, 1))
  op <- eigencore:::as_operator(A)

  res <- eigencore:::reference_arnoldi_general(
    op, k = 2L, target = largest_real(), tol = 1e-10, vectors = FALSE
  )
  expect_null(res$vectors)
})
