# Epic D kernel upgrades: mid-sweep convergence checks (check_stride) for the
# native block thick-restart Hermitian Lanczos, and the matrix-free (R-callback)
# entry to that same kernel. Correctness is pinned against a dense eigen()
# oracle; the mid-sweep speedup is measured as a strict matvec reduction on a
# problem whose leading pairs converge well inside a generous subspace budget.

# A well-separated leading spectrum over a fast-decaying bulk. The top pairs
# converge in a few blocks, so a legacy single-sweep solve fills the whole m_max
# budget before its one convergence check while mid-sweep checks stop early.
clustered_decay_sym <- function(n, seed = 1L) {
  set.seed(seed)
  lam <- c(100, 60, 35, 20, 12, rep(1, n - 5L) * exp(-seq_len(n - 5L) / 5))
  Q <- qr.Q(qr(matrix(stats::rnorm(n * n), n)))
  A <- Q %*% diag(lam) %*% t(Q)
  (A + t(A)) / 2
}

matrix_free_hermitian_op <- function(A) {
  n <- nrow(A)
  linear_operator(
    dim = c(n, n),
    apply = function(X, alpha = 1, beta = 0, Y = NULL) {
      Z <- alpha * (A %*% X)
      if (is.null(Y) || beta == 0) Z else Z + beta * Y
    },
    structure = hermitian(),
    name = "matrix-free hermitian test operator"
  )
}

# `mult` near-degenerate leading eigenvalues within 1e-9 of 100, then a
# well-separated tail over a decaying bulk. A small mid-sweep subspace resolves
# only part of the cluster, so a naive early lock would return a set that misses
# an unresolved cluster member while keeping clean per-pair residuals.
degenerate_cluster_sym <- function(n, mult, seed = 1L) {
  set.seed(seed)
  top <- 100 + (seq_len(mult) - 1L) * 1e-9
  tail <- c(60, 35, 20)
  bulk_len <- n - mult - length(tail)
  lam <- c(top, tail, rep(1, bulk_len) * exp(-seq_len(bulk_len) / 5))
  Q <- qr.Q(qr(matrix(stats::rnorm(n * n), n)))
  A <- Q %*% diag(lam) %*% t(Q)
  (A + t(A)) / 2
}

# Build an n x n symmetric operator with a prescribed leading spectrum over a
# decaying positive bulk, in a seed-dependent random basis.
spectrum_sym <- function(lam_top, n, seed = 1L) {
  set.seed(seed)
  bulk_len <- n - length(lam_top)
  lam <- c(lam_top, rep(1, bulk_len) * exp(-seq_len(bulk_len) / 5))
  Q <- qr.Q(qr(matrix(stats::rnorm(n * n), n)))
  A <- Q %*% diag(lam) %*% t(Q)
  (A + t(A)) / 2
}

test_that("mid-sweep checks match the oracle and cut matvecs on dense and CSC", {
  n <- 200L
  k <- 4L
  A <- clustered_decay_sym(n, seed = 11L)
  truth <- sort(eigen(A, symmetric = TRUE)$values, decreasing = TRUE)[seq_len(k)]
  m_max <- 80L

  for (storage in c("dense", "csc")) {
    op <- if (identical(storage, "dense")) A else methods::as(A, "dgCMatrix")
    base <- eig_partial(
      op, k = k, target = largest(),
      method = lanczos(max_subspace = m_max, block = 2L), seed = 5
    )
    mid <- eig_partial(
      op, k = k, target = largest(),
      method = lanczos(max_subspace = m_max, block = 2L, check_stride = 1L),
      seed = 5
    )
    expect_equal(sort(values(mid), decreasing = TRUE), truth, tolerance = 1e-7)
    expect_certificate_clean(mid, tol = 1e-7)
    # Converges mid-sweep: substantially fewer operator applications than the
    # legacy single-check-per-sweep solve. The strong bound guards against a
    # regression where the cluster-clearance guard over-defers a well-separated
    # spectrum back to a full sweep (base 40 -> mid < 20, observed ~7-8).
    expect_lt(mid$matvecs, base$matvecs %/% 2L)
    expect_identical(mid$restart$check_stride, 1L)
  }
})

test_that("check_stride = 0 is regression-identical to the unqualified solve", {
  n <- 160L
  k <- 3L
  A <- clustered_decay_sym(n, seed = 4L)
  m_max <- 70L
  ref <- eig_partial(
    A, k = k, target = largest(),
    method = lanczos(max_subspace = m_max, block = 2L), seed = 8
  )
  zero <- eig_partial(
    A, k = k, target = largest(),
    method = lanczos(max_subspace = m_max, block = 2L, check_stride = 0L),
    seed = 8
  )
  expect_identical(values(zero), values(ref))
  expect_identical(zero$matvecs, ref$matvecs)
  expect_identical(vectors(zero), vectors(ref))
  expect_identical(zero$restart$restarts_used, ref$restart$restarts_used)
})

test_that("matrix-free block thick-restart matches the oracle", {
  n <- 120L
  k <- 4L
  A <- clustered_decay_sym(n, seed = 6L)
  truth <- sort(eigen(A, symmetric = TRUE)$values, decreasing = TRUE)[seq_len(k)]
  op <- matrix_free_hermitian_op(A)

  fit <- eig_partial(
    op, k = k, target = largest(),
    method = lanczos(max_subspace = 60L, block = 2L, check_stride = 1L), seed = 5
  )
  expect_identical(
    fit$method,
    "native block Hermitian Lanczos (matrix-free callback, thick restart, locking)"
  )
  expect_equal(sort(values(fit), decreasing = TRUE), truth, tolerance = 1e-6)
  # Matrix-free certificates withhold `passed` (stochastic norm estimate);
  # assert the underlying evidence directly.
  cert <- certificate(fit)
  expect_true(all(cert$converged))
  expect_lte(cert$max_backward_error, 1e-6)
  expect_true(isTRUE(cert$orthogonality_passed))
})

test_that("matrix-free scalar (block == 1) stays on the reference path", {
  n <- 80L
  k <- 3L
  A <- clustered_decay_sym(n, seed = 9L)
  op <- matrix_free_hermitian_op(A)
  fit <- eig_partial(
    op, k = k, target = largest(), method = lanczos(), maxit = n, seed = 2
  )
  expect_identical(
    fit$method, "reference Hermitian Lanczos (prototype/oracle fallback)"
  )
})

test_that("matrix-free block warm start enters at full width and beats cold", {
  n <- 120L
  k <- 4L
  A <- clustered_decay_sym(n, seed = 12L)
  truth <- sort(eigen(A, symmetric = TRUE)$values, decreasing = TRUE)[seq_len(k)]
  op <- matrix_free_hermitian_op(A)
  meth <- lanczos(max_subspace = 60L, block = k, check_stride = 1L)

  cold <- eig_partial(op, k = k, target = largest(), method = meth, seed = 3)

  set.seed(21)
  warm_basis <- eigen(A, symmetric = TRUE)$vectors[, seq_len(k), drop = FALSE] +
    1e-3 * matrix(stats::rnorm(n * k), n, k)
  warm <- eig_partial(
    op, k = k, target = largest(), method = meth, seed = 3,
    initial_subspace = warm_basis
  )

  expect_equal(sort(values(warm), decreasing = TRUE), truth, tolerance = 1e-7)
  # An initial_subspace of rank k with block >= k enters at full width: no
  # rotation-compression to a single column.
  prov <- warm$initial_subspace
  expect_identical(prov$start_source, "user_supplied")
  expect_equal(prov$accepted, k)
  expect_false(prov$compressed)

  cert <- certificate(warm)
  expect_true(all(cert$converged))
  expect_lte(cert$max_backward_error, 1e-6)
  # Warm converges in fewer operator applications than cold under mid-sweep.
  expect_lt(warm$matvecs, cold$matvecs)
})

test_that("mid-sweep checks return the correct VALUE SET on a degenerate cluster", {
  # A multiplicity-3 cluster (strictly greater than block for block in {1,2})
  # with k = 4 spanning the whole cluster. A small mid-sweep subspace resolves
  # only two of the three cluster members, so a naive early lock returns the
  # WRONG set {100,100,60,35} with genuine per-pair residuals and a passing
  # certificate. The cluster-clearance guard must defer such a window to the
  # full-subspace sweep boundary, reproducing the legacy (check_stride = 0) set.
  n <- 200L
  mult <- 3L
  k <- 4L
  A <- degenerate_cluster_sym(n, mult, seed = 3L)
  oracle <- sort(eigen(A, symmetric = TRUE)$values, decreasing = TRUE)[seq_len(k)]

  for (blk in c(1L, 2L)) {
    m_max <- if (blk == 1L) 60L else 80L
    for (cs in c(0L, 1L, 2L)) {
      fit <- eig_partial(
        A, k = k, target = largest(),
        method = lanczos(max_subspace = m_max, block = blk, check_stride = cs),
        seed = 7, tol = 1e-8
      )
      # Assert the VALUE SET, not merely per-pair residual cleanliness.
      expect_equal(
        sort(values(fit), decreasing = TRUE), oracle, tolerance = 1e-6,
        info = sprintf("block=%d check_stride=%d", blk, cs)
      )
      expect_certificate_clean(fit, tol = 1e-6)
    }
  }
})

test_that("mid-sweep returns the correct set for a one-copy-resolved miss (block < mult)", {
  # An isolated top eigenvalue followed by a multiplicity-2 near-degenerate pair
  # at the SECOND position, k spanning the pair, block=1 (< 2), tight m_max=40.
  # A small mid-sweep subspace resolves only one copy of the pair, so the window
  # {200,100,60,35} looks complete with NO tight pair anywhere -- neither the
  # cluster-clearance cut nor a workspace-complement probe can see the missing
  # copy (it is weakly present, not absent). The deflated-complement check finds
  # the second 100 (it survives deflation of the window) and defers. Verify the
  # returned VALUE SET matches the dense oracle across seeds, for every stride.
  n <- 200L
  k <- 4L
  lam_top <- c(200, 100, 100 + 1e-10, 60, 35)
  for (seed in c(1L, 2L, 3L, 7L, 11L)) {
    A <- spectrum_sym(lam_top, n, seed = seed)
    oracle <- sort(eigen(A, symmetric = TRUE)$values, decreasing = TRUE)[seq_len(k)]
    for (cs in c(0L, 1L, 2L)) {
      fit <- eig_partial(
        A, k = k, target = largest(),
        method = lanczos(max_subspace = 40L, block = 1L, check_stride = cs),
        seed = seed, tol = 1e-8
      )
      expect_equal(
        sort(values(fit), decreasing = TRUE), oracle, tolerance = 1e-6,
        info = sprintf("seed=%d check_stride=%d", seed, cs)
      )
    }
  }
})

test_that("mid-sweep returns the correct set for a multiplicity-2 pair at the top (block=1)", {
  # The seed-dependent variant: the near-degenerate pair sits at the very top.
  n <- 200L
  k <- 4L
  lam_top <- c(100, 100 + 1e-10, 60, 35, 20)
  for (seed in c(1L, 2L, 3L, 7L, 11L)) {
    A <- spectrum_sym(lam_top, n, seed = seed)
    oracle <- sort(eigen(A, symmetric = TRUE)$values, decreasing = TRUE)[seq_len(k)]
    for (cs in c(0L, 1L, 2L)) {
      fit <- eig_partial(
        A, k = k, target = largest(),
        method = lanczos(max_subspace = 40L, block = 1L, check_stride = cs),
        seed = seed, tol = 1e-8
      )
      expect_equal(
        sort(values(fit), decreasing = TRUE), oracle, tolerance = 1e-6,
        info = sprintf("seed=%d check_stride=%d", seed, cs)
      )
    }
  }
})

test_that("complement probe re-seeds past a step-1 breakdown (axis-aligned isolated eigenvalue)", {
  # The first deflated-complement probe seed is the deterministic coordinate
  # e_17 (0-based) = column 18 (1-based). Place an isolated eigenvalue (3) exactly
  # on that coordinate (decoupled), with the hidden {100,100+1e-10} near-degenerate
  # pair elsewhere. A single-run probe seeded there breaks down at step 1, sees
  # only "3", and would wrongly certify the window {200,100,60,35}. The probe must
  # deflate that exhausted coordinate, re-seed, find the second 100, and defer.
  n <- 200L
  k <- 4L
  axis1 <- ((1L * 17L + 0L) %% n) + 1L        # first probe coordinate, 1-based
  build <- function(seed) {
    set.seed(seed)
    rest <- setdiff(seq_len(n), axis1)
    lam_rest <- c(
      200, 100, 100 + 1e-10, 60, 35,
      rep(1, length(rest) - 5L) * exp(-seq_len(length(rest) - 5L) / 5)
    )
    Qr <- qr.Q(qr(matrix(stats::rnorm((n - 1L) * (n - 1L)), n - 1L)))
    M <- Qr %*% diag(lam_rest) %*% t(Qr)
    M <- (M + t(M)) / 2
    A <- matrix(0, n, n)
    A[rest, rest] <- M
    A[axis1, axis1] <- 3                       # eigenvector e_axis1, eigenvalue 3
    (A + t(A)) / 2
  }
  for (seed in 1:8) {
    A <- build(seed)
    oracle <- sort(eigen(A, symmetric = TRUE)$values, decreasing = TRUE)[seq_len(k)]
    for (cs in c(0L, 1L, 2L)) {
      fit <- eig_partial(
        A, k = k, target = largest(),
        method = lanczos(max_subspace = 40L, block = 1L, check_stride = cs),
        seed = seed, tol = 1e-8
      )
      expect_equal(
        sort(values(fit), decreasing = TRUE), oracle, tolerance = 1e-6,
        info = sprintf("seed=%d check_stride=%d", seed, cs)
      )
    }
  }
})

test_that("diagonal operator with a tight pair at the window edge defers to the boundary", {
  # Every eigenvector of a diagonal operator is a coordinate axis, so every probe
  # seed breaks down at step one and no segment is ever mixed. With the complement
  # never certifiable, the probe must return inconclusive (defer) and the full
  # sweep boundary reproduces the legacy set -- including a tight pair straddling
  # the k-th / (k+1)-th boundary.
  n <- 200L
  k <- 4L
  set.seed(101)
  diagvals <- c(
    100, 80, 60, 40, 40 + 1e-10, 20,
    sort(stats::runif(n - 6L, 0.1, 10), decreasing = TRUE)
  )
  A <- diag(diagvals)
  oracle <- sort(diagvals, decreasing = TRUE)[seq_len(k)]
  for (seed in 1:5) {
    for (cs in c(0L, 1L)) {
      fit <- eig_partial(
        A, k = k, target = largest(),
        method = lanczos(max_subspace = 40L, block = 1L, check_stride = cs),
        seed = seed, tol = 1e-8
      )
      expect_equal(
        sort(values(fit), decreasing = TRUE), oracle, tolerance = 1e-6,
        info = sprintf("seed=%d check_stride=%d", seed, cs)
      )
    }
  }
})
