# Tranche 4b SVD regressions: C37 (wide smallest-target Golub-Kahan),
# C42 (retained IRLBA native attempt robustness) and the P4 Golub-Kahan norm
# switch.

t4b_svd_fixture <- function(m, n, dmin, seed) {
  set.seed(seed)
  k <- min(m, n)
  U <- qr.Q(qr(matrix(stats::rnorm(m * k), m, k)))
  V <- qr.Q(qr(matrix(stats::rnorm(n * k), n, k)))
  d <- 10^seq(0, log10(dmin), length.out = k)
  list(A = U %*% (d * t(V)), d = d)
}

t4b_true_backward_error <- function(A, d, u, v) {
  left <- as.matrix(A %*% v) - sweep(u, 2L, d, "*")
  right <- as.matrix(Matrix::crossprod(A, u)) - sweep(v, 2L, d, "*")
  max(sqrt(colSums(left^2) + colSums(right^2))) / max(d, svd(as.matrix(A), 0, 0)$d[[1L]])
}

test_that("C37: wide smallest-target Golub-Kahan certifies like the tall case", {
  x <- t4b_svd_fixture(120L, 900L, 1e-6, 3701L)
  for (sparse in c(FALSE, TRUE)) {
    for (orientation in c("wide", "tall")) {
      A <- if (identical(orientation, "wide")) x$A else t(x$A)
      if (sparse) A <- methods::as(A, "dgCMatrix")
      fit <- svd_partial(
        A,
        rank = 3L,
        target = smallest(),
        method = golub_kahan(),
        tol = 1e-8,
        seed = 37L
      )
      expect_true(fit$certificate$passed, info = paste(orientation, sparse))
      expect_equal(sort(fit$d), sort(tail(x$d, 3L)), tolerance = 1e-8,
                   info = paste(orientation, sparse))
      expect_identical(
        isTRUE(fit$restart$internal_transposed),
        identical(orientation, "wide"),
        info = paste(orientation, sparse)
      )
      # Certificate is on the original A, and agrees with a direct check.
      expect_lt(t4b_true_backward_error(A, fit$d, fit$u, fit$v), 1e-8)
      expect_equal(dim(fit$u), c(nrow(A), 3L))
      expect_equal(dim(fit$v), c(ncol(A), 3L))
    }
  }
})

test_that("C37: auto sparse wide smallest SVD certifies through the fallback", {
  x <- t4b_svd_fixture(150L, 1200L, 1e-6, 3702L)
  A <- methods::as(x$A, "dgCMatrix")
  fit <- svd_partial(A, rank = 3L, target = smallest(), tol = 1e-8, seed = 38L)
  expect_true(fit$certificate$passed)
  expect_equal(sort(fit$d), sort(tail(x$d, 3L)), tolerance = 1e-8)
  expect_true(inherits(A, "dgCMatrix"))
  expect_lt(t4b_true_backward_error(A, fit$d, fit$u, fit$v), 1e-8)
})

test_that("C37: matrix-free wide smallest SVD runs on the adjoint view", {
  x <- t4b_svd_fixture(80L, 500L, 1e-5, 3703L)
  A <- x$A
  op <- linear_operator(
    dim = dim(A),
    apply = function(X, alpha = 1, beta = 0, Y = NULL) {
      out <- alpha * (A %*% X)
      if (!is.null(Y) && beta != 0) out <- out + beta * Y
      out
    },
    apply_adjoint = function(X, alpha = 1, beta = 0, Y = NULL) {
      out <- alpha * crossprod(A, X)
      if (!is.null(Y) && beta != 0) out <- out + beta * Y
      out
    }
  )
  fit <- eigencore:::native_golub_kahan_svd(op, rank = 2L, target = smallest(), tol = 1e-8)
  expect_true(fit$certificate$passed)
  expect_true(fit$restart$internal_transposed)
  expect_true(fit$restart$native_callback)
  expect_equal(sort(fit$d), sort(tail(x$d, 2L)), tolerance = 1e-8)
  expect_equal(dim(fit$u), c(80L, 2L))
  expect_lt(t4b_true_backward_error(A, fit$d, fit$u, fit$v), 1e-8)
})

test_that("C37: original-domain warm start maps onto the transposed problem", {
  x <- t4b_svd_fixture(60L, 300L, 1e-5, 3704L)
  set.seed(3705)
  fit <- eigencore:::native_golub_kahan_svd(
    x$A, rank = 2L, target = smallest(), tol = 1e-8,
    internal_start = stats::rnorm(300L)
  )
  expect_true(fit$restart$internal_transposed)
  expect_true(fit$restart$warm_started)
  expect_true(fit$certificate$passed)
})

test_that("C37: largest targets on wide inputs keep the as-given orientation", {
  x <- t4b_svd_fixture(60L, 300L, 1e-3, 3706L)
  fit <- eigencore:::native_golub_kahan_svd(x$A, rank = 2L, tol = 1e-8)
  expect_false(fit$restart$internal_transposed)
  expect_true(fit$certificate$passed)
})

test_that("C42: retained IRLBA native attempt certifies across seeds without fallback", {
  outcomes <- vapply(700:707, function(seed) {
    set.seed(seed)
    wide <- Matrix::t(Matrix::rsparsematrix(1500L, 160L, density = 0.02))
    set.seed(seed)
    fit <- eigencore:::native_irlba_lbd_retained_svd(
      wide, rank = 5L, work = 12L, retained = 7L, max_restarts = 7L,
      tol = 1e-8, vectors = "both"
    )
    c(
      certified = isTRUE(fit$certificate$passed),
      native = isTRUE(fit$restart$irlba_lbd_retained_native_attempted) &&
        !isTRUE(fit$restart$fallback_used),
      bounded = fit$matvecs <= 600L
    )
  }, logical(3))
  expect_true(all(outcomes["certified", ]))
  expect_true(all(outcomes["native", ]))
  expect_true(all(outcomes["bounded", ]))
})

test_that("C42: thick restart is what rescues the hard retained IRLBA case", {
  set.seed(701)
  wide <- Matrix::t(Matrix::rsparsematrix(2000L, 200L, density = 0.02))
  run <- function(thick_restarts) {
    set.seed(701)
    eigencore:::native_irlba_lbd_retained_svd(
      wide, rank = 5L, work = 12L, retained = 7L, max_restarts = 7L,
      tol = 1e-8, vectors = "both", thick_restarts = thick_restarts
    )
  }
  restarted <- run(NULL)
  expect_true(restarted$certificate$passed)
  expect_false(restarted$restart$fallback_used)
  expect_gte(restarted$restart$irlba_lbd_thick_restarts, 1L)
  expect_gt(restarted$restart$irlba_lbd_thick_restart_keep, 5L)
  oracle <- svd(as.matrix(wide), 0, 0)$d[1:5]
  expect_equal(restarted$d, oracle, tolerance = 1e-8)

  unrestarted <- run(0L)
  expect_equal(unrestarted$restart$irlba_lbd_thick_restarts, 0L)
  expect_true(unrestarted$restart$fallback_used)
  # The fallback itself must still certify (it used to return 9.3, 0, 0, 0, 0
  # after breaking down on its converged warm start).
  expect_true(unrestarted$certificate$passed)
  expect_equal(unrestarted$d, oracle, tolerance = 1e-8)
})

test_that("C42: Golub-Kahan restarts from a fresh vector after an invariant-subspace breakdown", {
  set.seed(702)
  wide <- Matrix::t(Matrix::rsparsematrix(800L, 120L, density = 0.03))
  s <- svd(as.matrix(wide))
  for (reorthogonalize in c(FALSE, TRUE)) {
    # Active domain: one-sided runs transpose the wide operator.
    start <- if (reorthogonalize) s$v[, 1L] else s$u[, 1L]
    fit <- eigencore:::native_golub_kahan_svd(
      wide, rank = 4L, tol = 1e-8,
      reorthogonalize = reorthogonalize, internal_start = start
    )
    expect_true(fit$certificate$passed, info = reorthogonalize)
    expect_gte(fit$restart$breakdown_restarts, 1L)
    expect_equal(fit$d, s$d[1:4], tolerance = 1e-8)
  }
})

test_that("P4: Golub-Kahan with BLAS norms matches the dense oracle", {
  set.seed(704)
  A <- Matrix::rsparsematrix(700L, 140L, density = 0.03)
  oracle <- svd(as.matrix(A), 0, 0)$d
  fit <- svd_partial(A, rank = 6L, method = golub_kahan(), tol = 1e-10, seed = 4L)
  expect_true(fit$certificate$passed)
  expect_equal(fit$d, oracle[1:6], tolerance = 1e-10)
  small <- svd_partial(A, rank = 3L, target = smallest(), method = golub_kahan(),
                       tol = 1e-10, seed = 5L)
  expect_true(small$certificate$passed)
  expect_equal(sort(small$d), sort(tail(oracle, 3L)), tolerance = 1e-9)
})
