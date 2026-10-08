# Tranche 2 Lanczos regressions: restart reuse (P2), relative convergence and
# breakdown tests (C16), pivoted tridiagonal shift-invert (C21) and the
# selected-eigenpair projected solves (P5).

t2_sparse_symmetric <- function(n, density, seed) {
  set.seed(seed)
  A <- Matrix::rsparsematrix(n, n, density = density)
  as(as(Matrix::forceSymmetric(A + Matrix::t(A)), "generalMatrix"), "CsparseMatrix")
}

t2_subspace_gap <- function(V, W) {
  # Largest principal-angle sine between two orthonormal column spaces.
  s <- svd(crossprod(V, W), nu = 0L, nv = 0L)$d
  sqrt(max(0, 1 - min(s)^2))
}

test_that("thick restart reuses kept Ritz images and agrees with eigen()", {
  S <- t2_sparse_symmetric(400L, 0.02, seed = 21)
  ref <- eigen(as.matrix(S), symmetric = TRUE)
  k <- 6L
  cases <- list(
    list(target = largest(), order = order(ref$values, decreasing = TRUE)),
    list(target = smallest(), order = order(ref$values)),
    list(target = largest_magnitude(),
         order = order(abs(ref$values), decreasing = TRUE))
  )
  for (case in cases) {
    fit <- eig_partial(S, k = k, target = case$target,
                       method = lanczos(max_subspace = 16L, max_restarts = 200L),
                       seed = 4, tol = 1e-10)
    expect_true(certificate(fit)$passed)
    expect_gt(fit$restart$restarts_used, 0L)
    idx <- case$order[seq_len(k)]
    expect_equal(sort(values(fit)), sort(ref$values[idx]), tolerance = 1e-9)
    gap <- t2_subspace_gap(vectors(fit), ref$vectors[, idx, drop = FALSE])
    expect_lt(gap, 1e-6)
    # P2: a restart applies the operator only to the continuation tail, never
    # again to the kept Ritz vectors, so every operator column is a Lanczos
    # step, a restart tail, the start block or a certification column.
    r <- fit$restart
    expect_equal(r$operator_columns - r$certification_operator_columns,
                 fit$iterations + r$restarts_used + 1L)
  }
})

test_that("block thick restart (block = 2) agrees with eigen() after restarts", {
  S <- t2_sparse_symmetric(300L, 0.03, seed = 22)
  ref <- eigen(as.matrix(S), symmetric = TRUE, only.values = TRUE)$values
  fit <- eig_partial(S, k = 4L, target = largest(),
                     method = lanczos(block = 2L, max_subspace = 14L,
                                      max_restarts = 200L),
                     seed = 3, tol = 1e-10)
  expect_equal(fit$restart$block, 2L)
  expect_gt(fit$restart$restarts_used, 0L)
  expect_true(certificate(fit)$passed)
  expect_equal(values(fit), sort(ref, decreasing = TRUE)[1:4], tolerance = 1e-9)
})

test_that("dense thick restart with repeated eigenvalues agrees with eigen()", {
  A <- symmetric_with_spectrum(c(9, 9, 7, 5, 5, 5, seq(4, -4, length.out = 114)),
                               seed = 23)
  ref <- eigen(A, symmetric = TRUE, only.values = TRUE)$values
  fit <- eig_partial(A, k = 5L, target = largest(),
                     method = lanczos(max_subspace = 12L, max_restarts = 300L),
                     seed = 5, tol = 1e-10)
  expect_true(certificate(fit)$passed)
  # Values come back in locking order (a later-resolved degenerate copy locks
  # after the next distinct value; unchanged pre-T2 behaviour), so compare sets.
  expect_equal(sort(values(fit), decreasing = TRUE), ref[1:5], tolerance = 1e-8)
})

test_that("eigs_sym is invariant under scaling the matrix by 1e-12 and 1e12", {
  S <- t2_sparse_symmetric(500L, 0.01, seed = 11)
  set.seed(3)
  D <- crossprod(matrix(rnorm(200 * 200), 200)) / 200
  for (M in list(S, D)) {
    for (which in c("LA", "SA", "LM")) {
      fits <- lapply(c(1, 1e-12, 1e12), function(s) {
        set.seed(5)
        eigencore::eigs_sym(M * s, 5L, which)
      })
      base <- fits[[1L]]
      expect_equal(base$nconv, 5L)
      scale_ref <- max(abs(base$values))
      for (i in 2:3) {
        s <- c(1, 1e-12, 1e12)[i]
        expect_equal(fits[[i]]$nconv, base$nconv)
        expect_lt(max(abs(fits[[i]]$values / s - base$values)) / scale_ref, 1e-9)
      }
    }
  }
})

test_that("scaled operators keep relative breakdown and convergence decisions", {
  A <- diag(c(10, 8, 6, rep(1, 57)))
  for (s in c(1e-12, 1, 1e12)) {
    fit <- eig_partial(A * s, k = 3L, target = largest(), seed = 27, tol = 1e-10)
    expect_true(certificate(fit)$passed)
    expect_equal(values(fit) / s, c(10, 8, 6), tolerance = 1e-10)
  }
})

test_that("native scalar Lanczos convergence flags are scale invariant", {
  set.seed(8)
  A <- symmetric_with_spectrum(seq(-3, 12, length.out = 60), seed = 8)
  start <- rnorm(60)
  run <- function(s) {
    .Call("eigencore_lanczos_dense", A * s, 40L, start, 4L, 1L, 1e-10,
          PACKAGE = "eigencore")
  }
  base <- run(1)
  for (s in c(1e-12, 1e12)) {
    scaled <- run(s)
    expect_identical(scaled$iterations, base$iterations)
    expect_identical(scaled$history_nconv, base$history_nconv)
    expect_equal(scaled$alpha / s, base$alpha, tolerance = 1e-12)
  }
})

test_that("tridiagonal shift-invert pivots through a zero leading pivot", {
  # Zero diagonal after the shift: the unpivoted Thomas recurrence hits a zero
  # first pivot, while the shifted matrix itself is well conditioned (n even:
  # eigenvalues 2 cos(j pi / (n + 1)) stay away from 0).
  n <- 40L
  sigma <- 0.3
  A <- Matrix::bandSparse(
    n, k = c(-1, 0, 1),
    diagonals = list(rep(1, n - 1L), rep(sigma, n), rep(1, n - 1L))
  )
  A <- as(A, "CsparseMatrix")
  expected <- eigen(as.matrix(A), symmetric = TRUE, only.values = TRUE)$values
  expected <- expected[order(abs(expected - sigma))][1:3]

  fit <- eig_partial(A, k = 3L, target = nearest(sigma),
                     method = shift_invert(sigma = sigma), tol = 1e-10,
                     allow_dense_fallback = "never")
  expect_identical(fit$method,
                   eigencore:::native_tridiagonal_shift_invert_label())
  expect_false(isTRUE(fit$transform$sigma_perturbed))
  expect_false(any(grepl("perturbed", fit$warnings)))
  expect_equal(fit$transform$factorization_cache$factorization,
               "LAPACK dgttrf/dgttrs")
  expect_equal(sort(values(fit)), sort(expected), tolerance = 1e-9)
  expect_certificate_clean(fit)
})

test_that("tridiagonal shift-invert is accurate at a pivot-hostile tiny pivot", {
  # sigma 1e-10 below the diagonal: Thomas pivots are 1e-10 and -1e10 (growth
  # 1e20), pivoted LU swaps rows and stays backward stable.
  n <- 40L
  d <- 0.3
  sigma <- d - 1e-10
  lower <- rep(1, n - 1L)
  diag_shifted <- rep(d - sigma, n)
  start <- rnorm(n)
  native <- .Call("eigencore_shift_invert_lanczos_tridiagonal",
                  lower, diag_shifted, lower, as.integer(n), start, 3L, 3L,
                  1e-12, PACKAGE = "eigencore")
  T_shift <- diag(diag_shifted)
  T_shift[cbind(2:n, 1:(n - 1L))] <- lower
  T_shift[cbind(1:(n - 1L), 2:n)] <- lower
  mu_ref <- 1 / eigen(T_shift, symmetric = TRUE, only.values = TRUE)$values
  mu_ref <- mu_ref[order(abs(mu_ref), decreasing = TRUE)][1:3]
  expect_equal(sort(native$ritz_values), sort(mu_ref), tolerance = 1e-10)
  expect_gt(native$factorization_cache$condition_estimate, sqrt(.Machine$double.eps))
  res <- T_shift %*% native$ritz_vectors -
    native$ritz_vectors %*% diag(1 / native$ritz_values)
  expect_lt(max(abs(res)), 1e-10)

  # Through the public API: no sigma perturbation, clean certificate.
  A <- Matrix::bandSparse(
    n, k = c(-1, 0, 1),
    diagonals = list(lower, rep(d, n), lower)
  )
  fit <- eig_partial(as(A, "CsparseMatrix"), k = 3L, target = nearest(sigma),
                     method = shift_invert(sigma = sigma), tol = 1e-10,
                     allow_dense_fallback = "never")
  expect_false(isTRUE(fit$transform$sigma_perturbed))
  expect_certificate_clean(fit)
})

test_that("tridiagonal shift-invert still rejects a singular shift", {
  lower <- c(1, 1)
  diag_shifted <- c(0, 0, 0)   # tridiag(1, 0, 1) of order 3 is singular
  expect_error(
    .Call("eigencore_shift_invert_lanczos_tridiagonal",
          lower, diag_shifted, lower, 3L, c(1, 2, 3), 1L, 3L, 1e-10,
          PACKAGE = "eigencore"),
    "zero pivot|near-singular"
  )
})

test_that("dense shift-invert uses Bunch-Kaufman and stays accurate for interior shifts", {
  vals <- c(-6, -2.5, -1, 0.2, 0.9, 2, 3.5, seq(5, 20, length.out = 33))
  A <- symmetric_with_spectrum(vals, seed = 31)
  sigma <- 0.5
  fit <- eig_partial(A, k = 3L, target = nearest(sigma),
                     method = shift_invert(sigma = sigma), tol = 1e-10)
  expected <- vals[order(abs(vals - sigma))][1:3]
  expect_identical(fit$transform$factorization_cache$factorization,
                   "LAPACK dsytrf/dsytrs")
  expect_equal(sort(values(fit)), sort(expected), tolerance = 1e-9)
  expect_certificate_clean(fit)
  # Exactly singular shift is still reported.
  expect_error(
    .Call("eigencore_shift_invert_lanczos_dense", diag(c(1, 2, 3)), 2, 3L,
          c(1, 1, 1), 1L, 3L, 1e-10, PACKAGE = "eigencore"),
    "singular"
  )
})
