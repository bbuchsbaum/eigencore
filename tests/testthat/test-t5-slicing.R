# Tranche 5 phase 2: interval(a, b) targets and spectrum slicing, the C60
# smallest-target LDL' route, reuse of the shift-invert factor for counts,
# and the native LDL' solve kernel.

# Option scoping without a withr dependency.
sl_local_options <- function(opts, .env = parent.frame()) {
  old <- options(opts)
  do.call(on.exit, list(call("options", old), add = TRUE), envir = .env)
  invisible(old)
}

sl_with_options <- function(opts, code) {
  old <- options(opts)
  on.exit(options(old), add = TRUE)
  force(code)
}

sl_orthogonal <- function(n, seed) {
  set.seed(seed)
  qr.Q(qr(matrix(stats::rnorm(n * n), n)))
}

sl_dense_with_spectrum <- function(d, seed) {
  Q <- sl_orthogonal(length(d), seed)
  A <- Q %*% (d * t(Q))
  (A + t(A)) / 2
}

# Sparse symmetric matrix with exactly known spectrum `d` (random 2 x 2
# rotations on disjoint index pairs of diag(d)).
sl_rotated_sparse <- function(d, seed) {
  n <- length(d)
  set.seed(seed)
  P <- sample(n)
  half <- n %/% 2L
  a <- P[2L * seq_len(half) - 1L]
  b <- P[2L * seq_len(half)]
  th <- stats::runif(half, 0, 2 * pi)
  ii <- c(a, a, b, b)
  jj <- c(a, b, a, b)
  xx <- c(cos(th), -sin(th), sin(th), cos(th))
  if (n %% 2L) {
    ii <- c(ii, P[n]); jj <- c(jj, P[n]); xx <- c(xx, 1)
  }
  R <- Matrix::sparseMatrix(i = ii, j = jj, x = xx, dims = c(n, n))
  S <- R %*% Matrix::Diagonal(n, d) %*% Matrix::t(R)
  Matrix::forceSymmetric(methods::as((S + Matrix::t(S)) / 2, "CsparseMatrix"))
}

sl_laplacian_2d <- function(m) {
  T1 <- Matrix::bandSparse(m, k = c(0, 1),
                           diagonals = list(rep(2, m), rep(-1, m - 1L)),
                           symmetric = TRUE)
  I <- Matrix::Diagonal(m)
  methods::as(kronecker(T1, I) + kronecker(I, T1), "generalMatrix")
}

# Exact spectra of the Dirichlet (and Neumann) 2-D grid Laplacians.
sl_laplacian_2d_values <- function(m, neumann = FALSE) {
  t1 <- if (neumann) 4 * sin(pi * (0:(m - 1L)) / (2 * m))^2 else
    4 * sin(pi * seq_len(m) / (2 * (m + 1)))^2
  sort(as.vector(outer(t1, t1, `+`)))
}

sl_expect_interval <- function(fit, expected, tol = 1e-8, info = "",
                               status = c("inertia_verified", "exact")) {
  v <- values(fit)
  expect_length(v, length(expected))
  if (length(expected)) {
    expect_equal(sort(v), sort(expected), tolerance = tol, info = info)
  }
  cert <- certificate(fit)
  expect_true(cert$passed, info = info)
  expect_true(cert$target_completeness %in% status, info = info)
}

test_that("interval() validates its end points and labels the target", {
  tg <- interval(1, 2)
  expect_s3_class(tg, "eigencore_target")
  expect_identical(tg$kind, "interval")
  expect_identical(tg$value, list(lower = 1, upper = 2))
  expect_match(eigencore:::target_label(tg), "^interval\\(1, 2\\)$")
  expect_error(interval(2, 1), "a < b")
  expect_error(interval(1, 1), "a < b")
  expect_error(interval(NA, 1), "single number")
  expect_error(interval(c(1, 2), 3), "single number")
  expect_error(interval(-Inf, Inf), "eig_full")
  expect_identical(interval(-Inf, 0)$value$lower, -Inf)
  expect_identical(interval(0, Inf)$value$upper, Inf)
})

test_that("k is optional only for interval targets and is a cap there", {
  A <- diag(c(1, 2, 3, 4, 5))
  expect_error(eig_partial(A), "k \\(the number of eigenpairs\\) is required")
  fit <- eig_partial(A, target = interval(1.5, 4.5))
  expect_equal(values(fit), c(2, 3, 4))
  expect_identical(fit$requested, 3L)
  fit <- eig_partial(A, k = 4, target = interval(1.5, 4.5))
  expect_equal(values(fit), c(2, 3, 4))
  expect_error(eig_partial(A, k = 2, target = interval(1.5, 4.5)),
               "holds 3 eigenvalues, more than k = 2")
  expect_error(eig_partial(A, target = interval(1, 2), method = lanczos()),
               "method = auto")
  expect_error(eig_partial(matrix(c(1, 2, 0, 1), 2), target = interval(0, 3)),
               "Hermitian")
})

test_that("dense interval targets match eigen() including edge cases", {
  d <- c(-3, -1, 0.5, 1, 1, 1, 2, 4, 4, 7, 10, 12)
  A <- sl_dense_with_spectrum(d, 11)
  cases <- list(
    c(0.75, 5), c(1, 4), c(4, 4.5), c(5, 6), c(-Inf, 1), c(4, Inf),
    c(-100, 100), c(9.99, 10.01)
  )
  for (ab in cases) {
    fit <- eig_partial(A, target = interval(ab[[1L]], ab[[2L]]))
    expect_identical(fit$method, eigencore:::interval_dense_label())
    expected <- d[d >= ab[[1L]] - 1e-9 & d <= ab[[2L]] + 1e-9]
    sl_expect_interval(fit, expected, info = paste(ab, collapse = ","),
                       status = "exact")
    if (length(expected)) {
      V <- fit$vectors
      expect_lt(max(abs(A %*% V - V %*% diag(values(fit), length(expected)))), 1e-10)
    }
  }
  # Repeated eigenvalue at both end points: all copies are inside.
  fit <- eig_partial(A, target = interval(1, 4))
  expect_equal(values(fit), c(1, 1, 1, 2, 4, 4), tolerance = 1e-10)
  fit <- eig_partial(A, target = interval(5, 6))
  expect_length(values(fit), 0L)
  expect_true(certificate(fit)$passed)
  expect_equal(dim(fit$vectors), c(12L, 0L))
  fit <- eig_partial(A, target = interval(0, 3), vectors = FALSE)
  expect_null(fit$vectors)
  expect_equal(values(fit), c(0.5, 1, 1, 1, 2), tolerance = 1e-10)
})

test_that("generalized and complex dense interval targets", {
  set.seed(4)
  A <- crossprod(matrix(rnorm(60 * 60), 60)) / 60
  B <- crossprod(matrix(rnorm(60 * 60), 60)) / 60 + diag(60)
  ev <- Re(eigen(solve(B, A), only.values = TRUE)$values)
  fit <- eig_partial(A, B = B, target = interval(0.2, 0.6))
  expect_identical(fit$method, eigencore:::interval_dense_generalized_label())
  sl_expect_interval(fit, ev[ev >= 0.2 & ev <= 0.6], status = "exact")
  X <- fit$vectors
  expect_lt(max(abs(crossprod(X, B %*% X) - diag(ncol(X)))), 1e-8)
  Bd <- Matrix::Diagonal(60, seq(1, 2, length.out = 60))
  evd <- Re(eigen(solve(as.matrix(Bd), A), only.values = TRUE)$values)
  fit <- eig_partial(A, B = Bd, target = interval(0.5, 1.5))
  sl_expect_interval(fit, evd[evd >= 0.5 & evd <= 1.5], status = "exact")

  set.seed(5)
  Z <- matrix(complex(real = rnorm(400), imaginary = rnorm(400)), 20)
  H <- (Z + Conj(t(Z))) / 2
  evh <- eigen(H, only.values = TRUE)$values
  fit <- eig_partial(H, target = interval(-2, 2))
  sl_expect_interval(fit, evh[evh >= -2 & evh <= 2], status = "exact")
})

test_that("sparse interval slicing matches the known spectrum", {
  sl_local_options(list(eigencore.interval_dense_limit = 0L,
                            eigencore.threads = 1L))
  d <- c(rep(1:12, each = 25), seq(0.05, 0.95, length.out = 40))
  S <- sl_rotated_sparse(d, 7)
  cases <- list(
    c(3, 6),            # end points exactly on 25-fold eigenvalues
    c(3.2, 3.8),        # empty
    c(4.5, 5.5),        # one eigenvalue of multiplicity 25
    c(0.4, 0.6),        # simple eigenvalues only
    c(-Inf, 1.5),       # one-sided
    c(11, Inf)
  )
  for (ab in cases) {
    fit <- eig_partial(S, target = interval(ab[[1L]], ab[[2L]]))
    expect_identical(fit$method, eigencore:::interval_slicing_label())
    expected <- d[d >= ab[[1L]] - 1e-9 & d <= ab[[2L]] + 1e-9]
    sl_expect_interval(fit, expected, tol = 1e-7, info = paste(ab, collapse = ","),
                       status = "inertia_verified")
    expect_identical(fit$interval$count, length(expected))
    if (length(expected)) {
      V <- fit$vectors
      expect_lt(max(abs(crossprod(V) - diag(ncol(V)))), 1e-8)
    }
  }
  # End points on eigenvalues: those values sit within the residual bound of
  # the boundary and are resolved by the widened recount, not dropped.
  fit <- eig_partial(S, target = interval(3, 6))
  rec <- certificate(fit)$completeness
  expect_true(length(rec$boundary_ambiguous) > 0L)
  expect_identical(rec$widened_count, 100)
})

test_that("a wide sparse interval is sliced and merged without duplicates or gaps", {
  sl_local_options(list(eigencore.interval_dense_limit = 0L,
                            eigencore.interval_slice_size = 60L,
                            eigencore.threads = 1L))
  S <- sl_laplacian_2d(55)  # n = 3025, double eigenvalues
  ev <- sl_laplacian_2d_values(55)
  a <- 1.2
  b <- 2.6
  expected <- ev[ev >= a & ev <= b]
  expect_gt(length(expected), 350L)
  fit <- eig_partial(S, target = interval(a, b))
  sl_expect_interval(fit, expected, tol = 1e-8, status = "inertia_verified")
  rec <- fit$interval
  expect_gt(rec$slices, 4L)
  expect_equal(sum(rec$slice_table$count), length(expected))
  expect_length(rec$centre_count_mismatch, 0L)
  V <- fit$vectors
  expect_lt(max(abs(crossprod(V) - diag(ncol(V)))), 1e-8)
  expect_true(all(rec$slice_table$native_solve[rec$slice_table$count > 0]))
  # The per-slice LDL' factors share one symbolic analysis.
  expect_gt(rec$factorizations, rec$slices)
})

test_that("generalized sparse interval slicing", {
  sl_local_options(list(eigencore.interval_dense_limit = 0L,
                            eigencore.threads = 1L))
  n <- 400
  set.seed(9)
  A <- Matrix::rsparsematrix(n, n, density = 0.01, symmetric = TRUE) +
    Matrix::Diagonal(n, 3)
  B <- Matrix::forceSymmetric(methods::as(Matrix::bandSparse(
    n, k = c(0, 1), diagonals = list(rep(4, n), rep(1, n - 1L)), symmetric = TRUE),
    "CsparseMatrix"))
  ev <- Re(eigen(solve(as.matrix(B), as.matrix(A)), only.values = TRUE)$values)
  fit <- eig_partial(A, B = B, target = interval(0.6, 0.9))
  expect_identical(fit$method, eigencore:::interval_slicing_label())
  sl_expect_interval(fit, ev[ev >= 0.6 & ev <= 0.9], tol = 1e-7)
  X <- fit$vectors
  expect_lt(max(abs(as.matrix(crossprod(X, B %*% X)) - diag(ncol(X)))), 1e-8)
  Bd <- Matrix::Diagonal(n, seq(1, 3, length.out = n))
  evd <- Re(eigen(solve(as.matrix(Bd), as.matrix(A)), only.values = TRUE)$values)
  fit <- eig_partial(A, B = Bd, target = interval(1, 2))
  sl_expect_interval(fit, evd[evd >= 1 & evd <= 2], tol = 1e-7)
})

test_that("interval slicing survives unreliable factorisations", {
  sl_local_options(list(eigencore.interval_dense_limit = 0L,
                            eigencore.threads = 1L))
  d <- c(rep(c(1, 2, 3), each = 10), seq(4, 6, length.out = 50))
  S <- sl_rotated_sparse(d, 3)
  expected <- d[d >= 1.5 & d <= 4.5]
  # Every LDL' probe solve rejected: the slices fall back to sparse LU.
  sl_with_options(list(eigencore.shift_invert_ldl_tol = 1e-300), {
    fit <- eig_partial(S, target = interval(1.5, 4.5))
  })
  sl_expect_interval(fit, expected, tol = 1e-7)
  expect_true(any(fit$interval$slice_table$ldl_fallback))
  expect_match(paste(fit$warnings, collapse = " "), "sparse LU")
  # Every sparse count declared unreliable: counts fall back to dense
  # Bunch-Kaufman (n <= eigencore.inertia_dense_fallback_limit).
  sl_with_options(list(eigencore.inertia_pivot_tol = 1), {
    fit <- eig_partial(S, target = interval(1.5, 4.5))
  })
  sl_expect_interval(fit, expected, tol = 1e-7)
  expect_match(fit$interval$count_method, "sparse_cholmod_ldl")
})

test_that("interval plans are executable and replannable", {
  A <- diag(c(1, 2, 3, 4))
  plan <- plan_solver(eigen_problem(A, target = interval(1.5, 3.5)))
  expect_identical(plan$method, eigencore:::interval_dense_label())
  expect_identical(plan$requested, 4L)
  expect_true(any(grepl("k inferred", plan$reasons)))
  fit <- solve(plan)
  expect_equal(values(fit), c(2, 3))
  fit <- solve(plan, replan = TRUE)
  expect_equal(values(fit), c(2, 3))
  expect_error(eig_partial(A, target = interval(1, 2), initial_subspace = diag(4)[, 1]),
               "initial_subspace")
})

test_that("nearest() inertia counts reuse the shift-invert factor", {
  sl_local_options(list(eigencore.threads = 1L))
  S <- sl_laplacian_2d(40)
  fit <- eig_partial(S, k = 5, target = nearest(2.3),
                     method = shift_invert(2.3, max_subspace = 60))
  expect_identical(certificate(fit)$target_completeness, "not_checked")
  sl_local_options(list(eigencore.target_completeness = "inertia"))
  fit <- eig_partial(S, k = 5, target = nearest(2.3))
  cert <- certificate(fit)
  expect_true(cert$target_completeness %in% c("inertia_verified", "inertia_inconclusive"))
  expect_true(cert$completeness$reused_symbolic)
  expect_null(fit$transform$inertia_seed)
  ev <- sl_laplacian_2d_values(40)
  expect_equal(sort(values(fit)), sort(ev[order(abs(ev - 2.3))][1:5]), tolerance = 1e-8)
  # A count at the factored shift itself is answered without refactoring.
  ctx <- eigencore:::inertia_context(S)
  seed <- list(sigma = 2.3, factor = Matrix::Cholesky(
    Matrix::forceSymmetric(S, uplo = "U"), LDL = TRUE, super = FALSE, Imult = -2.3),
    tally = list(ok = TRUE, sigma = 2.3, neg = 7, zero = 0, pos = 1593),
    positive_definite = FALSE)
  eigencore:::inertia_seed_context(ctx, seed)
  before <- ctx$factorizations
  expect_identical(eigencore:::inertia_at(ctx, 2.3)$neg, 7)
  expect_identical(ctx$factorizations, before)
  expect_identical(ctx$seed_hits, 1L)
})

test_that("C60: smallest targets route to LDL' shift-invert only when cheap and PD", {
  skip_if_not(eigencore:::cholmod_bridge_available())
  sl_local_options(list(eigencore.smallest_ldl_min_n = 1000,
                            eigencore.threads = 1L))
  S <- sl_laplacian_2d(45)  # n = 2025
  plan <- plan_solver(eigen_problem(S, target = smallest()), k = 6)
  expect_identical(plan$method, eigencore:::shift_invert_sparse_label())
  expect_true(any(grepl("C60", plan$reasons)))
  fit <- eig_partial(S, k = 6, target = smallest())
  ev <- sl_laplacian_2d_values(45)
  expect_equal(sort(values(fit)), ev[1:6], tolerance = 1e-9)
  expect_true(certificate(fit)$passed)
  expect_identical(fit$sigma, 0)
  expect_true(fit$transform$factorization_cache$positive_definite_factor)
  expect_true(fit$restart$native_solve)

  # High-fill random sparse matrix keeps the Lanczos route.
  set.seed(1)
  n <- 3000
  R <- methods::as(Matrix::rsparsematrix(n, n, density = 6 / n, symmetric = TRUE) +
                     Matrix::Diagonal(n, 10), "generalMatrix")
  plan <- plan_solver(eigen_problem(R, target = smallest()), k = 5)
  expect_false(grepl("shift-invert", plan$method))

  # Singular PSD (graph Laplacian with Neumann ends): shift slightly below 0.
  m <- 45
  T0 <- Matrix::bandSparse(m, k = c(0, 1),
                           diagonals = list(c(1, rep(2, m - 2L), 1), rep(-1, m - 1L)),
                           symmetric = TRUE)
  I <- Matrix::Diagonal(m)
  G <- methods::as(kronecker(T0, I) + kronecker(I, T0), "generalMatrix")
  fit <- eig_partial(G, k = 4, target = smallest())
  evg <- sl_laplacian_2d_values(m, neumann = TRUE)
  expect_equal(sort(values(fit)), evg[1:4], tolerance = 1e-8)
  expect_lt(fit$sigma, 0)
  expect_true(certificate(fit)$passed)

  # Indefinite despite a positive diagonal: the solve proves no PD shift and
  # falls back to the Lanczos route.
  Sx <- S - Matrix::Diagonal(nrow(S), 0.05)
  plan <- plan_solver(eigen_problem(Sx, target = smallest()), k = 4)
  expect_identical(plan$method, eigencore:::shift_invert_sparse_label())
  fit <- eig_partial(Sx, k = 4, target = smallest())
  expect_true(fit$fallback_used)
  expect_identical(fit$fallback_reason$code, "spd_shift_rejected")
  evx <- sl_laplacian_2d_values(45) - 0.05
  expect_equal(sort(values(fit)), evx[1:4], tolerance = 1e-7)
  # The route can be switched off.
  sl_local_options(list(eigencore.smallest_ldl_route = FALSE))
  plan <- plan_solver(eigen_problem(S, target = smallest()), k = 6)
  expect_false(grepl("shift-invert", plan$method))
})

test_that("native LDL' solve kernels match Matrix::solve", {
  S <- sl_laplacian_2d(30)
  As <- methods::as(Matrix::forceSymmetric(S, uplo = "U"), "CsparseMatrix")
  F <- Matrix::Cholesky(As, LDL = TRUE, super = FALSE, perm = TRUE, Imult = -1.37)
  set.seed(2)
  X <- matrix(rnorm(nrow(S) * 3), nrow(S), 3)
  Y0 <- matrix(rnorm(nrow(S) * 3), nrow(S), 3)
  ref <- as.matrix(Matrix::solve(F, X, system = "A"))
  modes <- c("eigencore", if (eigencore:::cholmod_bridge_available()) "cholmod")
  for (mode in modes) {
    sl_with_options(list(eigencore.native_ldl_solve = mode), {
      spec <- eigencore:::ldl_native_solve_spec(F)
      expect_identical(spec$type, if (mode == "cholmod") "cholmod_solve" else "ldl_solve")
      kernel <- eigencore:::new_native_composite_kernel(spec)
      expect_false(is.null(kernel))
      out <- eigencore:::native_composite_block_apply(kernel, X)
      expect_equal(out, ref, tolerance = 1e-12, info = mode)
      out <- eigencore:::native_composite_block_apply(kernel, X, alpha = 2,
                                                      beta = -0.5, Y = Y0)
      expect_equal(out, 2 * ref - 0.5 * Y0, tolerance = 1e-12, info = mode)
      out <- eigencore:::native_composite_block_apply(kernel, X[, 1, drop = FALSE],
                                                      adjoint = TRUE)
      expect_equal(out, ref[, 1, drop = FALSE], tolerance = 1e-12, info = mode)
    })
  }
  sl_with_options(list(eigencore.native_ldl_solve = FALSE), {
    expect_null(eigencore:::ldl_native_solve_spec(F))
  })
  # A supernodal (non-simplicial) factor gets no native leaf.
  Fs <- Matrix::Cholesky(As, LDL = FALSE, super = TRUE, Imult = 0)
  expect_null(eigencore:::ldl_native_solve_spec(Fs))
  expect_null(eigencore:::new_native_composite_kernel(list(
    type = "ldl_solve", p = 0:1, i = 0L, x = 0, nz = 1L, perm = NULL)))

  # Shift-invert solves through the native kernel agree with the R callback.
  fit_native <- eig_partial(S, k = 4, target = nearest(1.37),
                            method = shift_invert(1.37), seed = 1)
  expect_true(fit_native$restart$native_solve)
  sl_with_options(list(eigencore.native_ldl_solve = FALSE), {
    fit_r <- eig_partial(S, k = 4, target = nearest(1.37),
                         method = shift_invert(1.37), seed = 1)
  })
  expect_false(fit_r$restart$native_solve)
  expect_equal(sort(values(fit_native)), sort(values(fit_r)), tolerance = 1e-10)
  # Generalized: R (A - sigma B)^{-1} R' as a native product.
  B <- Matrix::forceSymmetric(methods::as(Matrix::bandSparse(
    nrow(S), k = c(0, 1), diagonals = list(rep(4, nrow(S)), rep(1, nrow(S) - 1L)),
    symmetric = TRUE), "CsparseMatrix"))
  fit_g <- eig_partial(S, B = B, k = 3, target = nearest(0.5),
                       method = shift_invert(0.5))
  expect_true(fit_g$restart$native_solve)
  expect_true(certificate(fit_g)$passed)
  evg <- Re(eigen(solve(as.matrix(B), as.matrix(S)), only.values = TRUE)$values)
  expect_equal(sort(values(fit_g)), sort(evg[order(abs(evg - 0.5))][1:3]),
               tolerance = 1e-8)
})
