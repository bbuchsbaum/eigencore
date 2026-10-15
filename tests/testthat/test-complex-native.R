# Native complex kernels for complex Hermitian problems: the zhetrf inertia
# count (R/complex_inertia.R, src/complex_inertia.cpp) and the scalar
# thick-restart Lanczos for matrix-free complex operators
# (R/complex_lanczos.R).

complex_native_fixture <- function(n, lambda = NULL, seed = 1L) {
  set.seed(seed)
  Z <- matrix(complex(real = rnorm(n * n), imaginary = rnorm(n * n)), n)
  if (is.null(lambda)) {
    return((Z + Conj(t(Z))) / 2)
  }
  Q <- qr.Q(qr(Z))
  H <- Q %*% (lambda * Conj(t(Q)))
  (H + Conj(t(H))) / 2
}

complex_native_expect_counts <- function(r, ev, info = "") {
  expect_true(r$reliable, info = info)
  lo <- r$zero_band[[1L]]
  hi <- r$zero_band[[2L]]
  expect_identical(r$below, sum(ev < lo), info = info)
  expect_identical(r$zero, sum(ev >= lo & ev <= hi), info = info)
  expect_identical(r$above, sum(ev > hi), info = info)
}

test_that("native complex counts match eigen() on random, repeated and near-singular shifts", {
  for (n in c(1L, 2L, 3L, 17L, 60L)) {
    H <- complex_native_fixture(n, seed = n)
    ev <- eigen(H, only.values = TRUE)$values
    shifts <- c(min(ev) - 1, max(ev) + 1, mean(ev), ev[[1L]], ev[[n]],
                ev[[ceiling(n / 2)]] + 1e-10)
    for (s in shifts) {
      r <- eigen_count(H, s)
      expect_identical(r$method, "dense_bunch_kaufman_hermitian")
      complex_native_expect_counts(r, ev, info = paste("random n", n, "sigma", s))
    }
  }
  lambda <- rep(c(-2, 0.5, 3, 7), times = c(5L, 9L, 8L, 3L))
  H <- complex_native_fixture(length(lambda), lambda, seed = 4L)
  for (s in c(-2, 0.5, 3, 7, 0, 1, 0.5 + 1e-9, 3 - 1e-12, 8)) {
    complex_native_expect_counts(eigen_count(H, s), lambda,
                                 info = paste("repeated sigma", s))
  }
  # Graded spectrum: shifts at or next to tiny eigenvalues.
  lambda <- c(10^seq(-14, 0, length.out = 15), -10^seq(-12, 2, length.out = 15))
  H <- complex_native_fixture(length(lambda), lambda, seed = 5L)
  ev <- eigen(H, only.values = TRUE)$values
  for (s in c(0, 1e-15, -1e-13, 1e-8, lambda[[3L]], lambda[[20L]])) {
    complex_native_expect_counts(eigen_count(H, s), ev, info = paste("graded sigma", s))
  }
})

test_that("native complex counts handle pencils and agree with the embedding count", {
  n <- 20L
  H <- complex_native_fixture(n, seed = 11L)
  set.seed(12)
  Y <- matrix(complex(real = rnorm(n * n), imaginary = rnorm(n * n)), n)
  B <- Y %*% Conj(t(Y)) + diag(n)
  B <- (B + Conj(t(B))) / 2
  emb <- function(M) rbind(cbind(Re(M), -Im(M)), cbind(Im(M), Re(M)))
  L <- chol(emb(B))
  C <- backsolve(L, t(backsolve(L, emb(H), transpose = TRUE)), transpose = TRUE)
  gev <- sort(eigen((C + t(C)) / 2, only.values = TRUE)$values)[c(TRUE, FALSE)]
  Breal <- crossprod(matrix(rnorm(n * n), n)) + diag(n)
  for (s in c(-1, 0, 0.3, 2)) {
    r <- eigen_count(H, s, B = B)
    expect_identical(r$method, "dense_bunch_kaufman_hermitian")
    expect_true(r$reliable)
    expect_identical(r$below, sum(gev < s))
    r2 <- withr::with_options(list(eigencore.complex_inertia_native = FALSE),
                              eigen_count(H, s, B = Breal))
    r1 <- eigen_count(H, s, B = Breal)
    expect_identical(c(r1$below, r1$above), c(r2$below, r2$above))
  }
  expect_error(eigen_count(H, 0, B = -B), "positive definite")
})

test_that("realified dense solves count natively and verify completeness", {
  n <- 160L
  H <- complex_native_fixture(n, seed = 21L)
  ev <- eigen(H, only.values = TRUE)$values
  fit <- eig_partial(H, k = 4, target = largest(), method = lanczos(block = 2))
  cert <- certificate(fit)
  expect_true(cert$passed)
  expect_equal(sort(values(fit)), sort(ev[1:4]), tolerance = 1e-9)
  expect_identical(cert$target_completeness, "inertia_verified")
  expect_identical(cert$completeness$inertia_method, "dense_bunch_kaufman_hermitian")
  # The context counts the real 2n embedding (counts doubled).
  P2 <- eigencore:::complex_hermitian_realified_problem(eigencore:::eigen_problem(H, target = largest()))
  ctx <- eigencore:::inertia_context(P2$A)
  expect_identical(ctx$n, 2L * n)
  expect_identical(ctx$method, "dense_bunch_kaufman_hermitian")
  t0 <- eigencore:::inertia_at(ctx, 0)
  expect_identical(t0$neg, 2 * sum(ev < 0))
  # Dense complex pencil with a real SPD metric.
  set.seed(22)
  Bm <- crossprod(matrix(rnorm(n * n), n)) / n + diag(n)
  fitg <- eig_partial(H, B = Bm, k = 3, target = smallest())
  certg <- certificate(fitg)
  expect_true(certg$passed)
  gev <- sort(Re(eigen(solve(Bm, H), only.values = TRUE)$values))
  expect_equal(sort(values(fitg)), gev[1:3], tolerance = 1e-8)
})

complex_native_linop <- function(H, counter = NULL) {
  force(H)
  linear_operator(
    dim = dim(H),
    apply = function(X, alpha = 1, beta = 0, Y = NULL) {
      if (!is.null(counter)) counter$cols <- counter$cols + ncol(as.matrix(X))
      out <- alpha * (H %*% X)
      if (is.null(Y) || beta == 0) out else out + beta * Y
    },
    dtype = "complex",
    structure = hermitian()
  )
}

test_that("matrix-free complex operators run the native complex Lanczos", {
  n <- 200L
  H <- complex_native_fixture(n, seed = 31L)
  ev <- eigen(H, only.values = TRUE)$values
  for (tg in list(largest(), smallest(), largest_magnitude())) {
    counter <- new.env()
    counter$cols <- 0
    fit <- eig_partial(complex_native_linop(H, counter), k = 5, target = tg)
    want <- switch(tg$kind, largest = ev[1:5], smallest = rev(ev)[1:5],
                   largest_magnitude = ev[order(abs(ev), decreasing = TRUE)][1:5])
    cert <- certificate(fit)
    expect_identical(fit$method,
                     "native complex Hermitian thick-restart Lanczos (matrix-free)")
    expect_true(cert$passed, info = tg$kind)
    expect_equal(sort(values(fit)), sort(want), tolerance = 1e-9, info = tg$kind)
    expect_true(cert$target_completeness %in% c("inertia_verified", "probed", "repaired"),
                info = tg$kind)
    expect_true(is.complex(vectors(fit)))
    # Fewer operator columns than the realified route.
    counter2 <- new.env()
    counter2$cols <- 0
    fit2 <- withr::with_options(
      list(eigencore.complex_native_lanczos = FALSE),
      eig_partial(complex_native_linop(H, counter2), k = 5, target = tg))
    expect_identical(fit2$method, "native realified complex Hermitian route (real 2n embedding)")
    expect_lt(counter$cols, counter2$cols)
  }
})

test_that("native complex Lanczos finds repeated eigenvalues and honours options", {
  lambda <- c(10, 10, 10, 9, 9, seq(8, 0, length.out = 115))
  H <- complex_native_fixture(length(lambda), lambda, seed = 41L)
  op <- complex_native_linop(H)
  fit <- eig_partial(op, k = 6, target = largest())
  expect_true(certificate(fit)$passed)
  expect_equal(sort(values(fit)), c(8, 9, 9, 10, 10, 10), tolerance = 1e-9)
  # values only, a warm start, lanczos() and a tight tolerance
  fit <- eig_partial(op, k = 3, target = smallest(), vectors = FALSE,
                     initial_subspace = matrix(1 + 0i, length(lambda), 1L))
  expect_null(fit$vectors)
  expect_equal(sort(values(fit)), sort(lambda)[1:3], tolerance = 1e-9)
  fit <- eig_partial(op, k = 2, target = largest(), method = lanczos(), tol = 1e-12)
  expect_true(certificate(fit)$passed)
  # Unsupported targets and block methods keep the realified route.
  fit <- eig_partial(op, k = 2, target = largest(), method = lanczos(block = 2))
  expect_identical(fit$method, "native realified complex Hermitian route (real 2n embedding)")
})
