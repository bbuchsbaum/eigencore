# Tranche 5 (capability gap 4): Sylvester inertia counting (eigen_count),
# CHOLMOD LDL' sparse shift-invert, and the deterministic inertia
# completeness certificate.

t5_orthogonal <- function(n, seed) {
  set.seed(seed)
  qr.Q(qr(matrix(stats::rnorm(n * n), n)))
}

t5_dense_with_spectrum <- function(d, seed) {
  Q <- t5_orthogonal(length(d), seed)
  A <- Q %*% (d * t(Q))
  (A + t(A)) / 2
}

# Sparse (not diagonal) symmetric matrix with exactly known spectrum: random
# 2x2 rotations on disjoint index pairs of diag(d).
t5_rotated_sparse <- function(d, seed) {
  n <- length(d)
  set.seed(seed)
  P <- sample(n)
  ii <- jj <- integer()
  xx <- numeric()
  for (t in seq_len(n %/% 2L)) {
    a <- P[2L * t - 1L]
    b <- P[2L * t]
    th <- stats::runif(1L, 0, 2 * pi)
    ii <- c(ii, a, a, b, b)
    jj <- c(jj, a, b, a, b)
    xx <- c(xx, cos(th), -sin(th), sin(th), cos(th))
  }
  if (n %% 2L) {
    ii <- c(ii, P[n]); jj <- c(jj, P[n]); xx <- c(xx, 1)
  }
  R <- Matrix::sparseMatrix(i = ii, j = jj, x = xx, dims = c(n, n))
  S <- R %*% Matrix::Diagonal(n, d) %*% Matrix::t(R)
  Matrix::forceSymmetric(methods::as((S + Matrix::t(S)) / 2, "CsparseMatrix"))
}

t5_laplacian_1d <- function(m) {
  Matrix::bandSparse(m, k = c(0, 1),
                     diagonals = list(rep(2, m), rep(-1, m - 1L)),
                     symmetric = TRUE)
}

t5_laplacian_2d <- function(m) {
  T1 <- t5_laplacian_1d(m)
  I <- Matrix::Diagonal(m)
  Matrix::forceSymmetric(methods::as(
    kronecker(T1, I) + kronecker(I, T1), "CsparseMatrix"))
}

# Counts must match the known spectrum `d` on the reported half-lines.
t5_expect_counts <- function(r, d, info = "") {
  expect_true(r$reliable, info = info)
  expect_identical(r$below, sum(d < r$zero_band[[1L]]), info = info)
  expect_identical(r$above, sum(d > r$zero_band[[2L]]), info = info)
  expect_identical(r$below + r$zero + r$above, length(d), info = info)
}

# Zero diagonal, couplings b_i between i and i + m (not tridiagonal):
# eigenvalues +/- b_i. The unpivoted LDL' of S - 0 I meets a zero pivot.
t5_zero_diagonal_pairs <- function(b) {
  m <- length(b)
  i <- seq_len(m)
  Matrix::sparseMatrix(i = c(i, i + m), j = c(i + m, i), x = c(b, b),
                       dims = c(2L * m, 2L * m))
}

with_t5_completeness <- function(mode, code) {
  old <- options(eigencore.target_completeness = mode)
  on.exit(options(old), add = TRUE)
  force(code)
}

test_that("eigen_count matches the spectrum of dense symmetric matrices", {
  d <- c(-3, -1, -1, 0.5, 2, 2, 2, 4, seq(5, 9, length.out = 32L))
  A <- t5_dense_with_spectrum(d, 11L)
  set.seed(1)
  sigmas <- c(stats::runif(25L, -4, 10), -3.5, 0, 1, 3, 4.5, 9.5)
  for (s in sigmas) {
    r <- eigen_count(A, s)
    expect_s3_class(r, "eigencore_inertia")
    expect_identical(r$method, "dense_bunch_kaufman")
    expect_false(r$perturbed)
    t5_expect_counts(r, d, info = paste("sigma", s))
    expect_identical(r$below, sum(d < s))
  }
})

test_that("shifts at (repeated) eigenvalues are perturbed into a zero band", {
  d <- c(-3, -1, -1, 0.5, 2, 2, 2, 4, seq(5, 9, length.out = 32L))
  A <- t5_dense_with_spectrum(d, 12L)
  for (s in c(-1, 2, 4, 9)) {
    r <- eigen_count(A, s)
    t5_expect_counts(r, d, info = paste("sigma", s))
    expect_true(r$perturbed)
    expect_identical(r$zero, sum(d == s))
    expect_lt(r$perturbation, 1e-6 * 20)
    expect_true(any(grepl("zero band", r$notes)))
  }
  # Without perturbation an at-eigenvalue shift is reported unreliable.
  r <- eigen_count(A, 2, perturb = FALSE)
  expect_false(r$reliable)
  expect_true(any(grepl("not certified", r$notes)))
  # Diagonal input: exact.
  r <- eigen_count(Matrix::Diagonal(x = d), 2)
  expect_identical(r$method, "diagonal")
  expect_true(r$reliable)
  expect_identical(c(r$below, r$zero, r$above), c(4L, 3L, 33L))
})

test_that("eigen_count matches sparse symmetric matrices through CHOLMOD LDL'", {
  set.seed(3)
  d <- sort(c(stats::rnorm(150L), 0.25, 0.25, 0.25, -2, -2))
  S <- t5_rotated_sparse(d, 31L)
  for (s in c(stats::runif(15L, -3, 3), -2, 0.25, min(d) - 1, max(d) + 1)) {
    r <- eigen_count(S, s)
    expect_identical(r$method, "sparse_cholmod_ldl")
    t5_expect_counts(r, d, info = paste("sigma", s))
    expect_gt(r$diagnostics$factor_nnz, 0)
  }
  # Random sparse with eigen() as the oracle; symmetric dgCMatrix input.
  set.seed(4)
  R <- Matrix::rsparsematrix(400L, 400L, density = 0.01, symmetric = TRUE,
                             rand.x = stats::rnorm)
  G <- methods::as(methods::as(R, "generalMatrix"), "CsparseMatrix")
  ev <- eigen(as.matrix(R), symmetric = TRUE, only.values = TRUE)$values
  for (s in c(-1.7, -0.3, 0.05, 0.9, 2.2)) {
    r <- eigen_count(G, s)
    expect_true(r$reliable)
    expect_identical(r$below, sum(ev < r$zero_band[[1L]]))
    expect_identical(r$above, sum(ev > r$zero_band[[2L]]))
  }
})

test_that("tridiagonal and diagonal inputs use the Sturm / exact fast paths", {
  m <- 300L
  T1 <- t5_laplacian_1d(m)
  lam <- 2 - 2 * cos(pi * seq_len(m) / (m + 1))
  for (s in c(0.001, 0.5, 1, 2, 3.3, 3.99)) {
    r <- eigen_count(T1, s)
    expect_identical(r$method, "tridiagonal_sturm")
    t5_expect_counts(r, lam, info = paste("sigma", s))
  }
  # Generalized with diagonal B keeps the Sturm path.
  b <- seq(1, 2, length.out = m)
  Bd <- Matrix::Diagonal(x = b)
  Ad <- as.matrix(T1)
  geo <- eigen(diag(1 / sqrt(b)) %*% Ad %*% diag(1 / sqrt(b)), symmetric = TRUE,
               only.values = TRUE)$values
  r <- eigen_count(T1, 1.1, B = Bd)
  expect_identical(r$method, "tridiagonal_sturm")
  expect_true(r$generalized)
  expect_identical(r$below, sum(geo < 1.1))
})

test_that("generalized SPD pencils: dense and sparse", {
  set.seed(5)
  n <- 120L
  A <- crossprod(matrix(stats::rnorm(n * n), n)) / n - 1
  A <- (A + t(A)) / 2
  B <- crossprod(matrix(stats::rnorm(n * n), n)) / n + diag(n)
  L <- chol(B)
  C <- backsolve(L, t(backsolve(L, A, transpose = TRUE)), transpose = TRUE)
  gev <- eigen((C + t(C)) / 2, symmetric = TRUE, only.values = TRUE)$values
  for (s in c(-0.8, -0.2, 0.1, 0.6)) {
    r <- eigen_count(A, s, B = B)
    expect_true(r$generalized)
    t5_expect_counts(r, gev, info = paste("dense sigma", s))
  }
  # Sparse A, sparse SPD (non-diagonal) B.
  d <- seq(-2, 2, length.out = 90L)
  As <- t5_rotated_sparse(d, 51L)
  Bs <- Matrix::forceSymmetric(methods::as(
    t5_laplacian_1d(90L) + 2 * Matrix::Diagonal(90L), "CsparseMatrix"))
  Ld <- chol(as.matrix(Bs))
  Cd <- backsolve(Ld, t(backsolve(Ld, as.matrix(As), transpose = TRUE)), transpose = TRUE)
  sev <- eigen((Cd + t(Cd)) / 2, symmetric = TRUE, only.values = TRUE)$values
  for (s in c(-0.4, -0.05, 0.02, 0.3)) {
    r <- eigen_count(As, s, B = Bs)
    expect_identical(r$method, "sparse_cholmod_ldl")
    t5_expect_counts(r, sev, info = paste("sparse sigma", s))
  }
  expect_error(eigen_count(A, 0, B = -B), "positive definite")
})

test_that("complex Hermitian counts use the real symmetric embedding", {
  set.seed(6)
  n <- 30L
  H <- matrix(complex(real = stats::rnorm(n * n), imaginary = stats::rnorm(n * n)), n)
  H <- H + Conj(t(H))
  ev <- eigen(H, only.values = TRUE)$values
  for (s in c(-3, 0, 1.5, 6)) {
    r <- eigen_count(H, s)
    expect_identical(r$method, "dense_bunch_kaufman_hermitian_embedding")
    t5_expect_counts(r, ev, info = paste("sigma", s))
  }
})

test_that("input validation", {
  expect_error(eigen_count(matrix(1:6, 2), 0), "square")
  expect_error(eigen_count(matrix(c(1, 2, 3, 4), 2), 0), "symmetric")
  expect_error(eigen_count(diag(3), NA), "sigma")
  expect_error(eigen_count(diag(3), c(1, 2)), "sigma")
  op <- linear_operator(c(3L, 3L), function(X, alpha = 1, beta = 0, Y = NULL) alpha * X,
                        structure = hermitian())
  expect_error(eigen_count(op, 0), "explicit matrix")
  # Operators with a matrix source are accepted.
  r <- eigen_count(as_operator(diag(c(1, 2, 3))), 2.5)
  expect_identical(r$below, 2L)
  expect_output(print(r), "below: 2")
})

test_that("unpivoted LDL' breakdown is perturbed, falls back, or is reported unreliable", {
  # Zero diagonal, eigenvalues +/- b: the unpivoted LDL' of S - 0 I meets an
  # exactly zero pivot (CHOLMOD aborts); nearby shifts recover the count.
  S <- t5_zero_diagonal_pairs(seq(1, 2, length.out = 40L))
  r <- eigen_count(S, 0, perturb = FALSE)
  expect_false(r$reliable)
  expect_true(is.na(r$below))
  expect_true(any(grepl("factorisation failed", r$notes)))
  r <- eigen_count(S, 0)
  expect_true(r$reliable)
  expect_true(r$perturbed || grepl("_fallback$", r$method))
  expect_identical(c(r$below, r$zero, r$above), c(40L, 0L, 40L))
  expect_false(r$diagnostics$attempts$reliable[[1L]])

  # A singular random sparse matrix (two-digit entries, empty rows) at
  # sigma = 0: no sparse LDL' attempt is trustworthy, so the count falls
  # back to dense Bunch-Kaufman (or, with the fallback disabled, is
  # reported unreliable rather than guessed).
  set.seed(2)
  R <- Matrix::rsparsematrix(600L, 600L, density = 0.003, symmetric = TRUE)
  ev <- eigen(as.matrix(R), symmetric = TRUE, only.values = TRUE)$values
  r <- eigen_count(R, 0)
  expect_true(r$reliable)
  expect_identical(r$below, sum(ev < r$zero_band[[1L]]))
  expect_identical(r$above, sum(ev > r$zero_band[[2L]]))
  if (grepl("_fallback$", r$method)) {
    expect_true(any(grepl("Bunch-Kaufman", r$notes)))
  }
  old <- options(eigencore.inertia_dense_fallback_limit = 0L)
  on.exit(options(old), add = TRUE)
  r2 <- eigen_count(R, 0)
  if (!r2$reliable) {
    expect_true(any(grepl("not certified", r2$notes)))
    expect_false(any(r2$diagnostics$attempts$reliable))
  } else {
    expect_identical(r2$below, sum(ev < r2$zero_band[[1L]]))
  }
})

test_that("shift-invert through CHOLMOD LDL' matches eigen()", {
  set.seed(7)
  d <- sort(stats::runif(300L, -5, 5))
  S <- t5_rotated_sparse(d, 71L)
  sigma <- 0.4321
  fit <- eig_partial(S, 6L, target = nearest(sigma),
                     method = shift_invert(sigma = sigma))
  expect_identical(fit$method,
                   "native thick-restart Hermitian Lanczos shift-invert (sparse LDL' solve callback)")
  expected <- d[order(abs(d - sigma))][1:6]
  expect_equal(sort(values(fit)), sort(expected), tolerance = 1e-8)
  expect_true(certificate(fit)$passed)
  cache <- fit$transform$factorization_cache
  expect_identical(cache$label_kind, "sparse_ldl")
  expect_true(cache$inertia_reliable)
  expect_identical(unname(cache$inertia[["below"]]), as.numeric(sum(d < sigma)))
  expect_gt(cache$factor_nnz, 0)
  # RSpectra-style front end.
  fit2 <- eigs_sym(S, 6L, sigma = sigma)
  expect_equal(sort(fit2$values), sort(expected), tolerance = 1e-8)
})

test_that("generalized shift-invert with sparse SPD B uses LDL'", {
  d <- seq(-2, 2, length.out = 90L)
  As <- t5_rotated_sparse(d, 52L)
  Bs <- methods::as(t5_laplacian_1d(90L) + 2 * Matrix::Diagonal(90L), "CsparseMatrix")
  Ld <- chol(as.matrix(Bs))
  Cd <- backsolve(Ld, t(backsolve(Ld, as.matrix(As), transpose = TRUE)), transpose = TRUE)
  gev <- eigen((Cd + t(Cd)) / 2, symmetric = TRUE, only.values = TRUE)$values
  sigma <- 0.0123
  fit <- eig_partial(As, B = Bs, k = 4L, target = nearest(sigma),
                     method = shift_invert(sigma = sigma),
                     allow_dense_fallback = "never")
  expect_identical(fit$method,
                   "native thick-restart generalized SPD Lanczos shift-invert (sparse LDL' solve callback)")
  expect_equal(sort(values(fit)), sort(gev[order(abs(gev - sigma))][1:4]),
               tolerance = 1e-7)
  expect_true(certificate(fit)$passed)
  expect_identical(fit$transform$factorization_cache$label_kind,
                   "sparse_ldl_generalized")
})

test_that("shift-invert falls back to Matrix::lu when the LDL' factor is unreliable", {
  # Zero diagonal: the unpivoted LDL' of A - 0 I hits an exact zero pivot.
  b <- seq(1, 3, length.out = 30L)
  S <- t5_zero_diagonal_pairs(b)
  fit <- eig_partial(S, 4L, target = nearest(0), method = shift_invert(sigma = 0))
  expect_identical(fit$method,
                   "native thick-restart Hermitian Lanczos shift-invert (sparse LU solve callback)")
  expect_equal(sort(values(fit)), sort(c(-1, 1, -b[2], b[2])), tolerance = 1e-8)
  expect_true(certificate(fit)$passed)
  cache <- fit$transform$factorization_cache
  expect_identical(cache$label_kind, "sparse_lu")
  expect_true(cache$ldl_attempted)
  expect_match(cache$ldl_fallback_reason, "CHOLMOD LDL' failed")
  expect_identical(fit$fallback_reason$code, "factorization_unreliable")
  # A probe tolerance no factor can meet forces the fallback too.
  d <- seq(-3, 3, length.out = 80L)
  S2 <- t5_rotated_sparse(d, 81L)
  old <- options(eigencore.shift_invert_ldl_tol = 0)
  on.exit(options(old), add = TRUE)
  fit <- eig_partial(S2, 3L, target = nearest(0.05), method = shift_invert(sigma = 0.05))
  expect_identical(fit$transform$factorization_cache$label_kind, "sparse_lu")
  expect_match(fit$transform$factorization_cache$ldl_fallback_reason, "backward error")
  expect_equal(sort(values(fit)), sort(d[order(abs(d - 0.05))][1:3]), tolerance = 1e-8)
})

test_that("9,9,7,7,7: inertia verifies (after repair) on sparse and dense routes", {
  n <- 60L
  d <- c(9, 9, 7, 7, 7, seq(5, 0.1, length.out = n - 5L))
  for (s in 1:4) {
    # Sparse non-diagonal: CHOLMOD route.
    S <- t5_rotated_sparse(d, 900L + s)
    fit <- eig_partial(S, 5L, largest(), method = lanczos(max_subspace = 12L),
                       seed = s)
    cert <- certificate(fit)
    expect_identical(cert$target_completeness, "inertia_verified", info = paste("sparse", s))
    expect_true(cert$passed)
    expect_true(cert$target_passed)
    expect_identical(cert$completeness$inertia_method, "sparse_cholmod_ldl")
    expect_identical(cert$completeness$count_upper, 5)
    expect_equal(sort(values(fit)), c(7, 7, 7, 9, 9), tolerance = 1e-7)
    # Smallest of -S.
    fit <- eig_partial(-S, 5L, smallest(), method = lanczos(max_subspace = 12L),
                       seed = s)
    expect_identical(certificate(fit)$target_completeness, "inertia_verified")
    expect_equal(sort(values(fit)), -c(9, 9, 7, 7, 7), tolerance = 1e-7)
    # Dense route.
    A <- t5_dense_with_spectrum(d, 950L + s)
    fit <- eig_partial(A, 5L, largest(), method = lanczos(max_subspace = 12L),
                       seed = s)
    cert <- certificate(fit)
    expect_identical(cert$target_completeness, "inertia_verified", info = paste("dense", s))
    expect_identical(cert$completeness$inertia_method, "dense_bunch_kaufman")
    expect_equal(sort(values(fit)), c(7, 7, 7, 9, 9), tolerance = 1e-7)
  }
  # Diagonal sparse input (the probe test's construction) is repaired.
  set.seed(1001L)
  P <- sample(n)
  D <- Matrix::sparseMatrix(i = P, j = P, x = d, dims = c(n, n))
  fit <- eig_partial(D, 5L, largest(), method = lanczos(), seed = 1L)
  cert <- certificate(fit)
  expect_identical(cert$target_completeness, "inertia_verified")
  expect_true(cert$completeness$repaired)
  expect_true(any(grepl("repaired", fit$warnings)))
  expect_equal(sort(values(fit)), c(7, 7, 7, 9, 9), tolerance = 1e-7)
})

test_that("an unrepaired missing copy is inertia_failed and withholds the certificate", {
  n <- 60L
  d <- c(9, 9, 7, 7, 7, seq(5, 0.1, length.out = n - 5L))
  set.seed(1001L)
  P <- sample(n)
  D <- Matrix::sparseMatrix(i = P, j = P, x = d, dims = c(n, n))
  old <- options(eigencore.completeness_max_rounds = 0L)
  on.exit(options(old), add = TRUE)
  fit <- eig_partial(D, 5L, largest(), method = lanczos(), seed = 1L)
  cert <- certificate(fit)
  expect_identical(cert$target_completeness, "inertia_failed")
  expect_false(cert$passed)
  expect_false(cert$target_passed)
  expect_true(cert$residual_passed)
  expect_true(any(grepl("inertia count", fit$warnings)))
})

test_that("a repeated eigenvalue straddling the target edge is verified as a tie", {
  # k = 3 of 9, 9, 7, 7, 7: every choice of one 7 is a correct answer. The
  # counts prove both 9s are returned and the third value lies in the edge
  # window (R/completeness_hermitian.R, hermitian_completeness_tie()).
  d <- c(9, 9, 7, 7, 7, seq(5, 0.1, length.out = 35L))
  A <- t5_dense_with_spectrum(d, 13L)
  fit <- eig_partial(A, 3L, largest(), method = lanczos(), seed = 1L)
  cert <- certificate(fit)
  expect_identical(cert$target_completeness, "inertia_verified")
  expect_true(cert$target_passed)
  expect_true(cert$passed)
  expect_true(cert$completeness$tie)
  expect_identical(as.integer(cert$completeness$tie_eigenvalues), 3L)
  expect_match(cert$completeness$reason, "tie")
  expect_equal(sort(values(fit)), c(7, 9, 9), tolerance = 1e-8)
})

test_that("an edge cluster the counts cannot resolve stays inertia_inconclusive", {
  # Returned 4.2 is within its residual bound of 5, but not far enough inside
  # the lower threshold to prove that 5 itself is returned: no claim.
  A <- diag(c(5, 2, 1, 0.5, 0.2))
  gate <- list(ctx = eigencore:::inertia_context(A), delta_A = 0, delta_B = 0,
               metric = NULL)
  problem <- list(A = as_operator(A))
  vals <- c(4.2, 2)
  parts <- eigencore:::hermitian_completeness_parts(largest(), vals)
  res <- eigencore:::hermitian_completeness_count(problem, gate, vals, c(0.8, 0.8),
                                                  0, parts)
  expect_identical(res$status, "inertia_inconclusive")
  expect_match(res$record$reason, "not separated")
  expect_null(res$record$tie)
  cert <- eigencore:::require_verified_completeness(list(certificate =
    eigencore:::certificate_with_completeness(list(passed = TRUE), res$status, res$record)))$certificate
  expect_false(cert$passed)
  expect_true(cert$residual_passed)
  expect_true(any(grepl("inconclusive", cert$notes)))
})

test_that("generalized SPD Lanczos is inertia-verified against A - t B", {
  n <- 60L
  bdiag <- seq(1, 3, length.out = n)
  d <- c(9, 9, 7, 7, 7, seq(5, 0.1, length.out = n - 5L))
  for (s in 1:3) {
    set.seed(7000L + s)
    P <- sample(n)
    C <- Matrix::sparseMatrix(i = P, j = P, x = d, dims = c(n, n))
    Bh <- Matrix::Diagonal(n, sqrt(bdiag))
    A <- as.matrix(Bh %*% C %*% Bh)
    A <- (A + t(A)) / 2
    fit <- eig_partial(A, 5L, largest(), B = Matrix::Diagonal(n, bdiag),
                       method = lanczos(max_subspace = 12L), seed = s)
    cert <- certificate(fit)
    expect_identical(cert$target_completeness, "inertia_verified")
    expect_true(cert$passed)
    expect_gt(cert$completeness$metric_lower_bound, 0)
    expect_equal(sort(values(fit)), c(7, 7, 7, 9, 9), tolerance = 1e-7)
  }
})

test_that("nearest and largest_magnitude targets are verified by interval counts", {
  set.seed(8)
  d <- sort(stats::runif(200L, -4, 4))
  S <- t5_rotated_sparse(d, 82L)
  # auto mode counts interior (nearest) targets when the cost gate passes ...
  fit <- eig_partial(S, 5L, target = nearest(0.7), method = shift_invert(sigma = 0.7))
  cert <- certificate(fit)
  expect_identical(cert$target_completeness, "inertia_verified")
  expect_true(cert$passed)
  expect_identical(cert$completeness$gate, "cost gate passed")
  # ... and always on request.
  fit <- with_t5_completeness(
    "inertia",
    eig_partial(S, 5L, target = nearest(0.7), method = shift_invert(sigma = 0.7))
  )
  cert <- certificate(fit)
  expect_identical(cert$target_completeness, "inertia_verified")
  expect_identical(cert$completeness$kind, "nearest")
  expect_equal(sort(values(fit)), sort(d[order(abs(d - 0.7))][1:5]), tolerance = 1e-8)
  fit <- eig_partial(S, 4L, target = largest_magnitude(), method = lanczos())
  cert <- certificate(fit)
  expect_identical(cert$target_completeness, "inertia_verified")
  expect_equal(sort(abs(values(fit))), sort(abs(d))[197:200], tolerance = 1e-7)
})

test_that("completeness modes and the cost gate", {
  d <- c(9, 8, 7, 6, seq(5, 0.1, length.out = 56L))
  S <- t5_rotated_sparse(d, 99L)
  fit <- eig_partial(S, 4L, largest(), method = lanczos(), seed = 1L)
  expect_identical(certificate(fit)$target_completeness, "inertia_verified")
  expect_identical(certificate(fit)$completeness$gate, "cost gate passed")
  fit <- eig_partial(S, 4L, largest(), method = lanczos(completeness = "probe"), seed = 1L)
  expect_identical(certificate(fit)$target_completeness, "probed")
  fit <- eig_partial(S, 4L, largest(), method = lanczos(completeness = "none"), seed = 1L)
  expect_identical(certificate(fit)$target_completeness, "not_checked")
  fit <- with_t5_completeness("inertia",
                              eig_partial(S, 4L, largest(), method = lanczos(), seed = 1L))
  expect_identical(certificate(fit)$target_completeness, "inertia_verified")
  expect_error(lanczos(completeness = "sometimes"), "auto")
  # A zero budget sends auto mode to the probe, and says why.
  old <- options(eigencore.completeness_inertia_seconds = 0,
                 eigencore.completeness_inertia_ratio = 0)
  on.exit(options(old), add = TRUE)
  fit <- eig_partial(S, 4L, largest(), method = lanczos(), seed = 1L)
  cert <- certificate(fit)
  expect_identical(cert$target_completeness, "probed")
  expect_match(cert$completeness$inertia_gate, "cost")
  # "inertia" mode ignores the gate.
  fit <- with_t5_completeness("inertia",
                              eig_partial(S, 4L, largest(), method = lanczos(), seed = 1L))
  expect_identical(certificate(fit)$target_completeness, "inertia_verified")
  # Small matrix-free operators are materialised (n applies) and counted ...
  op <- linear_operator(dim(S), function(X, alpha = 1, beta = 0, Y = NULL) {
    Z <- alpha * as.matrix(S %*% X)
    if (is.null(Y) || beta == 0) Z else Z + beta * Y
  }, structure = hermitian())
  fit <- with_t5_completeness("inertia",
                              eig_partial(op, 4L, largest(), method = lanczos(block = 2L),
                                          seed = 1L))
  expect_identical(certificate(fit)$target_completeness, "inertia_verified")
  expect_true(certificate(fit)$completeness$materialized)
  # ... larger ones keep the probe.
  old_limit <- options(eigencore.completeness_materialize_limit = 10L)
  on.exit(options(old_limit), add = TRUE)
  fit <- with_t5_completeness("inertia",
                              eig_partial(op, 4L, largest(), method = lanczos(block = 2L),
                                          seed = 1L))
  expect_identical(certificate(fit)$target_completeness, "probed")
  expect_match(certificate(fit)$completeness$inertia_gate, "matrix source")
})

test_that("internal interval counts reuse one context", {
  m <- 30L
  L2 <- t5_laplacian_2d(m)
  lam <- as.vector(outer(2 - 2 * cos(pi * seq_len(m) / (m + 1)),
                         2 - 2 * cos(pi * seq_len(m) / (m + 1)), "+"))
  ctx <- eigencore:::inertia_context(L2)
  expect_identical(ctx$kind, "sparse")
  iv <- eigencore:::inertia_interval_count(ctx, 1.01, 2.03)
  expect_true(iv$reliable)
  expect_equal(iv$count, sum(lam > iv$interval[[1L]] & lam < iv$interval[[2L]]))
  cost <- eigencore:::inertia_factor_cost(ctx)
  if (eigencore:::cholmod_bridge_available()) {
    expect_identical(cost$source, "cholmod_analyze")
    expect_gt(cost$flops, 0)
  }
})

test_that("a count-proven gap is repaired even when the probe sees nothing", {
  # Exact eigenpairs 9, 9, 7, 7, 5 of a spectrum 9, 9, 7, 7, 7, ...: every
  # residual is tiny, the count proves a 7 is missing, and a one-step,
  # one-column probe is too weak to show it, so the repair starts the
  # deflated complement solve from a deterministic block.
  n <- 60L
  d <- c(9, 9, 7, 7, 7, seq(5, 0.1, length.out = n - 5L))
  Q <- t5_orthogonal(n, 3L)
  A <- Q %*% (d * t(Q))
  A <- (A + t(A)) / 2
  V <- Q[, c(1, 2, 3, 4, 6)]
  vals <- d[c(1, 2, 3, 4, 6)]
  op <- as_operator(A)
  cert <- eigencore:::certify_eigen_operator(op, vals, V, tol = 1e-8)
  expect_true(cert$passed)
  problem <- list(A = op, metric = NULL, target = largest())
  ctx <- eigencore:::inertia_context(op)
  bare <- eigencore:::inertia_completeness_check(op, vals, cert$residuals,
                                                 cert$orthogonality, largest(),
                                                 ctx = ctx)
  expect_identical(bare$status, "inertia_failed")
  old <- options(eigencore.completeness_probe_steps = 1L,
                 eigencore.completeness_probe_block = 1L)
  on.exit(options(old), add = TRUE)
  chk <- eigencore:::inertia_completeness_run(problem, vals, V, cert, tol = 1e-8,
                                              ctx = ctx)
  expect_identical(chk$status, "inertia_verified")
  expect_true(chk$repaired)
  expect_equal(sort(chk$values), c(7, 7, 7, 9, 9), tolerance = 1e-8)
})
