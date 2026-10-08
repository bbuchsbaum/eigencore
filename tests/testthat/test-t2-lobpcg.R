# Tranche-2 LOBPCG regressions: Knyazev search directions with soft locking
# (C19), two-pass constraint projection re-applied after normalisation (C20),
# pivoted tridiagonal preconditioner (C21), block B-orthonormalisation and
# A X reuse (P3), norm-relative stopping (C16/P4), cached SPD checks (P13).

t2_lobpcg_rotated <- function(values, seed) {
  set.seed(seed)
  n <- length(values)
  Q <- qr.Q(qr(matrix(stats::rnorm(n * n), n, n)))
  A <- Q %*% diag(values) %*% t(Q)
  (A + t(A)) / 2
}

t2_lobpcg_spd <- function(n, seed, cond = 10) {
  set.seed(seed)
  Q <- qr.Q(qr(matrix(stats::rnorm(n * n), n, n)))
  B <- Q %*% diag(exp(seq(0, log(cond), length.out = n))) %*% t(Q)
  (B + t(B)) / 2
}

# Generalized symmetric-definite reference via Cholesky reduction
# (what geigen::geigen(A, B, symmetric = TRUE) computes).
t2_lobpcg_geigen <- function(A, B) {
  R <- chol(B)
  Ri <- backsolve(R, diag(nrow(B)))
  C <- t(Ri) %*% A %*% Ri
  ev <- eigen((C + t(C)) / 2, symmetric = TRUE)
  list(values = ev$values, vectors = Ri %*% ev$vectors)
}

# Largest principal-angle sine between span(U) and span(V) in the B metric.
t2_lobpcg_subspace_sine <- function(U, V, B = NULL) {
  if (is.null(B)) B <- diag(nrow(U))
  R <- chol(B)
  Uq <- qr.Q(qr(R %*% U))
  Vq <- qr.Q(qr(R %*% V))
  sqrt(max(0, 1 - min(svd(crossprod(Uq, Vq))$d)^2))
}

test_that("native LOBPCG matches eigen() for both ends of a dense problem", {
  values <- c(seq(1, 3, length.out = 6), seq(10, 60, length.out = 54))
  A <- t2_lobpcg_rotated(values, 701)
  ref <- eigen(A, symmetric = TRUE)
  op <- as_operator(A, structure = hermitian())

  small <- eigencore:::native_lobpcg_hermitian(op, 4, smallest(), tol = 1e-10,
                                               maxit = 400L, seed = 11)
  expect_true(small$certificate$passed)
  expect_equal(small$values, rev(ref$values)[1:4], tolerance = 1e-9)
  expect_lt(t2_lobpcg_subspace_sine(small$vectors, ref$vectors[, 60:57]), 1e-6)
  expect_lt(max(abs(crossprod(small$vectors) - diag(4))), 1e-12)

  large <- eigencore:::native_lobpcg_hermitian(op, 3, largest(), tol = 1e-10,
                                               maxit = 400L, seed = 12)
  expect_true(large$certificate$passed)
  expect_equal(large$values, ref$values[1:3], tolerance = 1e-9)
})

test_that("native generalized LOBPCG matches a Cholesky-reduced dense reference", {
  n <- 50
  A <- t2_lobpcg_rotated(c(0.5, 0.9, 1.4, seq(5, 40, length.out = n - 3)), 702)
  B <- t2_lobpcg_spd(n, 703, cond = 50)
  ref <- t2_lobpcg_geigen(A, B)
  opA <- as_operator(A, structure = hermitian())
  opB <- as_operator(B, structure = hermitian())

  fit <- eigencore:::native_generalized_lobpcg_hermitian(
    opA, opB, 3, smallest(), tol = 1e-10, maxit = 500L, seed = 13
  )
  expect_true(fit$certificate$passed)
  expect_equal(fit$values, rev(ref$values)[1:3], tolerance = 1e-8)
  expect_lt(max(abs(crossprod(fit$vectors, B %*% fit$vectors) - diag(3))), 1e-10)
  expect_lt(t2_lobpcg_subspace_sine(fit$vectors, ref$vectors[, n:(n - 2)], B), 1e-6)

  # sparse A with diagonal B goes through the CSC/diagonal native kernel
  As <- Matrix::bandSparse(200, k = c(-1, 0, 1),
                           diagonals = list(rep(-1, 199), rep(2, 200), rep(-1, 199)))
  As <- methods::as(methods::as(As, "generalMatrix"), "CsparseMatrix")
  set.seed(704)
  bd <- stats::runif(200, 1, 3)
  Bs <- Matrix::Diagonal(200, x = bd)
  ref_s <- t2_lobpcg_geigen(as.matrix(As), diag(bd))
  fit_s <- eigencore:::native_generalized_lobpcg_hermitian(
    as_operator(As, structure = hermitian()),
    as_operator(Bs, structure = hermitian()),
    3, smallest(), tol = 1e-9, maxit = 600L, seed = 14,
    preconditioner = shifted_tridiagonal_preconditioner(As, shift = 1e-3)
  )
  expect_true(fit_s$certificate$passed)
  expect_equal(fit_s$values, rev(ref_s$values)[1:3], tolerance = 1e-8)
})

test_that("clustered and repeated eigenvalues converge to the invariant subspace", {
  # Exactly repeated eigenvalue (rotation-invariant Ritz vectors) and a tight
  # cluster below a gap: P must follow the Ritz vectors' signs and rotation.
  values <- c(1, 1, 1, 1 + 1e-7, 1 + 2e-7, seq(3, 8, length.out = 55))
  A <- t2_lobpcg_rotated(values, 705)
  ref <- eigen(A, symmetric = TRUE)
  op <- as_operator(A, structure = hermitian())
  fit <- eigencore:::native_lobpcg_hermitian(op, 5, smallest(), tol = 1e-10,
                                             maxit = 300L, seed = 15)
  expect_true(fit$certificate$passed)
  expect_equal(sort(fit$values), sort(values)[1:5], tolerance = 1e-9)
  expect_lt(t2_lobpcg_subspace_sine(fit$vectors, ref$vectors[, 60:56]), 1e-6)
  expect_lt(max(abs(crossprod(fit$vectors) - diag(5))), 1e-12)

  # generalized pencil with a repeated eigenvalue
  B <- t2_lobpcg_spd(60, 706, cond = 20)
  R <- chol(B)
  Ag <- t(R) %*% A %*% R
  Ag <- (Ag + t(Ag)) / 2
  fit_g <- eigencore:::native_generalized_lobpcg_hermitian(
    as_operator(Ag, structure = hermitian()),
    as_operator(B, structure = hermitian()),
    5, smallest(), tol = 1e-10, maxit = 400L, seed = 16
  )
  expect_true(fit_g$certificate$passed)
  expect_equal(sort(fit_g$values), sort(values)[1:5], tolerance = 1e-9)
  expect_lt(max(abs(crossprod(fit_g$vectors, B %*% fit_g$vectors) - diag(5))), 1e-10)
  # R x spans the invariant subspace of A for the cluster
  expect_lt(t2_lobpcg_subspace_sine(R %*% fit_g$vectors, ref$vectors[, 60:56]), 1e-6)
})

test_that("constraints stay B-orthogonal to 1e-12 after the solve", {
  n <- 40
  A <- t2_lobpcg_rotated(seq_len(n), 707)
  ref <- eigen(A, symmetric = TRUE)
  op <- as_operator(A, structure = hermitian())

  # exact (invariant) constraints: the deflated pairs are true eigenpairs
  Ce <- ref$vectors[, n:(n - 1)] %*% matrix(c(1, 2, -1, 1), 2, 2)
  fit <- eigencore:::native_lobpcg_hermitian(op, 3, smallest(), tol = 1e-10,
                                             maxit = 400L, seed = 17,
                                             constraints = Ce)
  expect_equal(fit$constraints_rank, 2L)
  expect_true(fit$certificate$passed)
  expect_equal(fit$values, c(3, 4, 5), tolerance = 1e-9)
  expect_lt(max(abs(crossprod(qr.Q(qr(Ce)), fit$vectors))), 1e-12)

  # inexact constraints: Ritz values of the constrained problem
  set.seed(7071)
  C <- ref$vectors[, n:(n - 1)] + 1e-3 * matrix(stats::rnorm(n * 2), n, 2)
  C <- sweep(C, 2L, sqrt(colSums(C^2)), "/")
  fit_i <- eigencore:::native_lobpcg_hermitian(op, 3, smallest(), tol = 1e-10,
                                               maxit = 200L, seed = 17,
                                               constraints = C)
  expect_lt(max(abs(crossprod(C, fit_i$vectors))), 1e-12)
  Zs <- qr.Q(qr(C), complete = TRUE)[, -(1:2)]
  ref_s <- eigen(t(Zs) %*% A %*% Zs, symmetric = TRUE)
  expect_equal(fit_i$values, rev(ref_s$values)[1:3], tolerance = 1e-10)

  # generalized: constraints are B-orthogonal to the computed vectors and the
  # values match the pencil restricted to the B-complement of C
  B <- t2_lobpcg_spd(n, 708, cond = 30)
  fit_g <- eigencore:::native_generalized_lobpcg_hermitian(
    op, as_operator(B, structure = hermitian()), 3, smallest(),
    tol = 1e-10, maxit = 200L, seed = 18, constraints = C
  )
  expect_equal(fit_g$constraints_rank, 2L)
  expect_lt(max(abs(crossprod(C, B %*% fit_g$vectors))), 1e-12)
  Z <- qr.Q(qr(B %*% C), complete = TRUE)[, -(1:2)]    # C' B Z = 0
  ref_c <- t2_lobpcg_geigen(t(Z) %*% A %*% Z, t(Z) %*% B %*% Z)
  expect_equal(fit_g$values, rev(ref_c$values)[1:3], tolerance = 1e-10)

  # generalized with exact constraints (B-eigenvectors of the pencil) certifies
  ref_g <- t2_lobpcg_geigen(A, B)
  Cg <- ref_g$vectors[, n:(n - 1)]
  fit_ge <- eigencore:::native_generalized_lobpcg_hermitian(
    op, as_operator(B, structure = hermitian()), 2, smallest(),
    tol = 1e-10, maxit = 300L, seed = 25, constraints = Cg
  )
  expect_true(fit_ge$certificate$passed)
  expect_equal(fit_ge$values, rev(ref_g$values)[3:4], tolerance = 1e-9)
  expect_lt(max(abs(crossprod(Cg, B %*% fit_ge$vectors))), 1e-12)

  # rank-deficient constraint block: duplicated column counts once
  fit_d <- eigencore:::native_lobpcg_hermitian(op, 2, smallest(), tol = 1e-9,
                                               maxit = 300L, seed = 19,
                                               constraints = cbind(C, C[, 1]))
  expect_equal(fit_d$constraints_rank, 2L)
  expect_lt(max(abs(crossprod(C, fit_d$vectors))), 1e-12)
})

test_that("stopping rule is invariant to scaling A", {
  values <- c(seq(0.1, 0.5, length.out = 4), seq(2, 10, length.out = 46))
  A <- t2_lobpcg_rotated(values, 709)
  ref <- eigen(A, symmetric = TRUE)
  fit1 <- eigencore:::native_lobpcg_hermitian(
    as_operator(A, structure = hermitian()), 3, smallest(),
    tol = 1e-9, maxit = 300L, seed = 20
  )
  fit2 <- eigencore:::native_lobpcg_hermitian(
    as_operator(A * 1e-10, structure = hermitian()), 3, smallest(),
    tol = 1e-9, maxit = 300L, seed = 20
  )
  expect_true(fit1$certificate$passed)
  expect_true(fit2$certificate$passed)
  expect_lte(abs(fit1$iterations - fit2$iterations), 1L)
  expect_gt(fit2$iterations, 2L)
  expect_equal(fit2$values / 1e-10, rev(ref$values)[1:3], tolerance = 1e-8)
  expect_equal(fit1$values, rev(ref$values)[1:3], tolerance = 1e-8)

  B <- t2_lobpcg_spd(50, 710, cond = 10)
  g1 <- eigencore:::native_generalized_lobpcg_hermitian(
    as_operator(A, structure = hermitian()), as_operator(B, structure = hermitian()),
    2, smallest(), tol = 1e-9, maxit = 300L, seed = 21
  )
  g2 <- eigencore:::native_generalized_lobpcg_hermitian(
    as_operator(A * 1e-10, structure = hermitian()),
    as_operator(B * 1e8, structure = hermitian()),
    2, smallest(), tol = 1e-9, maxit = 300L, seed = 21
  )
  expect_true(g1$certificate$passed)
  expect_true(g2$certificate$passed)
  expect_lte(abs(g1$iterations - g2$iterations), 1L)
  expect_equal(g2$values / 1e-18, g1$values, tolerance = 1e-7)
})

test_that("one operator block application per iteration", {
  values <- c(seq(1, 2, length.out = 4), seq(5, 50, length.out = 76))
  A <- t2_lobpcg_rotated(values, 711)
  fit <- eigencore:::native_lobpcg_hermitian(
    as_operator(A, structure = hermitian()), 3, smallest(),
    tol = 1e-10, maxit = 300L, seed = 22
  )
  expect_true(fit$certificate$passed)
  expect_gt(fit$iterations, 5L)
  # previously 2 * iterations - 1 block applications (A X and A [X W P]);
  # now A [W P] once per iteration plus the initial and periodic/final
  # explicit refreshes of A X.
  expect_lte(fit$matvecs, fit$iterations + ceiling(fit$iterations / 16) + 2L)
})

test_that("tridiagonal preconditioner factors indefinite shifted systems with pivoting", {
  n <- 40
  L <- Matrix::bandSparse(n, k = c(-1, 0, 1),
                          diagonals = list(rep(-1, n - 1), rep(2, n), rep(-1, n - 1)))
  L <- methods::as(methods::as(L, "generalMatrix"), "CsparseMatrix")
  # L - I has a zero pivot at the second step of unpivoted elimination.
  T0 <- L - Matrix::Diagonal(n)
  pre <- shifted_tridiagonal_preconditioner(T0, shift = 0)
  fit <- eigencore:::native_lobpcg_hermitian(
    as_operator(L, structure = hermitian()), 2, smallest(),
    tol = 1e-8, maxit = 300L, seed = 23, preconditioner = pre
  )
  # The unpivoted Thomas sweep used to fail here (status -5); the pivoted
  # factorisation applies the (indefinite, so not convergence-accelerating)
  # preconditioner without breaking down.
  expect_gt(fit$preconditioner_calls, 0L)
  expect_true(all(is.finite(fit$values)))
  expect_true(all(is.finite(fit$residuals)))
  ref <- sort(eigen(as.matrix(L), symmetric = TRUE, only.values = TRUE)$values)
  expect_true(all(fit$values >= ref[1] - 1e-12))

  fit_t <- eigencore:::native_lobpcg_tridiagonal_hermitian(
    as_operator(L, structure = hermitian()), 2, smallest(),
    tol = 1e-9, maxit = 80L, shift = 1e-3, seed = 24
  )
  expect_true(fit_t$certificate$passed)
  expect_equal(fit_t$preconditioner$factorization, "tridiagonal_lu_dgttrf")
})

test_that("dense SPD metric checks share one cached Cholesky factorisation", {
  B <- t2_lobpcg_spd(30, 712, cond = 5)
  Bop <- as_operator(B, structure = hermitian())
  computed <- 0L
  count <- function(source) {
    computed <<- computed + 1L
    TRUE
  }
  key <- paste0("t2_test_", sample.int(1e6, 1L))
  for (i in 1:3) {
    expect_true(eigencore:::spd_metric_cache_get(B, count, kind = key))
  }
  expect_identical(computed, 1L)
  rm(list = key, envir = eigencore:::.spd_metric_cache)

  expect_true(eigencore:::generalized_spd_metric_known(Bop))
  F1 <- eigencore:::generalized_spd_metric_dense_factor(B)
  F2 <- eigencore:::generalized_spd_metric_dense_factor(eigencore:::source_or_null(Bop))
  expect_identical(F1, F2)
  expect_equal(crossprod(F1), (B + t(B)) / 2, tolerance = 1e-12)
  expect_false(eigencore:::generalized_spd_metric_known(
    as_operator(-B, structure = hermitian())
  ))
  expect_null(eigencore:::generalized_spd_metric_dense_factor(-B))
})
