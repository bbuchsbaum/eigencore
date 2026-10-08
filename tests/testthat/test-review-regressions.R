test_that("double centering removes the grand mean once", {
  set.seed(301)
  X <- matrix(rnorm(40 * 6, mean = 3), 40, 6)
  truth <- X - outer(rowMeans(X), rep(1, 6)) - outer(rep(1, 40), colMeans(X)) +
    mean(X)
  expect_equal(rowMeans(truth), rep(0, 40))
  expect_equal(colMeans(truth), rep(0, 6))

  dense <- center(X, rows = TRUE, columns = TRUE)
  expect_equal(svd_partial(dense, rank = 2)$d, svd(truth)$d[1:2],
               tolerance = 1e-8)

  sparse <- center(methods::as(Matrix::Matrix(X, sparse = TRUE), "CsparseMatrix"),
                   rows = TRUE, columns = TRUE)
  Z <- matrix(rnorm(12), 6, 2)
  expect_equal(sparse$apply(Z), truth %*% Z)
  W <- matrix(rnorm(80), 40, 2)
  expect_equal(sparse$apply_adjoint(W), t(truth) %*% W)
})

test_that("column-centered sparse operators certify with an exact norm", {
  set.seed(302)
  A <- Matrix::rsparsematrix(300, 20, density = 0.2)
  op <- center(A, columns = TRUE)
  dense <- sweep(as.matrix(A), 2L, Matrix::colMeans(A), `-`)
  expect_equal(op$metadata$frobenius_norm, sqrt(sum(dense^2)))
  fit <- svd_partial(op, rank = 3)
  expect_false(fit$certificate$scale_is_estimate)
  expect_equal(fit$d, svd(dense)$d[1:3], tolerance = 1e-8)
})

test_that("complex Frobenius norms keep imaginary parts", {
  Z <- matrix(complex(real = 0, imaginary = 1:4), 2, 2)
  expect_equal(eigencore:::matrix_norm(Z), sqrt(sum(Mod(Z)^2)))
  S <- Matrix::forceSymmetric(Matrix::Matrix(matrix(c(2, 1, 1, 3), 2), sparse = TRUE))
  expect_equal(eigencore:::matrix_norm(S), sqrt(4 + 1 + 1 + 9))
})

test_that("eigs_sym honours sigma, lower, ordering, and opts", {
  set.seed(303)
  n <- 80
  M <- crossprod(matrix(rnorm(n * n), n))
  ref <- eigen(M, symmetric = TRUE, only.values = TRUE)$values
  sigma <- ref[40] + 1e-3
  near <- eigs_sym(M, 3, sigma = sigma)
  expect_equal(sort(near$values), sort(ref[order(abs(ref - sigma))[1:3]]),
               tolerance = 1e-8)
  expect_false(is.unsorted(rev(near$values)))

  L <- M
  L[upper.tri(L)] <- 0
  expect_equal(eigs_sym(L, 3)$values, ref[1:3], tolerance = 1e-8)
  U <- M
  U[lower.tri(U)] <- 0
  expect_equal(eigs_sym(U, 3, lower = FALSE)$values, ref[1:3], tolerance = 1e-8)

  sa <- eigs_sym(M, 3, which = "SA")
  expect_equal(sa$values, rev(ref)[3:1], tolerance = 1e-8)

  no_vec <- eigs_sym(M, 2, opts = list(retvec = FALSE, tol = 1e-10))
  expect_null(no_vec$vectors)
  expect_warning(eigs_sym(M, 2, opts = list(maxitr = 10)), "not used")
  expect_warning(eigs_sym(M, 2, opts = list(bogus = 1)), "Unknown opts")
  expect_error(eigs_sym(M, 2, which = "XX"), "Unknown ARPACK selector")
  expect_error(eigs_sym(M, 2, sigma = 1, which = "LA"), "only which = 'LM'")
})

test_that("eigs_sym defaults to largest magnitude like RSpectra", {
  A <- diag(c(-10, 1, 2, 3, 5))
  expect_equal(eigs_sym(A, 2)$values, c(5, -10))
})

test_that("RSpectra function interfaces are accepted", {
  set.seed(304)
  n <- 50
  M <- crossprod(matrix(rnorm(n * n), n))
  f <- function(x, args) args %*% x
  fit <- eigs_sym(f, 2, n = n, args = M)
  expect_equal(fit$values, eigen(M, TRUE, only.values = TRUE)$values[1:2],
               tolerance = 1e-6)

  A <- matrix(rnorm(60 * 8), 60, 8)
  sv <- svds(function(x, args) args %*% x, 2,
             Atrans = function(x, args) crossprod(args, x),
             dim = dim(A), args = A)
  expect_equal(sv$d, svd(A)$d[1:2], tolerance = 1e-6)
})

test_that("svds honours nu, nv, center, and scale", {
  set.seed(305)
  A <- matrix(rnorm(50 * 6, mean = 2), 50, 6)
  sv <- svds(A, 4, nu = 1, nv = 2)
  expect_equal(ncol(sv$u), 1L)
  expect_equal(ncol(sv$v), 2L)

  pc <- svds(A, 2, opts = list(center = TRUE, scale = TRUE))
  expect_equal(pc$d, svd(scale(A))$d[1:2], tolerance = 1e-8)
  sparse <- methods::as(Matrix::Matrix(A, sparse = TRUE), "CsparseMatrix")
  pc_sparse <- svds(sparse, 2, opts = list(center = TRUE))
  expect_equal(pc_sparse$d, svd(scale(A, scale = FALSE))$d[1:2],
               tolerance = 1e-8)
})

test_that("seeded solves restore the global random stream", {
  set.seed(10)
  a <- runif(1)
  set.seed(10)
  invisible(eig_partial(diag(5:1), k = 2, seed = 1))
  expect_equal(runif(1), a)
})

test_that("complex Hermitian matrices with repeated eigenvalues keep orthonormal vectors", {
  set.seed(306)
  Q <- qr.Q(qr(matrix(complex(real = rnorm(36), imaginary = rnorm(36)), 6)))
  H <- Q %*% diag(c(1, 1, 1, 2, 3, 4)) %*% Conj(t(Q))
  H <- (H + Conj(t(H))) / 2
  fit <- eig_full(H)
  expect_lt(max(Mod(Conj(t(fit$vectors)) %*% fit$vectors - diag(6))), 1e-10)
  expect_true(fit$certificate$passed)

  B <- diag(c(1, 2, 3, 4, 5, 6)) + 0i
  pencil <- eig_full(H, B)
  expect_lt(max(Mod(Conj(t(pencil$vectors)) %*% B %*% pencil$vectors - diag(6))),
            1e-10)
})

test_that("eig_full falls back to QZ for symmetric indefinite B", {
  set.seed(307)
  A <- crossprod(matrix(rnorm(25), 5))
  B <- diag(c(1, -1, 2, 1, 3))
  fit <- eig_full(A, B)
  ref <- eigen(solve(B, A), only.values = TRUE)$values
  expect_equal(sort(Re(fit$values)), sort(Re(ref)), tolerance = 1e-8)
  expect_error(eig_full(A, B, structure = hermitian()), "potrf")
})

test_that("non-finite inputs are rejected up front", {
  M <- diag(3)
  M[1, 2] <- NA
  expect_error(eig_partial(M, k = 1), "NA, NaN, or Inf")
  expect_error(svd_partial(Matrix::Matrix(M, sparse = TRUE), rank = 1),
               "NA, NaN, or Inf")
  expect_error(eig_full(M), "NA, NaN, or Inf")
})

test_that("k and rank are validated centrally", {
  A <- diag(5:1)
  expect_error(eig_partial(A, k = 0), "between 1 and 5")
  expect_error(eig_partial(A, k = 6), "between 1 and 5")
  expect_error(eig_partial(A, k = 2.5), "whole number")
  expect_error(svd_partial(matrix(1, 4, 3), rank = 4), "between 1 and 3")
  expect_error(eig_partial(A, k = 2, target = both_ends(2, 2)),
               "k = k_low \\+ k_high")
})

test_that("SVD is scale invariant on the Gram path", {
  set.seed(308)
  m <- 400
  n <- 30
  U <- qr.Q(qr(matrix(rnorm(m * n), m)))
  V <- qr.Q(qr(matrix(rnorm(n * n), n)))
  d <- 10^seq(0, -3, length.out = n)
  A <- U %*% (d * t(V))
  for (s in c(1e-10, 1e10)) {
    fit <- svd_partial(A * s, rank = 3)
    expect_equal(fit$d / s, d[1:3], tolerance = 1e-8)
    sparse <- methods::as(Matrix::Matrix(A * s, sparse = TRUE), "CsparseMatrix")
    expect_equal(svd_partial(sparse, rank = 3)$d / s, d[1:3], tolerance = 1e-8)
  }
})

test_that("tridiagonal preconditioner accepts symmetric storage", {
  T3 <- matrix(c(2, -1, 0, -1, 2, -1, 0, -1, 2), 3)
  p_base <- shifted_tridiagonal_preconditioner(T3, shift = 0.5)
  p_sym <- shifted_tridiagonal_preconditioner(
    Matrix::forceSymmetric(Matrix::Matrix(T3, sparse = TRUE)), shift = 0.5
  )
  R <- matrix(c(1, 2, 3), 3, 1)
  expected <- solve(T3 + 0.5 * diag(3), R)
  expect_equal(unclass(p_base(R)), expected, ignore_attr = TRUE)
  expect_equal(unclass(p_sym(R)), expected, ignore_attr = TRUE)
})

test_that("rescaled nonsymmetric matrices are not classified as symmetric", {
  set.seed(309)
  A <- matrix(rnorm(400), 20)
  ref <- eigen(A, only.values = TRUE)$values
  for (s in c(1, 1e-9, 1e-12)) {
    B <- A * s
    S <- methods::as(Matrix::Matrix(B, sparse = TRUE), "generalMatrix")
    expect_identical(as_operator(B)$structure$kind, "general")
    expect_identical(as_operator(S)$structure$kind, "general")
    fit <- eig_partial(B, k = 2, target = largest_magnitude())
    expect_equal(sort(Mod(fit$values)) / s, sort(Mod(ref))[19:20],
                 tolerance = 1e-8)
  }
  expect_identical(as_operator(diag(3) * 1e-300)$structure$kind, "hermitian")
})
