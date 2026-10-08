# Tranche 2 dense drivers (review items P5, C21): dsyevr / dgesdd / dbdsdc /
# dstevr / dgttrf-dgttrs replace dsyev / dgesvd / dense-gesvd / dstev /
# unpivoted Thomas.

t2_subspace_gap <- function(U, V) {
  # Largest sine of the principal angles between span(U) and span(V)
  # (both orthonormal); invariant to column signs and rotations within the
  # subspace.
  P <- V - U %*% crossprod(U, V)
  if (!length(P)) {
    return(0)
  }
  max(svd(P, nu = 0, nv = 0)$d)
}

t2_rel <- function(x, y) {
  max(abs(x - y)) / max(1, max(abs(y)))
}

test_that("eig_full symmetric path matches base eigen with dsyevr", {
  set.seed(2101)
  n <- 60
  X <- matrix(rnorm(n * n), n)
  A <- crossprod(X) - 20 * diag(n)
  fit <- eig_full(A)
  ref <- eigen(A, symmetric = TRUE)

  expect_lt(t2_rel(sort(values(fit)), sort(ref$values)), 1e-10)
  expect_true(certificate(fit)$passed)
  V <- vectors(fit)
  expect_lt(max(abs(crossprod(V) - diag(n))), 1e-10)
  ord_fit <- order(values(fit))
  ord_ref <- order(ref$values)
  # Simple spectrum here: each eigenvector agrees up to sign.
  for (j in c(1L, 17L, n)) {
    expect_lt(t2_subspace_gap(V[, ord_fit[j], drop = FALSE],
                              ref$vectors[, ord_ref[j], drop = FALSE]), 1e-8)
  }

  raw <- eigencore:::native_dense_symmetric_eigen(A)
  expect_identical(raw$driver, "dsyevr")
  expect_false(is.unsorted(raw$values))
})

test_that("eig_full values-only path computes no vectors and stays uncertified", {
  set.seed(2102)
  n <- 40
  X <- matrix(rnorm(n * n), n)
  A <- X + t(X)
  fit <- eig_full(A, vectors = FALSE)
  ref <- eigen(A, symmetric = TRUE, only.values = TRUE)$values

  expect_null(fit$vectors)
  expect_lt(t2_rel(sort(values(fit)), sort(ref)), 1e-10)
  cert <- certificate(fit)
  expect_false(isTRUE(cert$passed))

  raw <- eigencore:::native_dense_symmetric_eigen(A, vectors = FALSE)
  expect_null(raw$vectors)
  expect_equal(raw$values, sort(ref), tolerance = 1e-12)

  Z <- matrix(complex(real = rnorm(n * n), imaginary = rnorm(n * n)), n)
  H <- Z + Conj(t(Z))
  cfit <- eig_full(H, vectors = FALSE)
  expect_null(cfit$vectors)
  expect_lt(t2_rel(sort(values(cfit)),
                   sort(eigen(H, symmetric = TRUE, only.values = TRUE)$values)),
            1e-10)
  craw <- eigencore:::native_dense_complex_hermitian_eigen(H, vectors = FALSE)
  expect_null(craw$vectors)
  cvec <- eig_full(H)
  expect_true(certificate(cvec)$passed)
})

test_that("dsyevr keeps an orthonormal basis for repeated eigenvalues", {
  set.seed(2103)
  n <- 50
  Q <- qr.Q(qr(matrix(rnorm(n * n), n)))
  lambda <- c(rep(3, 20), rep(-1, 15), 10 + seq_len(15))
  A <- Q %*% (lambda * t(Q))
  A <- (A + t(A)) / 2
  fit <- eig_full(A)
  V <- vectors(fit)

  expect_true(certificate(fit)$passed)
  expect_lt(max(abs(crossprod(V) - diag(n))), 1e-10)
  expect_lt(t2_rel(sort(values(fit)), sort(lambda)), 1e-10)
  # The 20-dimensional eigenspace of 3 is recovered as a subspace.
  idx3 <- which(abs(values(fit) - 3) < 1e-8)
  expect_length(idx3, 20L)
  expect_lt(t2_subspace_gap(V[, idx3], Q[, 1:20]), 1e-8)
})

test_that("eig_partial dense Hermitian uses dsyevr with subset selection", {
  set.seed(2104)
  n <- 40
  X <- matrix(rnorm(n * n), n)
  A <- X + t(X)
  ref <- eigen(A, symmetric = TRUE)

  top <- eig_partial(A, k = 15, target = largest())
  expect_identical(top$restart$eigensolver, "lapack_dsyevr_selected")
  expect_lt(t2_rel(values(top), ref$values[1:15]), 1e-10)
  expect_true(certificate(top)$passed)

  low <- eig_partial(A, k = 15, target = smallest())
  expect_lt(t2_rel(values(low), rev(ref$values)[1:15]), 1e-10)
  expect_lt(t2_subspace_gap(vectors(low), ref$vectors[, n:(n - 14)]), 1e-8)

  sel <- eigencore:::native_dense_symmetric_eigen_selected(A, 5, largest(),
                                                           vectors = FALSE)
  expect_null(sel$vectors)
  expect_equal(sel$values, ref$values[1:5], tolerance = 1e-12)

  mag <- eig_partial(A, k = 12, target = largest_magnitude())
  expect_identical(mag$restart$eigensolver, "lapack_dsyevr_full")
  expect_true(certificate(mag)$passed)
})

test_that("native dense SVD uses dgesdd and matches base svd", {
  set.seed(2105)
  for (dims in list(c(70L, 30L), c(30L, 70L), c(25L, 25L))) {
    A <- matrix(rnorm(prod(dims)), dims[1L], dims[2L])
    fit <- eigencore:::native_dense_svd(A)
    ref <- svd(A)
    expect_identical(fit$driver, "dgesdd")
    expect_lt(t2_rel(fit$d, ref$d), 1e-10)
    expect_lt(max(abs(fit$u %*% (fit$d * t(fit$v)) - A)), 1e-10 * max(ref$d))
    expect_lt(max(abs(crossprod(fit$u) - diag(min(dims)))), 1e-10)
    expect_lt(max(abs(crossprod(fit$v) - diag(min(dims)))), 1e-10)
    expect_lt(t2_subspace_gap(fit$u[, 1:3], ref$u[, 1:3]), 1e-8)
  }
  A <- matrix(rnorm(80 * 30), 80)
  sfit <- svd_partial(A, rank = 4, target = smallest())
  expect_lt(t2_rel(sfit$d, sort(svd(A)$d)[1:4]), 1e-10)
  expect_true(certificate(sfit)$passed)
})

test_that("bidiagonal SVD on (d, e) matches dense SVD including near-singular cases", {
  set.seed(2106)
  bidiag <- function(alpha, beta) {
    n <- length(alpha)
    B <- diag(alpha, n)
    if (n > 1) B[cbind(seq_len(n - 1), 2:n)] <- beta
    B
  }
  cases <- list(
    list(alpha = runif(30, 0.5, 2), beta = runif(29, 0.1, 1)),
    list(alpha = c(1, 1e-14, 2, 3, 1e-300, 4), beta = c(1, 1, 0, 1e-12, 2)),
    list(alpha = c(0, 0, 0), beta = c(1, 0)),
    list(alpha = c(-2, 3, -1), beta = c(0.5, -0.25)),
    list(alpha = 5, beta = numeric())
  )
  for (case in cases) {
    B <- bidiag(case$alpha, case$beta)
    fit <- eigencore:::native_bidiagonal_svd(case$alpha, case$beta)
    n <- length(case$alpha)
    expect_equal(dim(fit$u), c(n, n))
    expect_equal(dim(fit$v), c(n, n))
    expect_false(is.unsorted(rev(fit$d)))
    expect_true(all(fit$d >= 0))
    expect_lt(max(abs(fit$d - svd(B)$d)), 1e-12 * max(1, svd(B)$d[1]))
    expect_lt(max(abs(fit$u %*% (fit$d * t(fit$v)) - B)), 1e-12 * max(1, fit$d[1]))
    expect_lt(max(abs(crossprod(fit$u) - diag(n))), 1e-12)
    expect_lt(max(abs(crossprod(fit$v) - diag(n))), 1e-12)
  }
})

test_that("tridiagonal eigen via dstevr matches dense eigen", {
  set.seed(2107)
  for (n in c(1L, 2L, 7L, 120L)) {
    alpha <- rnorm(n)
    beta <- rnorm(max(n - 1L, 0L))
    T <- diag(alpha, n)
    if (n > 1) {
      T[cbind(1:(n - 1), 2:n)] <- beta
      T[cbind(2:n, 1:(n - 1))] <- beta
    }
    fit <- eigencore:::native_tridiagonal_eigen(alpha, beta)
    ref <- eigen(T, symmetric = TRUE)
    expect_false(is.unsorted(fit$values))
    expect_lt(t2_rel(fit$values, rev(ref$values)), 1e-10)
    expect_lt(max(abs(T %*% fit$vectors - fit$vectors %*% diag(fit$values, n))),
              1e-10 * max(1, abs(ref$values)))
    expect_lt(max(abs(crossprod(fit$vectors) - diag(n))), 1e-10)
  }
  # Glued Wilkinson-like matrix with tight clusters.
  n <- 41
  alpha <- abs(seq(-20, 20))
  beta <- rep(1, n - 1)
  fit <- eigencore:::native_tridiagonal_eigen(alpha, beta)
  expect_lt(max(abs(crossprod(fit$vectors) - diag(n))), 1e-10)
})

test_that("tridiagonal solve pivots on indefinite systems and rejects singular ones", {
  solve_tri <- function(lower, diag, upper, B) {
    .Call("eigencore_tridiagonal_solve", as.numeric(lower), as.numeric(diag),
          as.numeric(upper), as.matrix(B), PACKAGE = "eigencore")
  }
  dense_tri <- function(lower, diag, upper) {
    n <- length(diag)
    T <- diag(diag, n)
    if (n > 1) {
      T[cbind(2:n, 1:(n - 1))] <- lower
      T[cbind(1:(n - 1), 2:n)] <- upper
    }
    T
  }
  # Zero leading diagonal: the unpivoted Thomas algorithm breaks down here.
  lower <- c(1, 1, 1)
  diag <- c(0, 0, 0, 0)
  upper <- c(1, 1, 1)
  B <- cbind(c(1, 2, 3, 4), c(-1, 0, 1, 0))
  X <- solve_tri(lower, diag, upper, B)
  expect_lt(max(abs(dense_tri(lower, diag, upper) %*% X - B)), 1e-12)

  # Indefinite shifted Laplacian, shift between eigenvalues; tiny leading
  # pivot (1e-20) that Thomas would amplify catastrophically.
  set.seed(2108)
  n <- 200
  lower <- rep(-1, n - 1)
  upper <- rep(-1, n - 1)
  diag <- rep(2, n) - 1.37
  diag[1] <- 1e-20
  B <- matrix(rnorm(n * 3), n)
  T <- dense_tri(lower, diag, upper)
  X <- solve_tri(lower, diag, upper, B)
  expect_lt(max(abs(T %*% X - B)) / max(abs(B)), 1e-10)
  expect_lt(max(abs(X - solve(T, B))) / max(abs(X)), 1e-8)

  # Exactly singular: [1 1; 1 1].
  expect_error(solve_tri(1, c(1, 1), 1, diag(2)), "singular")
  # Singular to working precision.
  expect_error(solve_tri(1, c(1, 1 + 1e-17), 1, diag(2)), "singular")
  expect_error(solve_tri(c(1, NA), c(1, 2, 3), c(1, 1), diag(3)), "non-finite")
})
