# C12: certificates report the normwise 2-norm backward error with a
# denominator that never exceeds ||A||_2, so the reported value is an upper
# bound on the true backward error and `passed` is sound.

true_eigen_backward <- function(A, values, vectors, B = NULL) {
  nA <- norm(as.matrix(A), "2")
  nB <- if (is.null(B)) 1 else norm(as.matrix(B), "2")
  Bv <- if (is.null(B)) vectors else as.matrix(B) %*% vectors
  r <- sqrt(colSums(Mod(as.matrix(A) %*% vectors - sweep(Bv, 2L, values, `*`))^2))
  r / ((nA + Mod(values) * nB) * sqrt(colSums(Mod(vectors)^2)))
}

true_svd_backward <- function(A, d, u, v) {
  A <- as.matrix(A)
  left <- sqrt(colSums((A %*% v - sweep(u, 2L, d, `*`))^2))
  right <- sqrt(colSums((crossprod(A, u) - sweep(v, 2L, d, `*`))^2))
  sqrt(left^2 + right^2) / norm(A, "2")
}

# The reported bound may undercut a residual recomputed in R only by the
# rounding error of evaluating that residual (about n * eps after scaling by
# ||A||_2): converged pairs sit at ~1e-16, where different BLAS kernels round
# differently. Perturbed pairs below keep the comparison sharp.
expect_bounds_truth <- function(reported, truth, n) {
  slack <- n * .Machine$double.eps
  expect_true(all(reported >= truth * (1 - 1e-10) - slack),
              info = paste(signif(reported, 3), signif(truth, 3), collapse = " | "))
}

test_that("C12: reported eigen backward error bounds the true normwise one", {
  set.seed(101)
  for (trial in 1:4) {
    n <- 60L
    X <- matrix(rnorm(n * n), n)
    A <- crossprod(X) / n - diag(n) * (trial - 2)
    for (target in list(largest(), smallest())) {
      fit <- eig_partial(A, k = 4, target = target, tol = 1e-10)
      cert <- fit$certificate
      truth <- true_eigen_backward(A, fit$values, fit$vectors)
      expect_bounds_truth(cert$backward_error, truth, n)
      expect_true(cert$norm_values[["A"]] <= norm(A, "2") * (1 + 1e-12))
      expect_false(cert$scale_is_estimate)
      expect_match(cert$norm_bound_type, "^two_norm_(exact|lower_bound)\\+identity_exact$")
    }
  }
  # Perturbed (unconverged) pairs, nonsymmetric matrix, right certificate.
  G <- matrix(rnorm(40 * 40), 40)
  e <- eigen(G)
  vec <- e$vectors[, 1:3] + 1e-4 * complex(real = rnorm(120), imaginary = rnorm(120))
  cert <- eigencore:::certify_dense_general_eigen(G, e$values[1:3], vec, tol = 1e-8)
  truth <- true_eigen_backward(G, e$values[1:3], vec)
  expect_bounds_truth(cert$backward_error, truth, 40)
})

test_that("C12: reported SVD backward error bounds the true normwise one", {
  set.seed(102)
  for (trial in 1:3) {
    A <- matrix(rnorm(80 * 30), 80, 30)
    for (target in list(largest(), smallest())) {
      fit <- svd_partial(A, rank = 3, target = target, tol = 1e-10)
      cert <- fit$certificate
      truth <- true_svd_backward(A, fit$d, fit$u, fit$v)
      expect_bounds_truth(cert$backward_error, truth, 80)
      expect_true(cert$norm_values[["A"]] <= norm(A, "2") * (1 + 1e-12))
    }
    # Deliberately inaccurate triplets.
    s <- svd(A)
    u <- s$u[, 1:3] + 1e-5 * matrix(rnorm(80 * 3), 80)
    v <- s$v[, 1:3]
    cert <- eigencore:::certify_svd(A, s$d[1:3], u, v, tol = 1e-8)
    truth <- true_svd_backward(A, s$d[1:3], u, v)
    expect_bounds_truth(cert$backward_error, truth, 80)
    S <- Matrix::rsparsematrix(200, 50, density = 0.1)
    cert_op <- eigencore:::certify_svd_operator(
      as_operator(S), s$d[1:3], u[1:3, , drop = FALSE][rep(1:3, length.out = 200), ],
      v[rep(1:30, length.out = 50), ], tol = 1e-8
    )
    truth_op <- true_svd_backward(
      S, s$d[1:3], u[1:3, , drop = FALSE][rep(1:3, length.out = 200), ],
      v[rep(1:30, length.out = 50), ]
    )
    expect_bounds_truth(cert_op$backward_error, truth_op, 200)
  }
})

test_that("C12: a flat spectrum residual that passed under Frobenius now fails", {
  set.seed(103)
  n <- 400L
  Q <- qr.Q(qr(matrix(rnorm(n * n), n)))
  d <- seq(1.9, 2, length.out = n)
  A <- Q %*% (d * t(Q))
  A <- (A + t(A)) / 2
  eps <- 1e-6
  x <- Q[, 1] + eps * Q[, n]
  lambda <- d[[1]]
  r <- sqrt(sum((A %*% x - lambda * x)^2))
  xn <- sqrt(sum(x^2))
  tol <- 1e-8
  frobenius_eta <- r / ((norm(A, "F") + abs(lambda)) * xn)
  two_norm_eta <- r / ((norm(A, "2") + abs(lambda)) * xn)
  # The residual sits between tol * ||A||_2 and tol * ||A||_F.
  expect_lt(frobenius_eta, tol)
  expect_gt(two_norm_eta, tol)

  cert <- eigencore:::certify_eigen(A, lambda, matrix(x), tol = tol)
  expect_false(cert$passed)
  expect_gte(cert$backward_error, two_norm_eta * (1 - 1e-10))
  cert_op <- eigencore:::certify_eigen_operator(as_operator(A), lambda, matrix(x),
                                                tol = tol)
  expect_false(cert_op$passed)
  expect_gte(cert_op$backward_error, two_norm_eta * (1 - 1e-10))
})

test_that("C12: matrix-free operators certify with lower-bound norms", {
  set.seed(104)
  n <- 80L
  Q <- qr.Q(qr(matrix(rnorm(n * n), n)))
  A <- Q %*% (c(10, 8, 6, seq(2, 0.1, length.out = n - 3)) * t(Q))
  A <- (A + t(A)) / 2
  op <- linear_operator(
    dim = dim(A),
    apply = function(X, alpha = 1, beta = 0, Y = NULL) {
      Z <- alpha * (A %*% X)
      if (is.null(Y) || beta == 0) Z else Z + beta * Y
    },
    structure = hermitian()
  )
  fit <- eig_partial(op, k = 3, target = largest(), tol = 1e-9)
  expect_true(fit$certificate$passed)
  expect_false(fit$certificate$scale_is_estimate)
  expect_identical(fit$certificate$norm_bound_type,
                   "two_norm_lower_bound+identity_exact")
  expect_false(any(grepl("stochastic", fit$certificate$notes)))
  truth <- true_eigen_backward(A, fit$values, fit$vectors)
  expect_true(all(fit$certificate$backward_error >= truth * (1 - 1e-10)))

  B <- matrix(rnorm(60 * 25), 60, 25)
  svd_op <- linear_operator(
    dim = dim(B),
    apply = function(X, alpha = 1, beta = 0, Y = NULL) alpha * (B %*% X),
    apply_adjoint = function(X, alpha = 1, beta = 0, Y = NULL) alpha * crossprod(B, X)
  )
  sfit <- svd_partial(svd_op, rank = 2, target = largest(), tol = 1e-9)
  expect_true(sfit$certificate$passed)
  expect_identical(sfit$certificate$norm_bound_type, "two_norm_lower_bound")
  expect_true(all(sfit$certificate$backward_error >=
                    true_svd_backward(B, sfit$d, sfit$u, sfit$v) * (1 - 1e-10)))
})

test_that("C12: diagonal operators use the exact two-norm", {
  d <- c(5, -7, 3, 1, 0.5)
  D <- Matrix::Diagonal(x = d)
  fit <- eig_partial(D, k = 2, target = largest_magnitude())
  cert <- fit$certificate
  expect_true(cert$passed)
  expect_identical(cert$norm_bound_type, "two_norm_exact+identity_exact")
  expect_identical(cert$norm_source, "diagonal+identity")
  expect_equal(unname(cert$norm_values[["A"]]), 7)
  expect_equal(cert$scale,
               (7 + abs(fit$values)) * sqrt(colSums(fit$vectors^2)))

  # Exact scale even when the certified pair is far from the dominant end.
  x <- c(0, 0, 0, 0, 1)
  cert_small <- eigencore:::certify_eigen_operator(as_operator(D), 0.5, matrix(x))
  expect_identical(cert_small$norm_bound_type, "two_norm_exact+identity_exact")
  expect_equal(unname(cert_small$norm_values[["A"]]), 7)

  sfit <- svd_partial(D, rank = 2, target = largest())
  expect_identical(sfit$certificate$norm_bound_type, "two_norm_exact")
  expect_equal(unname(sfit$certificate$norm_values[["A"]]), 7)
})

test_that("C12: the Krylov norm bound refines smallest targets without using the RNG", {
  set.seed(105)
  n <- 150L
  Q <- qr.Q(qr(matrix(rnorm(n * n), n)))
  lam <- seq(1, 200, length.out = n)
  A <- Q %*% (lam * t(Q))
  A <- (A + t(A)) / 2
  calls <- 0L
  op <- linear_operator(
    dim = dim(A),
    apply = function(X, alpha = 1, beta = 0, Y = NULL) {
      calls <<- calls + 1L
      Z <- alpha * (A %*% X)
      if (is.null(Y) || beta == 0) Z else Z + beta * Y
    },
    structure = hermitian()
  )
  # Smallest pairs with residuals of ~1e-6: failing against the free bound
  # ||A x|| / ||x|| ~ lambda_min, passing against ||A||_2 = 200.
  values <- lam[1:2]
  vectors <- Q[, 1:2] + 5e-9 * Q[, n:(n - 1)]
  tol <- 1e-8
  truth <- true_eigen_backward(A, values, vectors)
  expect_true(all(truth < tol))

  set.seed(99)
  seed_before <- .Random.seed
  cert <- eigencore:::certify_eigen_operator(op, values, vectors, tol = tol)
  expect_identical(.Random.seed, seed_before)
  expect_true(cert$passed)
  expect_identical(cert$norm_source, "lanczos+identity")
  expect_lte(cert$norm_values[["A"]], 200 * (1 + 1e-12))
  expect_gt(cert$norm_values[["A"]], 0.9 * 200)
  expect_true(all(cert$backward_error >= truth * (1 - 1e-10)))

  # Memoised per operator: a second certificate adds only its own apply.
  used <- calls
  cert2 <- eigencore:::certify_eigen_operator(op, values, vectors, tol = tol)
  expect_identical(calls, used + 1L)
  expect_identical(cert2$norm_values, cert$norm_values)

  # Pairs that fail even against an upper bound do not trigger refinement.
  dense_op <- as_operator(A)
  bad <- Q[, 1:2] + 1e-2 * Q[, n:(n - 1)]
  cert_bad <- eigencore:::certify_eigen_operator(dense_op, values, bad, tol = tol)
  expect_false(cert_bad$passed)
  expect_false(grepl("lanczos", cert_bad$norm_source))
})

test_that("C39: complex eigenvector certificates apply real operators in two real passes", {
  set.seed(106)
  n <- 30L
  G <- matrix(rnorm(n * n), n)
  real_only <- linear_operator(
    dim = dim(G),
    apply = function(X, alpha = 1, beta = 0, Y = NULL) {
      if (is.complex(X)) stop("complex block reached a real apply")
      alpha * (G %*% X)
    },
    apply_adjoint = function(X, alpha = 1, beta = 0, Y = NULL) {
      if (is.complex(X)) stop("complex block reached a real apply")
      alpha * crossprod(G, X)
    }
  )
  e <- eigen(G)
  idx <- which(Im(e$values) != 0)[1:2]
  right <- eigencore:::certify_general_eigen_operator(real_only, e$values[idx],
                                                      e$vectors[, idx], tol = 1e-8)
  dense <- eigencore:::certify_dense_general_eigen(G, e$values[idx],
                                                   e$vectors[, idx], tol = 1e-8)
  expect_true(right$passed)
  expect_equal(right$residuals, dense$residuals, tolerance = 1e-8)

  left <- eigen(t(G))
  lidx <- vapply(e$values[idx], function(z) which.min(Mod(left$values - z)), 1L)
  W <- left$vectors[, lidx]
  cert_left <- eigencore:::certify_left_eigen_operator(real_only, e$values[idx], W,
                                                       tol = 1e-8)
  expect_true(all(cert_left$converged))

  S <- Matrix::rsparsematrix(n, n, density = 0.2) + Matrix::Diagonal(n)
  S <- methods::as(S, "generalMatrix")
  es <- eigen(as.matrix(S))
  sidx <- which(Im(es$values) != 0)[1:2]
  sparse_cert <- eigencore:::certify_general_eigen_operator(
    as_operator(S), es$values[sidx], es$vectors[, sidx], tol = 1e-8
  )
  expect_true(sparse_cert$passed)
})
