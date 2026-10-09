# Minimised reproducers for invariant violations found by the oracle sweep
# (helper-oracle.R; docs/test-assurance.md). Each block names the sweep case
# that found it. Tracker ids O1..On are listed in docs/test-assurance.md.

oracle_true_backward <- function(A, values, vectors) {
  nA <- max(abs(eigen(A, symmetric = TRUE, only.values = TRUE)$values))
  R <- A %*% vectors - sweep(vectors, 2L, values, `*`)
  sqrt(colSums(R^2)) / ((nA + abs(values)) * sqrt(colSums(vectors^2)))
}

test_that("O1: backward errors of operators with ||A||_2 < eps are not floored (case 17, 342)", {
  # The certificate divided by max(||A|| + |lambda|, eps): for a 1e-24-norm
  # operator that made the reported backward error ~1e8 times too small and
  # certified pairs whose true backward error was 1e-5.
  set.seed(1)
  X <- matrix(rnorm(12 * 8), 12, 8) * 1e-12
  A <- crossprod(X)
  fit <- eig_partial(crossprod_operator(X), k = 2, tol = 1e-8)
  true_be <- oracle_true_backward(A, fit$values, fit$vectors)
  expect_true(all(fit$certificate$backward_error >= 0.99 * true_be))
  if (isTRUE(fit$certificate$passed)) {
    expect_lte(max(true_be), 1.05e-8)
  }

  D <- diag(c(5, 4, 3, 2, 1) * 1e-20)
  sfit <- svd_partial(D, rank = 2)
  r <- sqrt(colSums((D %*% sfit$v - sweep(sfit$u, 2L, sfit$d, `*`))^2) +
              colSums((crossprod(D, sfit$u) - sweep(sfit$v, 2L, sfit$d, `*`))^2))
  expect_true(all(sfit$certificate$backward_error >= 0.99 * r / 5e-20))

  # The scale itself is never floored at eps.
  op <- as_operator(diag(c(3, 2, 1) * 1e-30))
  cert <- eigencore:::certify_eigen_operator(op, c(3e-30, 2e-30),
                                            diag(3)[, 1:2] + 1e-3, tol = 1e-8)
  expect_lt(max(cert$scale), 1e-20)
  expect_false(cert$passed)
})

test_that("O2: a certificate never passes with fewer than k pairs (case 188)", {
  # Shift-invert Krylov on a spectrum with repeated eigenvalues exhausts its
  # Krylov space at the number of distinct eigenvalues and returned
  # 1, 7, 9 for k = 6 with certificate$passed = TRUE.
  set.seed(3)
  Q <- qr.Q(qr(matrix(rnorm(64), 8)))
  A <- Q %*% diag(c(1, 1, 7, 7, 7, 9, 9, 9)) %*% t(Q)
  A <- (A + t(A)) / 2
  fit <- eig_partial(A, k = 6, target = nearest(2), method = shift_invert(2))
  if (length(fit$values) < 6L) {
    expect_false(fit$certificate$passed)
    expect_true(any(grepl("requested pairs were returned", fit$certificate$notes)))
  } else {
    expect_equal(sort(fit$values), c(1, 1, 7, 7, 7, 9), tolerance = 1e-8)
  }
})

test_that("O3: shift_invert(sigma) refuses targets other than nearest(sigma) (case 191)", {
  # The shift-invert routes computed the eigenvalues nearest sigma but
  # labelled and ordered them by the problem's target: target = smallest()
  # returned a certified "smallest" set that was not the smallest.
  A <- diag(c(-6.7, -2.9, -0.8, 1.7, 4.9, 7.2))
  expect_error(
    eig_partial(A, k = 3, target = smallest_magnitude(), method = shift_invert(2.9)),
    "nearest sigma"
  )
  expect_error(
    eig_partial(A, k = 3, target = smallest(), method = shift_invert(2.9)),
    "nearest sigma"
  )
  expect_error(
    eig_partial(A, k = 3, target = nearest(1), method = shift_invert(2.9)),
    "nearest sigma"
  )
  # Without an explicit target, shift_invert(sigma) means nearest(sigma).
  fit <- eig_partial(A, k = 3, method = shift_invert(2.9))
  expect_equal(sort(fit$values), c(-0.8, 1.7, 4.9), tolerance = 1e-8)
  expect_identical(fit$target, "nearest(2.9)")
  # smallest_magnitude with sigma = 0 is the same set as nearest(0).
  fit0 <- eig_partial(A, k = 2, target = smallest_magnitude(), method = shift_invert(0))
  expect_equal(sort(fit0$values), c(-0.8, 1.7), tolerance = 1e-8)
})

test_that("O4: tiny-scale complex general matrices are not routed as Hermitian (case 129, 309)", {
  # isSymmetric.matrix() compares absolutely when mean(Mod(x)) < tol, so a
  # 1e-12-scale general complex matrix was solved with zheev.
  set.seed(4)
  n <- 8
  Z <- matrix(complex(real = rnorm(n * n), imaginary = rnorm(n * n)), n) * 1e-12
  expect_identical(as_operator(Z)$structure$kind, "general")
  fit <- eig_partial(Z, k = 1, target = largest_real())
  ev <- eigen(Z, only.values = TRUE)$values
  expect_equal(Re(fit$values), max(Re(ev)), tolerance = 1e-6)
  H <- Z + Conj(t(Z))
  expect_identical(as_operator(H)$structure$kind, "hermitian")
})

test_that("O5: reference Golub-Kahan breakdown is relative; no DSYRK error on tiny input (case 282)", {
  # An absolute breakdown threshold ended the 1e-12-scale run before its
  # first step; the empty U/V then reached dsyrk with ldc = 0.
  D <- Matrix::Diagonal(x = c(3, 2, 1.5, 1, 0.5) * 1e-12)
  fit <- svd_partial(D, rank = 1, method = golub_kahan())
  expect_equal(fit$d, 3e-12, tolerance = 1e-8)
  fit2 <- svd_partial(D, rank = 1, target = nearest(1.2e-12), method = golub_kahan())
  expect_length(fit2$d, 1L)
  expect_equal(fit2$d, 1e-12, tolerance = 1e-6)
  # the zero operator still breaks down cleanly
  Z <- Matrix::Diagonal(x = rep(0, 4))
  expect_no_error(svd_partial(Z, rank = 1, method = golub_kahan()))
})

# ---------------------------------------------------------------------------
# Known, not fixed here (see docs/test-assurance.md). These are certified
# results whose SET is wrong while the certificate reports
# target_completeness = "not_checked" (or has no completeness field): the
# per-pair residual certificate is sound, target identity is not checked on
# these routes.
# ---------------------------------------------------------------------------

repeated_spectrum_matrix <- function(values, seed = 5) {
  set.seed(seed)
  n <- length(values)
  Q <- qr.Q(qr(matrix(rnorm(n * n), n)))
  A <- Q %*% diag(values) %*% t(Q)
  (A + t(A)) / 2
}

test_that("O6 (known): nearest(sigma) misses copies of repeated eigenvalues (case 194)", {
  skip("known: O6 / C50 (completeness probe does not cover nearest/shift-invert)")
  A <- repeated_spectrum_matrix(c(-1, 1, 7, 7, 7, 7, 9, 9))
  fit <- eig_partial(A, k = 6, target = nearest(7.5))
  if (isTRUE(fit$certificate$passed)) {
    expect_equal(sort(fit$values), c(7, 7, 7, 7, 9, 9), tolerance = 1e-8)
  }
})

test_that("O7 (known): both_ends / nearest on the reference Lanczos route miss copies (case 341)", {
  skip("known: O7 / C50 (reference Lanczos route is not probed)")
  A <- repeated_spectrum_matrix(rep(c(-1, 1, 3, 5, 7, 9), length.out = 50))
  fit <- eig_partial(A, k = 5, target = both_ends(2, 3), method = lanczos(block = 2))
  if (isTRUE(fit$certificate$passed)) {
    expect_equal(sort(fit$values), c(-1, -1, 9, 9, 9), tolerance = 1e-8)
  }
})

test_that("O8 (known): SVD routes miss copies of repeated singular values (case 215)", {
  skip("known: O8 (no target-completeness check for SVD)")
  set.seed(215)
  d <- sample(c(9, 7, 5, 3, 1), 20, replace = TRUE) + 1e-10 * rnorm(20)
  fit <- svd_partial(Matrix::Diagonal(x = d), rank = 6)
  if (isTRUE(fit$certificate$passed)) {
    expect_equal(fit$d, sort(d, decreasing = TRUE)[1:6], tolerance = 1e-8)
  }
})

test_that("O9 (known): matrix-free smallest_magnitude misses a multiple zero eigenvalue (case 133)", {
  skip("known: O9 / C50 (smallest_magnitude not probed)")
  A <- repeated_spectrum_matrix(c(rep(0, 10), seq(-3, 3, length.out = 20)))
  f <- function(x, args) as.numeric(A %*% x)
  fit <- eigs_sym(f, k = 6, which = "SM", n = 30, opts = list(ncv = 20))
  if (isTRUE(fit$certificate$passed)) {
    expect_lt(max(abs(fit$values)), 1e-8)
  }
})

test_that("O10 (known): sparse both_ends takes the unrestarted reference Lanczos and does not certify", {
  skip("known: O10 (quality: sparse/matrix-free both_ends has no native route)")
  set.seed(1)
  n <- 40
  vals <- sample(c(-1, 1), n, replace = TRUE) * (seq_len(n) + stats::runif(n, 0, 0.5))
  Q <- qr.Q(qr(matrix(rnorm(n * n), n)))
  A <- Q %*% (vals * t(Q))
  A <- (A + t(A)) / 2
  S <- methods::as(methods::as(Matrix::Matrix(A, sparse = TRUE), "generalMatrix"),
                   "CsparseMatrix")
  fit <- eig_partial(S, k = 1, target = both_ends(0, 1), tol = 1e-10)
  expect_true(fit$certificate$passed)
})

test_that("O11: svds(nu = 0, nv = 0) is still certified (shim solves with both sides)", {
  set.seed(9)
  S <- Matrix::rsparsematrix(80, 30, 0.2)
  for (nunv in list(c(0L, 0L), c(1L, 0L), c(0L, 3L))) {
    fit <- expect_no_warning(svds(S, 3, nu = nunv[1], nv = nunv[2],
                                  opts = list(tol = 1e-10, center = TRUE)))
    expect_true(fit$certificate$passed)
    expect_identical(fit$nconv, 3L)
    expect_equal(NCOL(fit$u) * (nunv[1] > 0), nunv[1])
    expect_equal(NCOL(fit$v) * (nunv[2] > 0), nunv[2])
  }
})

test_that("O12: randomized SVD ranks nearest(sigma) by singular value, not its square (case 502)", {
  # The wide-core Gram path ordered the Gram eigenvalues d^2 by |d^2 - sigma|.
  set.seed(12)
  m <- 12
  n <- 80
  d <- c(10, 8, 6, 5, 4, 3, 2.5, 2, 1.5, 1, 0.5, 0.25)
  U <- qr.Q(qr(matrix(rnorm(m * m), m)))
  V <- qr.Q(qr(matrix(rnorm(n * m), n)))
  A <- U %*% diag(d) %*% t(V)
  fit <- suppressWarnings(svd_partial(A, rank = 3, target = nearest(4.2),
                                      method = randomized(), seed = 1))
  expect_equal(sort(fit$d), c(3, 4, 5), tolerance = 1e-6)
  expect_identical(order(abs(fit$d - 4.2)), seq_along(fit$d))
})

test_that("O13: implicit smallest_magnitude on a singular sparse nonsymmetric matrix (case 6743)", {
  # The sparse LU of A - 0 I factored, but Matrix::solve() refused it at apply
  # time, surfacing as "Krylov-Schur Arnoldi operator apply failed with
  # status=-8"; singularity is now detected at factor time and the implicit
  # route perturbs sigma.
  # The generator's matrix: V J V^-1 with a 3x3 Jordan block at 0, sparse V.
  S <- oracle_build(oracle_case(6743L))$obj
  expect_s4_class(S, "dgCMatrix")
  expect_no_error(fit <- eig_partial(S, k = 3, target = smallest_magnitude()))
  expect_error(eig_partial(S, k = 2, target = nearest(0), method = shift_invert(0)),
               "singular")
})

test_that("O14: eigs_sym keeps per-pair certificate fields aligned with re-sorted values", {
  set.seed(14)
  S <- Matrix::rsparsematrix(300, 300, 0.02, symmetric = TRUE)
  A <- as.matrix(S)
  fit <- eigs_sym(S, 4, which = "SA")
  r <- sqrt(colSums((A %*% fit$vectors - sweep(fit$vectors, 2L, fit$values, `*`))^2))
  # residuals span orders of magnitude; compare on the log scale
  expect_lt(max(abs(log10(fit$certificate$residuals) - log10(r))), 0.5)
  expect_true(all(diff(fit$values) <= 0))
})

test_that("O15: unsupported method/size combinations fail with a clear message", {
  X <- matrix(rnorm(30), 6)
  expect_error(svd_partial(X, 1, method = lanczos()), "eigensolver method")
  expect_error(svd_partial(X, 1, method = lobpcg()), "eigensolver method")
  set.seed(15)
  S <- Matrix::rsparsematrix(20, 20, 0.3, symmetric = TRUE)
  expect_error(eig_partial(S, 20), "operator dimension n = 20")
  expect_error(eig_partial(S, 19, method = lanczos(block = 3)), "smaller block")
})

# Certified results whose set is not the target set on non-repeated spectra
# (completeness "not_checked"); reproduced from sweep case ids.
expect_certified_set_ok <- function(id) {
  rec <- oracle_run_case(id)
  expect_true(!isTRUE(rec$certified) || isTRUE(rec$set_ok),
              label = sprintf("case %d [%s] certified set", id, rec$describe))
}

test_that("O16 (known): nonsymmetric Arnoldi certifies pairs outside the LI/SI/LM target set", {
  skip("known: O16 (nonsymmetric target identity is not checked)")
  for (id in c(2692L, 2783L, 4543L, 7072L, 13198L)) expect_certified_set_ok(id)
})

test_that("O17 (known): LOBPCG magnitude targets on indefinite problems miss the other end", {
  skip("known: O17 (LOBPCG magnitude targets are not completeness-checked)")
  for (id in c(8847L, 10647L, 13154L)) expect_certified_set_ok(id)
})

test_that("O18 (known): matrix-free nonsymmetric smallest_magnitude certifies non-smallest pairs", {
  skip("known: O18 (no shift-invert for matrix-free nonsymmetric SM; not checked)")
  set.seed(518)
  M <- matrix(rnorm(2500), 50) / sqrt(50)
  op <- linear_operator(
    dim = c(50, 50),
    apply = function(X, alpha = 1, beta = 0, Y = NULL) {
      Z <- alpha * (M %*% X)
      if (is.null(Y) || beta == 0) Z else Z + beta * Y
    },
    apply_adjoint = function(X, alpha = 1, beta = 0, Y = NULL) {
      Z <- alpha * crossprod(M, X)
      if (is.null(Y) || beta == 0) Z else Z + beta * Y
    })
  fit <- eig_partial(op, 1, target = smallest_magnitude(), tol = 1e-10)
  if (isTRUE(fit$certificate$passed)) {
    expect_equal(Mod(fit$values), min(Mod(eigen(M, only.values = TRUE)$values)),
                 tolerance = 1e-6)
  }
})
