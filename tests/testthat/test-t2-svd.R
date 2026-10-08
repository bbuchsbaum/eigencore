# Tranche 2 SVD regressions: certificate soundness (C13), Ritz extraction
# (P12), Gram build (P17) and CSC structure validation (C28).

t2_two_sided <- function(A, d, u, v) {
  A <- as.matrix(A)
  left <- sqrt(colSums((A %*% v - sweep(u, 2L, d, `*`))^2))
  right <- sqrt(colSums((crossprod(A, u) - sweep(v, 2L, d, `*`))^2))
  k <- length(d)
  list(
    left = left,
    right = right,
    combined = sqrt(left^2 + right^2),
    backward = sqrt(left^2 + right^2) / sqrt(sum(A^2)),
    orth_u = max(abs(crossprod(u) - diag(k))),
    orth_v = max(abs(crossprod(v) - diag(k)))
  )
}

t2_graded_matrix <- function(m, n, floor_exponent, seed) {
  set.seed(seed)
  sv <- 10^seq(0, floor_exponent, length.out = min(m, n))
  Q1 <- qr.Q(qr(matrix(rnorm(m * min(m, n)), m, min(m, n))))
  Q2 <- qr.Q(qr(matrix(rnorm(n * min(m, n)), n, min(m, n))))
  Q1 %*% (sv * t(Q2))
}

# Residuals reported by a certificate must agree with residuals recomputed
# from the original matrix: relative agreement where the true residual is
# above roundoff, and never an underestimate beyond roundoff.
t2_expect_residuals_match <- function(reported, true, floor) {
  reported <- as.numeric(reported)
  true <- as.numeric(true)
  expect_true(
    all(abs(reported - true) <= 1e-2 * true + floor),
    info = paste0(
      "reported=", paste(signif(reported, 3), collapse = ","),
      " true=", paste(signif(true, 3), collapse = ",")
    )
  )
}

test_that("C13a: materialized Gram SVD certificate tests both sides against A", {
  skip_if_not_installed("Matrix")
  for (wide in c(FALSE, TRUE)) {
    A <- t2_graded_matrix(200L, 30L, -4, seed = 1301)
    if (wide) A <- t(A)
    S <- methods::as(Matrix::Matrix(A, sparse = TRUE), "CsparseMatrix")
    fit <- eigencore:::native_gram_svd(S, 5L, target = smallest(), tol = 1e-12)
    expect_true(isTRUE(fit$restart$materialized_gram))
    truth <- t2_two_sided(A, fit$d, fit$u, fit$v)
    floor <- 64 * .Machine$double.eps * sqrt(sum(A^2))
    t2_expect_residuals_match(fit$certificate$residuals$left, truth$left, floor)
    t2_expect_residuals_match(fit$certificate$residuals$right, truth$right, floor)
    if (isTRUE(fit$certificate$passed)) {
      expect_lte(max(truth$backward), 1e-12)
    }
  }
})

test_that("C13a: public svd_partial never passes a Gram certificate the true residual fails", {
  skip_if_not_installed("Matrix")
  for (wide in c(FALSE, TRUE)) {
    for (floor_exp in c(-4, -4.5)) {
      A <- t2_graded_matrix(200L, 30L, floor_exp, seed = 1302)
      if (wide) A <- t(A)
      S <- methods::as(Matrix::Matrix(A, sparse = TRUE), "CsparseMatrix")
      fit <- svd_partial(S, 5L, target = smallest(), tol = 1e-12)
      truth <- t2_two_sided(A, fit$d, fit$u, fit$v)
      if (isTRUE(fit$certificate$passed)) {
        expect_lte(max(truth$backward), 1e-12)
        expect_lte(max(truth$orth_u, truth$orth_v), sqrt(.Machine$double.eps))
      }
    }
  }
})

test_that("C13b: native CSC Gram kernels report residuals and orthogonality of the returned U, V", {
  skip_if_not_installed("Matrix")
  for (side in c("right", "left")) {
    for (floor_exp in c(-4, -5)) {
      A <- t2_graded_matrix(200L, 30L, floor_exp, seed = 1303)
      if (side == "left") A <- t(A)
      S <- methods::as(Matrix::Matrix(A, sparse = TRUE), "CsparseMatrix")
      entry <- if (side == "right") "eigencore_csc_right_gram_svd" else "eigencore_csc_left_gram_svd"
      native <- .Call(entry, S@i, S@p, S@x, S@Dim, 30L, 1e-12, PACKAGE = "eigencore")
      truth <- t2_two_sided(A, native$d, native$u, native$v)
      floor <- 64 * .Machine$double.eps * sqrt(sum(A^2))
      t2_expect_residuals_match(native$diagnostics$left, truth$left, floor)
      t2_expect_residuals_match(native$diagnostics$right, truth$right, floor)
      orth <- native$diagnostics$orthogonality
      t2_expect_residuals_match(orth[["U"]], truth$orth_u, 64 * .Machine$double.eps)
      t2_expect_residuals_match(orth[["V"]], truth$orth_v, 64 * .Machine$double.eps)
    }
  }
})

test_that("C13b: default svd_partial on a wide sparse matrix does not certify a non-orthogonal V", {
  skip_if_not_installed("Matrix")
  # Before the fix the left-Gram fast path inferred V'V from U'GU and reported
  # an orthogonality loss of 1.3e-8 (< sqrt(eps)) while the returned V had a
  # true loss of 3.2e-8, so the certificate passed at the default tolerance.
  set.seed(5)
  A <- t2_graded_matrix(400L, 60L, -10, seed = 5)
  A <- t(A)
  S <- methods::as(Matrix::Matrix(A, sparse = TRUE), "CsparseMatrix")
  fit <- svd_partial(S, 28L, tol = 1e-8)
  truth <- t2_two_sided(A, fit$d, fit$u, fit$v)
  expect_true(fit$certificate$passed)
  expect_lte(max(truth$orth_u, truth$orth_v), sqrt(.Machine$double.eps))
  expect_lte(max(truth$backward), 1e-8)
  t2_expect_residuals_match(fit$certificate$orthogonality[["V"]], truth$orth_v,
                            64 * .Machine$double.eps)
})

test_that("C13b: Gram fast path and left-Gram path certify in original coordinates", {
  skip_if_not_installed("Matrix")
  set.seed(1304)
  for (wide in c(FALSE, TRUE)) {
    A <- t2_graded_matrix(300L, 40L, -2, seed = 1304)
    if (wide) A <- t(A)
    S <- methods::as(Matrix::Matrix(A, sparse = TRUE), "CsparseMatrix")
    fit <- svd_partial(S, 6L)
    expect_true(fit$certificate$passed)
    truth <- t2_two_sided(A, fit$d, fit$u, fit$v)
    expect_lte(max(truth$backward), 1e-8)
    floor <- 64 * .Machine$double.eps * sqrt(sum(A^2))
    t2_expect_residuals_match(fit$certificate$residuals$left, truth$left, floor)
    t2_expect_residuals_match(fit$certificate$residuals$right, truth$right, floor)
    expect_equal(fit$d, svd(A)$d[1:6], tolerance = 1e-10)
  }
})

test_that("C13c: cached-Av certificates recompute A v instead of trusting the cache", {
  skip_if_not_installed("Matrix")
  set.seed(1305)
  m <- 12L
  n <- 5L
  A <- matrix(rnorm(m * n), m, n)
  s <- svd(A)
  k <- 3L
  # Columns orthogonal to range(A): A^T z = 0. Adding them to U leaves
  # A^T U unchanged, so the right residual cannot see the corruption.
  Z <- qr.Q(qr(cbind(s$u, matrix(rnorm(m * k), m, k))))[, n + seq_len(k)]
  scale <- sqrt(1 + 0.3^2)
  u_bad <- (s$u[, 1:k] + 0.3 * Z) / scale
  d_bad <- s$d[1:k] / scale
  v <- s$v[, 1:k]
  stale_av <- sweep(u_bad, 2L, d_bad, `*`)
  truth <- t2_two_sided(A, d_bad, u_bad, v)
  expect_lt(max(truth$right), 1e-12)
  expect_gt(min(truth$left), 0.1)

  dense <- eigencore:::native_dense_svd_certificate_cached_av(A, d_bad, u_bad, v, stale_av)
  expect_equal(dense$left, truth$left, tolerance = 1e-10)
  expect_false(any(dense$converged))

  S <- Matrix::Matrix(A, sparse = TRUE)
  op <- eigencore:::as_operator(S)
  cert <- eigencore:::certify_svd_operator_cached_av(op, d_bad, u_bad, v, stale_av)
  expect_equal(cert$residuals$left, truth$left, tolerance = 1e-10)
  expect_false(cert$passed)
})

test_that("P12: block Golub-Kahan Ritz extraction matches dense SVD of AV", {
  set.seed(1306)
  n <- 30L
  m <- 45L
  p <- 9L
  A <- matrix(rnorm(m * n), m, n)
  V <- qr.Q(qr(matrix(rnorm(n * p), n, p)))
  AV <- A %*% V
  ref <- svd(AV)
  for (target in list(largest(), smallest())) {
    rk <- 4L
    ritz <- eigencore:::native_block_golub_kahan_ritz(V, AV, rank = rk, target = target,
                                                      active_cols = p)
    idx <- if (identical(target$kind, "largest")) seq_len(rk) else rev(seq_len(p))[seq_len(rk)]
    expect_equal(ritz$d, ref$d[idx], tolerance = 1e-12)
    signs <- sign(colSums(ritz$u * ref$u[, idx]))
    expect_equal(sweep(ritz$u, 2L, signs, `*`), ref$u[, idx], tolerance = 1e-10)
    expect_equal(sweep(ritz$v, 2L, signs, `*`), V %*% ref$v[, idx], tolerance = 1e-10)
    expect_equal(ritz$Avectors, A %*% ritz$v, tolerance = 1e-10)
    expect_equal(ritz$coefficients, sweep(ref$v[, idx], 2L, signs, `*`), tolerance = 1e-10)
  }
  # Wide active block (more columns than rows) and a rank-deficient AV.
  AV2 <- AV[1:6, ]
  ritz2 <- eigencore:::native_block_golub_kahan_ritz(V, AV2, rank = 6L, target = largest(),
                                                     active_cols = p)
  expect_equal(ritz2$d, svd(AV2)$d, tolerance = 1e-12)
  AV3 <- cbind(AV[, 1:4], AV[, 1:4] %*% matrix(rnorm(20), 4, 5))
  ritz3 <- eigencore:::native_block_golub_kahan_ritz(V, AV3, rank = 6L, target = largest(),
                                                     active_cols = p)
  expect_equal(ritz3$d, svd(AV3)$d[1:6], tolerance = 1e-10)
  expect_equal(crossprod(ritz3$u[, 1:4]), diag(4), tolerance = 1e-10)
})

test_that("P17: explicit Gram SVD kernels match dense SVD on sparse input", {
  skip_if_not_installed("Matrix")
  set.seed(1307)
  S <- Matrix::rsparsematrix(400L, 60L, density = 0.08)
  A <- as.matrix(S)
  ref <- svd(A)
  right <- .Call("eigencore_csc_right_gram_svd", S@i, S@p, S@x, S@Dim, 8L, 1e-8,
                 PACKAGE = "eigencore")
  expect_equal(right$d, ref$d[1:8], tolerance = 1e-10)
  expect_true(all(right$diagnostics$converged))
  St <- Matrix::t(S)
  left <- .Call("eigencore_csc_left_gram_svd", St@i, St@p, St@x, St@Dim, 8L, 1e-8,
                PACKAGE = "eigencore")
  expect_equal(left$d, ref$d[1:8], tolerance = 1e-10)
  expect_true(all(left$diagnostics$converged))
})

test_that("C28: CSC kernels reject malformed column pointers and row indices", {
  skip_if_not_installed("Matrix")
  S <- Matrix::rsparsematrix(20L, 8L, density = 0.3)
  i <- S@i
  p <- S@p
  x <- S@x
  dim <- S@Dim
  X <- matrix(1, 8L, 2L)
  Y <- matrix(0, 20L, 2L)
  call_apply <- function(i, p, x) {
    .Call("eigencore_csc_block_apply", i, p, x, dim, X, 1.0, 0.0, Y, FALSE,
          PACKAGE = "eigencore")
  }
  call_gram <- function(i, p, x) {
    .Call("eigencore_csc_right_gram_svd", i, p, x, dim, 2L, 1e-8, PACKAGE = "eigencore")
  }
  call_cert <- function(i, p, x) {
    .Call("eigencore_csc_svd_certificate", i, p, x, dim, c(1, 1),
          diag(1, 20L, 2L), diag(1, 8L, 2L), 1, 1e-8, PACKAGE = "eigencore")
  }
  for (fn in list(call_apply, call_gram, call_cert)) {
    expect_no_error(fn(i, p, x))
    bad_p_len <- p[-1L]
    expect_error(fn(i, bad_p_len, x), "CSC")
    bad_p_mono <- p
    bad_p_mono[3L] <- bad_p_mono[5L] + 1L
    expect_error(fn(i, bad_p_mono, x), "CSC")
    bad_p_end <- p
    bad_p_end[length(p)] <- length(x) + 5L
    expect_error(fn(i, bad_p_end, x), "CSC")
    bad_i <- i
    bad_i[1L] <- 20L
    expect_error(fn(bad_i, p, x), "CSC")
    bad_i[1L] <- -1L
    expect_error(fn(bad_i, p, x), "CSC")
  }
})
