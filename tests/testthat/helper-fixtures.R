orthogonal_fixture <- function(n, seed = 1L) {
  set.seed(seed)
  qr.Q(qr(matrix(rnorm(n * n), nrow = n)))
}

symmetric_with_spectrum <- function(values, seed = 1L) {
  Q <- orthogonal_fixture(length(values), seed = seed)
  Q %*% diag(values, nrow = length(values)) %*% t(Q)
}

rectangular_with_singular_values <- function(values, m, n, seed = 1L) {
  r <- length(values)
  U <- qr.Q(qr(matrix(rnorm(m * m), nrow = m)))[, seq_len(r), drop = FALSE]
  V <- qr.Q(qr(matrix(rnorm(n * n), nrow = n)))[, seq_len(r), drop = FALSE]
  U %*% diag(values, nrow = r) %*% t(V)
}

subspace_distance <- function(X, Y) {
  X <- qr.Q(qr(as.matrix(X)))
  Y <- qr.Q(qr(as.matrix(Y)))
  s <- svd(crossprod(X, Y), nu = 0, nv = 0)$d
  sqrt(max(0, 1 - min(pmin(s, 1))^2))
}

expect_certificate_clean <- function(fit, tol = 1e-8) {
  cert <- certificate(fit)
  expect_true(cert$passed)
  expect_true(all(is.finite(unlist(cert$residuals, use.names = FALSE))))
  expect_true(all(is.finite(cert$backward_error)))
  expect_lte(cert$max_backward_error, tol)
}

#' @noRd
spectrum_pattern <- function(pattern = c("uniform", "clustered", "exponential",
                                         "geometric", "two_cluster"),
                             n) {
  pattern <- match.arg(pattern)
  switch(
    pattern,
    uniform     = seq.int(1, n, length.out = n),
    clustered   = c(rep(1, n - 3L) + seq.int(0, by = 1e-3, length.out = n - 3L),
                    n / 2, n / 2 + 1, n),
    exponential = exp(seq.int(-1, -n, length.out = n)),
    geometric   = 0.7 ^ (seq.int(0, n - 1L)),
    two_cluster = c(seq.int(1, 1 + 1e-3, length.out = floor(n / 2)),
                    seq.int(n - 1, n, length.out = n - floor(n / 2)))
  )
}

#' @noRd
random_symmetric_with_spectrum <- function(n, pattern = "uniform", seed = 1L) {
  vals <- spectrum_pattern(pattern, n)
  symmetric_with_spectrum(vals, seed = seed)
}

#' @noRd
random_rectangular_with_pattern <- function(m, n, rank = NULL,
                                            pattern = "uniform", seed = 1L) {
  rank <- rank %||% min(m, n)
  vals <- spectrum_pattern(pattern, rank)
  rectangular_with_singular_values(vals, m, n, seed = seed)
}
