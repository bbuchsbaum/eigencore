#' @keywords internal
reference_lanczos_hermitian <- function(op, k, target = largest(), tol = 1e-8,
                                        maxit = NULL, vectors = TRUE,
                                        reorthogonalize = TRUE) {
  op <- as_operator(op)
  if (op$dim[1L] != op$dim[2L]) {
    stop("Hermitian Lanczos requires a square operator.", call. = FALSE)
  }
  if (!identical(op$structure$kind, "hermitian")) {
    stop("Hermitian Lanczos requires an operator with hermitian() structure.", call. = FALSE)
  }

  n <- op$dim[1L]
  if (is.null(maxit)) {
    maxit <- min(n, max(20L, 4L * k + 20L))
  } else {
    maxit <- min(n, as.integer(maxit))
  }
  if (maxit < k) {
    stop("maxit/max_subspace must be at least k.", call. = FALSE)
  }

  Q <- matrix(0, n, maxit)
  alpha <- numeric(maxit)
  beta <- numeric(maxit)
  q_prev <- numeric(n)
  q <- stats::rnorm(n)
  q_norm <- sqrt(sum(q^2))
  if (q_norm == 0) {
    q[1L] <- 1
    q_norm <- 1
  }
  q <- q / q_norm

  nops <- 0L
  final <- NULL
  iterations <- 0L

  for (j in seq_len(maxit)) {
    iterations <- j
    Q[, j] <- q
    z <- apply_operator(op, matrix(q, n, 1L))[, 1L]
    nops <- nops + 1L

    if (j > 1L) {
      z <- z - beta[[j - 1L]] * q_prev
    }
    alpha[[j]] <- sum(q * z)
    z <- z - alpha[[j]] * q

    if (isTRUE(reorthogonalize)) {
      Qj <- Q[, seq_len(j), drop = FALSE]
      z <- z - Qj %*% crossprod(Qj, z)
      z <- z - Qj %*% crossprod(Qj, z)
    }

    beta[[j]] <- sqrt(sum(z^2))
    if (j >= k) {
      final <- reference_lanczos_ritz(op, Q[, seq_len(j), drop = FALSE], alpha[seq_len(j)],
                                      beta[seq_len(j)], k, target, tol)
      if (all(final$certificate$converged)) {
        break
      }
    }
    if (beta[[j]] <= max(100 * .Machine$double.eps, tol * 1e-3)) {
      break
    }

    q_prev <- q
    q <- z / beta[[j]]
  }

  if (is.null(final)) {
    final <- reference_lanczos_ritz(op, Q[, seq_len(iterations), drop = FALSE],
                                    alpha[seq_len(iterations)], beta[seq_len(iterations)],
                                    k, target, tol)
  }
  if (!isTRUE(vectors)) {
    final$vectors <- NULL
  }
  final$iterations <- iterations
  final$matvecs <- nops
  final
}

#' @keywords internal
reference_lanczos_ritz <- function(op, Q, alpha, beta, k, target, tol) {
  T <- tridiagonal_matrix(alpha, beta)
  eig <- eigen(T, symmetric = TRUE)
  idx <- order_indices(eig$values, target)
  idx <- idx[seq_len(min(k, length(idx)))]
  values <- eig$values[idx]
  vectors <- Q %*% eig$vectors[, idx, drop = FALSE]
  vectors <- mgs2(vectors)$Q
  # Recompute Rayleigh quotients after orthogonalizing Ritz vectors.
  AV <- apply_operator(op, vectors)
  values <- as.numeric(diag(crossprod(vectors, AV)))
  idx2 <- order_indices(values, target)
  values <- values[idx2]
  vectors <- vectors[, idx2, drop = FALSE]
  cert <- certify_eigen_operator(op, values, vectors, tol = tol)

  list(
    values = values,
    vectors = vectors,
    residuals = cert$residuals,
    backward_error = cert$backward_error,
    orthogonality = cert$orthogonality,
    certificate = cert
  )
}

#' @keywords internal
tridiagonal_matrix <- function(alpha, beta) {
  m <- length(alpha)
  T <- diag(alpha, m)
  if (m > 1L) {
    off <- beta[seq_len(m - 1L)]
    T[cbind(seq_len(m - 1L), 2:m)] <- off
    T[cbind(2:m, seq_len(m - 1L))] <- off
  }
  T
}
