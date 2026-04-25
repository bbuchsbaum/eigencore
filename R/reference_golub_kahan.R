#' @keywords internal
reference_golub_kahan_svd <- function(op, rank, target = largest(), tol = 1e-8,
                                      maxit = NULL, vectors = c("both", "left", "right", "none"),
                                      reorthogonalize = TRUE) {
  vectors <- match.arg(vectors)
  op <- as_operator(op)
  if (is.null(op$apply_adjoint)) {
    stop("Golub-Kahan SVD requires an adjoint operator.", call. = FALSE)
  }

  m <- op$dim[1L]
  n <- op$dim[2L]
  limit <- min(m, n)
  if (is.null(maxit)) {
    maxit <- min(limit, max(20L, 4L * rank + 20L))
  } else {
    maxit <- min(limit, as.integer(maxit))
  }
  if (maxit < rank) {
    stop("maxit/max_subspace must be at least rank.", call. = FALSE)
  }

  U <- matrix(0, m, maxit)
  V <- matrix(0, n, maxit)
  alpha <- numeric(maxit)
  beta <- numeric(maxit)

  v <- stats::rnorm(n)
  v_norm <- sqrt(sum(v^2))
  if (v_norm == 0) {
    v[1L] <- 1
    v_norm <- 1
  }
  v <- v / v_norm
  u_prev <- numeric(m)
  beta_prev <- 0

  nops <- 0L
  final <- NULL
  iterations <- 0L

  for (j in seq_len(maxit)) {
    iterations <- j
    V[, j] <- v

    u <- apply_operator(op, matrix(v, n, 1L))[, 1L] - beta_prev * u_prev
    nops <- nops + 1L
    if (isTRUE(reorthogonalize) && j > 1L) {
      Uprev <- U[, seq_len(j - 1L), drop = FALSE]
      u <- u - Uprev %*% crossprod(Uprev, u)
      u <- u - Uprev %*% crossprod(Uprev, u)
    }
    alpha[[j]] <- sqrt(sum(u^2))
    if (alpha[[j]] <= max(100 * .Machine$double.eps, tol * 1e-3)) {
      break
    }
    u <- u / alpha[[j]]
    U[, j] <- u

    z <- apply_adjoint_operator(op, matrix(u, m, 1L))[, 1L] - alpha[[j]] * v
    nops <- nops + 1L
    if (isTRUE(reorthogonalize)) {
      Vj <- V[, seq_len(j), drop = FALSE]
      z <- z - Vj %*% crossprod(Vj, z)
      z <- z - Vj %*% crossprod(Vj, z)
    }
    beta[[j]] <- sqrt(sum(z^2))

    if (j >= rank) {
      final <- reference_golub_kahan_ritz(
        op,
        U[, seq_len(j), drop = FALSE],
        V[, seq_len(j), drop = FALSE],
        alpha[seq_len(j)],
        beta[seq_len(j)],
        rank,
        target,
        tol
      )
      if (all(final$certificate$converged)) {
        break
      }
    }
    if (beta[[j]] <= max(100 * .Machine$double.eps, tol * 1e-3)) {
      break
    }

    u_prev <- u
    beta_prev <- beta[[j]]
    v <- z / beta[[j]]
  }

  if (is.null(final)) {
    final <- reference_golub_kahan_ritz(
      op,
      U[, seq_len(iterations), drop = FALSE],
      V[, seq_len(iterations), drop = FALSE],
      alpha[seq_len(iterations)],
      beta[seq_len(iterations)],
      rank,
      target,
      tol
    )
  }

  if (vectors == "left") {
    final$v <- NULL
  } else if (vectors == "right") {
    final$u <- NULL
  } else if (vectors == "none") {
    final$u <- NULL
    final$v <- NULL
  }
  final$iterations <- iterations
  final$matvecs <- nops
  final
}

#' @keywords internal
reference_golub_kahan_ritz <- function(op, U, V, alpha, beta, rank, target, tol) {
  B <- bidiagonal_matrix(alpha, beta)
  bd <- svd(B, nu = nrow(B), nv = ncol(B))
  idx <- order_indices(bd$d, target)
  idx <- idx[seq_len(min(rank, length(idx)))]
  d <- bd$d[idx]
  u <- U %*% bd$u[, idx, drop = FALSE]
  v <- V %*% bd$v[, idx, drop = FALSE]

  # Stabilize returned vectors before certification.
  u <- mgs2(u)$Q
  v <- mgs2(v)$Q
  AV <- apply_operator(op, v)
  d <- pmax(0, as.numeric(diag(crossprod(u, AV))))
  idx2 <- order_indices(d, target)
  d <- d[idx2]
  u <- u[, idx2, drop = FALSE]
  v <- v[, idx2, drop = FALSE]

  cert <- certify_svd_operator(op, d, u, v, tol = tol)
  list(
    d = d,
    u = u,
    v = v,
    values = d,
    residuals = cert$residuals,
    backward_error = cert$backward_error,
    orthogonality = cert$orthogonality,
    certificate = cert
  )
}

#' @keywords internal
bidiagonal_matrix <- function(alpha, beta) {
  p <- length(alpha)
  B <- matrix(0, p, p)
  diag(B) <- alpha
  if (p > 1L) {
    off <- beta[seq_len(p - 1L)]
    B[cbind(seq_len(p - 1L), 2:p)] <- off
  }
  B
}
