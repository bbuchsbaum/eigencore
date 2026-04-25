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
  u_reorth_workspace <- basis_workspace(m, maxit, 1L)
  v_reorth_workspace <- basis_workspace(n, maxit, 1L)

  for (j in seq_len(maxit)) {
    iterations <- j
    V[, j] <- v

    u <- apply_operator(op, matrix(v, n, 1L))[, 1L] - beta_prev * u_prev
    nops <- nops + 1L
    if (isTRUE(reorthogonalize) && j > 1L) {
      Uprev <- U[, seq_len(j - 1L), drop = FALSE]
      u <- reorthogonalize_against(matrix(u, m, 1L), Uprev, passes = 2L, workspace = u_reorth_workspace)[, 1L]
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
      z <- reorthogonalize_against(matrix(z, n, 1L), Vj, passes = 2L, workspace = v_reorth_workspace)[, 1L]
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
native_golub_kahan_svd <- function(op, rank, target = largest(), tol = 1e-8,
                                   maxit = NULL,
                                   vectors = c("both", "left", "right", "none")) {
  vectors <- match.arg(vectors)
  op <- as_operator(op)
  if (is.null(op$apply_adjoint)) {
    stop("Native Golub-Kahan SVD requires an adjoint operator.", call. = FALSE)
  }

  m <- op$dim[1L]
  n <- op$dim[2L]
  limit <- min(m, n)
  fixed_maxit <- !is.null(maxit)
  if (is.null(maxit)) {
    maxit <- min(limit, max(20L, 4L * rank + 20L))
  } else {
    maxit <- min(limit, as.integer(maxit))
  }
  if (maxit < rank) {
    stop("maxit/max_subspace must be at least rank.", call. = FALSE)
  }

  start <- stats::rnorm(n)
  storage <- op$metadata$storage %||% NULL
  source <- source_or_null(op)
  run_native <- function(active_maxit) {
    if (identical(storage, "dgCMatrix")) {
      A <- op$metadata$matrix
      .Call(
        "eigencore_golub_kahan_csc",
        methods::slot(A, "i"),
        methods::slot(A, "p"),
        methods::slot(A, "x"),
        methods::slot(A, "Dim"),
        as.integer(active_maxit),
        as.numeric(start),
        PACKAGE = "eigencore"
      )
    } else if (is.matrix(source) && is.double(source)) {
      .Call(
        "eigencore_golub_kahan_dense",
        source,
        as.integer(active_maxit),
        as.numeric(start),
        PACKAGE = "eigencore"
      )
    } else {
      stop("Native Golub-Kahan currently supports dense double matrices and dgCMatrix operators only.", call. = FALSE)
    }
  }

  final <- NULL
  iter <- NULL
  active_maxit <- maxit
  retries <- 0L
  repeat {
    iter <- run_native(active_maxit)
    used <- seq_len(iter$iterations)
    final <- native_golub_kahan_ritz(
      op,
      iter$U[, used, drop = FALSE],
      iter$V[, used, drop = FALSE],
      iter$alpha[used],
      iter$beta[used],
      rank,
      target,
      tol
    )

    if (all(final$certificate$converged) || fixed_maxit || active_maxit >= limit) {
      break
    }
    retries <- retries + 1L
    active_maxit <- min(
      limit,
      max(active_maxit + max(10L, 2L * rank), as.integer(ceiling(1.5 * active_maxit)))
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
  final$iterations <- iter$iterations
  final$matvecs <- iter$matvecs
  final$restart <- list(
    kind = "adaptive_subspace_growth",
    implemented = TRUE,
    ritz_native = TRUE,
    retries = retries,
    final_max_subspace = active_maxit,
    fixed_max_subspace = fixed_maxit
  )
  final
}

#' @keywords internal
reference_randomized_svd <- function(op, rank, target = largest(), tol = 1e-8,
                                     oversample = 10L, n_iter = 2L,
                                     vectors = c("both", "left", "right", "none")) {
  vectors <- match.arg(vectors)
  op <- as_operator(op)
  if (is.null(op$apply_adjoint)) {
    stop("Randomized SVD requires an adjoint operator.", call. = FALSE)
  }
  m <- op$dim[1L]
  n <- op$dim[2L]
  limit <- min(m, n)
  rank <- min(as.integer(rank), limit)
  oversample <- max(0L, as.integer(oversample))
  n_iter <- max(0L, as.integer(n_iter))
  l <- min(limit, rank + oversample)

  Omega <- matrix(stats::rnorm(n * l), nrow = n, ncol = l)
  Y <- apply_operator(op, Omega)
  matvecs <- 1L
  Q <- qr.Q(qr(Y))
  for (iter in seq_len(n_iter)) {
    Z <- apply_adjoint_operator(op, Q)
    Y <- apply_operator(op, Z)
    matvecs <- matvecs + 2L
    Q <- qr.Q(qr(Y))
  }

  B_t <- apply_adjoint_operator(op, Q)
  matvecs <- matvecs + 1L
  core <- t(B_t)
  small <- svd(core, nu = min(nrow(core), rank), nv = min(ncol(core), rank))
  idx <- order_indices(small$d, target)
  idx <- idx[seq_len(min(rank, length(idx)))]
  d <- small$d[idx]
  u <- Q %*% small$u[, idx, drop = FALSE]
  v <- small$v[, idx, drop = FALSE]
  cert <- certify_svd_operator(op, d, u, v, tol = tol)

  if (vectors == "left") {
    v <- NULL
  } else if (vectors == "right") {
    u <- NULL
  } else if (vectors == "none") {
    u <- NULL
    v <- NULL
  }
  list(
    d = d,
    u = u,
    v = v,
    values = d,
    residuals = cert$residuals,
    backward_error = cert$backward_error,
    orthogonality = cert$orthogonality,
    certificate = cert,
    iterations = n_iter + 1L,
    matvecs = matvecs,
    restart = list(
      kind = "randomized_range_finder",
      implemented = TRUE,
      native = FALSE,
      oversample = oversample,
      n_iter = n_iter,
      sample_dimension = l
    )
  )
}

#' @keywords internal
reference_golub_kahan_ritz <- function(op, U, V, alpha, beta, rank, target, tol) {
  bd <- native_bidiagonal_svd(alpha, beta)
  idx <- order_indices(bd$d, target)
  idx <- idx[seq_len(min(rank, length(idx)))]
  d <- bd$d[idx]
  u <- U %*% bd$u[, idx, drop = FALSE]
  v <- V %*% bd$v[, idx, drop = FALSE]

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
native_bidiagonal_svd <- function(alpha, beta) {
  .Call(
    "eigencore_bidiagonal_svd",
    as.numeric(alpha),
    as.numeric(beta),
    PACKAGE = "eigencore"
  )
}

#' @keywords internal
native_svd_target_kind <- function(target) {
  kind <- if (inherits(target, "eigencore_target")) target$kind else "largest"
  switch(
    kind,
    largest = 1L,
    smallest = 2L,
    largest_magnitude = 3L,
    smallest_magnitude = 4L,
    stop(
      "Native Golub-Kahan SVD does not support target ", target_label(target),
      "; use a dense oracle fallback until shift-invert/refined interior SVD is implemented.",
      call. = FALSE
    )
  )
}

#' @keywords internal
native_golub_kahan_ritz <- function(op, U, V, alpha, beta, rank, target, tol) {
  ritz <- .Call(
    "eigencore_golub_kahan_ritz",
    U,
    V,
    as.numeric(alpha),
    as.numeric(beta),
    as.integer(rank),
    as.integer(native_svd_target_kind(target)),
    PACKAGE = "eigencore"
  )
  cert <- certify_svd_operator(op, ritz$d, ritz$u, ritz$v, tol = tol)
  list(
    d = ritz$d,
    u = ritz$u,
    v = ritz$v,
    values = ritz$d,
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
