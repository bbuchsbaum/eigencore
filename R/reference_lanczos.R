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
  reorth_workspace <- basis_workspace(n, maxit, 1L)

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
      z <- reorthogonalize_against(matrix(z, n, 1L), Qj, passes = 2L, workspace = reorth_workspace)[, 1L]
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
native_lanczos_hermitian <- function(op, k, target = largest(), tol = 1e-8,
                                     maxit = NULL, max_restarts = NULL,
                                     vectors = TRUE) {
  op <- as_operator(op)
  if (op$dim[1L] != op$dim[2L]) {
    stop("Native Hermitian Lanczos requires a square operator.", call. = FALSE)
  }
  if (!identical(op$structure$kind, "hermitian")) {
    stop("Native Hermitian Lanczos requires an operator with hermitian() structure.", call. = FALSE)
  }

  n <- op$dim[1L]
  m_max <- if (is.null(maxit)) {
    min(n, max(as.integer(k) + 1L, 3L * as.integer(k) + 20L))
  } else {
    min(n, as.integer(maxit))
  }
  if (m_max < as.integer(k) + 1L) {
    stop("max_subspace must be at least k + 1.", call. = FALSE)
  }
  if (is.null(max_restarts)) {
    max_restarts <- 100L
  } else {
    max_restarts <- as.integer(max_restarts)
    if (max_restarts < 0L) {
      stop("max_restarts must be non-negative.", call. = FALSE)
    }
  }

  start <- stats::rnorm(n)
  storage <- op$metadata$storage %||% NULL
  source <- source_or_null(op)
  iter <- if (identical(storage, "dgCMatrix")) {
    A <- op$metadata$matrix
    .Call(
      "eigencore_thick_restart_lanczos_csc",
      methods::slot(A, "i"),
      methods::slot(A, "p"),
      methods::slot(A, "x"),
      methods::slot(A, "Dim"),
      as.integer(k),
      as.integer(m_max),
      as.integer(lanczos_target_kind(target)),
      as.numeric(tol),
      as.integer(max_restarts),
      as.numeric(start),
      PACKAGE = "eigencore"
    )
  } else if (is.matrix(source) && is.double(source)) {
    .Call(
      "eigencore_thick_restart_lanczos_dense",
      source,
      as.integer(k),
      as.integer(m_max),
      as.integer(lanczos_target_kind(target)),
      as.numeric(tol),
      as.integer(max_restarts),
      as.numeric(start),
      PACKAGE = "eigencore"
    )
  } else {
    stop("Native Hermitian Lanczos currently supports dense double matrices and dgCMatrix operators only.", call. = FALSE)
  }

  values <- iter$values
  vec_matrix <- iter$vectors
  cert <- certify_eigen_operator_residuals(op, values, vec_matrix, iter$residuals, tol = tol)
  if (!isTRUE(vectors)) {
    vec_matrix <- NULL
  }

  list(
    values = values,
    vectors = vec_matrix,
    residuals = cert$residuals,
    backward_error = cert$backward_error,
    orthogonality = cert$orthogonality,
    certificate = cert,
    iterations = iter$iterations,
    matvecs = iter$matvecs,
    convergence_history = data.frame(
      restart = seq_len(iter$restarts + 1L) - 1L,
      iteration = c(rep(NA_integer_, iter$restarts), iter$iterations),
      n_locked = c(rep(NA_integer_, iter$restarts), iter$n_locked)
    ),
    locked = seq_len(iter$n_locked),
    restart = list(
      kind = "thick_restart",
      implemented = TRUE,
      locking = "in_native_loop",
      locked = seq_len(iter$n_locked),
      locked_count = iter$n_locked,
      restarts_used = iter$restarts,
      max_restarts = max_restarts,
      max_subspace = m_max,
      final_active_subspace = iter$m_active_final %||% NA_integer_
    )
  )
}

#' @keywords internal
native_block_lanczos_hermitian <- function(op, k, target = largest(), tol = 1e-8,
                                           maxit = NULL, block = 4L,
                                           vectors = TRUE) {
  op <- as_operator(op)
  if (op$dim[1L] != op$dim[2L]) {
    stop("Native block Hermitian Lanczos requires a square operator.", call. = FALSE)
  }
  if (!identical(op$structure$kind, "hermitian")) {
    stop("Native block Hermitian Lanczos requires an operator with hermitian() structure.", call. = FALSE)
  }

  n <- op$dim[1L]
  block <- as.integer(block)
  if (length(block) != 1L || is.na(block) || block < 2L) {
    stop("native block Hermitian Lanczos requires block >= 2.", call. = FALSE)
  }
  m_max <- if (is.null(maxit)) {
    min(n, max(as.integer(k) + block, 8L * as.integer(k) + 4L * block))
  } else {
    min(n, as.integer(maxit))
  }
  if (m_max < as.integer(k)) {
    stop("max_subspace must be at least k.", call. = FALSE)
  }

  start <- matrix(stats::rnorm(n * block), nrow = n, ncol = block)
  storage <- op$metadata$storage %||% NULL
  source <- source_or_null(op)
  iter <- if (identical(storage, "dgCMatrix")) {
    A <- op$metadata$matrix
    .Call(
      "eigencore_block_lanczos_csc",
      methods::slot(A, "i"),
      methods::slot(A, "p"),
      methods::slot(A, "x"),
      methods::slot(A, "Dim"),
      as.integer(k),
      as.integer(m_max),
      as.integer(block),
      as.integer(lanczos_target_kind(target)),
      as.numeric(tol),
      start,
      PACKAGE = "eigencore"
    )
  } else if (is.matrix(source) && is.double(source)) {
    .Call(
      "eigencore_block_lanczos_dense",
      source,
      as.integer(k),
      as.integer(m_max),
      as.integer(block),
      as.integer(lanczos_target_kind(target)),
      as.numeric(tol),
      start,
      PACKAGE = "eigencore"
    )
  } else {
    stop("Native block Hermitian Lanczos currently supports dense double matrices and dgCMatrix operators only.", call. = FALSE)
  }

  values <- iter$values
  vec_matrix <- iter$vectors
  cert <- certify_eigen_operator_residuals(op, values, vec_matrix, iter$residuals, tol = tol)
  if (!isTRUE(vectors)) {
    vec_matrix <- NULL
  }

  list(
    values = values,
    vectors = vec_matrix,
    residuals = cert$residuals,
    backward_error = cert$backward_error,
    orthogonality = cert$orthogonality,
    certificate = cert,
    iterations = iter$iterations,
    matvecs = iter$matvecs,
    convergence_history = data.frame(
      restart = 0L,
      iteration = iter$iterations,
      n_locked = sum(iter$converged)
    ),
    locked = which(iter$converged),
    restart = list(
      kind = "block_krylov_rayleigh_ritz",
      implemented = TRUE,
      locking = "none",
      locked = which(iter$converged),
      locked_count = sum(iter$converged),
      restarts_used = 0L,
      max_restarts = 0L,
      max_subspace = m_max,
      final_active_subspace = iter$m_active_final %||% NA_integer_,
      block = block
    )
  )
}

#' @keywords internal
lanczos_target_kind <- function(target) {
  kind <- if (inherits(target, "eigencore_target")) target$kind else "largest"
  switch(
    kind,
    largest = 1L,
    smallest = 2L,
    largest_magnitude = 3L,
    smallest_magnitude = 4L,
    stop(
      "Native Hermitian Lanczos does not support target ", target_label(target),
      "; use a dense oracle fallback or shift_invert() when that path is implemented.",
      call. = FALSE
    )
  )
}

#' @keywords internal
native_lanczos_target_supported <- function(target) {
  kind <- if (inherits(target, "eigencore_target")) target$kind else "largest"
  kind %in% c("largest", "smallest", "largest_magnitude", "smallest_magnitude")
}

#' @keywords internal
reference_lanczos_ritz <- function(op, Q, alpha, beta, k, target, tol) {
  eig <- native_tridiagonal_eigen(alpha, beta)
  idx <- order_indices(eig$values, target)
  idx <- idx[seq_len(min(k, length(idx)))]
  values <- eig$values[idx]
  vectors <- Q %*% eig$vectors[, idx, drop = FALSE]
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
native_tridiagonal_eigen <- function(alpha, beta) {
  .Call(
    "eigencore_tridiagonal_eigen",
    as.numeric(alpha),
    as.numeric(beta),
    PACKAGE = "eigencore"
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
