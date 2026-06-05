# Shared scaffolding for R reference block solvers. These helpers keep
# validation, restart bookkeeping, and lock accounting consistent while the
# solver-specific files continue to own their recurrence and Ritz extraction.

#' @keywords internal
reference_block_positive_int <- function(x, name) {
  x <- as.integer(x)
  if (length(x) != 1L || is.na(x) || x < 1L) {
    stop(name, " must be a positive integer.", call. = FALSE)
  }
  x
}

#' @keywords internal
reference_block_nonnegative_int <- function(x, name) {
  x <- as.integer(x)
  if (length(x) != 1L || is.na(x) || x < 0L) {
    stop(name, " must be a non-negative integer.", call. = FALSE)
  }
  x
}

#' @keywords internal
reference_block_restart_controls <- function(requested,
                                             requested_name,
                                             limit,
                                             block,
                                             default_block,
                                             max_subspace,
                                             default_max_subspace,
                                             max_subspace_min = NULL,
                                             max_restarts = 100L,
                                             pad = NULL,
                                             clamp_requested = FALSE,
                                             clamp_block = FALSE) {
  requested <- reference_block_positive_int(requested, requested_name)
  limit <- as.integer(limit)
  if (length(limit) != 1L || is.na(limit) || limit < 1L) {
    stop("limit must be a positive integer.", call. = FALSE)
  }
  if (isTRUE(clamp_requested)) {
    requested <- min(requested, limit)
  }

  if (is.function(default_block)) {
    default_block <- default_block(requested)
  }
  block <- if (is.null(block)) default_block else as.integer(block)
  block <- reference_block_positive_int(block, "block")
  if (isTRUE(clamp_block)) {
    block <- min(block, limit)
  }

  max_restarts <- reference_block_nonnegative_int(max_restarts, "max_restarts")
  if (is.function(default_max_subspace)) {
    default_max_subspace <- default_max_subspace(requested, block)
  }
  if (is.null(max_subspace_min)) {
    max_subspace_min <- requested
  } else if (is.function(max_subspace_min)) {
    max_subspace_min <- max_subspace_min(requested, block)
  }
  max_subspace <- if (is.null(max_subspace)) {
    default_max_subspace
  } else {
    as.integer(max_subspace)
  }
  max_subspace <- min(limit, reference_block_positive_int(max_subspace, "max_subspace"))
  if (max_subspace < max_subspace_min) {
    stop("max_subspace must be at least ", requested_name, ".", call. = FALSE)
  }

  pad <- if (is.null(pad)) {
    min(requested, max(block, 5L))
  } else {
    reference_block_nonnegative_int(pad, "pad")
  }

  list(
    requested = requested,
    block = block,
    max_subspace = max_subspace,
    max_restarts = max_restarts,
    pad = pad
  )
}

#' @keywords internal
reference_scalar_subspace_controls <- function(requested,
                                               requested_name,
                                               limit,
                                               maxit,
                                               default_maxit,
                                               maxit_name = "maxit") {
  requested <- reference_block_positive_int(requested, requested_name)
  limit <- as.integer(limit)
  if (length(limit) != 1L || is.na(limit) || limit < 1L) {
    stop("limit must be a positive integer.", call. = FALSE)
  }

  if (is.function(default_maxit)) {
    default_maxit <- default_maxit(requested)
  }
  maxit <- if (is.null(maxit)) {
    default_maxit
  } else {
    maxit
  }
  maxit <- reference_block_positive_int(maxit, maxit_name)
  maxit <- min(limit, maxit)
  if (maxit < requested) {
    stop(maxit_name, "/max_subspace must be at least ", requested_name, ".",
         call. = FALSE)
  }

  list(
    requested = requested,
    maxit = maxit
  )
}

#' @keywords internal
reference_block_accept <- function(X, Q_lock, V, max_cols,
                                   tol = sqrt(.Machine$double.eps)) {
  if (max_cols <= 0L || ncol(X) == 0L) {
    return(list(Q = matrix(0, nrow(X), 0L), ortho_passes = 0L))
  }
  against <- cbind(Q_lock, V)
  Y <- if (ncol(against) > 0L) {
    reorthogonalize_against(X, against, passes = 2L)
  } else {
    X
  }
  decomp <- native_mgs2(Y, tol = tol)
  Q <- decomp$Q
  if (ncol(Q) > max_cols) {
    Q <- Q[, seq_len(max_cols), drop = FALSE]
  }
  list(Q = Q, ortho_passes = 2L)
}

#' @keywords internal
reference_block_initial_basis <- function(n, block, Q_lock, max_subspace) {
  accepted <- reference_block_accept(
    matrix(stats::rnorm(n * block), nrow = n, ncol = block),
    Q_lock,
    matrix(0, n, 0L),
    max_cols = max_subspace
  )
  Q <- accepted$Q
  if (ncol(Q) == 0L) {
    Q <- diag(1, n, min(block, n))
  }
  list(Q = Q, ortho_passes = accepted$ortho_passes)
}

#' @keywords internal
reference_block_lock_count <- function(converged, remaining) {
  lock_now <- 0L
  for (i in seq_len(min(remaining, length(converged)))) {
    if (!isTRUE(converged[[i]])) {
      break
    }
    lock_now <- lock_now + 1L
  }
  lock_now
}

#' @keywords internal
reference_block_unlocked_indices <- function(lock_now, n_values) {
  if (lock_now < n_values) {
    seq.int(lock_now + 1L, n_values)
  } else {
    integer()
  }
}

#' @keywords internal
reference_block_take_count <- function(n_values, remaining, pad) {
  min(n_values, max(remaining + pad, remaining))
}

#' @keywords internal
reference_block_keep_count <- function(n_unlocked, remaining, pad, max_subspace,
                                       block) {
  min(n_unlocked, max(remaining + pad, remaining), max_subspace - block)
}

#' @keywords internal
reference_block_history_frame <- function(restart, iterations, n_locked, cert) {
  data.frame(
    restart = restart,
    iteration = iterations,
    n_locked = n_locked,
    max_residual = cert$max_residual,
    max_backward_error = cert$max_backward_error %||% NA_real_
  )
}
