#' @keywords internal
#' @details
#' Despite the historical name `reference_block_lanczos_*`, this routine
#' implements **block-Krylov subspace iteration with reorthogonalization
#' and thick restart**, not the three-term block Lanczos recurrence. New
#' Krylov columns are produced from `AV[, last_block]` directly and then
#' reorthogonalized against the current basis (lines 73-89), without the
#' explicit alpha/beta block subtractions that distinguish a Lanczos
#' recurrence. Use this only as the block-subspace-iteration oracle; do
#' not rely on it as a parity reference for the native block Lanczos
#' three-term path. The honest name is exposed below as
#' `reference_block_subspace_iteration_thick_restart_hermitian`.
reference_block_subspace_iteration_thick_restart_hermitian <- function(op, k, target = largest(),
                                                            tol = 1e-8,
                                                            max_subspace = NULL,
                                                            max_restarts = 100L,
                                                            block = NULL,
                                                            pad = NULL,
                                                            vectors = TRUE) {
  op <- as_operator(op)
  if (op$dim[1L] != op$dim[2L]) {
    stop("Block Hermitian Lanczos requires a square operator.", call. = FALSE)
  }
  if (!identical(op$structure$kind, "hermitian")) {
    stop("Block Hermitian Lanczos requires an operator with hermitian() structure.", call. = FALSE)
  }

  n <- op$dim[1L]
  k <- as.integer(k)
  if (length(k) != 1L || is.na(k) || k < 1L) {
    stop("k must be a positive integer.", call. = FALSE)
  }
  block <- if (is.null(block)) {
    2L
  } else {
    as.integer(block)
  }
  if (length(block) != 1L || is.na(block) || block < 1L) {
    stop("block must be a positive integer.", call. = FALSE)
  }
  max_restarts <- as.integer(max_restarts)
  if (length(max_restarts) != 1L || is.na(max_restarts) || max_restarts < 0L) {
    stop("max_restarts must be a non-negative integer.", call. = FALSE)
  }
  if (is.null(max_subspace)) {
    max_subspace <- default_block_lanczos_max_subspace(k, block)
  }
  max_subspace <- min(n, as.integer(max_subspace))
  if (length(max_subspace) != 1L || is.na(max_subspace) || max_subspace < k) {
    stop("max_subspace must be at least k.", call. = FALSE)
  }
  pad <- if (is.null(pad)) min(k, max(block, 5L)) else as.integer(pad)
  if (length(pad) != 1L || is.na(pad) || pad < 0L) {
    stop("pad must be a non-negative integer.", call. = FALSE)
  }

  q_lock <- matrix(0, n, 0L)
  locked_values <- numeric()
  locked_vectors <- matrix(0, n, 0L)
  locked_residuals <- numeric()
  iterations <- 0L
  matvecs <- 0L
  ortho_passes <- 0L
  locking_events <- 0L
  history <- list()

  accepted <- reference_block_accept(
    matrix(stats::rnorm(n * block), nrow = n, ncol = block),
    q_lock,
    matrix(0, n, 0L),
    max_cols = max_subspace
  )
  V <- accepted$Q
  ortho_passes <- ortho_passes + accepted$ortho_passes
  if (ncol(V) == 0L) {
    V <- diag(1, n, min(block, n))
  }
  AV <- apply_operator(op, V)
  matvecs <- matvecs + 1L

  last_candidates <- NULL
  restarts_used <- 0L
  for (restart in 0:max_restarts) {
    while (ncol(V) < max_subspace) {
      source_cols <- seq.int(max(1L, ncol(V) - block + 1L), ncol(V))
      accepted <- reference_block_accept(
        AV[, source_cols, drop = FALSE],
        q_lock,
        V,
        max_cols = max_subspace - ncol(V)
      )
      ortho_passes <- ortho_passes + accepted$ortho_passes
      if (ncol(accepted$Q) == 0L) {
        break
      }
      V <- cbind(V, accepted$Q)
      AV <- cbind(AV, apply_operator(op, accepted$Q))
      matvecs <- matvecs + 1L
      iterations <- iterations + 1L
    }

    rr <- reference_block_lanczos_rr(op, V, AV, target = target, tol = tol)
    remaining <- k - length(locked_values)
    take <- min(length(rr$values), max(remaining + pad, remaining))
    if (take < 1L) {
      break
    }
    rr <- reference_block_lanczos_slice(rr, seq_len(take))

    lock_now <- 0L
    for (i in seq_len(min(remaining, length(rr$values)))) {
      if (!isTRUE(rr$certificate$converged[[i]])) {
        break
      }
      lock_now <- lock_now + 1L
    }
    if (lock_now > 0L) {
      locked_values <- c(locked_values, rr$values[seq_len(lock_now)])
      locked_vectors <- cbind(locked_vectors, rr$vectors[, seq_len(lock_now), drop = FALSE])
      locked_residuals <- c(locked_residuals, rr$residuals[seq_len(lock_now)])
      q_lock <- locked_vectors
      locking_events <- locking_events + 1L
    }

    history[[length(history) + 1L]] <- data.frame(
      restart = restart,
      iteration = iterations,
      n_locked = length(locked_values),
      max_residual = rr$certificate$max_residual
    )

    unlocked_idx <- if (lock_now < length(rr$values)) {
      seq.int(lock_now + 1L, length(rr$values))
    } else {
      integer()
    }
    last_candidates <- reference_block_lanczos_slice(rr, unlocked_idx)
    if (length(locked_values) >= k || restart == max_restarts) {
      restarts_used <- restart
      break
    }

    k_remaining <- k - length(locked_values)
    k_keep <- min(length(unlocked_idx), max(k_remaining + pad, k_remaining), max_subspace - block)
    if (k_keep > 0L) {
      keep_idx <- seq_len(k_keep)
      V <- last_candidates$vectors[, keep_idx, drop = FALSE]
      AV <- last_candidates$Avectors[, keep_idx, drop = FALSE]
    } else {
      V <- matrix(0, n, 0L)
      AV <- matrix(0, n, 0L)
    }

    tail_count <- min(block, length(unlocked_idx))
    tail <- if (tail_count > 0L) {
      last_candidates$residual_vectors[, seq_len(tail_count), drop = FALSE]
    } else {
      matrix(stats::rnorm(n * block), nrow = n, ncol = block)
    }
    accepted <- reference_block_accept(tail, q_lock, V, max_cols = max_subspace - ncol(V))
    ortho_passes <- ortho_passes + accepted$ortho_passes
    if (ncol(accepted$Q) == 0L) {
      accepted <- reference_block_accept(
        matrix(stats::rnorm(n * block), nrow = n, ncol = block),
        q_lock,
        V,
        max_cols = max_subspace - ncol(V)
      )
      ortho_passes <- ortho_passes + accepted$ortho_passes
    }
    if (ncol(accepted$Q) == 0L) {
      restarts_used <- restart
      break
    }
    V <- cbind(V, accepted$Q)
    AV <- cbind(AV, apply_operator(op, accepted$Q))
    matvecs <- matvecs + 1L
    restarts_used <- restart + 1L
  }

  values <- locked_values
  vecs <- locked_vectors
  residuals <- locked_residuals
  if (length(values) < k && !is.null(last_candidates) && length(last_candidates$values)) {
    need <- min(k - length(values), length(last_candidates$values))
    values <- c(values, last_candidates$values[seq_len(need)])
    vecs <- cbind(vecs, last_candidates$vectors[, seq_len(need), drop = FALSE])
    residuals <- c(residuals, last_candidates$residuals[seq_len(need)])
  }
  if (length(values) > k) {
    values <- values[seq_len(k)]
    vecs <- vecs[, seq_len(k), drop = FALSE]
    residuals <- residuals[seq_len(k)]
  }

  # Recompute residuals against the freshly-cbinded vecs after pad/lock.
  # Reusing residuals from earlier basis snapshots can declare convergence on
  # pairs whose residuals no longer correspond to the final returned vectors.
  cert <- certify_eigen_operator(op, values, vecs, tol = tol)
  if (!isTRUE(vectors)) {
    vecs <- NULL
  }

  list(
    values = values,
    vectors = vecs,
    residuals = cert$residuals,
    backward_error = cert$backward_error,
    orthogonality = cert$orthogonality,
    certificate = cert,
    iterations = iterations,
    matvecs = matvecs,
    restarts = restarts_used,
    ortho_passes = ortho_passes,
    locking_events = locking_events,
    block = block,
    convergence_history = if (length(history)) do.call(rbind, history) else data.frame(),
    locked = seq_len(length(locked_values)),
    restart = list(
      kind = "block_subspace_iteration_thick_restart_reference",
      implemented = TRUE,
      locking = "reference_loop",
      locked_count = length(locked_values),
      restarts_used = restarts_used,
      max_restarts = max_restarts,
      max_subspace = max_subspace,
      final_active_subspace = ncol(V),
      block = block,
      ortho_passes = ortho_passes,
      locking_events = locking_events
    )
  )
}

#' @keywords internal
#' Deprecated alias kept until external test sites migrate. The implementation
#' is block-Krylov subspace iteration with thick restart, not the three-term
#' Lanczos recurrence; new callers must use the honestly-named function.
reference_block_lanczos_thick_restart_hermitian <-
  reference_block_subspace_iteration_thick_restart_hermitian

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
reference_block_lanczos_rr <- function(op, V, AV, target, tol) {
  H <- crossprod(V, AV)
  H <- (H + t(H)) / 2
  eig <- eigen(H, symmetric = TRUE)
  idx <- order_indices(eig$values, target)
  S <- eig$vectors[, idx, drop = FALSE]
  values <- eig$values[idx]
  vectors <- V %*% S
  Avectors <- AV %*% S
  residual_vectors <- Avectors - sweep(vectors, 2L, values, `*`)
  residuals <- col_norms(residual_vectors)
  cert <- certify_eigen_operator_residuals(op, values, vectors, residuals, tol = tol)
  list(
    values = values,
    vectors = vectors,
    Avectors = Avectors,
    residual_vectors = residual_vectors,
    residuals = residuals,
    certificate = cert
  )
}

#' @keywords internal
reference_block_lanczos_slice <- function(rr, idx) {
  if (!length(idx)) {
    n <- if (is.null(rr$vectors)) 0L else nrow(rr$vectors)
    return(list(
      values = numeric(),
      vectors = matrix(0, n, 0L),
      Avectors = matrix(0, n, 0L),
      residual_vectors = matrix(0, n, 0L),
      residuals = numeric(),
      certificate = rr$certificate
    ))
  }
  rr$values <- rr$values[idx]
  rr$vectors <- rr$vectors[, idx, drop = FALSE]
  rr$Avectors <- rr$Avectors[, idx, drop = FALSE]
  rr$residual_vectors <- rr$residual_vectors[, idx, drop = FALSE]
  rr$residuals <- rr$residuals[idx]
  rr
}
