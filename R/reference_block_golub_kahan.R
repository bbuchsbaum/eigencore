#' @keywords internal
reference_block_golub_kahan_thick_restart_svd <- function(op, rank, target = largest(),
                                                          tol = 1e-8,
                                                          max_subspace = NULL,
                                                          max_restarts = 100L,
                                                          block = NULL,
                                                          pad = NULL,
                                                          vectors = c("both", "left", "right", "none")) {
  vectors <- match.arg(vectors)
  op <- as_operator(op)
  if (is.null(op$apply_adjoint)) {
    stop("Block Golub-Kahan SVD requires an adjoint operator.", call. = FALSE)
  }

  m <- op$dim[[1L]]
  n <- op$dim[[2L]]
  limit <- min(m, n)
  controls <- reference_block_restart_controls(
    requested = rank,
    requested_name = "rank",
    limit = limit,
    block = block,
    default_block = function(requested) {
      min(8L, max(2L, as.integer(ceiling(requested / 4))))
    },
    max_subspace = max_subspace,
    default_max_subspace = function(requested, block) {
      max(3L * requested + 20L, 6L * block + 20L)
    },
    max_restarts = max_restarts,
    pad = pad,
    clamp_requested = TRUE,
    clamp_block = TRUE
  )
  rank <- controls$requested
  block <- controls$block
  max_subspace <- controls$max_subspace
  max_restarts <- controls$max_restarts
  pad <- controls$pad

  qv_lock <- matrix(0, n, 0L)
  qu_lock <- matrix(0, m, 0L)
  locked_d <- numeric()
  locked_u <- matrix(0, m, 0L)
  locked_v <- matrix(0, n, 0L)
  locked_residuals <- numeric()
  iterations <- 0L
  matvecs <- 0L
  adjoint_matvecs <- 0L
  ortho_passes <- 0L
  locking_events <- 0L
  history <- list()

  accepted <- reference_block_initial_basis(n, block, qv_lock, max_subspace)
  V <- accepted$Q
  ortho_passes <- ortho_passes + accepted$ortho_passes
  AV <- apply_operator(op, V)
  matvecs <- matvecs + 1L
  U_basis <- matrix(0, m, 0L)

  last_candidates <- NULL
  restarts_used <- 0L
  for (restart in 0:max_restarts) {
    while (ncol(V) < max_subspace) {
      source_cols <- seq.int(max(1L, ncol(V) - block + 1L), ncol(V))
      left_seed <- AV[, source_cols, drop = FALSE]
      accepted_u <- reference_block_accept(
        left_seed,
        qu_lock,
        U_basis,
        max_cols = block
      )
      ortho_passes <- ortho_passes + accepted_u$ortho_passes
      if (ncol(accepted_u$Q) == 0L) {
        break
      }
      U_basis <- cbind(U_basis, accepted_u$Q)

      right_seed <- apply_adjoint_operator(op, accepted_u$Q)
      matvecs <- matvecs + 1L
      adjoint_matvecs <- adjoint_matvecs + 1L
      accepted_v <- reference_block_accept(
        right_seed,
        qv_lock,
        V,
        max_cols = max_subspace - ncol(V)
      )
      ortho_passes <- ortho_passes + accepted_v$ortho_passes
      if (ncol(accepted_v$Q) == 0L) {
        break
      }
      V <- cbind(V, accepted_v$Q)
      AV <- cbind(AV, apply_operator(op, accepted_v$Q))
      matvecs <- matvecs + 1L
      iterations <- iterations + 1L
    }

    rr <- reference_block_golub_kahan_rr(op, V, AV, target = target, tol = tol)
    matvecs <- matvecs + rr$matvecs
    adjoint_matvecs <- adjoint_matvecs + rr$adjoint_matvecs
    remaining <- rank - length(locked_d)
    take <- reference_block_take_count(length(rr$d), remaining, pad)
    if (take < 1L) {
      break
    }
    rr <- reference_block_golub_kahan_slice(rr, seq_len(take))

    lock_now <- reference_block_lock_count(rr$certificate$converged, remaining)
    if (lock_now > 0L) {
      locked_d <- c(locked_d, rr$d[seq_len(lock_now)])
      locked_u <- cbind(locked_u, rr$u[, seq_len(lock_now), drop = FALSE])
      locked_v <- cbind(locked_v, rr$v[, seq_len(lock_now), drop = FALSE])
      locked_residuals <- c(locked_residuals, rr$residuals[seq_len(lock_now)])
      qu_lock <- locked_u
      qv_lock <- locked_v
      locking_events <- locking_events + 1L
    }

    history[[length(history) + 1L]] <- reference_block_history_frame(
      restart,
      iterations,
      length(locked_d),
      rr$certificate
    )

    unlocked_idx <- reference_block_unlocked_indices(lock_now, length(rr$d))
    last_candidates <- reference_block_golub_kahan_slice(rr, unlocked_idx)
    if (length(locked_d) >= rank || restart == max_restarts) {
      restarts_used <- restart
      break
    }

    remaining <- rank - length(locked_d)
    k_keep <- reference_block_keep_count(
      length(unlocked_idx),
      remaining,
      pad,
      max_subspace,
      block
    )
    if (k_keep > 0L) {
      keep_idx <- seq_len(k_keep)
      V <- last_candidates$v[, keep_idx, drop = FALSE]
      AV <- last_candidates$Avectors[, keep_idx, drop = FALSE]
      U_basis <- last_candidates$u[, keep_idx, drop = FALSE]
    } else {
      V <- matrix(0, n, 0L)
      AV <- matrix(0, m, 0L)
      U_basis <- matrix(0, m, 0L)
    }

    tail_count <- min(block, length(unlocked_idx))
    tail <- if (tail_count > 0L) {
      last_candidates$right_residual_vectors[, seq_len(tail_count), drop = FALSE]
    } else {
      matrix(stats::rnorm(n * block), nrow = n, ncol = block)
    }
    accepted <- reference_block_accept(tail, qv_lock, V, max_cols = max_subspace - ncol(V))
    ortho_passes <- ortho_passes + accepted$ortho_passes
    if (ncol(accepted$Q) == 0L) {
      accepted <- reference_block_accept(
        matrix(stats::rnorm(n * block), nrow = n, ncol = block),
        qv_lock,
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

  d <- locked_d
  u <- locked_u
  v <- locked_v
  residuals <- locked_residuals
  if (length(d) < rank && !is.null(last_candidates) && length(last_candidates$d)) {
    need <- min(rank - length(d), length(last_candidates$d))
    d <- c(d, last_candidates$d[seq_len(need)])
    u <- cbind(u, last_candidates$u[, seq_len(need), drop = FALSE])
    v <- cbind(v, last_candidates$v[, seq_len(need), drop = FALSE])
    residuals <- c(residuals, last_candidates$residuals[seq_len(need)])
  }
  if (length(d) > rank) {
    d <- d[seq_len(rank)]
    u <- u[, seq_len(rank), drop = FALSE]
    v <- v[, seq_len(rank), drop = FALSE]
    residuals <- residuals[seq_len(rank)]
  }

  cert <- certify_svd_operator(op, d, u, v, tol = tol)
  final <- list(
    d = d,
    u = u,
    v = v,
    values = d,
    residuals = cert$residuals,
    backward_error = cert$backward_error,
    orthogonality = cert$orthogonality,
    certificate = cert,
    iterations = iterations,
    matvecs = matvecs,
    adjoint_matvecs = adjoint_matvecs,
    restarts = restarts_used,
    ortho_passes = ortho_passes,
    locking_events = locking_events,
    block = block,
    convergence_history = if (length(history)) do.call(rbind, history) else data.frame(),
    locked = seq_len(length(locked_d)),
    restart = list(
      kind = "block_golub_kahan_thick_restart_reference",
      implemented = TRUE,
      native = FALSE,
      locking = "reference_loop",
      locked_count = length(locked_d),
      restarts_used = restarts_used,
      max_restarts = max_restarts,
      max_subspace = max_subspace,
      final_active_subspace = ncol(V),
      block = block,
      pad = pad,
      ortho_passes = ortho_passes,
      locking_events = locking_events,
      matvecs = matvecs,
      adjoint_matvecs = adjoint_matvecs
    )
  )
  final <- complete_zero_singular_triplets(op, final, rank, target, tol)

  if (vectors == "left") {
    final$v <- NULL
  } else if (vectors == "right") {
    final$u <- NULL
  } else if (vectors == "none") {
    final$u <- NULL
    final$v <- NULL
  }
  final
}

#' @keywords internal
reference_block_golub_kahan_rr <- function(op, V, AV, target, tol) {
  if (ncol(V) == 0L) {
    stop("Block Golub-Kahan Ritz extraction requires a non-empty right basis.",
         call. = FALSE)
  }
  small <- svd(AV, nu = min(nrow(AV), ncol(AV)), nv = ncol(V))
  idx <- order_indices(small$d, target)
  # Defensive clamp: order_indices ranks values within length(small$d), but
  # when the underlying SVD returns fewer numerically-significant singular
  # values than nu (rank-deficient blocks), an index from a wider sort could
  # exceed ncol(small$u). Drop any idx beyond the actual returned vector
  # count rather than letting it subscript out-of-bounds.
  idx <- idx[idx <= length(small$d) & idx <= ncol(small$u) & idx <= ncol(small$v)]
  d <- small$d[idx]
  u <- small$u[, idx, drop = FALSE]
  coeff <- small$v[, idx, drop = FALSE]
  v <- V %*% coeff
  Avectors <- AV %*% coeff
  left_residual_vectors <- Avectors - sweep(u, 2L, d, `*`)
  Atu <- apply_adjoint_operator(op, u)
  right_residual_vectors <- Atu - sweep(v, 2L, d, `*`)
  residuals <- sqrt(colSums(left_residual_vectors^2) + colSums(right_residual_vectors^2))
  cert <- certify_svd_operator(op, d, u, v, tol = tol)
  list(
    d = d,
    u = u,
    v = v,
    Avectors = Avectors,
    left_residual_vectors = left_residual_vectors,
    right_residual_vectors = right_residual_vectors,
    residuals = residuals,
    certificate = cert,
    matvecs = 1L,
    adjoint_matvecs = 1L
  )
}

#' @keywords internal
reference_block_golub_kahan_slice <- function(rr, idx) {
  if (!length(idx)) {
    m <- if (is.null(rr$u)) 0L else nrow(rr$u)
    n <- if (is.null(rr$v)) 0L else nrow(rr$v)
    return(list(
      d = numeric(),
      u = matrix(0, m, 0L),
      v = matrix(0, n, 0L),
      Avectors = matrix(0, m, 0L),
      left_residual_vectors = matrix(0, m, 0L),
      right_residual_vectors = matrix(0, n, 0L),
      residuals = numeric(),
      certificate = rr$certificate,
      matvecs = 0L,
      adjoint_matvecs = 0L
    ))
  }
  rr$d <- rr$d[idx]
  rr$u <- rr$u[, idx, drop = FALSE]
  rr$v <- rr$v[, idx, drop = FALSE]
  rr$Avectors <- rr$Avectors[, idx, drop = FALSE]
  rr$left_residual_vectors <- rr$left_residual_vectors[, idx, drop = FALSE]
  rr$right_residual_vectors <- rr$right_residual_vectors[, idx, drop = FALSE]
  rr$residuals <- rr$residuals[idx]
  rr$matvecs <- rr$matvecs %||% 0L
  rr$adjoint_matvecs <- rr$adjoint_matvecs %||% 0L
  rr
}
