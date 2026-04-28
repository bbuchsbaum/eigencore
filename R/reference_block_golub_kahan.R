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
  rank <- as.integer(rank)
  if (length(rank) != 1L || is.na(rank) || rank < 1L) {
    stop("rank must be a positive integer.", call. = FALSE)
  }
  rank <- min(rank, limit)

  block <- if (is.null(block)) {
    min(8L, max(2L, as.integer(ceiling(rank / 4))))
  } else {
    as.integer(block)
  }
  if (length(block) != 1L || is.na(block) || block < 1L) {
    stop("block must be a positive integer.", call. = FALSE)
  }
  block <- min(block, limit)

  if (is.null(max_subspace)) {
    max_subspace <- max(3L * rank + 20L, 6L * block + 20L)
  }
  max_subspace <- min(n, as.integer(max_subspace))
  if (length(max_subspace) != 1L || is.na(max_subspace) || max_subspace < rank) {
    stop("max_subspace must be at least rank.", call. = FALSE)
  }

  max_restarts <- as.integer(max_restarts)
  if (length(max_restarts) != 1L || is.na(max_restarts) || max_restarts < 0L) {
    stop("max_restarts must be a non-negative integer.", call. = FALSE)
  }
  pad <- if (is.null(pad)) min(rank, max(block, 5L)) else as.integer(pad)
  if (length(pad) != 1L || is.na(pad) || pad < 0L) {
    stop("pad must be a non-negative integer.", call. = FALSE)
  }

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

  accepted <- reference_block_accept(
    matrix(stats::rnorm(n * block), nrow = n, ncol = block),
    qv_lock,
    matrix(0, n, 0L),
    max_cols = max_subspace
  )
  ortho_passes <- ortho_passes + accepted$ortho_passes
  V <- accepted$Q
  if (ncol(V) == 0L) {
    V <- diag(1, n, min(block, n))
  }
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
    take <- min(length(rr$d), max(remaining + pad, remaining))
    if (take < 1L) {
      break
    }
    rr <- reference_block_golub_kahan_slice(rr, seq_len(take))

    lock_now <- 0L
    for (i in seq_len(min(remaining, length(rr$d)))) {
      if (!isTRUE(rr$certificate$converged[[i]])) {
        break
      }
      lock_now <- lock_now + 1L
    }
    if (lock_now > 0L) {
      locked_d <- c(locked_d, rr$d[seq_len(lock_now)])
      locked_u <- cbind(locked_u, rr$u[, seq_len(lock_now), drop = FALSE])
      locked_v <- cbind(locked_v, rr$v[, seq_len(lock_now), drop = FALSE])
      locked_residuals <- c(locked_residuals, rr$residuals[seq_len(lock_now)])
      qu_lock <- locked_u
      qv_lock <- locked_v
      locking_events <- locking_events + 1L
    }

    history[[length(history) + 1L]] <- data.frame(
      restart = restart,
      iteration = iterations,
      n_locked = length(locked_d),
      max_residual = rr$certificate$max_residual,
      max_backward_error = rr$certificate$max_backward_error
    )

    unlocked_idx <- if (lock_now < length(rr$d)) {
      seq.int(lock_now + 1L, length(rr$d))
    } else {
      integer()
    }
    last_candidates <- reference_block_golub_kahan_slice(rr, unlocked_idx)
    if (length(locked_d) >= rank || restart == max_restarts) {
      restarts_used <- restart
      break
    }

    remaining <- rank - length(locked_d)
    k_keep <- min(length(unlocked_idx), max(remaining + pad, remaining), max_subspace - block)
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
  idx <- idx[seq_len(length(idx))]
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

#' @keywords internal
native_block_golub_kahan_ritz <- function(V, AV, rank, target = largest(),
                                          active_cols = NULL) {
  V <- as.matrix(V)
  AV <- as.matrix(AV)
  active_cols <- as.integer(active_cols %||% ncol(V))
  .Call(
    "eigencore_block_golub_kahan_ritz",
    V,
    AV,
    as.integer(rank),
    as.integer(native_svd_target_kind(target)),
    active_cols,
    PACKAGE = "eigencore"
  )
}

#' @keywords internal
native_block_golub_kahan_basis <- function(op, max_subspace, block = NULL,
                                           start = NULL) {
  op <- as_operator(op)
  if (is.null(op$apply_adjoint)) {
    stop("Native block Golub-Kahan basis requires an adjoint operator.", call. = FALSE)
  }
  n <- op$dim[[2L]]
  block <- as.integer(block %||% 2L)
  if (length(block) != 1L || is.na(block) || block < 1L) {
    stop("block must be a positive integer.", call. = FALSE)
  }
  max_subspace <- min(n, as.integer(max_subspace))
  if (length(max_subspace) != 1L || is.na(max_subspace) || max_subspace < 1L) {
    stop("max_subspace must be a positive integer.", call. = FALSE)
  }
  start <- if (is.null(start)) {
    matrix(stats::rnorm(n * block), nrow = n, ncol = block)
  } else {
    as.matrix(start)
  }
  if (nrow(start) != n || ncol(start) < 1L) {
    stop("start must have nrow equal to ncol(op) and at least one column.", call. = FALSE)
  }

  storage <- op$metadata$storage %||% NULL
  source <- source_or_null(op)
  if (identical(storage, "dgCMatrix")) {
    A <- op$metadata$matrix
    .Call(
      "eigencore_block_golub_kahan_csc_basis",
      methods::slot(A, "i"),
      methods::slot(A, "p"),
      methods::slot(A, "x"),
      methods::slot(A, "Dim"),
      as.integer(max_subspace),
      start,
      PACKAGE = "eigencore"
    )
  } else if (is.matrix(source) && is.double(source)) {
    .Call(
      "eigencore_block_golub_kahan_dense_basis",
      source,
      as.integer(max_subspace),
      start,
      PACKAGE = "eigencore"
    )
  } else {
    stop("Native block Golub-Kahan basis currently supports dense double matrices and dgCMatrix operators only.",
         call. = FALSE)
  }
}

#' @keywords internal
native_block_golub_kahan_cycle_svd <- function(op, rank, target = largest(),
                                               tol = 1e-8,
                                               max_subspace = NULL,
                                               block = NULL,
                                               start = NULL,
                                               adaptive = is.null(max_subspace),
                                               max_attempts = NULL,
                                               adaptive_start = c("ritz", "ritz_lean", "initial"),
                                               vectors = c("both", "left", "right", "none")) {
  vectors <- match.arg(vectors)
  adaptive_start <- match.arg(adaptive_start)
  op <- as_operator(op)
  rank <- min(as.integer(rank), min(op$dim))
  n <- op$dim[[2L]]
  block <- as.integer(block %||% min(8L, max(2L, ceiling(rank / 4))))
  initial_max_subspace <- min(
    n,
    as.integer(max_subspace %||% max(3L * rank + 20L, 6L * block + 20L))
  )
  if (length(initial_max_subspace) != 1L || is.na(initial_max_subspace) || initial_max_subspace < rank) {
    stop("max_subspace must be at least rank.", call. = FALSE)
  }
  adaptive <- isTRUE(adaptive)
  max_attempts <- as.integer(max_attempts %||% if (adaptive) 4L else 1L)
  if (length(max_attempts) != 1L || is.na(max_attempts) || max_attempts < 1L) {
    stop("max_attempts must be a positive integer.", call. = FALSE)
  }
  start <- if (is.null(start)) {
    matrix(stats::rnorm(n * block), nrow = n, ncol = block)
  } else {
    as.matrix(start)
  }
  current_start <- start

  attempt_rows <- list()
  total_iterations <- 0L
  total_matvecs <- 0L
  total_ortho_passes <- 0L
  current_max_subspace <- initial_max_subspace
  out <- NULL
  attempt <- 0L
  repeat {
    attempt <- attempt + 1L
    basis <- native_block_golub_kahan_basis(
      op,
      max_subspace = current_max_subspace,
      block = block,
      start = current_start
    )
    ritz <- native_block_golub_kahan_ritz(
      basis$V,
      basis$AV,
      rank = rank,
      target = target,
      active_cols = basis$active_cols
    )
    cert <- certify_svd_operator(op, ritz$d, ritz$u, ritz$v, tol = tol)
    out <- list(
      d = ritz$d,
      u = ritz$u,
      v = ritz$v,
      values = ritz$d,
      residuals = cert$residuals,
      backward_error = cert$backward_error,
      orthogonality = cert$orthogonality,
      certificate = cert,
      iterations = basis$iterations,
      matvecs = basis$matvecs,
      block = block,
      restart = list(
        kind = "block_golub_kahan_native_basis_cycle",
        implemented = TRUE,
        native = TRUE,
        thick_restart = FALSE,
        adaptive = adaptive,
        adaptive_start = adaptive_start,
        attempt = attempt,
        active_cols = basis$active_cols,
        active_left_cols = basis$active_left_cols,
        max_subspace = current_max_subspace,
        block = block,
        ortho_passes = basis$ortho_passes,
        matvecs = basis$matvecs
      )
    )
    out <- complete_zero_singular_triplets(op, out, rank, target, tol)

    total_iterations <- total_iterations + basis$iterations
    total_matvecs <- total_matvecs + basis$matvecs
    total_ortho_passes <- total_ortho_passes + basis$ortho_passes
    attempt_rows[[attempt]] <- data.frame(
      attempt = attempt,
      max_subspace = current_max_subspace,
      active_cols = basis$active_cols,
      start_cols = ncol(current_start),
      warm_started = attempt > 1L && adaptive_start %in% c("ritz", "ritz_lean"),
      iterations = basis$iterations,
      matvecs = basis$matvecs,
      ortho_passes = basis$ortho_passes,
      certificate_passed = isTRUE(out$certificate$passed),
      max_backward_error = out$certificate$max_backward_error,
      max_residual = out$certificate$max_residual
    )

    if (isTRUE(out$certificate$passed) || !adaptive ||
        attempt >= max_attempts || current_max_subspace >= n) {
      break
    }
    next_max_subspace <- min(
      n,
      max(
        current_max_subspace + max(2L * rank, 4L * block, 10L),
        as.integer(ceiling(1.5 * current_max_subspace))
      )
    )
    if (next_max_subspace <= current_max_subspace) {
      break
    }
    current_max_subspace <- next_max_subspace
    if (identical(adaptive_start, "ritz")) {
      keep_cols <- min(ncol(ritz$v), max(rank, block))
      current_start <- cbind(
        ritz$v[, seq_len(keep_cols), drop = FALSE],
        matrix(stats::rnorm(n * block), nrow = n, ncol = block)
      )
    } else if (identical(adaptive_start, "ritz_lean")) {
      keep_cols <- min(ncol(ritz$v), max(rank, block))
      current_start <- ritz$v[, seq_len(keep_cols), drop = FALSE]
    } else {
      current_start <- start
    }
  }

  out$iterations <- total_iterations
  out$matvecs <- total_matvecs
  out$restart$attempts <- attempt
  out$restart$initial_max_subspace <- initial_max_subspace
  out$restart$final_max_subspace <- out$restart$max_subspace
  out$restart$attempted_subspaces <- vapply(attempt_rows, `[[`, integer(1), "max_subspace")
  out$restart$attempt_history <- do.call(rbind, attempt_rows)
  out$restart$total_iterations <- total_iterations
  out$restart$total_matvecs <- total_matvecs
  out$restart$total_ortho_passes <- total_ortho_passes
  out$restart$final_attempt_iterations <- attempt_rows[[attempt]]$iterations
  out$restart$final_attempt_matvecs <- attempt_rows[[attempt]]$matvecs
  out$restart$final_attempt_ortho_passes <- attempt_rows[[attempt]]$ortho_passes
  out$restart$final_iterations <- out$restart$final_attempt_iterations
  out$restart$final_matvecs <- out$restart$final_attempt_matvecs
  if (vectors == "left") {
    out$v <- NULL
  } else if (vectors == "right") {
    out$u <- NULL
  } else if (vectors == "none") {
    out$u <- NULL
    out$v <- NULL
  }
  out
}
