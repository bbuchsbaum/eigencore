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
                                           start = NULL, start_av = NULL) {
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
  if (!is.null(start_av)) {
    start_av <- as.matrix(start_av)
    if (nrow(start_av) != op$dim[[1L]] || ncol(start_av) > ncol(start)) {
      stop("start_av must have nrow equal to nrow(op) and no more columns than start.", call. = FALSE)
    }
  }

  storage <- op$metadata$storage %||% NULL
  source <- source_or_null(op)
  if (identical(storage, "dgCMatrix")) {
    A <- op$metadata$matrix
    if (is.null(start_av)) {
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
    } else {
      .Call(
        "eigencore_block_golub_kahan_csc_basis_cached",
        methods::slot(A, "i"),
        methods::slot(A, "p"),
        methods::slot(A, "x"),
        methods::slot(A, "Dim"),
        as.integer(max_subspace),
        start,
        start_av,
        PACKAGE = "eigencore"
      )
    }
  } else if (is.matrix(source) && is.double(source)) {
    if (is.null(start_av)) {
      .Call(
        "eigencore_block_golub_kahan_dense_basis",
        source,
        as.integer(max_subspace),
        start,
        PACKAGE = "eigencore"
      )
    } else {
      .Call(
        "eigencore_block_golub_kahan_dense_basis_cached",
        source,
        as.integer(max_subspace),
        start,
        start_av,
        PACKAGE = "eigencore"
      )
    }
  } else {
    stop("Native block Golub-Kahan basis currently supports dense double matrices and dgCMatrix operators only.",
         call. = FALSE)
  }
}

#' @keywords internal
native_block_golub_kahan_fit <- function(op, max_subspace, rank,
                                         target = largest(), block = NULL,
                                         start = NULL, start_av = NULL) {
  op <- as_operator(op)
  if (is.null(op$apply_adjoint)) {
    stop("Native block Golub-Kahan fit requires an adjoint operator.", call. = FALSE)
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
  if (!is.null(start_av)) {
    start_av <- as.matrix(start_av)
    if (nrow(start_av) != op$dim[[1L]] || ncol(start_av) > ncol(start)) {
      stop("start_av must have nrow equal to nrow(op) and no more columns than start.", call. = FALSE)
    }
  }

  target_kind <- as.integer(native_svd_target_kind(target))
  storage <- op$metadata$storage %||% NULL
  source <- source_or_null(op)
  if (identical(storage, "dgCMatrix")) {
    A <- op$metadata$matrix
    if (is.null(start_av)) {
      .Call(
        "eigencore_block_golub_kahan_csc_fit",
        methods::slot(A, "i"),
        methods::slot(A, "p"),
        methods::slot(A, "x"),
        methods::slot(A, "Dim"),
        as.integer(max_subspace),
        start,
        as.integer(rank),
        target_kind,
        PACKAGE = "eigencore"
      )
    } else {
      .Call(
        "eigencore_block_golub_kahan_csc_fit_cached",
        methods::slot(A, "i"),
        methods::slot(A, "p"),
        methods::slot(A, "x"),
        methods::slot(A, "Dim"),
        as.integer(max_subspace),
        start,
        as.integer(rank),
        target_kind,
        start_av,
        PACKAGE = "eigencore"
      )
    }
  } else if (is.matrix(source) && is.double(source)) {
    if (is.null(start_av)) {
      .Call(
        "eigencore_block_golub_kahan_dense_fit",
        source,
        as.integer(max_subspace),
        start,
        as.integer(rank),
        target_kind,
        PACKAGE = "eigencore"
      )
    } else {
      .Call(
        "eigencore_block_golub_kahan_dense_fit_cached",
        source,
        as.integer(max_subspace),
        start,
        as.integer(rank),
        target_kind,
        start_av,
        PACKAGE = "eigencore"
      )
    }
  } else {
    stop("Native block Golub-Kahan fit currently supports dense double matrices and dgCMatrix operators only.",
         call. = FALSE)
  }
}

#' @keywords internal
native_block_golub_kahan_retained_restart_abi <- function(op, rank,
                                                          target = largest(),
                                                          block = NULL,
                                                          max_attempts = NULL,
                                                          max_subspace = NULL,
                                                          adaptive_start = c("ritz")) {
  adaptive_start <- match.arg(adaptive_start, choices = c("ritz"))
  op <- as_operator(op)
  if (is.null(op$apply_adjoint)) {
    stop("Retained block Golub-Kahan restart ABI requires an adjoint operator.",
         call. = FALSE)
  }
  storage <- op$metadata$storage %||% NULL
  source <- source_or_null(op)
  native_storage <- if (identical(storage, "dgCMatrix")) {
    "dgCMatrix"
  } else if (is.matrix(source) && is.double(source)) {
    "double_matrix"
  } else {
    NA_character_
  }
  if (is.na(native_storage)) {
    stop("Retained block Golub-Kahan restart ABI currently supports dense double matrices and dgCMatrix operators only.",
         call. = FALSE)
  }

  n <- op$dim[[2L]]
  rank <- min(as.integer(rank), min(op$dim))
  if (length(rank) != 1L || is.na(rank) || rank < 1L) {
    stop("rank must be a positive integer.", call. = FALSE)
  }
  block <- as.integer(block %||% min(8L, max(2L, ceiling(rank / 4))))
  if (length(block) != 1L || is.na(block) || block < 1L) {
    stop("block must be a positive integer.", call. = FALSE)
  }
  initial_max_subspace <- min(
    n,
    as.integer(max_subspace %||% max(3L * rank + 20L, 6L * block + 20L))
  )
  if (length(initial_max_subspace) != 1L || is.na(initial_max_subspace) ||
      initial_max_subspace < rank) {
    stop("max_subspace must be at least rank.", call. = FALSE)
  }
  max_attempts <- as.integer(max_attempts %||% 3L)
  if (length(max_attempts) != 1L || is.na(max_attempts) || max_attempts < 1L) {
    stop("max_attempts must be a positive integer.", call. = FALSE)
  }

  subspaces <- integer(max_attempts)
  subspaces[[1L]] <- initial_max_subspace
  if (max_attempts > 1L) {
    for (attempt in seq.int(2L, max_attempts)) {
      previous <- subspaces[[attempt - 1L]]
      subspaces[[attempt]] <- min(
        n,
        max(
          previous + max(2L * rank, 4L * block, 10L),
          as.integer(ceiling(1.5 * previous))
        )
      )
    }
  }
  if (any(diff(subspaces) <= 0L)) {
    subspaces <- subspaces[c(TRUE, diff(subspaces) > 0L)]
  }

  structure(
    list(
      version = 1L,
      implemented = TRUE,
      entry_points = c(
        dense = "eigencore_block_golub_kahan_dense_retained_cycle",
        csc = "eigencore_block_golub_kahan_csc_retained_cycle"
      ),
      native_storage = native_storage,
      dim = op$dim,
      rank = rank,
      target_kind = native_svd_target_kind(target),
      policy = adaptive_start,
      block = block,
      max_attempts = length(subspaces),
      max_subspace_sequence = subspaces,
      input_schema = list(
        initial_start = c(n, block),
        random_tails = c(n, block, max(0L, length(subspaces) - 1L)),
        tail_layout = "column-major n x (block * restart_count)",
        operator = native_storage
      ),
      retained_state = c(
        "right Ritz vectors V_keep",
        "QR-normalized cached operator images A V_keep",
        "native basis workspace V/AV/U",
        "orthogonalization scratch",
        "attempt history",
        "certification-gated cached A V_keep retention with uncached fallback"
      ),
      output_schema = c(
        "d", "u", "v", "Avectors", "certificate diagnostics",
        "iterations", "matvecs", "ortho_passes", "attempt_history",
        "stage_seconds"
      ),
      invariants = c(
        "no R-side restart block construction after entry",
        "retained V_keep is reorthogonalized by the native basis runner on each attempt",
        "retained V_keep and A V_keep are transformed together by native QR normalization",
        "certification uses original operator coordinates",
        "attempt history fields match native_block_golub_kahan_cycle_svd()",
        "default policy is Ritz-plus-random until a benchmark proves otherwise",
        "cached A V retention must certify or fall back to the uncached retained path"
      )
    ),
    class = "eigencore_block_golub_kahan_retained_restart_abi"
  )
}

#' @keywords internal
native_block_golub_kahan_retained_cycle_svd <- function(op, rank,
                                                        target = largest(),
                                                        tol = 1e-8,
                                                        max_subspace = NULL,
                                                        block = NULL,
                                                        start = NULL,
                                                        max_attempts = NULL,
                                                        vectors = c("both", "left", "right", "none"),
                                                        retained_av_cache = FALSE) {
  vectors <- match.arg(vectors)
  op <- as_operator(op)
  abi <- native_block_golub_kahan_retained_restart_abi(
    op,
    rank = rank,
    target = target,
    block = block,
    max_attempts = max_attempts,
    max_subspace = max_subspace
  )
  n <- op$dim[[2L]]
  block <- abi$block
  rank <- abi$rank
  start <- if (is.null(start)) {
    matrix(stats::rnorm(n * block), nrow = n, ncol = block)
  } else {
    as.matrix(start)
  }
  if (nrow(start) != n || ncol(start) != block) {
    stop("start must have dimensions ncol(op) x block for retained block Golub-Kahan.",
         call. = FALSE)
  }
  restart_count <- max(0L, abi$max_attempts - 1L)
  random_tails <- matrix(
    stats::rnorm(n * block * restart_count),
    nrow = n,
    ncol = block * restart_count
  )
  target_kind <- as.integer(abi$target_kind)
  norm_info <- operator_norm_for_certificate_info(op)
  source <- source_or_null(op)
  storage <- op$metadata$storage %||% NULL
  run_retained <- function(use_retained_av_cache) {
    if (identical(storage, "dgCMatrix")) {
      A <- op$metadata$matrix
      .Call(
        "eigencore_block_golub_kahan_csc_retained_cycle",
        methods::slot(A, "i"),
        methods::slot(A, "p"),
        methods::slot(A, "x"),
        methods::slot(A, "Dim"),
        as.integer(abi$max_subspace_sequence[[1L]]),
        start,
        random_tails,
        as.integer(abi$max_attempts),
        as.integer(rank),
        target_kind,
        as.numeric(norm_info$value),
        as.numeric(tol),
        as.logical(use_retained_av_cache),
        PACKAGE = "eigencore"
      )
    } else if (is.matrix(source) && is.double(source)) {
      .Call(
        "eigencore_block_golub_kahan_dense_retained_cycle",
        source,
        as.integer(abi$max_subspace_sequence[[1L]]),
        start,
        random_tails,
        as.integer(abi$max_attempts),
        as.integer(rank),
        target_kind,
        as.numeric(norm_info$value),
        as.numeric(tol),
        as.logical(use_retained_av_cache),
        PACKAGE = "eigencore"
      )
    } else {
      stop("Native retained block Golub-Kahan cycle currently supports dense double matrices and dgCMatrix operators only.",
           call. = FALSE)
    }
  }

  cache_attempted <- isTRUE(retained_av_cache)
  cache_failed <- FALSE
  cache_error <- NA_character_
  cache_max_backward_error <- NA_real_
  cache_max_residual <- NA_real_
  cache_iterations <- 0L
  cache_matvecs <- 0L
  cache_ortho_passes <- 0L
  cache_stage_seconds <- NULL
  ritz <- NULL
  cert <- NULL
  if (cache_attempted) {
    cached <- tryCatch(
      run_retained(TRUE),
      error = function(e) structure(
        list(error = conditionMessage(e)),
        class = "eigencore_retained_av_cache_error"
      )
    )
    if (inherits(cached, "eigencore_retained_av_cache_error")) {
      cache_failed <- TRUE
      cache_error <- cached$error
    } else {
      cache_iterations <- cached$iterations %||% 0L
      cache_matvecs <- cached$matvecs %||% 0L
      cache_ortho_passes <- cached$ortho_passes %||% 0L
      cache_stage_seconds <- cached$stage_seconds
      cached_cert <- certify_svd_operator_cached_av(
        op, cached$d, cached$u, cached$v, cached$Avectors, tol = tol
      )
      if (isTRUE(cached_cert$passed)) {
        ritz <- cached
        cert <- cached_cert
      } else {
        cache_failed <- TRUE
        cache_max_backward_error <- cached_cert$max_backward_error
        cache_max_residual <- cached_cert$max_residual
      }
    }
  }
  if (is.null(ritz)) {
    ritz <- run_retained(FALSE)
    cert <- certify_svd_operator_cached_av(
      op, ritz$d, ritz$u, ritz$v, ritz$Avectors, tol = tol
    )
  }
  attempt_history <- ritz$attempt_history
  if (is.data.frame(attempt_history) &&
      !"certificate_passed" %in% names(attempt_history)) {
    attempt_history$certificate_passed <- FALSE
    attempt_history$max_backward_error <- Inf
    attempt_history$max_residual <- Inf
    if (nrow(attempt_history)) {
      last <- nrow(attempt_history)
      attempt_history$certificate_passed[[last]] <- isTRUE(cert$passed)
      attempt_history$max_backward_error[[last]] <- cert$max_backward_error
      attempt_history$max_residual[[last]] <- cert$max_residual
    }
  }
  reported_iterations <- ritz$iterations + if (cache_failed) cache_iterations else 0L
  reported_matvecs <- ritz$matvecs + if (cache_failed) cache_matvecs else 0L
  reported_ortho_passes <- ritz$ortho_passes + if (cache_failed) cache_ortho_passes else 0L
  reported_stage_seconds <- ritz$stage_seconds
  if (cache_failed && !is.null(cache_stage_seconds)) {
    merged_names <- union(names(reported_stage_seconds), names(cache_stage_seconds))
    reported_stage_seconds <- stats::setNames(
      vapply(
        merged_names,
        function(nm) {
          (reported_stage_seconds[[nm]] %||% 0) + (cache_stage_seconds[[nm]] %||% 0)
        },
        numeric(1)
      ),
      merged_names
    )
  }
  out <- list(
    d = ritz$d,
    u = ritz$u,
    v = ritz$v,
    values = ritz$d,
    residuals = cert$residuals,
    backward_error = cert$backward_error,
    orthogonality = cert$orthogonality,
    certificate = cert,
    iterations = reported_iterations,
    matvecs = reported_matvecs,
    block = block,
    stage_seconds = reported_stage_seconds,
    restart = list(
      kind = "block_golub_kahan_native_retained_cycle",
      implemented = TRUE,
      native = TRUE,
      thick_restart = TRUE,
      retained_restart = TRUE,
      retained_restart_native = TRUE,
      retained_av_cache = any(attempt_history$cached_start_used %||% FALSE),
      retained_av_cache_attempted = cache_attempted,
      retained_av_cache_fallback = cache_failed,
      retained_av_cache_failed_backward_error = cache_max_backward_error,
      retained_av_cache_failed_residual = cache_max_residual,
      retained_av_cache_failed_iterations = if (cache_failed) cache_iterations else 0L,
      retained_av_cache_failed_matvecs = if (cache_failed) cache_matvecs else 0L,
      retained_av_cache_failed_ortho_passes = if (cache_failed) cache_ortho_passes else 0L,
      native_attempt_certification = TRUE,
      native_early_stop = is.data.frame(attempt_history) &&
        any(attempt_history$certificate_passed) &&
        nrow(attempt_history) < abi$max_attempts,
      certified_attempt = if (is.data.frame(attempt_history) &&
          any(attempt_history$certificate_passed)) {
        which(attempt_history$certificate_passed)[[1L]]
      } else {
        NA_integer_
      },
      retained_restart_abi_version = abi$version,
      native_workspace_bytes = ritz$native_workspace_bytes %||% NA_real_,
      basis_returned = FALSE,
      adaptive = TRUE,
      adaptive_start = abi$policy,
      attempts = if (is.data.frame(attempt_history)) nrow(attempt_history) else abi$max_attempts,
      attempt = if (is.data.frame(attempt_history)) nrow(attempt_history) else abi$max_attempts,
      active_cols = ritz$active_cols,
      active_left_cols = ritz$active_left_cols,
      initial_max_subspace = abi$max_subspace_sequence[[1L]],
      max_subspace = if (is.data.frame(attempt_history) && nrow(attempt_history)) {
        tail(attempt_history$max_subspace, 1L)[[1L]]
      } else {
        tail(abi$max_subspace_sequence, 1L)[[1L]]
      },
      final_max_subspace = if (is.data.frame(attempt_history) && nrow(attempt_history)) {
        tail(attempt_history$max_subspace, 1L)[[1L]]
      } else {
        tail(abi$max_subspace_sequence, 1L)[[1L]]
      },
      attempted_subspaces = if (is.data.frame(attempt_history) && nrow(attempt_history)) {
        attempt_history$max_subspace
      } else {
        abi$max_subspace_sequence
      },
      attempt_history = attempt_history,
      total_iterations = reported_iterations,
      total_matvecs = reported_matvecs,
      total_ortho_passes = reported_ortho_passes,
      fallback_attempted = cache_failed,
      fallback_used = cache_failed,
      fallback_method = if (cache_failed) "retained_uncached_after_cached_av_failure" else NA_character_,
      fallback_error = if (cache_failed && !is.na(cache_error)) cache_error else NA_character_,
      fallback_max_backward_error = if (cache_failed) cache_max_backward_error else NA_real_,
      final_attempt_iterations = if (is.data.frame(attempt_history) && nrow(attempt_history)) {
        tail(attempt_history$iterations, 1L)[[1L]]
      } else {
        NA_integer_
      },
      final_attempt_matvecs = if (is.data.frame(attempt_history) && nrow(attempt_history)) {
        tail(attempt_history$matvecs, 1L)[[1L]]
      } else {
        NA_integer_
      },
      final_attempt_ortho_passes = if (is.data.frame(attempt_history) && nrow(attempt_history)) {
        tail(attempt_history$ortho_passes, 1L)[[1L]]
      } else {
        NA_integer_
      },
      final_iterations = if (is.data.frame(attempt_history) && nrow(attempt_history)) {
        tail(attempt_history$iterations, 1L)[[1L]]
      } else {
        NA_integer_
      },
      final_matvecs = if (is.data.frame(attempt_history) && nrow(attempt_history)) {
        tail(attempt_history$matvecs, 1L)[[1L]]
      } else {
        NA_integer_
      },
      ortho_passes = reported_ortho_passes,
      matvecs = reported_matvecs,
      cached_start_used = isTRUE(ritz$cached_start_used),
      stage_seconds = reported_stage_seconds
    )
  )
  out <- complete_zero_singular_triplets(op, out, rank, target, tol)
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

block_golub_kahan_normalize_cached_start <- function(start, start_av) {
  if (is.null(start_av)) {
    return(list(start = start, start_av = NULL))
  }
  cached_cols <- ncol(start_av)
  cached_start <- start[, seq_len(cached_cols), drop = FALSE]
  qr_start <- qr(cached_start)
  q <- qr.Q(qr_start)
  r <- qr.R(qr_start)
  if (ncol(q) < cached_cols) {
    stop("cached block Golub-Kahan restart prefix is rank deficient.", call. = FALSE)
  }
  normalized_prefix <- q[, seq_len(cached_cols), drop = FALSE]
  normalized_start <- if (cached_cols < ncol(start)) {
    out <- matrix(0, nrow = nrow(start), ncol = ncol(start))
    out[, seq_len(cached_cols)] <- normalized_prefix
    tail_cols <- seq.int(cached_cols + 1L, ncol(start))
    out[, tail_cols] <- start[, tail_cols, drop = FALSE]
    out
  } else {
    normalized_prefix
  }
  list(
    start = normalized_start,
    start_av = start_av %*% backsolve(r, diag(cached_cols))
  )
}

block_golub_kahan_restart_with_random_tail <- function(prefix, n, block) {
  prefix_cols <- ncol(prefix)
  out <- matrix(0, nrow = n, ncol = prefix_cols + block)
  out[, seq_len(prefix_cols)] <- prefix
  out[, seq.int(prefix_cols + 1L, prefix_cols + block)] <- stats::rnorm(n * block)
  out
}

block_golub_kahan_restart_with_tail <- function(prefix, tail, n, block) {
  prefix_cols <- ncol(prefix)
  tail_cols <- min(ncol(tail), block)
  out <- matrix(0, nrow = n, ncol = prefix_cols + block)
  out[, seq_len(prefix_cols)] <- prefix
  if (tail_cols > 0L) {
    out[, seq.int(prefix_cols + 1L, prefix_cols + tail_cols)] <-
      tail[, seq_len(tail_cols), drop = FALSE]
  }
  if (tail_cols < block) {
    random_cols <- seq.int(prefix_cols + tail_cols + 1L, prefix_cols + block)
    out[, random_cols] <- stats::rnorm(n * length(random_cols))
  }
  out
}

#' @keywords internal
native_block_golub_kahan_cycle_svd <- function(op, rank, target = largest(),
                                               tol = 1e-8,
                                               max_subspace = NULL,
                                               block = NULL,
                                               start = NULL,
                                               adaptive = is.null(max_subspace),
                                               max_attempts = NULL,
                                               adaptive_start = c("ritz", "ritz_cached", "ritz_cached_random", "ritz_residual", "ritz_lean", "initial"),
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
  current_start_av <- NULL

  attempt_rows <- list()
  total_iterations <- 0L
  total_matvecs <- 0L
  total_ortho_passes <- 0L
  total_stage_seconds <- c(native_iteration = 0, ritz = 0)
  current_max_subspace <- initial_max_subspace
  out <- NULL
  attempt <- 0L
  repeat {
    attempt <- attempt + 1L
    cached_start <- block_golub_kahan_normalize_cached_start(
      current_start,
      current_start_av
    )
    ritz <- native_block_golub_kahan_fit(
      op,
      max_subspace = current_max_subspace,
      rank = rank,
      target = target,
      block = block,
      start = cached_start$start,
      start_av = cached_start$start_av
    )
    cert_info <- certify_svd_operator_cached_av(
      op, ritz$d, ritz$u, ritz$v, ritz$Avectors, tol = tol,
      return_residual_vectors = identical(adaptive_start, "ritz_residual")
    )
    wrapped_cert <- cert_info[["certificate", exact = TRUE]]
    if (!is.null(wrapped_cert)) {
      cert <- wrapped_cert
      right_residual_vectors <- cert_info[["right_residual_vectors", exact = TRUE]]
    } else {
      cert <- cert_info
      right_residual_vectors <- NULL
    }
    out <- list(
      d = ritz$d,
      u = ritz$u,
      v = ritz$v,
      values = ritz$d,
      residuals = cert$residuals,
      backward_error = cert$backward_error,
      orthogonality = cert$orthogonality,
      certificate = cert,
      iterations = ritz$iterations,
      matvecs = ritz$matvecs,
      block = block,
      restart = list(
        kind = "block_golub_kahan_native_basis_cycle",
        implemented = TRUE,
        native = TRUE,
        thick_restart = FALSE,
        basis_returned = FALSE,
        adaptive = adaptive,
        adaptive_start = adaptive_start,
        attempt = attempt,
        active_cols = ritz$active_cols,
        active_left_cols = ritz$active_left_cols,
        max_subspace = current_max_subspace,
        block = block,
        ortho_passes = ritz$ortho_passes,
        matvecs = ritz$matvecs,
        cached_start_used = isTRUE(ritz$cached_start_used)
      )
    )
    out <- complete_zero_singular_triplets(op, out, rank, target, tol)

    total_iterations <- total_iterations + ritz$iterations
    total_matvecs <- total_matvecs + ritz$matvecs
    total_ortho_passes <- total_ortho_passes + ritz$ortho_passes
    stage_seconds <- ritz$stage_seconds %||% numeric()
    for (stage_name in names(total_stage_seconds)) {
      total_stage_seconds[[stage_name]] <- total_stage_seconds[[stage_name]] +
        (stage_seconds[[stage_name]] %||% 0)
    }
    attempt_rows[[attempt]] <- data.frame(
      attempt = attempt,
      max_subspace = current_max_subspace,
      active_cols = ritz$active_cols,
      start_cols = ncol(current_start),
      cached_start_used = isTRUE(ritz$cached_start_used),
      warm_started = attempt > 1L && adaptive_start %in% c("ritz", "ritz_cached", "ritz_cached_random", "ritz_residual", "ritz_lean"),
      iterations = ritz$iterations,
      matvecs = ritz$matvecs,
      ortho_passes = ritz$ortho_passes,
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
      current_start <- block_golub_kahan_restart_with_random_tail(
        ritz$v[, seq_len(keep_cols), drop = FALSE],
        n,
        block
      )
      current_start_av <- NULL
    } else if (identical(adaptive_start, "ritz_cached")) {
      keep_cols <- min(ncol(ritz$v), max(rank, block))
      current_start <- ritz$v[, seq_len(keep_cols), drop = FALSE]
      current_start_av <- ritz$Avectors[, seq_len(keep_cols), drop = FALSE]
    } else if (identical(adaptive_start, "ritz_cached_random")) {
      keep_cols <- min(ncol(ritz$v), max(rank, block))
      current_start <- block_golub_kahan_restart_with_random_tail(
        ritz$v[, seq_len(keep_cols), drop = FALSE],
        n,
        block
      )
      current_start_av <- ritz$Avectors[, seq_len(keep_cols), drop = FALSE]
    } else if (identical(adaptive_start, "ritz_residual")) {
      keep_cols <- min(ncol(ritz$v), max(rank, block))
      current_start <- block_golub_kahan_restart_with_tail(
        ritz$v[, seq_len(keep_cols), drop = FALSE],
        right_residual_vectors %||% matrix(0, n, 0L),
        n,
        block
      )
      current_start_av <- NULL
    } else if (identical(adaptive_start, "ritz_lean")) {
      keep_cols <- min(ncol(ritz$v), max(rank, block))
      current_start <- ritz$v[, seq_len(keep_cols), drop = FALSE]
      current_start_av <- NULL
    } else {
      current_start <- start
      current_start_av <- NULL
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
  out$stage_seconds <- total_stage_seconds
  out$restart$stage_seconds <- total_stage_seconds
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
