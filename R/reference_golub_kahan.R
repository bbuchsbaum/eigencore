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
  controls <- reference_scalar_subspace_controls(
    requested = rank,
    requested_name = "rank",
    limit = limit,
    maxit = maxit,
    default_maxit = function(rank) max(20L, 4L * rank + 20L)
  )
  rank <- controls$requested
  maxit <- controls$maxit

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
      iterations <- j - 1L
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
native_irlba_lbd_restart_abi <- function(op, rank, target = largest(),
                                         work = NULL,
                                         max_restarts = NULL,
                                         retained = NULL,
                                         reorth_policy = c(
                                           "one_sided_small_side",
                                           "full_two_sided",
                                           "bpro_two_sided",
                                           "bpro_one_sided_guarded",
                                           "bpro_block_guarded"
                                         )) {
  reorth_policy <- match.arg(reorth_policy)
  op <- as_operator(op)
  if (is.null(op$apply_adjoint)) {
    stop("Retained one-sided IRLBA/LBD restart ABI requires an adjoint operator.",
         call. = FALSE)
  }
  if (!native_gram_svd_target_supported(target)) {
    stop("Retained one-sided IRLBA/LBD restart ABI currently supports only largest singular values.",
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
    stop("Retained one-sided IRLBA/LBD restart ABI currently supports dense double matrices and dgCMatrix operators only.",
         call. = FALSE)
  }

  dims <- as.integer(op$dim)
  limit <- min(dims)
  rank <- min(as.integer(rank), limit)
  if (length(rank) != 1L || is.na(rank) || rank < 1L) {
    stop("rank must be a positive integer.", call. = FALSE)
  }
  work <- as.integer(work %||% max(rank + 7L, 2L * rank + 1L))
  work <- min(limit, work)
  if (length(work) != 1L || is.na(work) || work < rank) {
    stop("work must be at least rank.", call. = FALSE)
  }
  retained <- as.integer(retained %||% min(work - 1L, rank + 2L))
  retained <- min(retained, work - 1L)
  if (length(retained) != 1L || is.na(retained) || retained < rank) {
    stop("retained must be at least rank and smaller than work.", call. = FALSE)
  }
  max_restarts <- as.integer(max_restarts %||% 5L)
  if (length(max_restarts) != 1L || is.na(max_restarts) || max_restarts < 0L) {
    stop("max_restarts must be a non-negative integer.", call. = FALSE)
  }

  transposed <- dims[[1L]] < dims[[2L]]
  active_dim <- if (transposed) rev(dims) else dims
  small_side <- min(active_dim)
  large_side <- max(active_dim)
  active_domain <- active_dim[[2L]]
  active_codomain <- active_dim[[1L]]
  restart_budget <- max(0L, max_restarts)

  structure(
    list(
      version = 1L,
      implemented = TRUE,
      entry_points = c(
        dense = "eigencore_irlba_lbd_dense_retained",
        csc = "eigencore_irlba_lbd_csc_retained"
      ),
      native_storage = native_storage,
      original_dim = dims,
      active_dim = active_dim,
      internal_transposed = transposed,
      internal_orientation = if (transposed) "transposed_wide_operator" else "as_given",
      rank = rank,
      target_kind = native_svd_target_kind(target),
      work = work,
      retained = retained,
      max_restarts = restart_budget,
      reorth_policy = reorth_policy,
      small_side_dimension = small_side,
      large_side_dimension = large_side,
      input_schema = list(
        initial_start = active_domain,
        retained_right_subspace = c(active_domain, retained),
        retained_left_subspace = c(active_codomain, retained),
        bidiagonal_alpha = work,
        bidiagonal_beta = work,
        restart_random_tail = c(active_domain, max(0L, work - retained)),
        operator = native_storage
      ),
      retained_state = c(
        "right Ritz subspace in the active operator orientation",
        "left Ritz subspace in the active operator orientation",
        "projected bidiagonal recurrence alpha/beta",
        "small SVD Ritz extraction state",
        "one-sided small-side orthogonalization state",
        "exact two-sided certificate in original coordinates"
      ),
      output_schema = c(
        "d", "u", "v", "certificate diagnostics", "iterations",
        "matvecs", "restart_count", "attempt_history", "stage_seconds",
        "orientation metadata"
      ),
      invariants = c(
        "wide operators run internally on A^T so the reorthogonalized side is the smaller dimension",
        "retained Ritz left and right subspaces are rotated together after each small SVD",
        "restart never discards the coupled bidiagonal recurrence without recording a fallback",
        "one-sided reorthogonalization is a speed policy only and is never a certificate",
        "if one-sided orthogonality checks fail, the implementation must switch to full or monitored reorthogonalization",
        "returned triplets are sorted and certified in the original operator coordinates",
        "failed small-work attempts are retained as restart state, not thrown away and rerun from scratch"
      )
    ),
    class = "eigencore_irlba_lbd_restart_abi"
  )
}

#' @keywords internal
native_irlba_lbd_deterministic_matrix <- function(n, cols, offset = 0L) {
  if (cols < 1L) {
    return(matrix(numeric(0), nrow = n, ncol = 0L))
  }
  idx <- seq_len(n * cols)
  matrix(
    sin((idx + offset) * 0.754877666) +
      cos((idx + 3L * offset + 11L) * 0.569840291),
    nrow = n,
    ncol = cols
  )
}

#' @keywords internal
native_irlba_lbd_pad_subspace <- function(X, n, cols, offset = 0L) {
  X <- as.matrix(X)
  if (!identical(nrow(X), as.integer(n))) {
    stop("retained scout subspace has the wrong active dimension.", call. = FALSE)
  }
  tol <- sqrt(.Machine$double.eps)
  Q <- matrix(numeric(0), nrow = n, ncol = 0L)
  append_candidate <- function(candidate) {
    z <- as.numeric(candidate)
    if (ncol(Q)) {
      z <- z - Q %*% as.numeric(crossprod(Q, z))
    }
    z_norm <- sqrt(sum(z^2))
    if (!is.finite(z_norm) || z_norm <= tol) {
      return(FALSE)
    }
    Q <<- cbind(Q, z / z_norm)
    TRUE
  }
  keep <- min(ncol(X), cols)
  if (keep > 0L) {
    for (col in seq_len(keep)) {
      append_candidate(X[, col])
    }
  }
  candidate <- 1L
  max_candidates <- max(10L * cols, n + 4L * cols)
  while (ncol(Q) < cols && candidate <= max_candidates) {
    if (candidate <= 4L * cols) {
      append_candidate(native_irlba_lbd_deterministic_matrix(
        n, 1L, offset = offset + candidate
      )[, 1L])
    } else {
      unit <- numeric(n)
      unit[((candidate - 4L * cols - 1L) %% n) + 1L] <- 1
      append_candidate(unit)
    }
    candidate <- candidate + 1L
  }
  if (ncol(Q) < cols) {
    stop("retained scout subspace could not be padded to full rank.", call. = FALSE)
  }
  Q[, seq_len(cols), drop = FALSE]
}

#' @keywords internal
native_irlba_lbd_retained_state_from_scout <- function(op, scout, abi = NULL,
                                                       target = largest(),
                                                       rank = NULL,
                                                       work = NULL,
                                                       retained = NULL,
                                                       max_restarts = NULL) {
  op <- as_operator(op)
  rank <- as.integer(rank %||% length(scout$d %||% scout$values))
  if (length(rank) != 1L || is.na(rank) || rank < 1L) {
    stop("rank must be supplied or inferable from the scout result.", call. = FALSE)
  }
  abi <- abi %||% native_irlba_lbd_restart_abi(
    op,
    rank = rank,
    target = target,
    work = work,
    retained = retained,
    max_restarts = max_restarts
  )
  if (is.null(scout$u) || is.null(scout$v)) {
    stop("retained IRLBA/LBD scout state requires both left and right scout vectors.",
         call. = FALSE)
  }

  active_right <- if (isTRUE(abi$internal_transposed)) scout$u else scout$v
  active_left <- if (isTRUE(abi$internal_transposed)) scout$v else scout$u
  retained_right <- native_irlba_lbd_pad_subspace(
    active_right,
    abi$input_schema$retained_right_subspace[[1L]],
    abi$retained,
    offset = 17L
  )
  retained_left <- native_irlba_lbd_pad_subspace(
    active_left,
    abi$input_schema$retained_left_subspace[[1L]],
    abi$retained,
    offset = 29L
  )
  tail_cols <- abi$input_schema$restart_random_tail[[2L]]
  restart_random_tail <- native_irlba_lbd_pad_subspace(
    matrix(numeric(0), abi$input_schema$restart_random_tail[[1L]], 0L),
    abi$input_schema$restart_random_tail[[1L]],
    tail_cols,
    offset = 41L
  )
  retained_from_scout <- min(ncol(active_right), ncol(active_left), abi$retained)
  structure(
    list(
      abi = abi,
      initial_start = retained_right[, 1L],
      retained_right_subspace = retained_right,
      retained_left_subspace = retained_left,
      alpha = numeric(abi$work),
      beta = numeric(abi$work),
      restart_random_tail = restart_random_tail,
      retained_from_scout = retained_from_scout,
      retained_padding = abi$retained - retained_from_scout,
      recurrence_available = FALSE,
      restart_state_kind = "ritz_subspace_only",
      internal_transposed = abi$internal_transposed,
      internal_orientation = abi$internal_orientation
    ),
    class = "eigencore_irlba_lbd_retained_state"
  )
}

#' @keywords internal
native_irlba_lbd_lock_diagnostics <- function(certificate, rank,
                                              source = "exact_final_certificate",
                                              fallback_reason = NA_character_) {
  converged <- certificate$converged %||% logical(0)
  hard_locked <- min(as.integer(rank), sum(isTRUE(certificate$passed) & converged))
  lock_complete <- hard_locked >= as.integer(rank)
  list(
    retained_locked_count = hard_locked,
    retained_deflation = FALSE,
    irlba_lbd_locking_policy =
      "hard locks require exact two-sided final SVD certificate",
    irlba_lbd_lock_source = source,
    irlba_lbd_soft_locked_count = 0L,
    irlba_lbd_hard_locked_count = hard_locked,
    irlba_lbd_locked_triplets_certified = isTRUE(lock_complete),
    irlba_lbd_locked_orthogonality_loss =
      certificate$max_orthogonality_loss %||% NA_real_,
    irlba_lbd_future_vectors_orthogonal_to_locks =
      isTRUE(certificate$orthogonality_passed),
    irlba_lbd_lock_fallback_reason = if (isTRUE(lock_complete)) {
      NA_character_
    } else {
      fallback_reason %||% "exact final certificate did not converge all requested triplets"
    }
  )
}

#' @keywords internal
native_irlba_lbd_attach_lock_diagnostics <- function(restart, certificate, rank,
                                                     source = "exact_final_certificate",
                                                     fallback_reason = NA_character_) {
  c(restart, native_irlba_lbd_lock_diagnostics(
    certificate = certificate,
    rank = rank,
    source = source,
    fallback_reason = fallback_reason
  ))
}

#' @keywords internal
native_irlba_lbd_bpro_guard_diagnostics <- function(abi, native, certificate,
                                                    fallback_used = FALSE,
                                                    fallback_reason = NA_character_) {
  mode <- abi$reorth_policy %||% NA_character_
  bpro_mode <- mode %in% c(
    "bpro_two_sided",
    "bpro_one_sided_guarded",
    "bpro_block_guarded"
  )
  reorthogonalize_u <- if (is.null(native)) {
    NA
  } else {
    as.logical(native$reorthogonalize_u)
  }
  reorthogonalize_v <- if (is.null(native)) {
    NA
  } else {
    as.logical(native$reorthogonalize_v)
  }
  one_sided_used <- isTRUE(bpro_mode) &&
    identical(mode, "bpro_one_sided_guarded") &&
    !isTRUE(fallback_used) &&
    xor(isTRUE(reorthogonalize_u), isTRUE(reorthogonalize_v))
  block_size <- if (!isTRUE(bpro_mode)) {
    NA_integer_
  } else if (identical(mode, "bpro_block_guarded")) {
    as.integer(min(abi$rank, abi$retained))
  } else {
    1L
  }
  guard_reason <- if (isTRUE(fallback_used)) {
    fallback_reason %||% "guarded BPRO native attempt fell back"
  } else if (isTRUE(bpro_mode) && !isTRUE(certificate$orthogonality_passed)) {
    "exact SVD orthogonality guard failed"
  } else {
    NA_character_
  }
  list(
    irlba_lbd_reorth_mode = mode,
    irlba_lbd_one_sided_reorth_used = isTRUE(one_sided_used),
    irlba_lbd_bpro_block_size = block_size,
    irlba_lbd_bpro_exact_orthogonality_loss =
      certificate$max_orthogonality_loss %||% NA_real_,
    irlba_lbd_bpro_exact_orthogonality_passed =
      isTRUE(certificate$orthogonality_passed),
    irlba_lbd_bpro_guard_fallback_reason = guard_reason
  )
}

#' @keywords internal
native_irlba_lbd_attach_bpro_guard_diagnostics <- function(restart, abi, native,
                                                           certificate,
                                                           fallback_used = FALSE,
                                                           fallback_reason = NA_character_) {
  c(restart, native_irlba_lbd_bpro_guard_diagnostics(
    abi = abi,
    native = native,
    certificate = certificate,
    fallback_used = fallback_used,
    fallback_reason = fallback_reason
  ))
}

#' @keywords internal
native_svd_certificate_from_diagnostics <- function(diagnostics, tol,
                                                    norm_info,
                                                    swap_sides = FALSE) {
  if (is.null(diagnostics)) {
    return(NULL)
  }
  left <- diagnostics$left
  right <- diagnostics$right
  orthogonality <- diagnostics$orthogonality
  if (isTRUE(swap_sides)) {
    left <- diagnostics$right
    right <- diagnostics$left
    if (length(orthogonality) >= 2L) {
      orthogonality <- c(
        U = unname(orthogonality[[2L]]),
        V = unname(orthogonality[[1L]])
      )
    }
  }
  new_certificate(
    tol = tol,
    residuals = list(
      left = left,
      right = right,
      combined = diagnostics$combined
    ),
    backward_error = diagnostics$backward_error,
    orthogonality = orthogonality,
    converged = diagnostics$converged,
    scale = diagnostics$scale,
    norm_bound_type = norm_info$norm_bound_type %||% "unspecified",
    scale_is_estimate = isTRUE(norm_info$scale_is_estimate)
  )
}

#' @keywords internal
native_irlba_lbd_retained_svd <- function(op, rank, target = largest(),
                                          tol = 1e-8,
                                          work = NULL,
                                          max_restarts = NULL,
                                          retained = NULL,
                                          vectors = c("both", "left", "right", "none"),
                                          reorth_policy = c(
                                            "one_sided_small_side",
                                            "full_two_sided",
                                            "bpro_two_sided",
                                            "bpro_one_sided_guarded",
                                            "bpro_block_guarded"
                                          )) {
  vectors <- match.arg(vectors)
  reorth_policy <- match.arg(reorth_policy)
  abi <- native_irlba_lbd_restart_abi(
    op,
    rank = rank,
    target = target,
    work = work,
    max_restarts = max_restarts,
    retained = retained,
    reorth_policy = reorth_policy
  )
  abi$tolerance <- as.numeric(tol)
  abi$vectors <- vectors

  original_op <- as_operator(op)
  active_source <- if (isTRUE(abi$internal_transposed)) {
    native_golub_kahan_transpose_source(original_op)
  } else {
    source_or_null(original_op) %||% original_op$metadata$matrix %||% NULL
  }
  if (is.null(active_source)) {
    stop("Native retained one-sided IRLBA/LBD requires a native dense or dgCMatrix source.",
         call. = FALSE)
  }
  active_op <- as_operator(active_source)
  active_dim <- as.integer(active_op$dim)
  active_domain <- active_dim[[2L]]
  active_codomain <- active_dim[[1L]]

  initial_start <- stats::rnorm(active_domain)
  initial_start <- initial_start / sqrt(sum(initial_start^2))
  small <- native_golub_kahan_svd(
    active_op,
    rank = rank,
    target = target,
    tol = tol,
    maxit = abi$work,
    vectors = "both",
    reorthogonalize = identical(reorth_policy, "full_two_sided"),
    internal_start = initial_start
  )
  if (isTRUE(small$certificate$passed)) {
    final <- small
    if (isTRUE(abi$internal_transposed)) {
      final <- native_golub_kahan_swap_transposed_result(original_op, final, tol)
    }
    final$restart$irlba_lbd_policy <-
      "retained one-sided LBD scout certified before native restart"
    final$restart$irlba_lbd_retained_native_attempted <- FALSE
    final$restart$irlba_lbd_retained_native_fallback_reason <- NA_character_
    final$restart$retained_restart <- FALSE
    final$restart$retained_restart_native <- FALSE
    final$restart$retained_restart_abi_version <- abi$version
    final$restart$native_attempt_certification <- FALSE
    final$restart$irlba_lbd_restart_state_kind <- "scout_certified_no_restart"
    final$restart$irlba_lbd_recurrence_available <- FALSE
    final$restart$irlba_lbd_augmented_recurrence <- FALSE
    final$restart$irlba_lbd_retained_seed_strategy <- NA_character_
    final$restart$irlba_lbd_retained_from_scout <- rank
    final$restart$irlba_lbd_retained_padding <- 0L
    final$restart$irlba_lbd_retained_fixed_work_attempts <- 0L
    final$restart$work <- abi$work
    final$restart$retained <- abi$retained
    final$restart$internal_orientation <- abi$internal_orientation
    final$restart$internal_transposed <- abi$internal_transposed
    final$restart$irlba_lbd_scout_matvecs <- small$matvecs
    final$restart$irlba_lbd_scout_accounted_seconds <-
      sum(small$stage_seconds %||% NA_real_, na.rm = TRUE)
    final$restart$irlba_lbd_scout_certificate_passed <- TRUE
    final$restart$attempted_subspaces <- abi$work
    final$restart$certified_attempt <- 1L
    final$restart$fallback_attempted <- FALSE
    final$restart$fallback_used <- FALSE
    final$restart$fallback_method <- NA_character_
    final$restart$attempt_history <- data.frame(
      attempt = 1L,
      max_subspace = abi$work,
      iterations = small$iterations,
      matvecs = small$matvecs,
      accounted_seconds = final$restart$irlba_lbd_scout_accounted_seconds,
      warm_started = FALSE,
      certificate_passed = isTRUE(final$certificate$passed),
      max_backward_error = final$certificate$max_backward_error,
      max_residual = final$certificate$max_residual,
      stringsAsFactors = FALSE
    )
    final$restart <- native_irlba_lbd_attach_lock_diagnostics(
      final$restart,
      final$certificate,
      rank = rank,
      source = "exact_scout_certificate"
    )
    final$restart <- native_irlba_lbd_attach_bpro_guard_diagnostics(
      final$restart,
      abi = abi,
      native = NULL,
      certificate = final$certificate
    )
    return(native_irlba_lbd_select_vectors(final, vectors))
  }
  retained_right <- native_irlba_lbd_pad_basis(
    small$v,
    n_rows = active_domain,
    cols = abi$retained,
    tol = tol
  )
  retained_left <- native_irlba_lbd_pad_basis(
    small$u,
    n_rows = active_codomain,
    cols = abi$retained,
    tol = tol
  )
  retained_from_scout <- min(ncol(small$v), ncol(small$u), abi$retained)
  retained_padding <- abi$retained - retained_from_scout
  alpha <- numeric(abi$work)
  beta <- numeric(abi$work)
  random_tails <- matrix(
    stats::rnorm(active_domain * max(0L, abi$work - abi$retained)),
    nrow = active_domain,
    ncol = max(0L, abi$work - abi$retained)
  )
  reorth_code <- match(
    reorth_policy,
    c(
      "one_sided_small_side",
      "full_two_sided",
      "bpro_two_sided",
      "bpro_one_sided_guarded",
      "bpro_block_guarded"
    )
  )

  native <- tryCatch(
    native_irlba_lbd_retained_call(
      active_source = active_source,
      initial_start = initial_start,
      retained_right = retained_right,
      retained_left = retained_left,
      alpha = alpha,
      beta = beta,
      random_tails = random_tails,
      work = abi$work,
      retained = abi$retained,
      max_restarts = abi$max_restarts,
      rank = abi$rank,
      target_kind = abi$target_kind,
      tol = tol,
      reorth_policy = reorth_code
    ),
    error = function(e) {
      structure(list(error = conditionMessage(e)), class = "eigencore_irlba_lbd_native_error")
    }
  )

  fallback_reason <- NA_character_
  fallback <- NULL
  if (inherits(native, "eigencore_irlba_lbd_native_error")) {
    fallback_reason <- native$error
  } else {
    native_certificate_reused <- FALSE
    native_certificate_swapped <- FALSE
    native_diag <- native[["certificate_diagnostics", exact = TRUE]]
    cert <- if (is.null(native_diag)) {
      NULL
    } else {
      native_svd_certificate_from_diagnostics(
        native_diag,
        tol = tol,
        norm_info = operator_norm_for_certificate_info(active_op)
      )
    }
    if (is.null(cert)) {
      cert <- certify_svd_operator(active_op, native$d, native$u, native$v, tol = tol)
    } else {
      native_certificate_reused <- TRUE
    }
    final <- list(
      d = native$d,
      u = native$u,
      v = native$v,
      values = native$d,
      residuals = cert$residuals,
      backward_error = cert$backward_error,
      orthogonality = cert$orthogonality,
      certificate = cert
    )
    final <- complete_zero_singular_triplets(active_op, final, rank, target, tol)
    if (isTRUE(final$zero_singular_completion)) {
      native_diag <- NULL
      native_certificate_reused <- FALSE
    }
    if (isTRUE(abi$internal_transposed)) {
      native_certificate_swapped <- !is.null(native_diag)
      final <- native_golub_kahan_swap_transposed_result(
        original_op,
        final,
        tol,
        certificate_diagnostics = native_diag
      )
    }
    final$iterations <- native$iterations + small$iterations
    final$matvecs <- native$matvecs + small$matvecs
    final$restarts <- native$restart_count
    final$stage_seconds <- c(
      scout = sum(small$stage_seconds %||% NA_real_, na.rm = TRUE),
      native_iteration = sum(
        native$stage_apply_seconds,
        native$stage_recurrence_seconds,
        native$stage_reorthogonalization_seconds,
        native$stage_projected_solve_seconds,
        na.rm = TRUE
      ),
      apply = native$stage_apply_seconds,
      recurrence = native$stage_recurrence_seconds,
      reorthogonalization = native$stage_reorthogonalization_seconds,
      projected_solve = native$stage_projected_solve_seconds,
      ritz = 0,
      retry_overhead = 0
    )
    final$restart <- native_irlba_lbd_restart_diagnostics(
      abi = abi,
      native = native,
      small = small,
      final = final,
      fallback_attempted = FALSE,
      fallback_used = FALSE,
      fallback_reason = NA_character_,
      native_certificate_reused = native_certificate_reused,
      native_certificate_swapped = native_certificate_swapped
    )
    if (isTRUE(final$certificate$passed)) {
      return(native_irlba_lbd_select_vectors(final, vectors))
    }
    fallback_reason <- paste0(
      "retained native IRLBA/LBD certificate failed: max backward error ",
      signif(final$certificate$max_backward_error, 4)
    )
  }

  fallback_reorthogonalize <- !identical(reorth_policy, "one_sided_small_side")
  fallback_expects_original_domain <- isTRUE(fallback_reorthogonalize)
  warm_start <- if (!inherits(native, "eigencore_irlba_lbd_native_error")) {
    if (isTRUE(abi$internal_transposed) && isTRUE(fallback_expects_original_domain) &&
        !is.null(native$u) && ncol(native$u) > 0L) {
      native$u[, 1L]
    } else if ((!isTRUE(abi$internal_transposed) || !isTRUE(fallback_expects_original_domain)) &&
        !is.null(native$v) && ncol(native$v) > 0L) {
      native$v[, 1L]
    } else {
      NULL
    }
  } else {
    NULL
  }
  if (is.null(warm_start)) {
    warm_start <- if (isTRUE(abi$internal_transposed) && isTRUE(fallback_expects_original_domain)) {
      retained_left[, 1L]
    } else {
      retained_right[, 1L]
    }
  }
  fallback <- native_golub_kahan_svd(
    original_op,
    rank = rank,
    target = target,
    tol = tol,
    vectors = "both",
    reorthogonalize = fallback_reorthogonalize,
    internal_start = warm_start
  )
  scout_matvecs <- small$matvecs %||% 0L
  scout_iterations <- small$iterations %||% 0L
  retained_matvecs <- if (!inherits(native, "eigencore_irlba_lbd_native_error")) {
      native$matvecs %||% 0L
  } else {
    0L
  }
  retained_iterations <- if (!inherits(native, "eigencore_irlba_lbd_native_error")) {
    native$iterations %||% 0L
  } else {
    0L
  }
  fallback_matvecs <- fallback$matvecs %||% 0L
  fallback_iterations <- fallback$iterations %||% 0L
  total_matvecs <- scout_matvecs + retained_matvecs + fallback_matvecs
  total_iterations <- scout_iterations + retained_iterations + fallback_iterations
  fallback$matvecs <- total_matvecs
  fallback$iterations <- total_iterations
  fallback$restart$irlba_lbd_policy <-
    "retained one-sided LBD native core with certified adaptive fallback"
  fallback$restart$irlba_lbd_retained_native_attempted <- TRUE
  fallback$restart$irlba_lbd_retained_native_fallback_reason <- fallback_reason
  fallback$restart$fallback_attempted <- TRUE
  fallback$restart$fallback_used <- TRUE
  fallback$restart$fallback_method <- "adaptive one-sided Golub-Kahan"
  fallback$restart$retained_restart <- TRUE
  fallback$restart$retained_restart_native <- !inherits(native, "eigencore_irlba_lbd_native_error")
  fallback$restart$retained_restart_abi_version <- abi$version
  fallback$restart$native_attempt_certification <- !inherits(native, "eigencore_irlba_lbd_native_error")
  fallback$restart$irlba_lbd_restart_state_kind <- if (!inherits(native, "eigencore_irlba_lbd_native_error")) {
    native$restart_state_kind %||% "ritz_subspace_only"
  } else {
    "ritz_subspace_only"
  }
  fallback$restart$irlba_lbd_recurrence_available <- if (!inherits(native, "eigencore_irlba_lbd_native_error")) {
    isTRUE(native$recurrence_available)
  } else {
    FALSE
  }
  fallback$restart$irlba_lbd_augmented_recurrence <- if (!inherits(native, "eigencore_irlba_lbd_native_error")) {
    isTRUE(native$augmented_recurrence)
  } else {
    FALSE
  }
  fallback$restart$irlba_lbd_retained_seed_strategy <- if (identical(
    fallback$restart$irlba_lbd_restart_state_kind,
    "residual_augmented_projection"
  )) {
    "ritz_residual_augmented_krylov_projection"
  } else {
    "ritz_subspace_seeded_fixed_work"
  }
  fallback$restart$irlba_lbd_retained_from_scout <- retained_from_scout
  fallback$restart$irlba_lbd_retained_padding <- retained_padding
  fallback$restart$irlba_lbd_retained_fixed_work_attempts <-
    if (!inherits(native, "eigencore_irlba_lbd_native_error")) {
      if (identical(
        fallback$restart$irlba_lbd_restart_state_kind,
        "residual_augmented_projection"
      )) {
        0L
      } else {
        nrow(native$attempt_history)
      }
    } else {
      0L
    }
  fallback$restart$irlba_lbd_residual_augmented_cols <-
    if (!inherits(native, "eigencore_irlba_lbd_native_error")) {
      native$residual_augmented_cols %||% NA_integer_
    } else {
      NA_integer_
    }
  fallback$restart$irlba_lbd_augmented_tail_steps <-
    if (!inherits(native, "eigencore_irlba_lbd_native_error")) {
      native$augmented_tail_steps %||% NA_integer_
    } else {
      NA_integer_
    }
  fallback$restart$irlba_lbd_augmented_basis_cols <-
    if (!inherits(native, "eigencore_irlba_lbd_native_error")) {
      native$augmented_basis_cols %||% NA_integer_
    } else {
      NA_integer_
    }
  fallback$restart$irlba_lbd_augmented_restart_cycles <-
    if (!inherits(native, "eigencore_irlba_lbd_native_error")) {
      native$augmented_restart_cycles %||% NA_integer_
    } else {
      NA_integer_
    }
  fallback$restart$irlba_lbd_augmented_kept_vectors <-
    if (!inherits(native, "eigencore_irlba_lbd_native_error")) {
      native$augmented_kept_vectors %||% NA_integer_
    } else {
      NA_integer_
    }
  fallback$restart$irlba_lbd_augmented_small_svds <-
    if (!inherits(native, "eigencore_irlba_lbd_native_error")) {
      native$augmented_small_svds %||% NA_integer_
    } else {
      NA_integer_
    }
  fallback$restart$irlba_lbd_augmented_cached_aq_cols <-
    if (!inherits(native, "eigencore_irlba_lbd_native_error")) {
      native$augmented_cached_aq_cols %||% NA_integer_
    } else {
      NA_integer_
    }
  fallback$restart$irlba_lbd_augmented_from_scratch_matvecs <-
    if (!inherits(native, "eigencore_irlba_lbd_native_error")) {
      native$augmented_from_scratch_matvecs %||% NA_integer_
    } else {
      NA_integer_
    }
  fallback$restart$irlba_lbd_augmented_matvec_savings <-
    if (!inherits(native, "eigencore_irlba_lbd_native_error")) {
      native$augmented_matvec_savings %||% NA_integer_
    } else {
      NA_integer_
    }
  fallback$restart$irlba_lbd_augmented_min_cheap_residual <-
    if (!inherits(native, "eigencore_irlba_lbd_native_error")) {
      native$augmented_min_cheap_residual %||% NA_real_
    } else {
      NA_real_
    }
  fallback$restart$irlba_lbd_augmented_final_cheap_residual <-
    if (!inherits(native, "eigencore_irlba_lbd_native_error")) {
      native$augmented_final_cheap_residual %||% NA_real_
    } else {
      NA_real_
    }
  fallback$restart$irlba_lbd_augmented_reduces_from_scratch_work <-
    if (!inherits(native, "eigencore_irlba_lbd_native_error")) {
      isTRUE(native$augmented_reduces_from_scratch_work)
    } else {
      NA
    }
  fallback$restart$irlba_lbd_bpro_policy <-
    if (!inherits(native, "eigencore_irlba_lbd_native_error")) {
      isTRUE(native$bpro_policy)
    } else {
      identical(reorth_policy, "bpro_two_sided")
    }
  fallback$restart$irlba_lbd_bpro_passes_per_append <-
    if (!inherits(native, "eigencore_irlba_lbd_native_error")) {
      native$bpro_reorthogonalization_passes_per_append %||% NA_integer_
    } else {
      NA_integer_
    }
  fallback$restart$irlba_lbd_bpro_monitoring_threshold <-
    if (!inherits(native, "eigencore_irlba_lbd_native_error")) {
      native$bpro_monitoring_threshold %||% NA_real_
    } else {
      NA_real_
    }
  fallback$restart$irlba_lbd_bpro_monitored_appends <-
    if (!inherits(native, "eigencore_irlba_lbd_native_error")) {
      native$bpro_monitored_appends %||% NA_integer_
    } else {
      NA_integer_
    }
  fallback$restart$irlba_lbd_bpro_threshold_reorthogonalizations <-
    if (!inherits(native, "eigencore_irlba_lbd_native_error")) {
      native$bpro_threshold_reorthogonalizations %||% NA_integer_
    } else {
      NA_integer_
    }
  fallback$restart$irlba_lbd_bpro_max_estimated_orthogonality_loss <-
    if (!inherits(native, "eigencore_irlba_lbd_native_error")) {
      native$bpro_max_estimated_orthogonality_loss %||% NA_real_
    } else {
      NA_real_
    }
  fallback$restart$irlba_lbd_bpro_max_post_append_orthogonality_loss <-
    if (!inherits(native, "eigencore_irlba_lbd_native_error")) {
      native$bpro_max_post_append_orthogonality_loss %||% NA_real_
    } else {
      NA_real_
    }
  fallback$restart$irlba_lbd_bpro_basis_orthogonality_loss <-
    if (!inherits(native, "eigencore_irlba_lbd_native_error")) {
      native$bpro_augmented_basis_orthogonality_loss %||% NA_real_
    } else {
      NA_real_
    }
  fallback$restart$irlba_lbd_bpro_escalation_recommended <-
    if (!inherits(native, "eigencore_irlba_lbd_native_error")) {
      isTRUE(native$bpro_escalation_recommended) ||
        (!is.na(fallback$certificate$max_orthogonality_loss) &&
          fallback$certificate$max_orthogonality_loss >
            fallback$certificate$orthogonality_tolerance)
    } else {
      NA
    }
  fallback$restart$work <- abi$work
  fallback$restart$retained <- abi$retained
  fallback$restart$internal_orientation <- abi$internal_orientation
  fallback$restart$internal_transposed <- abi$internal_transposed
  fallback$restart$irlba_lbd_scout_matvecs <- scout_matvecs
  fallback$restart$irlba_lbd_scout_accounted_seconds <-
    sum(small$stage_seconds %||% NA_real_, na.rm = TRUE)
  fallback$restart$irlba_lbd_scout_certificate_passed <- FALSE
  fallback$restart$irlba_lbd_fallback_matvecs <- fallback_matvecs
  fallback$restart$irlba_lbd_fallback_iterations <- fallback_iterations
  fallback$restart$irlba_lbd_total_matvecs <- total_matvecs
  fallback$restart$irlba_lbd_total_iterations <- total_iterations
  if (!inherits(native, "eigencore_irlba_lbd_native_error")) {
    fallback$restart$irlba_lbd_retained_attempt_history <- native$attempt_history
    fallback$restart$irlba_lbd_retained_matvecs <- retained_matvecs
    fallback$restart$irlba_lbd_retained_iterations <- retained_iterations
  }
  fallback$restart <- native_irlba_lbd_attach_lock_diagnostics(
    fallback$restart,
    fallback$certificate,
    rank = rank,
    source = "exact_fallback_certificate",
    fallback_reason = fallback_reason
  )
  fallback$restart <- native_irlba_lbd_attach_bpro_guard_diagnostics(
    fallback$restart,
    abi = abi,
    native = if (!inherits(native, "eigencore_irlba_lbd_native_error")) native else NULL,
    certificate = fallback$certificate,
    fallback_used = TRUE,
    fallback_reason = fallback_reason
  )
  native_irlba_lbd_select_vectors(fallback, vectors)
}

#' @keywords internal
native_irlba_lbd_pad_basis <- function(x, n_rows, cols,
                                       tol = sqrt(.Machine$double.eps)) {
  x <- if (is.null(x)) matrix(0, n_rows, 0L) else as.matrix(x)
  if (nrow(x) != n_rows) {
    stop("retained IRLBA/LBD basis has non-conformable row count.", call. = FALSE)
  }
  if (ncol(x) > cols) {
    x <- x[, seq_len(cols), drop = FALSE]
  }
  if (ncol(x) < cols) {
    x <- cbind(
      x,
      orthonormal_completion(
        x,
        n_rows = n_rows,
        needed = cols - ncol(x),
        tol = tol
      )
    )
  }
  x
}

#' @keywords internal
native_irlba_lbd_retained_call <- function(active_source, initial_start,
                                           retained_right, retained_left,
                                           alpha, beta, random_tails,
                                           work, retained, max_restarts,
                                           rank, target_kind, tol,
                                           reorth_policy) {
  if (inherits(active_source, "dgCMatrix")) {
    return(.Call(
      "eigencore_irlba_lbd_csc_retained",
      methods::slot(active_source, "i"),
      methods::slot(active_source, "p"),
      methods::slot(active_source, "x"),
      methods::slot(active_source, "Dim"),
      as.numeric(initial_start),
      retained_right,
      retained_left,
      as.numeric(alpha),
      as.numeric(beta),
      random_tails,
      as.integer(work),
      as.integer(retained),
      as.integer(max_restarts),
      as.integer(rank),
      as.integer(target_kind),
      as.numeric(tol),
      as.integer(reorth_policy),
      PACKAGE = "eigencore"
    ))
  }
  if (is.matrix(active_source) && is.double(active_source)) {
    return(.Call(
      "eigencore_irlba_lbd_dense_retained",
      active_source,
      as.numeric(initial_start),
      retained_right,
      retained_left,
      as.numeric(alpha),
      as.numeric(beta),
      random_tails,
      as.integer(work),
      as.integer(retained),
      as.integer(max_restarts),
      as.integer(rank),
      as.integer(target_kind),
      as.numeric(tol),
      as.integer(reorth_policy),
      PACKAGE = "eigencore"
    ))
  }
  stop("native retained IRLBA/LBD supports only dense double and dgCMatrix sources.",
       call. = FALSE)
}

#' @keywords internal
native_irlba_lbd_restart_diagnostics <- function(abi, native, small, final,
                                                 fallback_attempted,
                                                 fallback_used,
                                                 fallback_reason,
                                                 native_certificate_reused = FALSE,
                                                 native_certificate_swapped = FALSE) {
  restart <- list(
    kind = "irlba_lbd_native_retained_core",
    implemented = TRUE,
    irlba_lbd_policy = "retained one-sided LBD native core",
    irlba_lbd_retained_native_attempted = TRUE,
    irlba_lbd_retained_native_fallback_reason = NA_character_,
    retained_restart = TRUE,
    retained_restart_native = TRUE,
    retained_restart_abi_version = abi$version,
    work = abi$work,
    retained = abi$retained,
    max_restarts = abi$max_restarts,
    restart_count = native$restart_count,
    attempts = nrow(native$attempt_history),
    attempted_subspaces = native$attempt_history$max_subspace,
    attempt_history = native$attempt_history,
    certified_attempt = if (isTRUE(final$certificate$passed)) nrow(native$attempt_history) else NA_integer_,
    native_attempt_certification = TRUE,
    native_workspace_bytes = native$native_workspace_bytes,
    reorthogonalization_mode = abi$reorth_policy,
    reorthogonalize_u = as.logical(native$reorthogonalize_u),
    reorthogonalize_v = as.logical(native$reorthogonalize_v),
    reorthogonalization_passes = native$reorthogonalization_passes,
    irlba_lbd_restart_state_kind = native$restart_state_kind %||% "ritz_subspace_only",
    irlba_lbd_recurrence_available = isTRUE(native$recurrence_available),
    irlba_lbd_augmented_recurrence = isTRUE(native$augmented_recurrence),
    irlba_lbd_retained_seed_strategy = if (identical(
      native$restart_state_kind %||% "ritz_subspace_only",
      "residual_augmented_projection"
    )) {
      "ritz_residual_augmented_krylov_projection"
    } else {
      "ritz_subspace_seeded_fixed_work"
    },
    irlba_lbd_retained_from_scout = min(ncol(small$v), ncol(small$u), abi$retained),
    irlba_lbd_retained_padding = abi$retained - min(ncol(small$v), ncol(small$u), abi$retained),
    irlba_lbd_retained_fixed_work_attempts = if (identical(
      native$restart_state_kind %||% "ritz_subspace_only",
      "residual_augmented_projection"
    )) {
      0L
    } else {
      nrow(native$attempt_history)
    },
    irlba_lbd_residual_augmented_cols = native$residual_augmented_cols %||% NA_integer_,
    irlba_lbd_augmented_tail_steps = native$augmented_tail_steps %||% NA_integer_,
    irlba_lbd_augmented_basis_cols = native$augmented_basis_cols %||% NA_integer_,
    irlba_lbd_augmented_restart_cycles =
      native$augmented_restart_cycles %||% NA_integer_,
    irlba_lbd_augmented_kept_vectors =
      native$augmented_kept_vectors %||% NA_integer_,
    irlba_lbd_augmented_small_svds =
      native$augmented_small_svds %||% NA_integer_,
    irlba_lbd_augmented_cached_aq_cols =
      native$augmented_cached_aq_cols %||% NA_integer_,
    irlba_lbd_augmented_from_scratch_matvecs =
      native$augmented_from_scratch_matvecs %||% NA_integer_,
    irlba_lbd_augmented_matvec_savings =
      native$augmented_matvec_savings %||% NA_integer_,
    irlba_lbd_augmented_min_cheap_residual =
      native$augmented_min_cheap_residual %||% NA_real_,
    irlba_lbd_augmented_final_cheap_residual =
      native$augmented_final_cheap_residual %||% NA_real_,
    irlba_lbd_augmented_reduces_from_scratch_work =
      isTRUE(native$augmented_reduces_from_scratch_work),
    irlba_lbd_native_certificate_diagnostics_reused =
      isTRUE(native_certificate_reused),
    irlba_lbd_native_certificate_diagnostics_swapped =
      isTRUE(native_certificate_swapped),
    irlba_lbd_bpro_policy = isTRUE(native$bpro_policy),
    irlba_lbd_bpro_passes_per_append =
      native$bpro_reorthogonalization_passes_per_append %||% NA_integer_,
    irlba_lbd_bpro_monitoring_threshold =
      native$bpro_monitoring_threshold %||% NA_real_,
    irlba_lbd_bpro_monitored_appends =
      native$bpro_monitored_appends %||% NA_integer_,
    irlba_lbd_bpro_threshold_reorthogonalizations =
      native$bpro_threshold_reorthogonalizations %||% NA_integer_,
    irlba_lbd_bpro_max_estimated_orthogonality_loss =
      native$bpro_max_estimated_orthogonality_loss %||% NA_real_,
    irlba_lbd_bpro_max_post_append_orthogonality_loss =
      native$bpro_max_post_append_orthogonality_loss %||% NA_real_,
    irlba_lbd_bpro_basis_orthogonality_loss =
      native$bpro_augmented_basis_orthogonality_loss %||% NA_real_,
    irlba_lbd_bpro_escalation_recommended =
      isTRUE(native$bpro_escalation_recommended) ||
        (!is.na(final$certificate$max_orthogonality_loss) &&
          final$certificate$max_orthogonality_loss >
            final$certificate$orthogonality_tolerance),
    internal_orientation = abi$internal_orientation,
    internal_transposed = abi$internal_transposed,
    scout_matvecs = small$matvecs,
    scout_iterations = small$iterations,
    scout_certificate_passed = isTRUE(small$certificate$passed),
    irlba_lbd_scout_matvecs = small$matvecs,
    irlba_lbd_scout_accounted_seconds =
      sum(small$stage_seconds %||% NA_real_, na.rm = TRUE),
    irlba_lbd_scout_certificate_passed = isTRUE(small$certificate$passed),
    irlba_lbd_retained_matvecs = native$matvecs %||% NA_integer_,
    irlba_lbd_retained_iterations = native$iterations %||% NA_integer_,
    irlba_lbd_total_matvecs = (small$matvecs %||% 0L) + (native$matvecs %||% 0L),
    irlba_lbd_total_iterations =
      (small$iterations %||% 0L) + (native$iterations %||% 0L),
    irlba_lbd_fallback_matvecs = NA_integer_,
    irlba_lbd_fallback_iterations = NA_integer_,
    fallback_attempted = fallback_attempted,
    fallback_used = fallback_used,
    fallback_method = if (fallback_used) "adaptive one-sided Golub-Kahan" else NA_character_,
    fallback_reason = fallback_reason,
    converged = isTRUE(final$certificate$passed),
    nconv = sum(final$certificate$converged),
    max_backward_error = final$certificate$max_backward_error,
    stage_seconds = final$stage_seconds,
    zero_singular_completion = isTRUE(final$zero_singular_completion),
    zero_singular_threshold = final$zero_singular_threshold %||% NA_real_,
    certified_in_original_coordinates = TRUE
  )
  restart <- native_irlba_lbd_attach_lock_diagnostics(
    restart,
    final$certificate,
    rank = abi$rank,
    source = "exact_retained_restart_certificate",
    fallback_reason = fallback_reason
  )
  native_irlba_lbd_attach_bpro_guard_diagnostics(
    restart,
    abi = abi,
    native = native,
    certificate = final$certificate,
    fallback_used = fallback_used,
    fallback_reason = fallback_reason
  )
}

#' @keywords internal
native_irlba_lbd_select_vectors <- function(final, vectors) {
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
native_golub_kahan_svd <- function(op, rank, target = largest(), tol = 1e-8,
                                   maxit = NULL,
                                   vectors = c("both", "left", "right", "none"),
                                   reorthogonalize = TRUE,
                                   internal_start = NULL) {
  vectors <- match.arg(vectors)
  op <- as_operator(op)
  if (is.null(op$apply_adjoint)) {
    stop("Native Golub-Kahan SVD requires an adjoint operator.", call. = FALSE)
  }

  m <- op$dim[1L]
  n <- op$dim[2L]
  limit <- min(m, n)
  fixed_maxit <- !is.null(maxit)
  if (!is.null(maxit)) {
    maxit <- min(limit, as.integer(maxit))
  }

  external_op <- op
  internal_transposed <- FALSE
  internal_orientation <- "as_given"
  if (!isTRUE(reorthogonalize) && m < n) {
    transposed_source <- native_golub_kahan_transpose_source(op)
    if (!is.null(transposed_source)) {
      op <- as_operator(transposed_source)
      m <- op$dim[1L]
      n <- op$dim[2L]
      internal_transposed <- TRUE
      internal_orientation <- "transposed_wide_operator"
    }
  }
  limit <- min(m, n)
  if (is.null(maxit)) {
    maxit <- default_golub_kahan_initial_subspace(
      c(m, n),
      rank,
      reorthogonalize = reorthogonalize
    )
  } else {
    maxit <- min(limit, as.integer(maxit))
  }
  if (maxit < rank) {
    stop("maxit/max_subspace must be at least rank.", call. = FALSE)
  }

  start <- if (is.null(internal_start)) {
    stats::rnorm(n)
  } else {
    internal_start <- as.numeric(internal_start)
    if (length(internal_start) != n || any(!is.finite(internal_start))) {
      stop("internal_start must be a finite numeric vector matching the active operator domain.", call. = FALSE)
    }
    start_norm <- sqrt(sum(internal_start^2))
    if (!is.finite(start_norm) || start_norm <= 100 * .Machine$double.eps) {
      stop("internal_start must have nonzero norm.", call. = FALSE)
    }
    internal_start / start_norm
  }
  storage <- op$metadata$storage %||% NULL
  source <- source_or_null(op)
  projected_stop_requested <- isTRUE(getOption("eigencore.golub_kahan_projected_stop", FALSE))
  projected_stop_disable_reason <- NULL
  if (identical(storage, "dgCMatrix") && m >= 4L * n) {
    projected_stop_disable_reason <- "disabled for high-aspect tall sparse operators"
  }
  projected_stop_enabled <- projected_stop_requested &&
    !fixed_maxit &&
    is.null(projected_stop_disable_reason)
  prefix_diagnostics <- isTRUE(getOption("eigencore.golub_kahan_prefix_diagnostics", FALSE))
  reorthogonalization_mode <- if (isTRUE(reorthogonalize)) {
    "full_two_sided"
  } else {
    "one_sided_small_side"
  }
  reorthogonalize_u <- isTRUE(reorthogonalize) || m <= n
  reorthogonalize_v <- isTRUE(reorthogonalize) || n <= m
  run_native <- function(active_maxit) {
    if (identical(storage, "dgCMatrix")) {
      A <- op$metadata$matrix
      if (isTRUE(prefix_diagnostics)) {
        .Call(
          "eigencore_golub_kahan_csc",
          methods::slot(A, "i"),
          methods::slot(A, "p"),
          methods::slot(A, "x"),
          methods::slot(A, "Dim"),
          as.integer(active_maxit),
          as.numeric(start),
          as.integer(rank),
          as.integer(native_svd_target_kind(target)),
          as.numeric(tol),
          as.logical(projected_stop_enabled),
          PACKAGE = "eigencore"
        )
      } else {
        .Call(
          "eigencore_golub_kahan_csc_fit",
          methods::slot(A, "i"),
          methods::slot(A, "p"),
          methods::slot(A, "x"),
          methods::slot(A, "Dim"),
          as.integer(active_maxit),
          as.numeric(start),
          as.integer(rank),
          as.integer(native_svd_target_kind(target)),
          as.numeric(tol),
          as.logical(projected_stop_enabled),
          as.logical(reorthogonalize_u),
          as.logical(reorthogonalize_v),
          PACKAGE = "eigencore"
        )
      }
    } else if (is.matrix(source) && is.double(source)) {
      if (isTRUE(prefix_diagnostics)) {
        .Call(
          "eigencore_golub_kahan_dense",
          source,
          as.integer(active_maxit),
          as.numeric(start),
          as.integer(rank),
          as.integer(native_svd_target_kind(target)),
          as.numeric(tol),
          as.logical(projected_stop_enabled),
          PACKAGE = "eigencore"
        )
      } else {
        .Call(
          "eigencore_golub_kahan_dense_fit",
          source,
          as.integer(active_maxit),
          as.numeric(start),
          as.integer(rank),
          as.integer(native_svd_target_kind(target)),
          as.numeric(tol),
          as.logical(projected_stop_enabled),
          as.logical(reorthogonalize_u),
          as.logical(reorthogonalize_v),
          PACKAGE = "eigencore"
        )
      }
    } else {
      stop("Native Golub-Kahan currently supports dense double matrices and dgCMatrix operators only.", call. = FALSE)
    }
  }

  final <- NULL
  iter <- NULL
  active_maxit <- maxit
  retries <- 0L
  history <- list()
  total_iterations <- 0L
  total_matvecs <- 0L
  total_native_seconds <- 0
  total_ritz_seconds <- 0
  total_native_stage_seconds <- c(
    apply = 0,
    recurrence = 0,
    reorthogonalization = 0,
    projected_solve = 0
  )
  total_reorthogonalization_passes <- 0L
  repeat {
    native_started <- proc.time()[["elapsed"]]
    iter <- run_native(active_maxit)
    native_elapsed <- proc.time()[["elapsed"]] - native_started
    total_native_seconds <- total_native_seconds + native_elapsed
    native_stage_seconds <- c(
      apply = iter$stage_apply_seconds %||% NA_real_,
      recurrence = iter$stage_recurrence_seconds %||% NA_real_,
      reorthogonalization = iter$stage_reorthogonalization_seconds %||% NA_real_,
      projected_solve = iter$projected_seconds %||% NA_real_
    )
    total_native_stage_seconds <- total_native_stage_seconds +
      replace(native_stage_seconds, is.na(native_stage_seconds), 0)
    native_reorthogonalization_passes <- iter$reorthogonalization_passes %||% NA_integer_
    if (!is.na(native_reorthogonalization_passes)) {
      total_reorthogonalization_passes <- total_reorthogonalization_passes +
        as.integer(native_reorthogonalization_passes)
    }

    ritz_started <- proc.time()[["elapsed"]]
    final <- if (!is.null(iter$d) && !is.null(iter$u) && !is.null(iter$v)) {
      cert <- certify_svd_operator(op, iter$d, iter$u, iter$v, tol = tol)
      list(
        d = iter$d,
        u = iter$u,
        v = iter$v,
        values = iter$d,
        residuals = cert$residuals,
        backward_error = cert$backward_error,
        orthogonality = cert$orthogonality,
        certificate = cert
      )
    } else {
      native_golub_kahan_ritz(
        op,
        iter$U,
        iter$V,
        iter$alpha,
        iter$beta,
        rank,
        target,
        tol,
        active_iterations = iter$iterations
      )
    }
    final <- complete_zero_singular_triplets(op, final, rank, target, tol)
    ritz_elapsed <- proc.time()[["elapsed"]] - ritz_started
    total_ritz_seconds <- total_ritz_seconds + ritz_elapsed
    total_iterations <- total_iterations + iter$iterations
    total_matvecs <- total_matvecs + iter$matvecs
    history[[length(history) + 1L]] <- data.frame(
      retry = retries,
      max_subspace = active_maxit,
      iterations = iter$iterations,
      matvecs = iter$matvecs,
      cumulative_iterations = total_iterations,
      cumulative_matvecs = total_matvecs,
      native_seconds = native_elapsed,
      ritz_seconds = ritz_elapsed,
      nconv = sum(final$certificate$converged),
      certificate_passed = isTRUE(final$certificate$passed),
      max_residual = final$certificate$max_residual,
      max_backward_error = final$certificate$max_backward_error,
      stage_apply_seconds = native_stage_seconds[["apply"]],
      stage_recurrence_seconds = native_stage_seconds[["recurrence"]],
      stage_reorthogonalization_seconds = native_stage_seconds[["reorthogonalization"]],
      stage_projected_solve_seconds = native_stage_seconds[["projected_solve"]],
      reorthogonalization_passes = native_reorthogonalization_passes
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

  prefix_history <- if (isTRUE(prefix_diagnostics)) {
    native_golub_kahan_prefix_diagnostics(
      op,
      iter = iter,
      rank = rank,
      target = target,
      tol = tol
    )
  } else {
    data.frame(
      prefix_iterations = integer(0),
      prefix_matvecs = integer(0),
      nconv = integer(0),
      certificate_passed = logical(0),
      max_residual = numeric(0),
      max_backward_error = numeric(0),
      error = character(0),
      stringsAsFactors = FALSE
    )
  }
  certified_prefixes <- prefix_history$prefix_iterations[prefix_history$certificate_passed]
  first_certified_prefix <- if (length(certified_prefixes)) {
    min(certified_prefixes)
  } else {
    NA_integer_
  }
  final_prefix_overshoot <- if (is.na(first_certified_prefix)) {
    NA_integer_
  } else {
    max(0L, iter$iterations - first_certified_prefix)
  }

  if (isTRUE(internal_transposed)) {
    final <- native_golub_kahan_swap_transposed_result(external_op, final, tol)
  }

  if (vectors == "left") {
    final$v <- NULL
  } else if (vectors == "right") {
    final$u <- NULL
  } else if (vectors == "none") {
    final$u <- NULL
    final$v <- NULL
  }
  final$iterations <- total_iterations
  final$matvecs <- total_matvecs
  final$restarts <- retries
  final$stage_seconds <- c(
    native_iteration = total_native_seconds,
    apply = total_native_stage_seconds[["apply"]],
    recurrence = total_native_stage_seconds[["recurrence"]],
    reorthogonalization = total_native_stage_seconds[["reorthogonalization"]],
    projected_solve = total_native_stage_seconds[["projected_solve"]],
    ritz = total_ritz_seconds,
    retry_overhead = 0
  )
  final$convergence_history <- do.call(rbind, history)
  final$restart <- list(
    kind = "adaptive_subspace_growth",
    implemented = TRUE,
    ritz_native = TRUE,
    restart_policy = "grow subspace until certificate convergence or limit",
    retries = retries,
    attempts = retries + 1L,
    initial_max_subspace = maxit,
    final_max_subspace = active_maxit,
    fixed_max_subspace = fixed_maxit,
    converged = all(final$certificate$converged),
    nconv = sum(final$certificate$converged),
    max_backward_error = final$certificate$max_backward_error,
    total_iterations = total_iterations,
    total_matvecs = total_matvecs,
    final_iterations = iter$iterations,
    final_matvecs = iter$matvecs,
    projected_stop_enabled = projected_stop_enabled,
    projected_stop_requested = projected_stop_requested,
    projected_stop_disable_reason = projected_stop_disable_reason %||% NA_character_,
    projected_stop = isTRUE(iter$projected_stop),
    projected_nconv = iter$projected_nconv %||% NA_integer_,
    projected_max_residual = iter$projected_max_residual %||% NA_real_,
    projected_checks = iter$projected_checks %||% NA_integer_,
    projected_seconds = iter$projected_seconds %||% NA_real_,
    native_workspace_bytes = iter$native_workspace_bytes %||% NA_real_,
    basis_returned = isTRUE(iter$basis_returned %||% (!is.null(iter$U) && !is.null(iter$V))),
    reorthogonalization_passes = total_reorthogonalization_passes,
    reorthogonalization_mode = reorthogonalization_mode,
    reorthogonalize_u = reorthogonalize_u,
    reorthogonalize_v = reorthogonalize_v,
    warm_started = !is.null(internal_start),
    internal_orientation = internal_orientation,
    internal_transposed = internal_transposed,
    zero_singular_completion = isTRUE(final$zero_singular_completion),
    zero_singular_threshold = final$zero_singular_threshold %||% NA_real_,
    prefix_diagnostics = prefix_diagnostics,
    prefix_history = prefix_history,
    first_certified_prefix = first_certified_prefix,
    final_prefix_iteration_overshoot = final_prefix_overshoot,
    final_prefix_matvec_overshoot = if (is.na(final_prefix_overshoot)) NA_integer_ else 2L * final_prefix_overshoot,
    stage_seconds = final$stage_seconds,
    history = final$convergence_history
  )
  final
}

#' @keywords internal
native_golub_kahan_transpose_source <- function(op) {
  storage <- op$metadata$storage %||% NULL
  source <- source_or_null(op)
  if (identical(storage, "dgCMatrix")) {
    return(get("t", envir = asNamespace("Matrix"))(op$metadata$matrix))
  }
  if (is.matrix(source) && is.double(source)) {
    return(t(source))
  }
  NULL
}

#' @keywords internal
native_golub_kahan_swap_transposed_result <- function(original_op, final, tol,
                                                      certificate_diagnostics = NULL) {
  if (is.null(final$u) || is.null(final$v)) {
    return(final)
  }
  u <- final$v
  v <- final$u
  cert <- native_svd_certificate_from_diagnostics(
    certificate_diagnostics,
    tol = tol,
    norm_info = final$certificate,
    swap_sides = TRUE
  )
  if (is.null(cert)) {
    cert <- certify_svd_operator(original_op, final$d, u, v, tol = tol)
  }
  final$u <- u
  final$v <- v
  final$residuals <- cert$residuals
  final$backward_error <- cert$backward_error
  final$orthogonality <- cert$orthogonality
  final$certificate <- cert
  final
}

#' @keywords internal
native_golub_kahan_prefix_diagnostics <- function(op, iter, rank, target, tol) {
  iterations <- as.integer(iter$iterations %||% 0L)
  rank <- as.integer(rank)
  if (iterations < rank) {
    return(data.frame(
      prefix_iterations = integer(0),
      prefix_matvecs = integer(0),
      nconv = integer(0),
      certificate_passed = logical(0),
      max_residual = numeric(0),
      max_backward_error = numeric(0)
    ))
  }

  candidate_prefixes <- unique(c(
    rank,
    2L * rank,
    4L * rank,
    4L * rank + 20L,
    8L * rank + 20L,
    iterations
  ))
  prefixes <- sort(unique(as.integer(candidate_prefixes[
    candidate_prefixes >= rank & candidate_prefixes <= iterations
  ])))

  rows <- lapply(prefixes, function(prefix) {
    fit <- tryCatch(
      native_golub_kahan_ritz(
        op,
        iter$U,
        iter$V,
        iter$alpha,
        iter$beta,
        rank,
        target,
        tol,
        active_iterations = prefix
      ),
      error = function(e) {
        structure(list(error = conditionMessage(e)), class = "eigencore_prefix_error")
      }
    )
    if (inherits(fit, "eigencore_prefix_error")) {
      return(data.frame(
        prefix_iterations = prefix,
        prefix_matvecs = 2L * prefix,
        nconv = 0L,
        certificate_passed = FALSE,
        max_residual = Inf,
        max_backward_error = Inf,
        error = fit$error,
        stringsAsFactors = FALSE
      ))
    }
    data.frame(
      prefix_iterations = prefix,
      prefix_matvecs = 2L * prefix,
      nconv = sum(fit$certificate$converged),
      certificate_passed = isTRUE(fit$certificate$passed),
      max_residual = fit$certificate$max_residual,
      max_backward_error = fit$certificate$max_backward_error,
      error = "",
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

#' @keywords internal
native_gram_svd <- function(op, rank, target = largest(), tol = 1e-8,
                            vectors = c("both", "left", "right", "none")) {
  vectors <- match.arg(vectors)
  op <- as_operator(op)
  source <- source_or_null(op)
  storage <- op$metadata$storage %||% NULL
  if (!identical(storage, "dgCMatrix") && !(is.matrix(source) && is.double(source))) {
    stop("Native Gram SVD requires a dense double matrix or dgCMatrix operator.", call. = FALSE)
  }
  if (!native_gram_svd_target_supported(target)) {
    stop("Native Gram SVD special case supports only largest singular values.", call. = FALSE)
  }

  A <- if (identical(storage, "dgCMatrix")) op$metadata$matrix else source
  m <- op$dim[1L]
  n <- op$dim[2L]
  limit <- min(m, n)
  rank <- min(as.integer(rank), limit)
  if (rank < 1L) {
    stop("rank must be positive.", call. = FALSE)
  }

  if (identical(storage, "dgCMatrix") && m < n) {
    native <- .Call(
      "eigencore_csc_left_gram_svd",
      methods::slot(A, "i"),
      methods::slot(A, "p"),
      methods::slot(A, "x"),
      methods::slot(A, "Dim"),
      as.integer(rank),
      as.numeric(tol),
      PACKAGE = "eigencore"
    )
    zero_tol <- gram_svd_zero_tolerance(native$d, tol)
    cert <- if (any(native$d <= zero_tol)) {
      certify_svd_operator(op, native$d, native$u, native$v, tol = tol)
    } else {
      new_certificate(
        tol = tol,
        residuals = list(
          left = native$diagnostics$left,
          right = native$diagnostics$right,
          combined = native$diagnostics$combined
        ),
        backward_error = native$diagnostics$backward_error,
        orthogonality = native$diagnostics$orthogonality,
        converged = native$diagnostics$converged,
        scale = native$diagnostics$scale,
        norm_bound_type = "frobenius_exact"
      )
    }
    u <- native$u
    v <- native$v
    if (vectors == "left") {
      v <- NULL
    } else if (vectors == "right") {
      u <- NULL
    } else if (vectors == "none") {
      u <- NULL
      v <- NULL
    }
    return(list(
      d = native$d,
      u = u,
      v = v,
      values = native$d,
      residuals = cert$residuals,
      backward_error = cert$backward_error,
      orthogonality = cert$orthogonality,
      certificate = cert,
      iterations = 1L,
      matvecs = 1L,
      stage_seconds = native$stage_seconds,
      restart = list(
        kind = "gram_svd_special_case",
        implemented = TRUE,
        native = TRUE,
        gram_side = "left",
        gram_dimension = m,
        native_gram_kernel = "csc_left_gram",
        native_gram_eigensolver = native$eigensolver %||% "lapack_dsyevr",
        native_gram_subspace_max_backward_error =
          native$subspace_max_backward_error %||% NA_real_,
        native_implicit_normal_lanczos_max_backward_error =
          native$implicit_lanczos_max_backward_error %||% NA_real_,
        native_implicit_normal_lanczos_iterations =
          native$implicit_lanczos_iterations %||% 0L,
        native_gram_krylov_iterations =
          native$gram_krylov_iterations %||% 0L,
        normal_operator_implicit =
          identical(native$eigensolver %||% "", "implicit_normal_lanczos"),
        materialized_gram =
          !identical(native$eigensolver %||% "", "implicit_normal_lanczos"),
        stage_seconds = native$stage_seconds,
        zero_singular_completion = any(native$d <= zero_tol),
        zero_singular_threshold = zero_tol,
        certificate_reuses_gram_sides = !any(native$d <= zero_tol),
        certified_in_original_coordinates = TRUE
      )
    ))
  }

  if (n <= m) {
    # NN-3: refuse a silent O(n^2) densification of crossprod(A) on a sparse
    # operator when the right-Gram dimension exceeds the native gate. The
    # randomized-SVD refinement path calls native_gram_svd() outside the
    # planner gate, so we must enforce it here. Callers (e.g. randomized
    # SVD's tryCatch) will see this as a fallback signal and stay on the
    # honest sparse path.
    if (identical(storage, "dgCMatrix")) {
      gram_max <- as.integer(getOption("eigencore.gram_svd_max_dimension", 512L))
      if (n > gram_max) {
        stop(
          "native_gram_svd: refusing to densify A^T A for sparse dgCMatrix ",
          "operator with right-Gram dimension n = ", n, " > ",
          "eigencore.gram_svd_max_dimension = ", gram_max,
          ". Use the matrix-free SVD path instead.",
          call. = FALSE
        )
      }
      native <- .Call(
        "eigencore_csc_right_gram_svd",
        methods::slot(A, "i"),
        methods::slot(A, "p"),
        methods::slot(A, "x"),
        methods::slot(A, "Dim"),
        as.integer(rank),
        as.numeric(tol),
        PACKAGE = "eigencore"
      )
      zero_tol <- gram_svd_zero_tolerance(native$d, tol)
      if (!any(native$d <= zero_tol)) {
        cert <- new_certificate(
          tol = tol,
          residuals = list(
            left = native$diagnostics$left,
            right = native$diagnostics$right,
            combined = native$diagnostics$combined
          ),
          backward_error = native$diagnostics$backward_error,
          orthogonality = native$diagnostics$orthogonality,
          converged = native$diagnostics$converged,
          scale = native$diagnostics$scale,
          norm_bound_type = "frobenius_exact"
        )
        u <- native$u
        v <- native$v
        if (vectors == "left") {
          v <- NULL
        } else if (vectors == "right") {
          u <- NULL
        } else if (vectors == "none") {
          u <- NULL
          v <- NULL
        }
        return(list(
          d = native$d,
          u = u,
          v = v,
          values = native$d,
          residuals = cert$residuals,
          backward_error = cert$backward_error,
          orthogonality = cert$orthogonality,
          certificate = cert,
          iterations = 1L,
          matvecs = 1L,
          stage_seconds = native$stage_seconds,
          restart = list(
            kind = "gram_svd_special_case",
            implemented = TRUE,
            native = TRUE,
            gram_side = "right",
            gram_dimension = n,
            native_gram_kernel = "csc_right_gram",
            native_gram_eigensolver = native$eigensolver %||% "lapack_dsyevr",
            native_gram_subspace_max_backward_error =
              native$subspace_max_backward_error %||% NA_real_,
            native_implicit_normal_lanczos_max_backward_error =
              native$implicit_lanczos_max_backward_error %||% NA_real_,
            native_implicit_normal_lanczos_iterations =
              native$implicit_lanczos_iterations %||% 0L,
            native_gram_krylov_iterations =
              native$gram_krylov_iterations %||% 0L,
            normal_operator_implicit = FALSE,
            materialized_gram = TRUE,
            stage_seconds = native$stage_seconds,
            zero_singular_completion = FALSE,
            zero_singular_threshold = zero_tol,
            certificate_reuses_gram_sides = TRUE,
            certified_in_original_coordinates = TRUE
          )
        ))
      }
    }
    gram <- as.matrix(crossprod(A))
    small <- gram_svd_eigen_slice(gram, rank, target)
    d <- sqrt(pmax(small$values, 0))
    v_full <- small$vectors
    u_full <- as.matrix(A %*% v_full)
    av_full <- u_full
    atu_full <- sweep(v_full, 2L, d, `*`)
    zero_tol <- gram_svd_zero_tolerance(d, tol)
    nz <- d > zero_tol
    d[!nz] <- 0
    if (any(nz)) {
      u_full[, nz] <- sweep(u_full[, nz, drop = FALSE], 2L, d[nz], `/`)
    }
    if (any(!nz)) {
      v_full[, !nz] <- orthonormal_completion(
        v_full[, nz, drop = FALSE],
        n_rows = n,
        needed = sum(!nz),
        tol = zero_tol
      )
      u_full[, !nz] <- orthonormal_completion(
        u_full[, nz, drop = FALSE],
        n_rows = m,
        needed = sum(!nz),
        tol = zero_tol
      )
      av_full[, !nz] <- 0
      atu_full[, !nz] <- 0
    }
  } else {
    gram <- as.matrix(tcrossprod(A))
    small <- gram_svd_eigen_slice(gram, rank, target)
    d <- sqrt(pmax(small$values, 0))
    u_full <- small$vectors
    v_full <- as.matrix(crossprod(A, u_full))
    atu_full <- v_full
    av_full <- sweep(u_full, 2L, d, `*`)
    zero_tol <- gram_svd_zero_tolerance(d, tol)
    nz <- d > zero_tol
    d[!nz] <- 0
    if (any(nz)) {
      v_full[, nz] <- sweep(v_full[, nz, drop = FALSE], 2L, d[nz], `/`)
    }
    if (any(!nz)) {
      u_full[, !nz] <- orthonormal_completion(
        u_full[, nz, drop = FALSE],
        n_rows = m,
        needed = sum(!nz),
        tol = zero_tol
      )
      v_full[, !nz] <- orthonormal_completion(
        v_full[, nz, drop = FALSE],
        n_rows = n,
        needed = sum(!nz),
        tol = zero_tol
      )
      av_full[, !nz] <- 0
      atu_full[, !nz] <- 0
    }
  }

  cert <- if (any(!nz)) {
    certify_svd_operator(op, d, u_full, v_full, tol = tol)
  } else {
    certify_svd_operator_cached_sides(
      op, d, u_full, v_full, av_full, atu_full, tol = tol
    )
  }
  u <- u_full
  v <- v_full
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
    iterations = 1L,
    matvecs = 2L,
    restart = list(
      kind = "gram_svd_special_case",
      implemented = TRUE,
      native = TRUE,
      gram_side = if (n <= m) "right" else "left",
      gram_dimension = min(m, n),
      zero_singular_completion = any(!nz),
      zero_singular_threshold = zero_tol,
      certificate_reuses_gram_sides = !any(!nz),
      certified_in_original_coordinates = TRUE
    )
  )
}

#' @keywords internal
gram_svd_eigen_slice <- function(gram, rank, target) {
  rank <- min(as.integer(rank), nrow(gram))
  kind <- if (inherits(target, "eigencore_target")) target$kind else "largest"
  if (kind %in% c("largest", "largest_magnitude") && rank < nrow(gram)) {
    selected <- tryCatch(
      native_dense_symmetric_eigen_selected(gram, rank, largest()),
      error = function(e) NULL
    )
    if (!is.null(selected)) {
      return(list(
        values = selected$values,
        vectors = selected$vectors
      ))
    }
  }
  eig <- native_dense_symmetric_eigen(gram)
  idx <- order_indices(eig$values, target)
  idx <- idx[seq_len(min(rank, length(idx)))]
  list(
    values = eig$values[idx],
    vectors = eig$vectors[, idx, drop = FALSE]
  )
}

#' @keywords internal
gram_svd_zero_tolerance <- function(d, tol) {
  scale <- max(1, d, na.rm = TRUE)
  max(100 * .Machine$double.eps * scale, sqrt(.Machine$double.eps) * scale, tol * scale * 1e-3)
}

#' @keywords internal
orthonormal_completion <- function(Q, n_rows, needed, tol = sqrt(.Machine$double.eps)) {
  needed <- as.integer(needed)
  if (needed < 1L) {
    return(matrix(0, n_rows, 0L))
  }
  Q <- if (is.null(Q) || ncol(Q) == 0L) {
    matrix(0, n_rows, 0L)
  } else {
    as.matrix(Q)
  }
  out <- matrix(0, n_rows, needed)
  accepted <- 0L
  against <- Q
  for (j in seq_len(n_rows)) {
    z <- numeric(n_rows)
    z[[j]] <- 1
    if (ncol(against) > 0L) {
      for (pass in 1:2) {
        z <- z - as.vector(against %*% crossprod(against, z))
      }
    }
    z_norm <- sqrt(sum(z^2))
    if (z_norm > tol) {
      accepted <- accepted + 1L
      out[, accepted] <- z / z_norm
      against <- cbind(against, out[, accepted, drop = FALSE])
      if (accepted == needed) {
        break
      }
    }
  }
  if (accepted < needed) {
    stop("Could not complete an orthonormal nullspace basis for zero singular triplets.",
         call. = FALSE)
  }
  out
}

#' @keywords internal
complete_zero_singular_triplets <- function(op, final, rank, target, tol) {
  kind <- if (inherits(target, "eigencore_target")) target$kind else "largest"
  if (!kind %in% c("largest", "largest_magnitude") ||
      is.null(final$u) || is.null(final$v)) {
    final$zero_singular_completion <- FALSE
    final$zero_singular_threshold <- NA_real_
    return(final)
  }

  rank <- min(as.integer(rank), min(op$dim))
  if (length(final$d) >= rank && !any(final$d <= gram_svd_zero_tolerance(final$d, tol))) {
    final$zero_singular_completion <- FALSE
    final$zero_singular_threshold <- NA_real_
    return(final)
  }

  zero_tol <- gram_svd_zero_tolerance(final$d, tol)
  nz <- final$d > zero_tol
  nonzero_count <- min(sum(nz), rank)
  zero_count <- rank - nonzero_count
  if (zero_count <= 0L) {
    final$zero_singular_completion <- FALSE
    final$zero_singular_threshold <- zero_tol
    return(final)
  }

  u_nonzero <- final$u[, nz, drop = FALSE]
  v_nonzero <- final$v[, nz, drop = FALSE]
  if (ncol(u_nonzero) > nonzero_count) {
    u_nonzero <- u_nonzero[, seq_len(nonzero_count), drop = FALSE]
    v_nonzero <- v_nonzero[, seq_len(nonzero_count), drop = FALSE]
  }
  d_nonzero <- final$d[nz]
  if (length(d_nonzero) > nonzero_count) {
    d_nonzero <- d_nonzero[seq_len(nonzero_count)]
  }

  u_zero <- orthonormal_completion(
    u_nonzero,
    n_rows = op$dim[[1L]],
    needed = zero_count,
    tol = zero_tol
  )
  v_zero <- orthonormal_completion(
    v_nonzero,
    n_rows = op$dim[[2L]],
    needed = zero_count,
    tol = zero_tol
  )
  d <- c(d_nonzero, rep(0, zero_count))
  u <- cbind(u_nonzero, u_zero)
  v <- cbind(v_nonzero, v_zero)
  cert <- certify_svd_operator(op, d, u, v, tol = tol)

  final$d <- d
  final$u <- u
  final$v <- v
  final$values <- d
  final$residuals <- cert$residuals
  final$backward_error <- cert$backward_error
  final$orthogonality <- cert$orthogonality
  final$certificate <- cert
  final$zero_singular_completion <- TRUE
  final$zero_singular_threshold <- zero_tol
  final
}

#' @keywords internal
reference_randomized_svd <- function(op, rank, target = largest(), tol = 1e-8,
                                     oversample = 10L, n_iter = 2L,
                                     vectors = c("both", "left", "right", "none"),
                                     refine = TRUE,
                                     normalizer = c("qr", "lu", "none")) {
  record_stages <- isTRUE(getOption("eigencore.randomized_stage_timing", FALSE))
  stage_seconds <- c(
    random = 0,
    apply = 0,
    normalize = 0,
    small_svd = 0,
    vector_form = 0,
    certificate = 0,
    refinement = 0
  )
  tick <- function() if (record_stages) proc.time()[["elapsed"]] else 0
  add_stage <- function(name, start) {
    if (record_stages) {
      stage_seconds[[name]] <<- stage_seconds[[name]] + (tick() - start)
    }
  }
  vectors <- match.arg(vectors)
  normalizer <- match.arg(normalizer)
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
  apply_pair <- randomized_svd_apply_pair(op)
  adaptive_stop <- isTRUE(getOption("eigencore.randomized_adaptive_stop", TRUE))
  adaptive_candidate_allowed <- adaptive_stop && identical(normalizer, "qr")

  candidate_from_Q <- function(Q) {
    t0 <- tick()
    uses_direct_projection <- !is.null(apply_pair$project)
    B_t <- if (uses_direct_projection) {
      apply_pair$project(Q)
    } else {
      apply_pair$apply_adjoint(Q)
    }
    add_stage("apply", t0)
    projection_transposed <- isTRUE(attr(B_t, "transposed", exact = TRUE))
    core <- if (projection_transposed) {
      B_t
    } else {
      t(B_t)
    }
    t0 <- tick()
    small <- randomized_svd_core_decomposition(core, rank = rank, target = target)
    idx <- seq_along(small$d)
    idx <- idx[seq_len(min(rank, length(idx)))]
    d <- small$d[idx]
    add_stage("small_svd", t0)
    t0 <- tick()
    small_u <- small$u[, idx, drop = FALSE]
    u <- Q %*% small_u
    v <- small$v[, idx, drop = FALSE]
    add_stage("vector_form", t0)
    t0 <- tick()
    cert <- certify_randomized_svd_projection(
      op,
      apply_pair,
      core = core,
      small_u = small_u,
      d = d,
      u = u,
      v = v,
      tol = tol
    )
    add_stage("certificate", t0)
    list(
      d = d,
      u = u,
      v = v,
      cert = cert,
      core_solver = small$solver,
      projection_kind = if (uses_direct_projection) {
        apply_pair$project_kind %||% "direct_qt_a"
      } else {
        "adjoint_apply"
      },
      projection_transposed = projection_transposed
    )
  }

  t0 <- tick()
  Omega <- matrix(stats::rnorm(n * l), nrow = n, ncol = l)
  add_stage("random", t0)
  t0 <- tick()
  Y <- apply_pair$apply(Omega)
  add_stage("apply", t0)
  matvecs <- 1L
  t0 <- tick()
  Q <- randomized_svd_normalize(Y, normalizer, final = n_iter == 0L)
  add_stage("normalize", t0)
  candidate <- NULL
  initial_cert <- NULL
  early_stop_used <- FALSE
  iterations_used <- 1L
  if (n_iter == 0L || isTRUE(adaptive_candidate_allowed)) {
    candidate <- candidate_from_Q(Q)
    matvecs <- matvecs + 1L
    initial_cert <- candidate$cert
    early_stop_used <- n_iter > 0L && isTRUE(candidate$cert$passed)
  }

  if (!early_stop_used) {
    if (n_iter > 0L) {
      for (iter in seq_len(n_iter)) {
        t0 <- tick()
        Z <- apply_pair$apply_adjoint(Q)
        add_stage("apply", t0)
        t0 <- tick()
        Z <- randomized_svd_normalize(Z, normalizer, final = FALSE)
        add_stage("normalize", t0)
        t0 <- tick()
        Y <- apply_pair$apply(Z)
        add_stage("apply", t0)
        matvecs <- matvecs + 2L
        t0 <- tick()
        Q <- randomized_svd_normalize(Y, normalizer, final = iter == n_iter)
        add_stage("normalize", t0)
      }
      candidate <- candidate_from_Q(Q)
      matvecs <- matvecs + 1L
      iterations_used <- n_iter + 1L
    }
  }
  if (is.null(initial_cert)) {
    initial_cert <- candidate$cert
  }
  d <- candidate$d
  u <- candidate$u
  v <- candidate$v
  cert <- candidate$cert

  refinement <- NULL
  if (isTRUE(refine) && !isTRUE(cert$passed)) {
    t0 <- tick()
    refinement <- tryCatch(
      native_gram_svd(op, rank = rank, target = target, tol = tol, vectors = "both"),
      error = function(e) {
        structure(list(error = conditionMessage(e)), class = "eigencore_refinement_error")
      }
    )
    add_stage("refinement", t0)
    if (!inherits(refinement, "eigencore_refinement_error") &&
        isTRUE(refinement$certificate$passed)) {
      d <- refinement$d
      u <- refinement$u
      v <- refinement$v
      cert <- refinement$certificate
    }
  }

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
    stage_seconds = if (record_stages) stage_seconds else numeric(),
    iterations = iterations_used,
    matvecs = matvecs + if (!is.null(refinement) && !inherits(refinement, "eigencore_refinement_error")) {
      refinement$matvecs
    } else {
      0L
    },
    restart = list(
      kind = if (!is.null(refinement) && !inherits(refinement, "eigencore_refinement_error") &&
        isTRUE(refinement$certificate$passed)) {
        "randomized_range_finder_refined"
      } else {
        "randomized_range_finder"
      },
      implemented = TRUE,
      native = FALSE,
      oversample = oversample,
      n_iter = n_iter,
      normalizer = normalizer,
      apply_kind = apply_pair$kind,
      native_sketch = isTRUE(apply_pair$native_sketch),
      core_solver = candidate$core_solver %||% NA_character_,
      projection_kind = candidate$projection_kind %||% NA_character_,
      projection_transposed = isTRUE(candidate$projection_transposed),
      certificate_reuses_projection = TRUE,
      adaptive_stop = adaptive_stop,
      adaptive_stop_used = early_stop_used,
      iterations_used = iterations_used,
      stage_seconds = if (record_stages) stage_seconds else numeric(),
      sample_dimension = l,
      approximate = TRUE,
      certificate_policy = if (isTRUE(refine)) {
        "residual certificate, with deterministic refinement when needed"
      } else {
        "residual certificate only; stochastic sketch is not sufficient to pass"
      },
      refine = isTRUE(refine),
      initial_certificate_passed = isTRUE(initial_cert$passed),
      initial_max_backward_error = initial_cert$max_backward_error,
      refinement_attempted = !is.null(refinement),
      refinement_kind = if (!is.null(refinement) && !inherits(refinement, "eigencore_refinement_error")) {
        refinement$restart$kind
      } else {
        NA_character_
      },
      refinement_passed = !is.null(refinement) &&
        !inherits(refinement, "eigencore_refinement_error") &&
        isTRUE(refinement$certificate$passed),
      refinement_error = if (inherits(refinement, "eigencore_refinement_error")) {
        refinement$error
      } else {
        NA_character_
      }
    )
  )
}

#' @keywords internal
randomized_svd_core_decomposition <- function(core, rank, target = largest()) {
  rank <- min(as.integer(rank), min(dim(core)))
  if (rank < 1L) {
    return(list(d = numeric(), u = matrix(0, nrow(core), 0L), v = matrix(0, ncol(core), 0L),
                solver = "empty"))
  }
  if (nrow(core) <= 128L && ncol(core) >= 2L * nrow(core)) {
    gram <- tcrossprod(core)
    eig <- eigen(gram, symmetric = TRUE)
    idx <- order_indices(eig$values, target)
    idx <- idx[seq_len(min(rank, length(idx)))]
    values <- pmax(eig$values[idx], 0)
    d <- sqrt(values)
    u <- eig$vectors[, idx, drop = FALSE]
    v <- crossprod(core, u)
    for (col in seq_along(d)) {
      if (d[[col]] > 100 * .Machine$double.eps) {
        v[, col] <- v[, col] / d[[col]]
      } else {
        v[, col] <- 0
      }
    }
    return(list(d = d, u = u, v = v, solver = "left_gram_eigen"))
  }

  small <- svd(core, nu = min(nrow(core), rank), nv = min(ncol(core), rank))
  idx <- order_indices(small$d, target)
  idx <- idx[seq_len(min(rank, length(idx)))]
  list(
    d = small$d[idx],
    u = small$u[, idx, drop = FALSE],
    v = small$v[, idx, drop = FALSE],
    solver = "dense_svd"
  )
}

#' @keywords internal
certify_randomized_svd_projection <- function(op, apply_pair, core, small_u,
                                              d, u, v, tol = 1e-8) {
  left_residual_matrix <- apply_pair$apply(v) - sweep(u, 2L, d, `*`)
  right_applied <- crossprod(core, small_u)
  right_residual_matrix <- right_applied - sweep(v, 2L, d, `*`)
  left <- col_norms(left_residual_matrix)
  right <- col_norms(right_residual_matrix)
  combined <- sqrt(left^2 + right^2)
  norm_A <- operator_norm_for_certificate_info(op)
  scale <- svd_backward_scale(norm_A$value, d)
  backward <- combined / scale
  orth_u <- max(abs(crossprod(u) - diag(length(d))))
  orth_v <- max(abs(crossprod(v) - diag(length(d))))
  new_certificate(
    tol = tol,
    residuals = list(left = left, right = right, combined = combined),
    backward_error = backward,
    orthogonality = c(U = orth_u, V = orth_v),
    converged = backward <= tol,
    scale = scale,
    norm_bound_type = norm_A$norm_bound_type,
    scale_is_estimate = isTRUE(norm_A$scale_is_estimate)
  )
}

#' @keywords internal
randomized_svd_apply_pair <- function(op) {
  source <- source_or_null(op) %||% op$metadata$matrix %||% NULL
  if (is.matrix(source) && is.double(source)) {
    dense_randomized_apply <- function(X, transpose = FALSE) {
      .Call(
        "eigencore_dense_randomized_apply",
        source,
        X,
        as.logical(transpose),
        PACKAGE = "eigencore"
      )
    }
    return(list(
      kind = "dense_direct",
      native_sketch = TRUE,
      apply = function(X) dense_randomized_apply(X, transpose = FALSE),
      apply_adjoint = function(X) dense_randomized_apply(X, transpose = TRUE),
      project_kind = "native_direct_qt_a",
      project = function(Q) {
        .Call(
          "eigencore_dense_randomized_project_transposed",
          source,
          Q,
          PACKAGE = "eigencore"
        )
      }
    ))
  }
  if (inherits(source, "dgCMatrix")) {
    csc_randomized_apply <- function(X, transpose = FALSE) {
      .Call(
        "eigencore_csc_randomized_apply",
        methods::slot(source, "i"),
        methods::slot(source, "p"),
        methods::slot(source, "x"),
        methods::slot(source, "Dim"),
        X,
        as.logical(transpose),
        PACKAGE = "eigencore"
      )
    }
    return(list(
      kind = "csc_direct",
      native_sketch = TRUE,
      apply = function(X) csc_randomized_apply(X, transpose = FALSE),
      apply_adjoint = function(X) csc_randomized_apply(X, transpose = TRUE),
      project_kind = "native_direct_qt_a",
      project = function(Q) {
        .Call(
          "eigencore_csc_randomized_project_transposed",
          methods::slot(source, "i"),
          methods::slot(source, "p"),
          methods::slot(source, "x"),
          methods::slot(source, "Dim"),
          Q,
          PACKAGE = "eigencore"
        )
      }
    ))
  }
  list(
    kind = "operator",
    native_sketch = FALSE,
    apply = function(X) apply_operator(op, X),
    apply_adjoint = function(X) apply_adjoint_operator(op, X)
  )
}

#' @keywords internal
randomized_svd_normalize <- function(X, normalizer = c("qr", "lu", "none"),
                                     final = FALSE) {
  normalizer <- match.arg(normalizer)
  if (isTRUE(final) || identical(normalizer, "qr")) {
    return(qr.Q(qr(X)))
  }
  if (identical(normalizer, "none")) {
    return(X)
  }
  randomized_svd_lu_normalize(X)
}

#' @keywords internal
randomized_svd_lu_normalize <- function(X) {
  out <- tryCatch({
    lu <- Matrix::lu(Matrix::Matrix(X, sparse = FALSE))
    pieces <- Matrix::expand(lu)
    L <- as.matrix(pieces$L)
    if (ncol(L) > ncol(X)) {
      L <- L[, seq_len(ncol(X)), drop = FALSE]
    }
    L
  }, error = function(e) NULL)
  if (is.null(out) || !identical(nrow(out), nrow(X)) || ncol(out) < ncol(X)) {
    return(qr.Q(qr(X)))
  }
  out[, seq_len(ncol(X)), drop = FALSE]
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
native_gram_svd_target_supported <- function(target) {
  kind <- if (inherits(target, "eigencore_target")) target$kind else "largest"
  kind %in% c("largest", "largest_magnitude")
}

#' @keywords internal
native_golub_kahan_ritz <- function(op, U, V, alpha, beta, rank, target, tol,
                                    active_iterations = NULL) {
  active_iterations <- as.integer(active_iterations %||% length(alpha))
  ritz <- .Call(
    "eigencore_golub_kahan_ritz",
    U,
    V,
    as.numeric(alpha),
    as.numeric(beta),
    as.integer(rank),
    as.integer(native_svd_target_kind(target)),
    active_iterations,
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
