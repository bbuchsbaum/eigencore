#' @keywords internal
default_block_lanczos_max_subspace <- function(k, block) {
  k <- as.integer(k)
  block <- as.integer(block)
  k_term <- if (k >= 16L) {
    8L * k
  } else {
    6L * k + 10L
  }
  max(k_term, 6L * block + 20L)
}

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
  controls <- reference_scalar_subspace_controls(
    requested = k,
    requested_name = "k",
    limit = n,
    maxit = maxit,
    default_maxit = function(k) max(20L, 4L * k + 20L)
  )
  k <- controls$requested
  maxit <- controls$maxit

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

  out <- native_block_lanczos_hermitian(
    op,
    k = k,
    target = target,
    tol = tol,
    maxit = m_max,
    block = 1L,
    max_restarts = max_restarts,
    vectors = vectors,
    full_subspace = FALSE,
    certificate_fallback = FALSE
  )
  out$restart$kind <- "thick_restart"
  if (is.data.frame(out$convergence_history)) {
    out$convergence_history$iteration <- out$convergence_history$m_active
    out$convergence_history$n_locked <- out$convergence_history$locked_after
  }
  if (isTRUE(vectors) && !is.null(out$vectors)) {
    source <- source_or_null(op)
    Av <- if (is.matrix(source) && is.double(source)) {
      source %*% out$vectors
    } else {
      apply_operator(op, out$vectors)
    }
    residuals <- col_norms(Av - sweep(out$vectors, 2L, out$values, `*`))
    cert <- certify_eigen_operator_residuals(op, out$values, out$vectors,
                                             residuals, tol = tol)
    out$residuals <- cert$residuals
    out$backward_error <- cert$backward_error
    out$orthogonality <- cert$orthogonality
    out$certificate <- cert
  }
  out
}

#' @keywords internal
native_block_lanczos_hermitian <- function(op, k, target = largest(), tol = 1e-8,
                                           maxit = NULL, block = 2L,
                                           max_restarts = 100L,
                                           vectors = TRUE,
                                           full_subspace = TRUE,
                                           certificate_fallback = TRUE) {
  op <- as_operator(op)
  if (op$dim[1L] != op$dim[2L]) {
    stop("Native block Hermitian Lanczos requires a square operator.", call. = FALSE)
  }
  if (!identical(op$structure$kind, "hermitian")) {
    stop("Native block Hermitian Lanczos requires an operator with hermitian() structure.", call. = FALSE)
  }

  n <- op$dim[1L]
  block <- as.integer(block)
  if (length(block) != 1L || is.na(block) || block < 1L) {
    stop("native block Hermitian Lanczos requires block >= 1.", call. = FALSE)
  }
  m_max <- if (is.null(maxit)) {
    min(n, default_block_lanczos_max_subspace(k, block))
  } else {
    min(n, as.integer(maxit))
  }
  if (m_max < as.integer(k) + block) {
    stop("max_subspace must be at least k + block for block thick-restart Lanczos.", call. = FALSE)
  }
  max_restarts <- as.integer(max_restarts)
  if (length(max_restarts) != 1L || is.na(max_restarts) || max_restarts < 0L) {
    stop("max_restarts must be a non-negative integer.", call. = FALSE)
  }

  norm_A <- operator_norm_for_certificate_info(op)$value
  storage <- op$metadata$storage %||% NULL
  source <- source_or_null(op)
  target_kind <- lanczos_target_kind(target)
  if (isTRUE(full_subspace) && is.matrix(source) && is.double(source) && m_max >= n) {
    return(native_block_full_subspace_hermitian(
      op,
      k = k,
      target = target,
      tol = tol,
      block = block,
      max_subspace = m_max,
      max_restarts = max_restarts,
      vectors = vectors
    ))
  }

  start <- matrix(stats::rnorm(n * block), nrow = n, ncol = block)
  iter <- if (identical(storage, "dgCMatrix")) {
    A <- op$metadata$matrix
    .Call(
      "eigencore_block_thick_restart_lanczos_csc",
      methods::slot(A, "i"),
      methods::slot(A, "p"),
      methods::slot(A, "x"),
      methods::slot(A, "Dim"),
      as.integer(k),
      as.integer(m_max),
      as.integer(block),
      as.integer(target_kind),
      as.numeric(tol),
      as.integer(max_restarts),
      as.numeric(norm_A),
      start,
      PACKAGE = "eigencore"
    )
  } else if (is.matrix(source) && is.double(source)) {
    .Call(
      "eigencore_block_thick_restart_lanczos_dense",
      source,
      as.integer(k),
      as.integer(m_max),
      as.integer(block),
      as.integer(target_kind),
      as.numeric(tol),
      as.integer(max_restarts),
      as.numeric(norm_A),
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
  ortho_passes <- as.integer(iter$ortho_passes %||% NA_integer_)
  locking_events <- as.integer(iter$locking_events %||% 0L)
  restarts_used <- as.integer(iter$restarts %||% 0L)
  locked <- seq_len(iter$n_locked %||% 0L)
  operator_allocations <- as.numeric(iter$operator_allocations %||% NA_real_)
  operator_bytes_allocated <- as.numeric(iter$operator_bytes_allocated %||% NA_real_)
  stage_seconds <- iter$stage_seconds %||% numeric()
  restart_history <- iter$restart_history %||% NULL
  restart_history <- if (is.list(restart_history) && length(restart_history)) {
    data.frame(
      restart = as.integer(restart_history$restart),
      m_active = as.integer(restart_history$m_active),
      selected_count = as.integer(restart_history$selected_count),
      locked_before = as.integer(restart_history$locked_before),
      locked_after = as.integer(restart_history$locked_after),
      nconv_wanted = as.integer(restart_history$nconv_wanted),
      max_residual = as.numeric(restart_history$max_residual),
      max_backward_error = as.numeric(restart_history$max_backward_error)
    )
  } else {
    data.frame(
      restart = restarts_used,
      m_active = iter$m_active_final %||% NA_integer_,
      selected_count = NA_integer_,
      locked_before = NA_integer_,
      locked_after = iter$n_locked %||% sum(iter$converged),
      nconv_wanted = sum(iter$converged),
      max_residual = max(cert$residuals),
      max_backward_error = cert$max_backward_error
    )
  }

  block_restart <- list(
    kind = "block_thick_restart_candidate",
    implemented = TRUE,
    locking = "in_native_loop",
    locked = locked,
    locked_count = iter$n_locked %||% sum(iter$converged),
    restarts_used = restarts_used,
    max_restarts = max_restarts,
    max_subspace = m_max,
    final_active_subspace = iter$m_active_final %||% NA_integer_,
    block = block,
    ortho_passes = ortho_passes,
    locking_events = locking_events,
    operator_allocations = operator_allocations,
    operator_bytes_allocated = operator_bytes_allocated,
    stage_seconds = stage_seconds,
    history = restart_history
  )

  certificate_failed_orthogonality <- is.finite(cert$max_orthogonality_loss) &&
    is.finite(cert$orthogonality_tolerance) &&
    cert$max_orthogonality_loss > cert$orthogonality_tolerance
  if (isTRUE(certificate_fallback) &&
      !isTRUE(cert$passed) && isTRUE(certificate_failed_orthogonality) &&
      ((is.matrix(source) && is.double(source)) || identical(storage, "dgCMatrix"))) {
    fallback <- if (is.matrix(source) && is.double(source)) {
      native_block_full_subspace_hermitian(
        op,
        k = k,
        target = target,
        tol = tol,
        block = block,
        max_subspace = n,
        max_restarts = max_restarts,
        vectors = vectors
      )
    } else {
      native_block_lanczos_hermitian(
        op,
        k = k,
        target = target,
        tol = tol,
        maxit = m_max,
        block = 1L,
        max_restarts = max_restarts,
        vectors = vectors,
        full_subspace = FALSE,
        certificate_fallback = FALSE
      )
    }
    fallback$restarts <- restarts_used
    fallback$ortho_passes <- ortho_passes
    fallback$locking_events <- locking_events
    fallback$operator_allocations <- operator_allocations
    fallback$operator_bytes_allocated <- operator_bytes_allocated
    fallback$stage_seconds <- stage_seconds
    fallback$convergence_history <- restart_history
    fallback$locked <- seq_len(sum(fallback$certificate$converged))
    fallback$restart$kind <- if (is.matrix(source) && is.double(source)) {
      "block_dense_lapack_certificate_fallback"
    } else {
      "block_scalar_lanczos_certificate_fallback"
    }
    fallback$restart$fallback_used <- TRUE
    fallback$restart$fallback_reason <- "native block Lanczos certificate failed"
    fallback$restart$failed_block_certificate <- cert
    fallback$restart$failed_block_restart <- block_restart
    fallback$restart$history <- restart_history
    return(fallback)
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
    restarts = restarts_used,
    ortho_passes = ortho_passes,
    locking_events = locking_events,
    block = block,
    operator_allocations = operator_allocations,
    operator_bytes_allocated = operator_bytes_allocated,
    stage_seconds = stage_seconds,
    convergence_history = restart_history,
    locked = locked,
    restart = block_restart
  )
}

#' @keywords internal
native_block_full_subspace_hermitian <- function(op, k, target, tol, block,
                                                 max_subspace, max_restarts,
                                                 vectors = TRUE) {
  op <- as_operator(op)
  source <- source_or_null(op)
  if (!(is.matrix(source) && is.double(source))) {
    stop("Full-subspace block Hermitian fallback requires a dense double source.", call. = FALSE)
  }
  target_kind <- lanczos_target_kind(target)
  eig <- if (target_kind %in% c(1L, 2L)) {
    tryCatch(
      native_dense_symmetric_eigen_selected(source, k, target),
      error = function(e) native_dense_symmetric_eigen(source)
    )
  } else {
    native_dense_symmetric_eigen(source)
  }
  eig_values <- eig[[1L]]
  eig_vectors <- eig[[2L]]
  if (length(eig_values) > as.integer(k)) {
    idx <- order_indices(eig_values, target)
    idx <- idx[seq_len(min(as.integer(k), length(idx)))]
    values <- eig_values[idx]
    vec_matrix <- eig_vectors[, idx, drop = FALSE]
  } else {
    values <- eig_values
    vec_matrix <- eig_vectors
  }
  cert <- certify_eigen_operator(op, values, vec_matrix, tol = tol)
  if (!isTRUE(vectors)) {
    vec_matrix <- NULL
  }
  nconv <- sum(cert$converged)
  list(
    values = values,
    vectors = vec_matrix,
    residuals = cert$residuals,
    backward_error = cert$backward_error,
    orthogonality = cert$orthogonality,
    certificate = cert,
    iterations = 1L,
    matvecs = 0L,
    restarts = 0L,
    ortho_passes = 0L,
    locking_events = 0L,
    block = as.integer(block),
    restart = list(
      kind = "block_full_subspace_dense_lapack",
      locking = "not_required_full_subspace",
      locked_count = nconv,
      max_subspace = max_subspace,
      final_active_subspace = nrow(source),
      block = as.integer(block)
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
  out <- .Call(
    "eigencore_tridiagonal_eigen",
    as.numeric(alpha),
    as.numeric(beta),
    PACKAGE = "eigencore"
  )
  if (is.null(names(out)) && length(out) == 2L) {
    names(out) <- c("values", "vectors")
  }
  out
}

#' @keywords internal
native_tridiagonal_eigen_selected <- function(alpha, beta, k, target) {
  out <- .Call(
    "eigencore_tridiagonal_eigen_selected",
    as.numeric(alpha),
    as.numeric(beta),
    as.integer(k),
    as.integer(lanczos_target_kind(target)),
    PACKAGE = "eigencore"
  )
  if (is.null(names(out)) && length(out) == 2L) {
    names(out) <- c("values", "vectors")
  }
  out
}

#' @keywords internal
tridiagonal_matrix <- function(alpha, beta) {
  m <- length(alpha)
  tridiagonal <- diag(alpha, m)
  if (m > 1L) {
    off <- beta[seq_len(m - 1L)]
    tridiagonal[cbind(seq_len(m - 1L), 2:m)] <- off
    tridiagonal[cbind(2:m, seq_len(m - 1L))] <- off
  }
  tridiagonal
}
