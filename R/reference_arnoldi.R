#' @keywords internal
reference_arnoldi_label <- function() {
  "reference Arnoldi (prototype/oracle fallback)"
}

#' @keywords internal
native_arnoldi_label <- function() {
  "native Arnoldi cycle + native Ritz extraction (compatibility)"
}

#' @keywords internal
native_refined_arnoldi_label <- function() {
  "native Arnoldi cycle + native refined Ritz extraction (V2 tranche)"
}

#' @keywords internal
native_matrix_free_arnoldi_label <- function() {
  "native matrix-free Arnoldi callback cycle + native Ritz extraction"
}

#' @keywords internal
reference_arnoldi_target_supported <- function(target) {
  kind <- if (inherits(target, "eigencore_target")) target$kind else "largest"
  kind %in% c(
    "largest",
    "smallest",
    "largest_magnitude",
    "largest_real",
    "smallest_real",
    "largest_imaginary",
    "smallest_imaginary"
  )
}

#' @keywords internal
# For real operators A and A^T share the same eigenvalues, so the same target
# works for both the forward and adjoint solves.  For complex operators the
# adjoint eigenvalues are the conjugates of A's eigenvalues, which would
# require flipping imaginary-part targets (largest_imaginary <->
# smallest_imaginary).  That case does not arise today because every code path
# that reaches arnoldi_left_eigen_contract operates on real (dtype "double")
# operators.  If complex operator support is added, this function must be
# updated to conjugate imaginary-part targets.
adjoint_arnoldi_target <- function(target) {
  target
}

#' @keywords internal
match_left_eigenvectors <- function(left_values, left_vectors, values) {
  left_values <- as.vector(left_values)
  left_vectors <- as.matrix(left_vectors)
  wanted <- as.vector(values)
  used <- rep(FALSE, length(left_values))
  idx <- integer(length(wanted))
  distance <- rep(Inf, length(wanted))

  for (j in seq_along(wanted)) {
    available <- which(!used)
    if (!length(available)) break
    d <- abs(left_values[available] - wanted[[j]])
    pick <- available[[which.min(d)]]
    idx[[j]] <- pick
    distance[[j]] <- min(d)
    used[[pick]] <- TRUE
  }

  if (any(idx == 0L)) {
    return(NULL)
  }
  list(
    values = left_values[idx],
    vectors = left_vectors[, idx, drop = FALSE],
    match_distance = distance
  )
}

#' @keywords internal
normalize_left_eigenvectors <- function(left_vectors, right_vectors) {
  left_vectors <- as.matrix(left_vectors)
  right_vectors <- as.matrix(right_vectors)
  gram <- crossprod(left_vectors, right_vectors)
  diag_gram <- diag(gram)
  stable <- is.finite(Mod(diag_gram)) &
    Mod(diag_gram) > 100 * .Machine$double.eps
  if (any(stable)) {
    left_vectors[, stable] <- sweep(
      left_vectors[, stable, drop = FALSE],
      2L,
      diag_gram[stable],
      `/`
    )
  }
  left_vectors
}

#' @keywords internal
arnoldi_left_eigen_contract <- function(op, values, right_vectors, target,
                                        tol = 1e-8, maxit = NULL,
                                        max_restarts = 0L,
                                        extraction = "projected_ritz") {
  if (is.null(right_vectors)) {
    return(list(
      supported = FALSE,
      reason = "right eigenvectors were not returned",
      vectors = NULL,
      certificate = NULL,
      biorthogonality = NULL
    ))
  }

  adjoint_op <- tryCatch(
    adjoint(op),
    error = function(e) e
  )
  if (inherits(adjoint_op, "error")) {
    return(list(
      supported = FALSE,
      reason = conditionMessage(adjoint_op),
      vectors = NULL,
      certificate = NULL,
      biorthogonality = NULL
    ))
  }

  k <- length(values)
  left_solver <- if (native_arnoldi_available(adjoint_op) ||
      native_matrix_free_arnoldi_available(adjoint_op)) {
    native_arnoldi_general
  } else {
    reference_arnoldi_general
  }
  left_iter <- tryCatch(
    left_solver(
      adjoint_op,
      k = k,
      target = adjoint_arnoldi_target(target),
      tol = tol,
      maxit = maxit,
      max_restarts = max_restarts,
      vectors = TRUE,
      extraction = extraction
    ),
    error = function(e) e
  )
  if (inherits(left_iter, "error")) {
    return(list(
      supported = FALSE,
      reason = conditionMessage(left_iter),
      vectors = NULL,
      certificate = NULL,
      biorthogonality = NULL
    ))
  }

  matched <- match_left_eigenvectors(left_iter$values, left_iter$vectors, values)
  if (is.null(matched)) {
    return(list(
      supported = FALSE,
      reason = "left eigenvalue matching failed",
      vectors = NULL,
      certificate = NULL,
      biorthogonality = NULL
    ))
  }

  left_vectors <- normalize_left_eigenvectors(matched$vectors, right_vectors)
  cert <- certify_left_eigen_operator(
    op,
    values,
    left_vectors,
    right_vectors = right_vectors,
    tol = tol
  )
  biorthogonality <- crossprod(left_vectors, as.matrix(right_vectors))
  list(
    supported = TRUE,
    reason = "left eigenvectors computed from the adjoint operator",
    values = matched$values,
    match_distance = matched$match_distance,
    vectors = left_vectors,
    certificate = cert,
    biorthogonality = biorthogonality,
    method = left_iter$restart$kind %||% "adjoint_arnoldi"
  )
}

#' @keywords internal
native_arnoldi_available <- function(op) {
  op <- as_operator(op)
  identical(native_kernel_kind(op), "dense") ||
    identical(native_kernel_kind(op), "csc")
}

#' @keywords internal
native_refined_arnoldi_available <- function(op) {
  native_arnoldi_available(op)
}

#' @keywords internal
native_matrix_free_arnoldi_available <- function(op) {
  op <- as_operator(op)
  is.null(source_or_null(op)) &&
    is.null(op$metadata$matrix) &&
    is.function(op$apply) &&
    identical(op$dtype, "double") &&
    op$dim[[1L]] == op$dim[[2L]]
}

#' @keywords internal
native_arnoldi_default_max_subspace <- function(n, k) {
  n <- as.integer(n)
  k <- as.integer(k)
  min(n, max(k + 8L, 9L * k))
}

#' @keywords internal
native_arnoldi_cycle <- function(op, start, m) {
  op <- as_operator(op)
  kind <- native_kernel_kind(op)
  if (identical(kind, "dense")) {
    source <- source_or_null(op)
    return(.Call(
      "eigencore_arnoldi_dense_cycle",
      source,
      as.numeric(start),
      as.integer(m),
      PACKAGE = "eigencore"
    ))
  }
  if (identical(kind, "csc")) {
    source <- op$metadata$matrix
    return(.Call(
      "eigencore_arnoldi_csc_cycle",
      methods::slot(source, "i"),
      methods::slot(source, "p"),
      methods::slot(source, "x"),
      methods::slot(source, "Dim"),
      as.numeric(start),
      as.integer(m),
      PACKAGE = "eigencore"
    ))
  }
  if (native_matrix_free_arnoldi_available(op)) {
    return(.Call(
      "eigencore_arnoldi_r_operator_cycle",
      as.integer(op$dim),
      op$apply,
      as.numeric(start),
      as.integer(m),
      PACKAGE = "eigencore"
    ))
  }
  stop("native Arnoldi requires a dense double, dgCMatrix, or real matrix-free operator.", call. = FALSE)
}

#' @keywords internal
native_arnoldi_projected_ritz <- function(cycle) {
  .Call(
    "eigencore_arnoldi_ritz",
    cycle$V,
    cycle$H,
    as.integer(cycle$iterations),
    PACKAGE = "eigencore"
  )
}

#' @keywords internal
native_arnoldi_refined_ritz_vectors <- function(cycle, values) {
  .Call(
    "eigencore_arnoldi_refined_ritz",
    cycle$V,
    cycle$H,
    as.integer(cycle$iterations),
    as.complex(values),
    PACKAGE = "eigencore"
  )
}

#' @keywords internal
native_arnoldi_general <- function(op, k, target = largest(), tol = 1e-8,
                                   maxit = NULL, max_restarts = 0L,
                                   vectors = TRUE,
                                   extraction = "projected_ritz") {
  op <- as_operator(op)
  extraction <- match.arg(extraction, c("projected_ritz", "refined_ritz"))
  if (op$dim[[1L]] != op$dim[[2L]]) {
    stop("native Arnoldi requires a square operator.", call. = FALSE)
  }
  matrix_free_native <- native_matrix_free_arnoldi_available(op)
  if (!native_arnoldi_available(op) && !matrix_free_native) {
    stop("native Arnoldi requires a dense double, dgCMatrix, or real matrix-free operator.", call. = FALSE)
  }
  if (matrix_free_native && identical(extraction, "refined_ritz")) {
    stop(
      "native refined Ritz extraction currently supports dense double and dgCMatrix operators only; matrix-free refined Ritz extraction is future scope.",
      call. = FALSE
    )
  }
  if (!reference_arnoldi_target_supported(target)) {
    stop("native Arnoldi currently supports largest/smallest real-part and largest-magnitude targets.",
         call. = FALSE)
  }
  n <- op$dim[[1L]]
  k <- as.integer(k)
  m <- as.integer(maxit %||% native_arnoldi_default_max_subspace(n, k))
  m <- min(n, max(k + 1L, m))
  max_restarts <- as.integer(max_restarts %||% 0L)
  if (max_restarts < 0L) {
    stop("max_restarts must be non-negative.", call. = FALSE)
  }

  start <- stats::rnorm(n)
  start <- start / sqrt(sum(start^2))
  best <- NULL
  best_score <- NULL
  selected_attempt <- NA_integer_
  history <- vector("list", max_restarts + 1L)
  total_matvecs <- 0L
  total_iterations <- 0L
  total_reorthogonalization_passes <- 0L
  native_workspace_bytes <- 0L
  total_cycle_seconds <- 0
  total_ritz_extraction_seconds <- 0

  for (attempt in seq_len(max_restarts + 1L)) {
    cycle_start <- proc.time()[["elapsed"]]
    cycle <- native_arnoldi_cycle(op, start, m)
    cycle_seconds <- proc.time()[["elapsed"]] - cycle_start
    total_cycle_seconds <- total_cycle_seconds + cycle_seconds
    total_matvecs <- total_matvecs + cycle$matvecs
    total_iterations <- total_iterations + cycle$iterations
    total_reorthogonalization_passes <- total_reorthogonalization_passes +
      cycle$reorthogonalization_passes
    native_workspace_bytes <- native_workspace_bytes + cycle$native_workspace_bytes
    ritz_start <- proc.time()[["elapsed"]]
    ritz <- native_arnoldi_ritz(op, cycle, k, target, tol, extraction = extraction)
    ritz_seconds <- proc.time()[["elapsed"]] - ritz_start
    total_ritz_extraction_seconds <- total_ritz_extraction_seconds + ritz_seconds
    history[[attempt]] <- data.frame(
      attempt = attempt,
      extraction = ritz$extraction,
      max_subspace = m,
      iterations = cycle$iterations,
      matvecs = cycle$matvecs,
      cycle_seconds = cycle_seconds,
      ritz_extraction_seconds = ritz_seconds,
      certificate_passed = isTRUE(ritz$certificate$passed),
      nconv = sum(ritz$certificate$converged),
      max_backward_error = ritz$certificate$max_backward_error,
      max_residual = ritz$certificate$max_residual,
      min_refined_residual_estimate = if (length(ritz$refined_residual_estimates %||% numeric())) {
        min(ritz$refined_residual_estimates)
      } else {
        NA_real_
      },
      stringsAsFactors = FALSE
    )
    score <- reference_arnoldi_score(ritz$certificate)
    if (is.null(best_score) || reference_arnoldi_score_better(score, best_score)) {
      best <- ritz
      best_score <- score
      selected_attempt <- attempt
    }
    if (isTRUE(ritz$certificate$passed)) {
      break
    }
    start <- Re(ritz$vectors[, 1L])
    if (!all(is.finite(start)) || sum(start^2) <= 100 * .Machine$double.eps) {
      start <- stats::rnorm(n)
    }
    start <- start / sqrt(sum(start^2))
  }

  kept_history <- history[seq_len(length(Filter(Negate(is.null), history)))]
  attempt_history <- do.call(rbind, kept_history)
  list(
    values = best$values,
    vectors = if (isTRUE(vectors)) best$vectors else NULL,
    certificate = best$certificate,
    iterations = total_iterations,
    matvecs = total_matvecs,
    restart = list(
      kind = if (matrix_free_native) {
        "native_matrix_free_arnoldi_callback_cycle"
      } else {
        "native_arnoldi_cycle"
      },
      implemented = TRUE,
      native = TRUE,
      matrix_free = isTRUE(matrix_free_native),
      ritz_extraction_native = TRUE,
      extraction = extraction,
      refined_extraction_native = identical(extraction, "refined_ritz"),
      refined_residual_estimates = best$refined_residual_estimates %||% numeric(),
      krylov_schur = FALSE,
      krylov_schur_status = "not implemented; V2 tranche promotes native refined Ritz extraction only",
      v2_issue = "bd-01KTF6H41S9XDN286TR3V184P4",
      max_subspace = m,
      max_restarts = max_restarts,
      restart_count = nrow(attempt_history) - 1L,
      attempted_subspaces = attempt_history$max_subspace,
      attempt_history = attempt_history,
      selected_attempt = selected_attempt,
      target_supported = TRUE,
      certified_attempt = if (isTRUE(best$certificate$passed)) nrow(attempt_history) else NA_integer_,
      reorthogonalization_passes = total_reorthogonalization_passes,
      native_workspace_bytes = native_workspace_bytes,
      stage_seconds = c(
        cycle = total_cycle_seconds,
        ritz_extraction = total_ritz_extraction_seconds
      )
    )
  )
}

#' @keywords internal
reference_arnoldi_general <- function(op, k, target = largest(), tol = 1e-8,
                                      maxit = NULL, max_restarts = 0L,
                                      vectors = TRUE,
                                      extraction = "projected_ritz") {
  op <- as_operator(op)
  if (op$dim[[1L]] != op$dim[[2L]]) {
    stop("reference Arnoldi requires a square operator.", call. = FALSE)
  }
  if (!reference_arnoldi_target_supported(target)) {
    stop("reference Arnoldi currently supports largest/smallest real-part and largest-magnitude targets.",
         call. = FALSE)
  }
  n <- op$dim[[1L]]
  k <- as.integer(k)
  # Smaller default subspace than native_arnoldi_general (which uses 9*k).
  # The reference path is an oracle fallback: it is called when no native
  # kernel is available, so throughput is limited by R-level matrix-vector
  # products.  A subspace of roughly 2*k is large enough for certification on
  # well-conditioned problems and keeps memory and per-restart cost low.  The
  # native path can afford the larger 9*k subspace because the cycle runs in
  # compiled C++ and restarts are cheap.
  m <- as.integer(maxit %||% min(n, max(k + 8L, 2L * k + 4L)))
  m <- min(n, max(k + 1L, m))
  max_restarts <- as.integer(max_restarts %||% 0L)
  if (max_restarts < 0L) {
    stop("max_restarts must be non-negative.", call. = FALSE)
  }

  start <- stats::rnorm(n)
  start <- start / sqrt(sum(start^2))
  best <- NULL
  best_score <- NULL
  selected_attempt <- NA_integer_
  history <- vector("list", max_restarts + 1L)
  total_matvecs <- 0L
  total_iterations <- 0L
  total_cycle_seconds <- 0
  total_ritz_extraction_seconds <- 0

  for (attempt in seq_len(max_restarts + 1L)) {
    cycle_start <- proc.time()[["elapsed"]]
    cycle <- reference_arnoldi_cycle(op, start, m)
    cycle_seconds <- proc.time()[["elapsed"]] - cycle_start
    total_cycle_seconds <- total_cycle_seconds + cycle_seconds
    total_matvecs <- total_matvecs + cycle$matvecs
    total_iterations <- total_iterations + cycle$iterations
    ritz_start <- proc.time()[["elapsed"]]
    ritz <- reference_arnoldi_ritz(op, cycle, k, target, tol)
    ritz_seconds <- proc.time()[["elapsed"]] - ritz_start
    total_ritz_extraction_seconds <- total_ritz_extraction_seconds + ritz_seconds
    history[[attempt]] <- data.frame(
      attempt = attempt,
      extraction = "projected_ritz",
      max_subspace = m,
      iterations = cycle$iterations,
      matvecs = cycle$matvecs,
      cycle_seconds = cycle_seconds,
      ritz_extraction_seconds = ritz_seconds,
      certificate_passed = isTRUE(ritz$certificate$passed),
      nconv = sum(ritz$certificate$converged),
      max_backward_error = ritz$certificate$max_backward_error,
      max_residual = ritz$certificate$max_residual,
      stringsAsFactors = FALSE
    )
    score <- reference_arnoldi_score(ritz$certificate)
    if (is.null(best_score) || reference_arnoldi_score_better(score, best_score)) {
      best <- ritz
      best_score <- score
      selected_attempt <- attempt
    }
    if (isTRUE(ritz$certificate$passed)) {
      break
    }
    start <- Re(ritz$vectors[, 1L])
    if (!all(is.finite(start)) || sum(start^2) <= 100 * .Machine$double.eps) {
      start <- stats::rnorm(n)
    }
    start <- start / sqrt(sum(start^2))
  }

  kept_history <- history[seq_len(length(Filter(Negate(is.null), history)))]
  attempt_history <- do.call(rbind, kept_history)
  list(
    values = best$values,
    vectors = if (isTRUE(vectors)) best$vectors else NULL,
    certificate = best$certificate,
    iterations = total_iterations,
    matvecs = total_matvecs,
    restart = list(
      kind = "reference_arnoldi",
      implemented = TRUE,
      native = FALSE,
      max_subspace = m,
      max_restarts = max_restarts,
      restart_count = nrow(attempt_history) - 1L,
      attempted_subspaces = attempt_history$max_subspace,
      attempt_history = attempt_history,
      selected_attempt = selected_attempt,
      target_supported = TRUE,
      certified_attempt = if (isTRUE(best$certificate$passed)) nrow(attempt_history) else NA_integer_,
      stage_seconds = c(
        cycle = total_cycle_seconds,
        ritz_extraction = total_ritz_extraction_seconds
      )
    )
  )
}

#' @keywords internal
reference_arnoldi_cycle <- function(op, start, m) {
  n <- op$dim[[1L]]
  V <- matrix(0, n, m + 1L)
  H <- matrix(0, m + 1L, m)
  V[, 1L] <- start / sqrt(sum(Mod(start)^2))
  iterations <- 0L
  matvecs <- 0L
  for (j in seq_len(m)) {
    w <- apply_operator(op, matrix(V[, j], n, 1L))[, 1L]
    matvecs <- matvecs + 1L
    for (i in seq_len(j)) {
      H[i, j] <- sum(Conj(V[, i]) * w)
      w <- w - H[i, j] * V[, i]
    }
    for (i in seq_len(j)) {
      corr <- sum(Conj(V[, i]) * w)
      H[i, j] <- H[i, j] + corr
      w <- w - corr * V[, i]
    }
    beta <- sqrt(sum(Mod(w)^2))
    H[j + 1L, j] <- beta
    iterations <- j
    if (!is.finite(beta) || beta <= 100 * .Machine$double.eps || j == m) {
      break
    }
    V[, j + 1L] <- w / beta
  }
  list(V = V, H = H, iterations = iterations, matvecs = matvecs)
}

#' @keywords internal
reference_arnoldi_score <- function(certificate) {
  error <- certificate$max_backward_error %||% Inf
  if (length(error) != 1L || is.na(error) || !is.finite(error)) {
    error <- Inf
  }
  list(
    passed = isTRUE(certificate$passed),
    nconv = sum(certificate$converged %||% FALSE),
    max_backward_error = error
  )
}

#' @keywords internal
reference_arnoldi_score_better <- function(candidate, incumbent) {
  if (!identical(candidate$passed, incumbent$passed)) {
    return(isTRUE(candidate$passed))
  }
  if (!identical(candidate$nconv, incumbent$nconv)) {
    return(candidate$nconv > incumbent$nconv)
  }
  candidate$max_backward_error < incumbent$max_backward_error
}

#' @keywords internal
reference_arnoldi_ritz <- function(op, cycle, k, target, tol) {
  m <- cycle$iterations
  Hm <- cycle$H[seq_len(m), seq_len(m), drop = FALSE]
  eig <- eigen(Hm)
  arnoldi_ritz_from_eigen(
    op, eig$values, eig$vectors, cycle$V, m, k, target, tol,
    vectors_are_ritz = FALSE
  )
}

#' @keywords internal
native_arnoldi_ritz <- function(op, cycle, k, target, tol,
                                extraction = c("projected_ritz", "refined_ritz")) {
  extraction <- match.arg(extraction)
  m <- cycle$iterations
  eig <- native_arnoldi_projected_ritz(cycle)
  if (identical(extraction, "refined_ritz")) {
    idx <- order_indices(eig$values, target)
    idx <- idx[seq_len(min(k, length(idx)))]
    values <- eig$values[idx]
    refined <- native_arnoldi_refined_ritz_vectors(cycle, values)
    return(arnoldi_ritz_from_eigen(
      op,
      values,
      NULL,
      cycle$V,
      m,
      k,
      target,
      tol,
      vectors_are_ritz = TRUE,
      vectors_override = refined$vectors,
      extraction = "refined_ritz",
      refined_residual_estimates = refined$refined_residual_estimates
    ))
  }
  arnoldi_ritz_from_eigen(
    op, eig$values, eig$vectors, cycle$V, m, k, target, tol,
    vectors_are_ritz = TRUE,
    extraction = "projected_ritz"
  )
}

#' @keywords internal
arnoldi_ritz_from_eigen <- function(op, eigenvalues, eigenvectors, V, m, k, target, tol,
                                    vectors_are_ritz,
                                    vectors_override = NULL,
                                    extraction = "projected_ritz",
                                    refined_residual_estimates = NULL) {
  eig <- list(values = eigenvalues, vectors = eigenvectors)
  idx <- order_indices(eig$values, target)
  idx <- idx[seq_len(min(k, length(idx)))]
  values <- eig$values[idx]
  if (!is.null(vectors_override)) {
    vectors <- vectors_override[, seq_len(length(idx)), drop = FALSE]
  } else if (isTRUE(vectors_are_ritz)) {
    vectors <- eig$vectors[, idx, drop = FALSE]
  } else {
    vectors <- V[, seq_len(m), drop = FALSE] %*%
      eig$vectors[, idx, drop = FALSE]
  }
  norms <- sqrt(colSums(Mod(vectors)^2))
  vectors <- sweep(vectors, 2L, pmax(norms, .Machine$double.eps), `/`)

  if (max(abs(Im(values))) <= sqrt(tol) &&
      max(abs(Im(vectors))) <= sqrt(tol)) {
    values <- Re(values)
    vectors <- Re(vectors)
  }

  cert <- tryCatch(
    certify_general_eigen_operator(op, values, vectors, tol = tol),
    error = function(e) {
      empty_certificate(
        tol,
        note = paste(
          "reference Arnoldi could not certify this Ritz basis with the current operator apply path:",
          conditionMessage(e)
        )
      )
    }
  )
  list(
    values = values,
    vectors = vectors,
    certificate = cert,
    extraction = extraction,
    refined_residual_estimates = refined_residual_estimates
  )
}
