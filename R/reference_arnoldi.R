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
#' Targets the native Krylov-Schur kernel ranks directly: the reference set
#' plus smallest magnitude (C41; converges slowly when the wanted values are
#' interior in modulus -- `nearest(0)` via shift-invert is the fast route).
native_arnoldi_target_supported <- function(target) {
  kind <- if (inherits(target, "eigencore_target")) target$kind else "largest"
  reference_arnoldi_target_supported(target) || identical(kind, "smallest_magnitude")
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
                                        extraction = "projected_ritz",
                                        krylov_schur_maxit = native_krylov_schur_default_maxit()) {
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
  native_left <- native_arnoldi_available(adjoint_op) ||
    native_matrix_free_arnoldi_available(adjoint_op)
  left_iter <- tryCatch(
    if (native_left) {
      # The right solve already fixed which eigenvalues are wanted, so the
      # adjoint Krylov-Schur ranks Ritz values by distance to them instead of
      # re-deciding which ones are extremal (on clustered spectra the two
      # independent solves could otherwise settle on different sets).
      native_arnoldi_general(
        adjoint_op,
        k = k,
        target = adjoint_arnoldi_target(target),
        tol = tol,
        maxit = maxit,
        max_restarts = max_restarts,
        vectors = TRUE,
        extraction = extraction,
        target_values = if (identical(op$dtype, "complex")) Conj(values) else values,
        # Biorthogonality error scales like (left residual) / (eigenvalue
        # gap), so converge the adjoint Ritz pairs two digits beyond `tol`.
        krylov_schur_tol = max(tol * 1e-2, 100 * .Machine$double.eps),
        krylov_schur_maxit = krylov_schur_maxit
      )
    } else {
      reference_arnoldi_general(
        adjoint_op,
        k = k,
        target = adjoint_arnoldi_target(target),
        tol = tol,
        maxit = maxit,
        max_restarts = max_restarts,
        vectors = TRUE,
        extraction = extraction
      )
    },
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
    arnoldi_certificate_operator(op),
    values,
    left_vectors,
    right_vectors = right_vectors,
    tol = tol
  )
  biorthogonality <- crossprod(left_vectors, as.matrix(right_vectors))
  left_certification_columns <- as.integer(
    (left_iter$certification_operator_columns %||% k) + k
  )
  left_certification_block_calls <- as.integer(
    (left_iter$certification_operator_block_calls %||% 1L) + 1L
  )
  list(
    supported = TRUE,
    reason = "left eigenvectors computed from the adjoint operator",
    values = matched$values,
    match_distance = matched$match_distance,
    vectors = left_vectors,
    certificate = cert,
    biorthogonality = biorthogonality,
    method = left_iter$restart$kind %||% "adjoint_arnoldi",
    adjoint_block_calls = as.integer(left_iter$matvecs %||% 0L),
    adjoint_columns = as.integer(left_iter$matvecs %||% 0L),
    certification_adjoint_block_calls = left_certification_block_calls,
    certification_adjoint_columns = left_certification_columns
  )
}

#' @keywords internal
#' Resolve the native kernel an Arnoldi run can use for `op`: a dense double
#' source, a dgCMatrix, or the materialized transpose carried by the adjoint
#' of a dgCMatrix operator (so the left-eigenvector solve on A^T stays native
#' instead of falling back to the R-level reference Arnoldi).
native_arnoldi_kernel <- function(op) {
  op <- as_operator(op)
  kind <- native_kernel_kind(op)
  if (identical(kind, "dense")) {
    return(list(kind = "dense", matrix = source_or_null(op)))
  }
  if (identical(kind, "csc")) {
    return(list(kind = "csc", matrix = op$metadata$matrix, transpose = NULL))
  }
  if (identical(op$metadata$storage %||% NULL, "adjoint:dgCMatrix") &&
      identical(op$dtype, "double") &&
      inherits(op$metadata$matrix, "dgCMatrix")) {
    parent_matrix <- op$metadata$parent$metadata$matrix %||% NULL
    return(list(
      kind = "csc",
      matrix = op$metadata$matrix,
      transpose = if (inherits(parent_matrix, "dgCMatrix")) parent_matrix else NULL
    ))
  }
  NULL
}

#' @keywords internal
#' CSC storage of the transpose of the kernel's operator. The Krylov-Schur
#' kernel applies A x as (A^T)^T x, a row gather that is faster than the
#' scatter form of a forward CSC product; the adjoint of a dgCMatrix operator
#' already carries its parent, so only the forward case pays one transpose.
native_arnoldi_kernel_transpose <- function(kernel) {
  tr <- kernel$transpose %||% NULL
  if (inherits(tr, "dgCMatrix")) {
    return(tr)
  }
  methods::as(Matrix::t(kernel$matrix), "CsparseMatrix")
}

#' @keywords internal
native_arnoldi_available <- function(op) {
  !is.null(native_arnoldi_kernel(op))
}

#' @keywords internal
#' Certificates for complex Ritz pairs of a real sparse operator must not
#' densify the source matrix (`as.matrix(source) %*% vectors`). This view of
#' `op` drops the materialized matrix and applies the original operator to the
#' real and imaginary parts as one real block, so the residuals are still
#' computed against the original operator with real matvecs. Norm metadata is
#' kept, so the certificate scale is unchanged.
arnoldi_certificate_operator <- function(op) {
  op <- as_operator(op)
  if (!is.null(source_or_null(op)) || is.null(op$metadata$matrix) ||
      !identical(op$dtype, "double")) {
    return(op)
  }
  split_complex <- function(f) {
    force(f)
    function(X, alpha = 1, beta = 0, Y = NULL) {
      if (!is.complex(X)) {
        return(f(X, alpha = alpha, beta = beta, Y = Y))
      }
      X <- as.matrix(X)
      p <- ncol(X)
      out <- as.matrix(f(cbind(Re(X), Im(X))))
      res <- matrix(
        complex(
          real = out[, seq_len(p), drop = FALSE],
          imaginary = out[, p + seq_len(p), drop = FALSE]
        ),
        nrow(out),
        p
      )
      res <- alpha * res
      if (!is.null(Y) && beta != 0) {
        res <- res + beta * Y
      }
      res
    }
  }
  view <- op
  view$apply <- split_complex(op$apply)
  if (!is.null(op$apply_adjoint)) {
    view$apply_adjoint <- split_complex(op$apply_adjoint)
  }
  # The adjoint of a dgCMatrix operator carries the materialized transpose but
  # no norm metadata; its exact Frobenius norm is the same as the parent's, so
  # record it instead of letting the certificate fall back to a stochastic
  # estimate (which withholds `passed`).
  matrix <- op$metadata$matrix
  if (is.null(view$metadata$frobenius_norm) && inherits(matrix, "dgCMatrix")) {
    view$metadata$frobenius_norm <- sqrt(sum(methods::slot(matrix, "x")^2))
  }
  view$metadata$matrix <- NULL
  view
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
#' Sparse general-pencil transformed Arnoldi needs a larger Krylov budget than
#' standard nonsymmetric partial solves: refined Ritz extraction on B^{-1} A
#' with random starts was flaky at the default 9*k subspace on moderate n.
sparse_general_pencil_default_max_subspace <- function(n, k) {
  n <- as.integer(n)
  k <- as.integer(k)
  min(
    n,
    max(
      native_arnoldi_default_max_subspace(n, k),
      min(n, 30L + 5L * (k - 1L))
    )
  )
}

#' @keywords internal
#' Krylov-Schur subspace size (the ARPACK/RSpectra `ncv` default for the
#' nonsymmetric problem): max(2k + 1, 20), capped at n. Dense inputs use the
#' same restarted subspace instead of an n-dimensional basis.
native_krylov_schur_default_ncv <- function(n, k) {
  n <- as.integer(n)
  k <- as.integer(k)
  min(n, max(2L * k + 1L, 20L))
}

#' @keywords internal
#' Cap on Krylov-Schur outer iterations (restarts) per attempt.
native_krylov_schur_default_maxit <- function() {
  1000L
}

#' @keywords internal
arnoldi_target_code <- function(target) {
  kind <- if (inherits(target, "eigencore_target")) target$kind else "largest"
  code <- switch(
    kind,
    largest = 0L,
    largest_real = 0L,
    smallest = 1L,
    smallest_real = 1L,
    largest_magnitude = 2L,
    largest_imaginary = 3L,
    smallest_imaginary = 4L,
    smallest_magnitude = 5L,
    NA_integer_
  )
  if (is.na(code)) {
    stop("native Krylov-Schur Arnoldi does not support target '", kind, "'.",
         call. = FALSE)
  }
  code
}

#' @keywords internal
native_arnoldi_cycle <- function(op, start, m) {
  op <- as_operator(op)
  kernel <- native_arnoldi_kernel(op)
  if (identical(kernel$kind, "dense")) {
    return(.Call(
      "eigencore_arnoldi_dense_cycle",
      kernel$matrix,
      as.numeric(start),
      as.integer(m),
      PACKAGE = "eigencore"
    ))
  }
  if (identical(kernel$kind, "csc")) {
    source <- kernel$matrix
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
#' Native real Krylov-Schur restarted Arnoldi (src/arnoldi.cpp). Returns a
#' cycle-like list whose V (n x (p+1)) and H = [T_p; b^T] ((p+1) x p) satisfy
#' A V[, 1:p] = V H, so the projected and refined Ritz extractions apply
#' unchanged to the converged Krylov-Schur decomposition.
native_krylov_schur <- function(op, start, k, m, target, tol,
                                maxit = native_krylov_schur_default_maxit(),
                                target_values = NULL) {
  op <- as_operator(op)
  kernel <- native_arnoldi_kernel(op)
  code <- arnoldi_target_code(target)
  tol <- as.numeric(tol)
  # target_values (internal): rank Ritz values by distance to this set
  # instead of by `target` (used by the adjoint solve for left vectors).
  target_values <- if (length(target_values)) as.complex(target_values) else NULL
  if (identical(kernel$kind, "dense")) {
    return(.Call(
      "eigencore_arnoldi_ks_dense",
      kernel$matrix,
      as.numeric(start),
      as.integer(k),
      as.integer(m),
      code,
      target_values,
      tol,
      as.integer(maxit),
      PACKAGE = "eigencore"
    ))
  }
  if (identical(kernel$kind, "csc")) {
    source <- native_arnoldi_kernel_transpose(kernel)
    if (!inherits(source, "dgCMatrix")) {
      stop("native Krylov-Schur Arnoldi requires a dgCMatrix transpose.", call. = FALSE)
    }
    return(.Call(
      "eigencore_arnoldi_ks_csc",
      methods::slot(source, "i"),
      methods::slot(source, "p"),
      methods::slot(source, "x"),
      methods::slot(source, "Dim"),
      TRUE,
      as.numeric(start),
      as.integer(k),
      as.integer(m),
      code,
      target_values,
      tol,
      as.integer(maxit),
      PACKAGE = "eigencore"
    ))
  }
  if (native_matrix_free_arnoldi_available(op)) {
    return(.Call(
      "eigencore_arnoldi_ks_r_operator",
      as.integer(op$dim),
      op$apply,
      as.numeric(start),
      as.integer(k),
      as.integer(m),
      code,
      target_values,
      tol,
      as.integer(maxit),
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
native_arnoldi_ritz_coefficients <- function(cycle) {
  .Call(
    "eigencore_arnoldi_ritz_coefficients",
    cycle$H,
    as.integer(cycle$iterations),
    PACKAGE = "eigencore"
  )
}

#' @keywords internal
native_arnoldi_ritz_vectors <- function(cycle, coefficients) {
  .Call(
    "eigencore_arnoldi_ritz_vectors",
    cycle$V,
    as.integer(cycle$iterations),
    as.matrix(coefficients) + 0i,
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
#' Start vector for the next certification attempt: the real part of the sum
#' of the current wanted Ritz vectors, so the new Krylov space contains all of
#' them; random when that is degenerate.
arnoldi_next_start <- function(vectors, n) {
  start <- if (!is.null(vectors) && length(vectors)) {
    Re(rowSums(as.matrix(vectors))) + Im(rowSums(as.matrix(vectors)))
  } else {
    numeric()
  }
  norm <- if (length(start) == n && all(is.finite(start))) sqrt(sum(start^2)) else 0
  if (!is.finite(norm) || norm <= 0) {
    start <- stats::rnorm(n)
    norm <- sqrt(sum(start^2))
  }
  start / norm
}

#' @keywords internal
native_arnoldi_general <- function(op, k, target = largest(), tol = 1e-8,
                                   maxit = NULL, max_restarts = 0L,
                                   vectors = TRUE,
                                   extraction = "projected_ritz",
                                   krylov_schur_maxit = native_krylov_schur_default_maxit(),
                                   target_values = NULL,
                                   krylov_schur_tol = tol) {
  op <- as_operator(op)
  extraction <- match.arg(extraction, c("projected_ritz", "refined_ritz"))
  if (op$dim[[1L]] != op$dim[[2L]]) {
    stop("native Arnoldi requires a square operator.", call. = FALSE)
  }
  matrix_free_native <- !native_arnoldi_available(op) &&
    native_matrix_free_arnoldi_available(op)
  if (!native_arnoldi_available(op) && !matrix_free_native) {
    stop("native Arnoldi requires a dense double, dgCMatrix, or real matrix-free operator.", call. = FALSE)
  }
  if (matrix_free_native && identical(extraction, "refined_ritz")) {
    stop(
      "native refined Ritz extraction currently supports dense double and dgCMatrix operators only; matrix-free refined Ritz extraction is future scope.",
      call. = FALSE
    )
  }
  if (!native_arnoldi_target_supported(target)) {
    stop("native Arnoldi currently supports largest/smallest real-part and largest-magnitude targets.",
         call. = FALSE)
  }
  n <- op$dim[[1L]]
  k <- as.integer(k)
  if (length(k) != 1L || is.na(k) || k < 1L || k > n) {
    stop("k must be between 1 and the operator dimension.", call. = FALSE)
  }
  m <- as.integer(maxit %||% native_krylov_schur_default_ncv(n, k))
  # Krylov-Schur needs room for at least one new direction beyond the wanted
  # block (plus a possible conjugate partner) unless the space is all of R^n.
  m <- min(n, max(k + 2L, m))
  max_restarts <- as.integer(max_restarts %||% 0L)
  if (max_restarts < 0L) {
    stop("max_restarts must be non-negative.", call. = FALSE)
  }
  krylov_schur_maxit <- as.integer(krylov_schur_maxit)

  start <- stats::rnorm(n)
  start <- start / sqrt(sum(start^2))
  best <- NULL
  best_score <- NULL
  selected_attempt <- NA_integer_
  history <- vector("list", max_restarts + 1L)
  total_matvecs <- 0L
  total_iterations <- 0L
  total_ks_restarts <- 0L
  ks_stage_seconds <- c(apply = 0, orthogonalization = 0,
                        projected_schur = 0, restart = 0)
  total_reorthogonalization_passes <- 0L
  native_workspace_bytes <- 0L
  total_cycle_seconds <- 0
  total_ritz_extraction_seconds <- 0
  m_initial <- m

  for (attempt in seq_len(max_restarts + 1L)) {
    cycle_start <- proc.time()[["elapsed"]]
    cycle <- native_krylov_schur(
      op, start, k = k, m = m, target = target, tol = krylov_schur_tol,
      maxit = krylov_schur_maxit, target_values = target_values
    )
    cycle_seconds <- proc.time()[["elapsed"]] - cycle_start
    total_cycle_seconds <- total_cycle_seconds + cycle_seconds
    total_matvecs <- total_matvecs + cycle$matvecs
    total_iterations <- total_iterations + cycle$krylov_schur_restarts + 1L
    total_ks_restarts <- total_ks_restarts + cycle$krylov_schur_restarts
    if (length(cycle$stage_seconds)) {
      ks_stage_seconds <- ks_stage_seconds +
        unname(cycle$stage_seconds[names(ks_stage_seconds)])
    }
    total_reorthogonalization_passes <- total_reorthogonalization_passes +
      cycle$reorthogonalization_passes
    native_workspace_bytes <- native_workspace_bytes + cycle$native_workspace_bytes
    ritz_start <- proc.time()[["elapsed"]]
    ritz <- native_arnoldi_ritz(op, cycle, k, target, tol, extraction = extraction,
                                target_values = target_values)
    ritz_seconds <- proc.time()[["elapsed"]] - ritz_start
    total_ritz_extraction_seconds <- total_ritz_extraction_seconds + ritz_seconds
    history[[attempt]] <- data.frame(
      attempt = attempt,
      extraction = ritz$extraction,
      max_subspace = m,
      iterations = cycle$krylov_schur_restarts + 1L,
      matvecs = cycle$matvecs,
      krylov_schur_restarts = cycle$krylov_schur_restarts,
      krylov_schur_nconv = cycle$nconv,
      krylov_schur_converged = isTRUE(cycle$converged),
      retained_subspace = cycle$iterations,
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
    # Every residual already meets the tolerance, so `passed` was withheld for
    # a reason a larger subspace cannot fix (e.g. a stochastic norm estimate
    # for a matrix-free operator).
    converged_flags <- ritz$certificate$converged %||% logical()
    if (length(converged_flags) >= k && all(converged_flags)) {
      break
    }
    # Certification failed: retry from the current wanted Ritz vectors with a
    # larger Krylov-Schur subspace.
    start <- arnoldi_next_start(ritz$vectors, n)
    m <- min(n, max(m + k, 2L * m))
  }

  kept_history <- history[seq_along(Filter(Negate(is.null), history))]
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
      krylov_schur = TRUE,
      krylov_schur_status = native_krylov_schur_status(),
      krylov_schur_restarts = total_ks_restarts,
      krylov_schur_stage_seconds = ks_stage_seconds,
      krylov_schur_maxit = krylov_schur_maxit,
      v2_issue = "bd-01KTF6H41S9XDN286TR3V184P4",
      max_subspace = m_initial,
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
native_krylov_schur_status <- function() {
  paste(
    "native real Krylov-Schur restart (dgees + dtrsen reordering, CGS2",
    "orthogonalization, conjugate pairs kept together)"
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

  kept_history <- history[seq_along(Filter(Negate(is.null), history))]
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
  # Running ||A v_j|| scale so the breakdown test is relative (scale-invariant).
  scale <- 0
  for (j in seq_len(m)) {
    w <- apply_operator(op, matrix(V[, j], n, 1L))[, 1L]
    matvecs <- matvecs + 1L
    scale <- max(scale, sqrt(sum(Mod(w)^2)))
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
    breakdown <- is.finite(beta) && beta <= 100 * .Machine$double.eps * scale
    H[j + 1L, j] <- if (breakdown) 0 else beta
    iterations <- j
    if (!is.finite(beta) || breakdown || j == m) {
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
    vectors_are_ritz = FALSE,
    value_scale = arnoldi_projected_scale(cycle$H)
  )
}

#' @keywords internal
#' Ordering of Ritz values: by `target`, or by distance to `target_values`
#' when the caller already knows which eigenvalues it wants (adjoint solve).
arnoldi_order_indices <- function(values, target, target_values = NULL) {
  if (!length(target_values)) {
    return(order_indices(values, target))
  }
  distance <- vapply(
    values,
    function(v) min(Mod(v - target_values)),
    numeric(1)
  )
  order(distance)
}

#' @keywords internal
arnoldi_projected_scale <- function(H) {
  H <- as.matrix(H)
  if (!length(H)) {
    return(0)
  }
  value <- sqrt(sum(Mod(H)^2))
  if (is.finite(value)) value else 0
}

#' @keywords internal
#' Projected (or refined) Ritz pairs for the k wanted values of a native
#' Arnoldi/Krylov-Schur decomposition. Only the selected Ritz vectors are
#' formed (one dgemm against V), never all m of them.
native_arnoldi_ritz <- function(op, cycle, k, target, tol,
                                extraction = c("projected_ritz", "refined_ritz"),
                                target_values = NULL) {
  extraction <- match.arg(extraction)
  m <- cycle$iterations
  eig <- native_arnoldi_ritz_coefficients(cycle)
  idx <- arnoldi_order_indices(eig$values, target, target_values)
  idx <- idx[seq_len(min(k, length(idx)))]
  values <- eig$values[idx]
  value_scale <- arnoldi_projected_scale(cycle$H)
  if (identical(extraction, "refined_ritz")) {
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
      refined_residual_estimates = refined$refined_residual_estimates,
      value_scale = value_scale
    ))
  }
  vectors <- native_arnoldi_ritz_vectors(
    cycle,
    eig$coefficients[, idx, drop = FALSE]
  )
  arnoldi_ritz_from_eigen(
    op, values, NULL, cycle$V, m, k, target, tol,
    vectors_are_ritz = TRUE,
    vectors_override = vectors,
    extraction = "projected_ritz",
    value_scale = value_scale
  )
}

#' @keywords internal
#' Relative realness test for Ritz values of a real operator: an imaginary
#' part is treated as rounding when |Im(lambda)| <= 100 * eps * max(|lambda|,
#' scale), with `scale` the projected-matrix norm. LAPACK returns exact zeros
#' for real eigenvalues of a real matrix, so genuine conjugate pairs (however
#' close to the real axis) stay complex.
arnoldi_real_values <- function(values, scale = 0) {
  if (!is.complex(values)) {
    return(rep(TRUE, length(values)))
  }
  if (!length(values)) {
    return(logical())
  }
  scale <- max(c(scale, Mod(values)), na.rm = TRUE)
  abs(Im(values)) <= 100 * .Machine$double.eps * pmax(Mod(values), scale)
}

#' @keywords internal
#' Real representative of complex vectors attached to real eigenvalues of a
#' real operator: rotate each column so its largest entry is real and
#' positive, then keep the real part. For a real eigenvalue the real part of an
#' eigenvector is itself an eigenvector, so this only removes the arbitrary
#' complex phase LAPACK (e.g. zgesvd in refined extraction) may attach.
arnoldi_realify_vectors <- function(vectors) {
  vectors <- as.matrix(vectors)
  if (!is.complex(vectors)) {
    return(vectors)
  }
  out <- matrix(0, nrow(vectors), ncol(vectors))
  for (j in seq_len(ncol(vectors))) {
    v <- vectors[, j]
    pivot <- v[[which.max(Mod(v))]]
    if (Mod(pivot) > 0) {
      v <- v * (Conj(pivot) / Mod(pivot))
    }
    r <- Re(v)
    nr <- sqrt(sum(r^2))
    out[, j] <- if (is.finite(nr) && nr > 0) r / nr else r
  }
  out
}

#' @keywords internal
arnoldi_ritz_from_eigen <- function(op, eigenvalues, eigenvectors, V, m, k, target, tol,
                                    vectors_are_ritz,
                                    vectors_override = NULL,
                                    extraction = "projected_ritz",
                                    refined_residual_estimates = NULL,
                                    value_scale = NULL) {
  eig <- list(values = eigenvalues, vectors = eigenvectors)
  idx <- order_indices(eig$values, target)
  idx <- idx[seq_len(min(k, length(idx)))]
  values <- eig$values[idx]
  if (!is.null(vectors_override)) {
    # vectors_override columns follow `eigenvalues`, so reorder them together.
    vectors <- vectors_override[, idx, drop = FALSE]
    if (length(refined_residual_estimates) == length(eigenvalues)) {
      refined_residual_estimates <- refined_residual_estimates[idx]
    }
  } else if (isTRUE(vectors_are_ritz)) {
    vectors <- eig$vectors[, idx, drop = FALSE]
  } else {
    vectors <- V[, seq_len(m), drop = FALSE] %*%
      eig$vectors[, idx, drop = FALSE]
  }
  norms <- sqrt(colSums(Mod(vectors)^2))
  vectors <- sweep(vectors, 2L, pmax(norms, .Machine$double.eps), `/`)

  if (length(values) &&
      all(arnoldi_real_values(values, value_scale %||% 0))) {
    values <- Re(values)
    vectors <- arnoldi_realify_vectors(vectors)
  }

  cert <- tryCatch(
    certify_general_eigen_operator(
      arnoldi_certificate_operator(op), values, vectors, tol = tol
    ),
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
