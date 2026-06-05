#' @keywords internal
generalized_lanczos_label <- function() {
  "reference generalized SPD B-orthogonal Lanczos refinement"
}

#' @keywords internal
generalized_lanczos_target_supported <- function(target) {
  kind <- if (inherits(target, "eigencore_target")) target$kind else "largest"
  kind %in% c("largest", "smallest")
}

#' @keywords internal
generalized_lanczos_metric_solve_available <- function(Bop) {
  Bop <- as_operator(Bop)
  storage <- Bop$metadata$storage %||% NULL
  source <- source_or_null(Bop)
  (is.matrix(source) && is.double(source)) ||
    identical(storage, "ddiMatrix") ||
    identical(storage, "dgCMatrix")
}

#' @keywords internal
generalized_lanczos_supported <- function(op, Bop, target = smallest()) {
  op <- as_operator(op)
  Bop <- as_operator(Bop)
  identical(op$structure$kind, "hermitian") &&
    identical(Bop$structure$kind, "hermitian") &&
    generalized_lanczos_target_supported(target) &&
    generalized_spd_metric_known(Bop) &&
    generalized_lanczos_metric_solve_available(Bop)
}

#' @keywords internal
generalized_lanczos_plan_controls <- function(problem, k, method) {
  is_lanczos_method <- inherits(method, "eigencore_method") &&
    identical(method$kind, "lanczos")
  max_subspace <- if (is_lanczos_method) method$max_subspace else NULL
  n <- as.integer(problem$A$dim[1L])
  max_subspace <- as.integer(max_subspace %||% min(n, max(20L, 4L * as.integer(k) + 20L)))
  max_subspace <- min(n, max(as.integer(k), max_subspace))
  list(
    max_subspace = max_subspace,
    reorthogonalize = if (is_lanczos_method) isTRUE(method$reorthogonalize) else TRUE,
    generalized = TRUE,
    orthogonalization = "B-orthogonal Rayleigh-Ritz",
    metric_solve = generalized_lanczos_metric_solve_label(problem$metric),
    certification_policy = "residual certificate on original generalized SPD eigenproblem"
  )
}

#' @keywords internal
generalized_lanczos_metric_solve_label <- function(Bop) {
  Bop <- as_operator(Bop)
  storage <- Bop$metadata$storage %||% NULL
  source <- source_or_null(Bop)
  if (is.matrix(source) && is.double(source)) {
    "dense Cholesky solve for B"
  } else if (identical(storage, "ddiMatrix")) {
    "diagonal solve for B"
  } else if (identical(storage, "dgCMatrix")) {
    "sparse Cholesky solve for B"
  } else {
    "unavailable"
  }
}

#' @keywords internal
generalized_lanczos_prepare_metric_solve <- function(Bop) {
  Bop <- as_operator(Bop)
  if (!generalized_spd_metric_known(Bop)) {
    stop("generalized Lanczos metric B must be symmetric positive definite.", call. = FALSE)
  }
  source <- source_or_null(Bop)
  storage <- Bop$metadata$storage %||% NULL
  if (is.matrix(source) && is.double(source)) {
    B <- (source + t(source)) / 2
    factor <- chol(B)
    return(list(
      kind = "dense_cholesky",
      label = "dense Cholesky solve for B",
      solve = function(R) {
        backsolve(factor, forwardsolve(t(factor), as.matrix(R)))
      }
    ))
  }
  if (identical(storage, "ddiMatrix")) {
    B <- Bop$metadata$matrix
    values <- if (identical(methods::slot(B, "diag"), "U")) {
      rep(1, Bop$dim[1L])
    } else {
      methods::slot(B, "x")
    }
    if (length(values) != Bop$dim[1L] || any(!is.finite(values)) || any(values <= 0)) {
      stop("generalized Lanczos diagonal metric B must be positive finite.", call. = FALSE)
    }
    return(list(
      kind = "diagonal",
      label = "diagonal solve for B",
      solve = function(R) {
        sweep(as.matrix(R), 1L, values, `/`)
      }
    ))
  }
  if (identical(storage, "dgCMatrix")) {
    B <- Bop$metadata$matrix
    factor <- Matrix::Cholesky(B, LDL = FALSE)
    return(list(
      kind = "sparse_cholesky",
      label = "sparse Cholesky solve for B",
      solve = function(R) {
        as.matrix(Matrix::solve(factor, as.matrix(R)))
      }
    ))
  }
  stop("generalized Lanczos requires dense, diagonal, or sparse CSC SPD metric B.", call. = FALSE)
}

#' @keywords internal
reference_generalized_lanczos_hermitian <- function(op, Bop, k, target = smallest(),
                                                    tol = 1e-8, maxit = NULL,
                                                    vectors = TRUE,
                                                    reorthogonalize = TRUE) {
  op <- as_operator(op)
  Bop <- as_operator(Bop)
  if (!generalized_lanczos_supported(op, Bop, target = target)) {
    stop("generalized Lanczos requires a Hermitian problem with dense/diagonal/sparse SPD metric B.",
         call. = FALSE)
  }
  n <- as.integer(op$dim[1L])
  controls <- reference_scalar_subspace_controls(
    requested = k,
    requested_name = "k",
    limit = n,
    maxit = maxit,
    default_maxit = function(k) max(20L, 4L * k + 20L)
  )
  k <- controls$requested
  maxit <- controls$maxit
  metric <- generalized_lanczos_prepare_metric_solve(Bop)

  Q <- matrix(0, n, maxit)
  alpha <- numeric(maxit)
  beta <- numeric(maxit)
  q_prev <- numeric(n)
  q <- stats::rnorm(n)
  orth <- lobpcg_b_orthonormalize(matrix(q, ncol = 1L), Bop)
  q <- orth$Q[, 1L]
  orth_methods <- orth$method
  iterations <- 0L
  matvecs <- 0L
  metric_solves <- 0L
  final <- NULL
  history_max_relative_residual <- rep(NA_real_, maxit)
  history_nconv <- rep(NA_integer_, maxit)

  for (j in seq_len(maxit)) {
    iterations <- j
    Q[, j] <- q
    Aq <- apply_operator(op, matrix(q, ncol = 1L))[, 1L]
    matvecs <- matvecs + 1L
    z <- metric$solve(matrix(Aq, ncol = 1L))[, 1L]
    metric_solves <- metric_solves + 1L
    if (j > 1L) {
      z <- z - beta[[j - 1L]] * q_prev
    }
    alpha[[j]] <- sum(q * Aq)
    z <- z - alpha[[j]] * q

    if (isTRUE(reorthogonalize)) {
      Qj <- Q[, seq_len(j), drop = FALSE]
      for (pass in seq_len(2L)) {
        Bz <- apply_operator(Bop, matrix(z, ncol = 1L))
        z <- z - drop(Qj %*% crossprod(Qj, Bz))
      }
    }
    Bz <- apply_operator(Bop, matrix(z, ncol = 1L))[, 1L]
    beta[[j]] <- sqrt(max(0, sum(z * Bz)))

    if (j >= k) {
      final <- generalized_lanczos_ritz(
        op, Bop, Q[, seq_len(j), drop = FALSE], k, target, tol
      )
      history_max_relative_residual[[j]] <- final$certificate$max_backward_error %||% Inf
      history_nconv[[j]] <- sum(final$certificate$converged %||% FALSE)
      if (isTRUE(final$certificate$passed)) {
        break
      }
    }
    if (!is.finite(beta[[j]]) || beta[[j]] <= max(100 * .Machine$double.eps, tol * 1e-3)) {
      break
    }
    q_prev <- q
    q <- z / beta[[j]]
  }

  if (is.null(final)) {
    final <- generalized_lanczos_ritz(
      op, Bop, Q[, seq_len(iterations), drop = FALSE], k, target, tol
    )
    history_max_relative_residual[[iterations]] <- final$certificate$max_backward_error %||% Inf
    history_nconv[[iterations]] <- sum(final$certificate$converged %||% FALSE)
  }
  history <- data.frame(
    iteration = seq_len(iterations),
    max_relative_residual = history_max_relative_residual[seq_len(iterations)],
    nconv = history_nconv[seq_len(iterations)]
  )
  list(
    values = final$values,
    vectors = if (isTRUE(vectors)) final$vectors else NULL,
    residuals = final$certificate$residuals,
    backward_error = final$certificate$backward_error,
    orthogonality = final$certificate$orthogonality,
    certificate = final$certificate,
    iterations = iterations,
    matvecs = matvecs,
    metric_solves = metric_solves,
    convergence_history = history,
    generalized = TRUE,
    metric_solve = list(kind = metric$kind, label = metric$label),
    restart = list(
      kind = "generalized_b_orthogonal_lanczos",
      implemented = TRUE,
      native = FALSE,
      generalized = TRUE,
      max_subspace = maxit,
      metric_solve = metric$label,
      metric_solves = metric_solves,
      reorthogonalize = isTRUE(reorthogonalize),
      orthogonalization_methods = unique(c(orth_methods, "B_reorthogonalize")),
      alpha = alpha[seq_len(iterations)],
      beta = beta[seq_len(iterations)]
    )
  )
}

#' @keywords internal
generalized_lanczos_ritz <- function(op, Bop, Q, k, target, tol) {
  AQ <- apply_operator(op, Q)
  BQ <- apply_operator(Bop, Q)
  projected_A <- crossprod(Q, AQ)
  projected_A <- (projected_A + t(projected_A)) / 2
  projected_B <- crossprod(Q, BQ)
  projected_B <- (projected_B + t(projected_B)) / 2
  eig <- dense_generalized_spd_eigen(projected_A, projected_B, vectors = TRUE)
  idx <- order_indices(eig$values, target)
  idx <- idx[seq_len(min(k, length(idx)))]
  values <- eig$values[idx]
  vectors <- Q %*% eig$vectors[, idx, drop = FALSE]
  orth <- lobpcg_b_orthonormalize(vectors, Bop)
  vectors <- orth$Q
  cert <- certify_eigen_operator(op, values, vectors, Bop = Bop, tol = tol)
  list(values = values, vectors = vectors, certificate = cert)
}
