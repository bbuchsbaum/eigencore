#' @keywords internal
generalized_lanczos_label <- function() {
  "reference generalized SPD B-orthogonal Lanczos refinement"
}

#' @keywords internal
native_generalized_lanczos_label <- function() {
  "native transformed generalized SPD B-orthogonal Lanczos"
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
native_generalized_lanczos_supported <- function(op, Bop, target = smallest()) {
  op <- as_operator(op)
  Bop <- as_operator(Bop)
  if (!generalized_lanczos_supported(op, Bop, target = target)) {
    return(FALSE)
  }
  storage <- op$metadata$storage %||% NULL
  B_storage <- Bop$metadata$storage %||% NULL
  source <- source_or_null(op)
  B_source <- source_or_null(Bop)
  native_A <- identical(storage, "dgCMatrix") ||
    identical(storage, "ddiMatrix") ||
    (is.matrix(source) && is.double(source))
  diagonal_B <- identical(B_storage, "ddiMatrix")
  dense_B <- is.matrix(B_source) && is.double(B_source)
  native_A && (diagonal_B || (dense_B && is.matrix(source) && is.double(source)))
}

#' @keywords internal
generalized_lanczos_plan_controls <- function(problem, k, method) {
  is_lanczos_method <- inherits(method, "eigencore_method") &&
    identical(method$kind, "lanczos")
  max_subspace <- if (is_lanczos_method) method$max_subspace else NULL
  n <- as.integer(problem$A$dim[1L])
  max_subspace <- as.integer(max_subspace %||% min(n, max(20L, 4L * as.integer(k) + 20L)))
  max_subspace <- min(n, max(as.integer(k), max_subspace))
  native <- native_generalized_lanczos_supported(
    problem$A,
    problem$metric,
    target = problem$target
  )
  metric_solve <- generalized_lanczos_metric_solve_metadata(problem$metric)
  list(
    max_subspace = max_subspace,
    block = if (is_lanczos_method) as.integer(method$block %||% 1L) else 1L,
    max_restarts = if (is_lanczos_method) as.integer(method$max_restarts %||% 100L) else 100L,
    reorthogonalize = if (is_lanczos_method) isTRUE(method$reorthogonalize) else TRUE,
    generalized = TRUE,
    native = native,
    orthogonalization = if (native) {
      "B-orthogonal transformed native Lanczos"
    } else {
      "B-orthogonal Rayleigh-Ritz"
    },
    metric_solve = metric_solve$label,
    metric_solve_kind = metric_solve$kind,
    metric_solve_native = metric_solve$native,
    metric_transform = if (native) generalized_lanczos_metric_transform_label(problem$metric) else NULL,
    certification_policy = "residual certificate on original generalized SPD eigenproblem"
  )
}

#' @keywords internal
generalized_lanczos_sparse_tridiagonal_metric_parts <- function(Bop) {
  Bop <- as_operator(Bop)
  storage <- Bop$metadata$storage %||% NULL
  if (!identical(storage, "dgCMatrix")) {
    return(NULL)
  }
  B <- Bop$metadata$matrix %||% NULL
  if (!inherits(B, "CsparseMatrix")) {
    return(NULL)
  }
  parts <- shift_invert_tridiagonal_parts(B, shift = 0)
  if (is.null(parts) ||
      any(!is.finite(parts$diag)) ||
      any(!is.finite(parts$lower)) ||
      any(!is.finite(parts$upper))) {
    return(NULL)
  }
  parts
}

#' @keywords internal
generalized_lanczos_metric_solve_metadata <- function(Bop) {
  Bop <- as_operator(Bop)
  storage <- Bop$metadata$storage %||% NULL
  source <- source_or_null(Bop)
  if (is.matrix(source) && is.double(source)) {
    list(
      kind = "dense_cholesky",
      label = "dense Cholesky solve for B",
      native = FALSE,
      factorization = "dense_cholesky"
    )
  } else if (identical(storage, "ddiMatrix")) {
    list(
      kind = "diagonal",
      label = "diagonal solve for B",
      native = FALSE,
      factorization = "diagonal"
    )
  } else if (identical(storage, "dgCMatrix")) {
    native_parts <- generalized_lanczos_sparse_tridiagonal_metric_parts(Bop)
    if (!is.null(native_parts)) {
      list(
        kind = "native_sparse_tridiagonal_thomas",
        label = "native sparse tridiagonal Thomas solve for B",
        native = TRUE,
        factorization = "tridiagonal_thomas"
      )
    } else {
      list(
        kind = "sparse_cholesky",
        label = "sparse Cholesky solve for B",
        native = FALSE,
        factorization = "Matrix::Cholesky"
      )
    }
  } else {
    list(
      kind = "unavailable",
      label = "unavailable",
      native = FALSE,
      factorization = NA_character_
    )
  }
}

#' @keywords internal
generalized_lanczos_metric_solve_label <- function(Bop) {
  generalized_lanczos_metric_solve_metadata(Bop)$label
}

#' @keywords internal
generalized_lanczos_native_tridiagonal_solve <- function(parts, R) {
  .Call(
    "eigencore_tridiagonal_solve",
    as.numeric(parts$lower),
    as.numeric(parts$diag),
    as.numeric(parts$upper),
    as.matrix(R),
    PACKAGE = "eigencore"
  )
}

#' @keywords internal
generalized_lanczos_metric_transform_label <- function(Bop) {
  Bop <- as_operator(Bop)
  storage <- Bop$metadata$storage %||% NULL
  source <- source_or_null(Bop)
  if (is.matrix(source) && is.double(source)) {
    "dense Cholesky similarity transform for B"
  } else if (identical(storage, "ddiMatrix")) {
    "diagonal scaling similarity transform for B"
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
    metadata <- generalized_lanczos_metric_solve_metadata(Bop)
    return(list(
      kind = metadata$kind,
      label = metadata$label,
      native = metadata$native,
      factorization = metadata$factorization,
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
    metadata <- generalized_lanczos_metric_solve_metadata(Bop)
    return(list(
      kind = metadata$kind,
      label = metadata$label,
      native = metadata$native,
      factorization = metadata$factorization,
      solve = function(R) {
        sweep(as.matrix(R), 1L, values, `/`)
      }
    ))
  }
  if (identical(storage, "dgCMatrix")) {
    native_parts <- generalized_lanczos_sparse_tridiagonal_metric_parts(Bop)
    if (!is.null(native_parts)) {
      metadata <- generalized_lanczos_metric_solve_metadata(Bop)
      return(list(
        kind = metadata$kind,
        label = metadata$label,
        native = metadata$native,
        factorization = metadata$factorization,
        solve = function(R) {
          generalized_lanczos_native_tridiagonal_solve(native_parts, R)
        }
      ))
    }
    B <- Bop$metadata$matrix
    factor <- Matrix::Cholesky(B, LDL = FALSE)
    metadata <- generalized_lanczos_metric_solve_metadata(Bop)
    return(list(
      kind = metadata$kind,
      label = metadata$label,
      native = metadata$native,
      factorization = metadata$factorization,
      solve = function(R) {
        as.matrix(Matrix::solve(factor, as.matrix(R)))
      }
    ))
  }
  stop("generalized Lanczos requires dense, diagonal, or sparse CSC SPD metric B.", call. = FALSE)
}

#' @keywords internal
generalized_lanczos_diagonal_values <- function(Bop) {
  Bop <- as_operator(Bop)
  if (!identical(Bop$metadata$storage %||% NULL, "ddiMatrix")) {
    stop("diagonal metric transform requires a diagonal Matrix metric.", call. = FALSE)
  }
  B <- Bop$metadata$matrix
  values <- if (identical(methods::slot(B, "diag"), "U")) {
    rep(1, Bop$dim[1L])
  } else {
    methods::slot(B, "x")
  }
  if (length(values) != Bop$dim[1L] || any(!is.finite(values)) || any(values <= 0)) {
    stop("generalized Lanczos diagonal metric B must be positive finite.", call. = FALSE)
  }
  as.numeric(values)
}

#' @keywords internal
generalized_lanczos_diagonal_transform_operator <- function(op, Bop) {
  op <- as_operator(op)
  d <- generalized_lanczos_diagonal_values(Bop)
  inv_sqrt <- 1 / sqrt(d)
  storage <- op$metadata$storage %||% NULL
  source <- source_or_null(op)
  transformed <- if (identical(storage, "dgCMatrix")) {
    A <- methods::as(op$metadata$matrix, "dgCMatrix")
    p <- methods::slot(A, "p")
    rows <- methods::slot(A, "i") + 1L
    cols <- rep(seq_len(ncol(A)), diff(p))
    x <- methods::slot(A, "x") * inv_sqrt[rows] * inv_sqrt[cols]
    methods::slot(A, "x") <- x
    A
  } else if (identical(storage, "ddiMatrix")) {
    A <- op$metadata$matrix
    a_diag <- if (identical(methods::slot(A, "diag"), "U")) {
      rep(1, op$dim[1L])
    } else {
      methods::slot(A, "x")
    }
    methods::as(
      methods::as(Matrix::Diagonal(x = as.numeric(a_diag) / d), "generalMatrix"),
      "CsparseMatrix"
    )
  } else if (is.matrix(source) && is.double(source)) {
    C <- sweep(sweep((source + t(source)) / 2, 1L, inv_sqrt, `*`), 2L, inv_sqrt, `*`)
    (C + t(C)) / 2
  } else {
    stop("native generalized Lanczos with diagonal B requires dense, diagonal, or dgCMatrix A.",
         call. = FALSE)
  }
  Cop <- as_operator(transformed)
  Cop$structure <- hermitian()
  list(
    operator = Cop,
    back_transform = function(Y) sweep(as.matrix(Y), 1L, sqrt(d), `/`),
    kind = "diagonal_similarity",
    label = "diagonal scaling similarity transform for B",
    transformed_storage = Cop$metadata$storage %||% if (is.matrix(source_or_null(Cop))) "dense" else NA_character_
  )
}

#' @keywords internal
generalized_lanczos_dense_transform_operator <- function(op, Bop) {
  op <- as_operator(op)
  Bop <- as_operator(Bop)
  A <- source_or_null(op)
  B <- source_or_null(Bop)
  if (!(is.matrix(A) && is.double(A) && is.matrix(B) && is.double(B))) {
    stop("native generalized Lanczos dense metric transform requires dense A and dense B.",
         call. = FALSE)
  }
  A <- (A + t(A)) / 2
  B <- (B + t(B)) / 2
  factor <- chol(B)
  inv_factor <- backsolve(factor, diag(nrow(B)))
  C <- crossprod(inv_factor, A %*% inv_factor)
  C <- (C + t(C)) / 2
  Cop <- as_operator(C)
  Cop$structure <- hermitian()
  list(
    operator = Cop,
    back_transform = function(Y) backsolve(factor, as.matrix(Y)),
    kind = "dense_cholesky_similarity",
    label = "dense Cholesky similarity transform for B",
    transformed_storage = "dense"
  )
}

#' @keywords internal
native_generalized_lanczos_transform_operator <- function(op, Bop, target = smallest()) {
  op <- as_operator(op)
  Bop <- as_operator(Bop)
  if (!native_generalized_lanczos_supported(op, Bop, target = target)) {
    stop("native generalized Lanczos currently supports dense B with dense A or diagonal B with dense/diagonal/dgCMatrix A.",
         call. = FALSE)
  }
  B_source <- source_or_null(Bop)
  if (is.matrix(B_source) && is.double(B_source)) {
    generalized_lanczos_dense_transform_operator(op, Bop)
  } else {
    generalized_lanczos_diagonal_transform_operator(op, Bop)
  }
}

#' @keywords internal
native_generalized_lanczos_hermitian <- function(op, Bop, k, target = smallest(),
                                                 tol = 1e-8, maxit = NULL,
                                                 vectors = TRUE,
                                                 block = 1L,
                                                 max_restarts = 100L) {
  op <- as_operator(op)
  Bop <- as_operator(Bop)
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
  block <- as.integer(block %||% 1L)
  if (length(block) != 1L || is.na(block) || block < 1L) {
    block <- 1L
  }
  max_restarts <- as.integer(max_restarts %||% 100L)
  if (length(max_restarts) != 1L || is.na(max_restarts) || max_restarts < 0L) {
    max_restarts <- 100L
  }
  transformed <- native_generalized_lanczos_transform_operator(op, Bop, target = target)
  iter <- native_block_lanczos_hermitian(
    transformed$operator,
    k = k,
    target = target,
    tol = tol,
    maxit = maxit,
    block = block,
    max_restarts = max_restarts,
    vectors = TRUE,
    full_subspace = TRUE,
    certificate_fallback = TRUE
  )
  mapped_vectors <- transformed$back_transform(iter$vectors)
  orth <- lobpcg_b_orthonormalize(mapped_vectors, Bop)
  mapped_vectors <- orth$Q[, seq_len(min(k, ncol(orth$Q))), drop = FALSE]
  values <- iter$values[seq_len(ncol(mapped_vectors))]
  cert <- certify_eigen_operator(op, values, mapped_vectors, Bop = Bop, tol = tol)
  history <- iter$convergence_history %||% data.frame()
  list(
    values = values,
    vectors = if (isTRUE(vectors)) mapped_vectors else NULL,
    residuals = cert$residuals,
    backward_error = cert$backward_error,
    orthogonality = cert$orthogonality,
    certificate = cert,
    iterations = iter$iterations %||% 1L,
    matvecs = iter$matvecs %||% 0L,
    restarts = iter$restarts %||% 0L,
    ortho_passes = iter$ortho_passes %||% NA_integer_,
    locking_events = iter$locking_events %||% 0L,
    block = block,
    convergence_history = history,
    generalized = TRUE,
    metric_solve = list(kind = transformed$kind, label = transformed$label),
    restart = list(
      kind = "native_transformed_generalized_b_orthogonal_lanczos",
      implemented = TRUE,
      native = TRUE,
      native_kernels = TRUE,
      generalized = TRUE,
      transformed_standard_lanczos = TRUE,
      transformed_operator_storage = transformed$transformed_storage,
      max_subspace = maxit,
      max_restarts = max_restarts,
      restarts_used = iter$restarts %||% 0L,
      metric_transform = transformed$label,
      metric_solve = transformed$label,
      metric_solves = 0L,
      block = block,
      reorthogonalize = TRUE,
      orthogonalization_methods = unique(c(
        "native_transformed_standard_lanczos",
        orth$method
      )),
      native_lanczos_restart = iter$restart %||% NULL
    )
  )
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
      metric_solve_kind = metric$kind,
      metric_solve_native = isTRUE(metric$native),
      native_metric_solve = isTRUE(metric$native),
      metric_factorization = metric$factorization %||% NA_character_,
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
