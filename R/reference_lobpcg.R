#' @keywords internal
reference_lobpcg_hermitian <- function(op, k, target = smallest(), tol = 1e-8,
                                       maxit = 200L, preconditioner = NULL,
                                       seed = NULL, Bop = NULL) {
  op <- as_operator(op)
  Bop <- if (is.null(Bop)) NULL else as_operator(Bop)
  if (op$dim[1L] != op$dim[2L]) {
    stop("LOBPCG requires a square operator.", call. = FALSE)
  }
  if (!is.null(Bop) && (!identical(Bop$dim, op$dim))) {
    stop("LOBPCG metric B must have the same square dimension as A.", call. = FALSE)
  }
  if (!identical(op$structure$kind, "hermitian")) {
    stop("LOBPCG requires a Hermitian operator.", call. = FALSE)
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }

  n <- op$dim[1L]
  k <- as.integer(k)
  maxit <- as.integer(maxit)
  preconditioner_info <- eigencore_preconditioner_info(preconditioner, include_arrays = FALSE)
  X <- matrix(stats::rnorm(n * k), nrow = n, ncol = k)
  orth_native <- FALSE
  orth_all_native <- TRUE
  orth_methods <- character()
  record_orthogonalization <- function(orth) {
    orth_native <<- orth_native || isTRUE(orth$native)
    orth_all_native <<- orth_all_native && isTRUE(orth$native)
    orth_methods <<- unique(c(orth_methods, orth$method %||% "unknown"))
  }
  orth <- lobpcg_b_orthonormalize(X, Bop)
  record_orthogonalization(orth)
  X <- orth$Q
  P <- NULL
  values <- rep(NA_real_, k)
  residual_norms <- rep(Inf, k)
  iterations <- 0L
  matvecs <- 0L
  preconditioner_calls <- 0L
  history_max_relative_residual <- rep(NA_real_, maxit)
  history_nconv <- rep(NA_integer_, maxit)

  for (iter in seq_len(maxit)) {
    iterations <- iter
    AX <- apply_operator(op, X)
    BX <- if (is.null(Bop)) X else apply_operator(Bop, X)
    matvecs <- matvecs + 1L
    values <- colSums(X * AX)
    R <- AX - sweep(BX, 2L, values, `*`)
    residual_norms <- col_norms(R)
    relative_residuals <- residual_norms / pmax(abs(values), 1)
    history_max_relative_residual[iter] <- max(relative_residuals)
    history_nconv[iter] <- sum(relative_residuals <= tol)
    if (max(relative_residuals) <= tol) {
      break
    }

    W <- if (is.null(preconditioner)) {
      R
    } else {
      preconditioner_calls <- preconditioner_calls + 1L
      as.matrix(preconditioner(R))
    }
    if (!identical(dim(W), dim(R))) {
      stop("LOBPCG preconditioner must return a block with the same dimensions as its input.", call. = FALSE)
    }
    S <- if (is.null(P)) cbind(X, W) else cbind(X, W, P)
    orth <- lobpcg_b_orthonormalize(S, Bop)
    record_orthogonalization(orth)
    Q <- orth$Q
    BQ <- orth$BQ
    AQ <- apply_operator(op, Q)
    matvecs <- matvecs + 1L
    H <- crossprod(Q, AQ)
    H <- (H + t(H)) / 2
    small <- eigen(H, symmetric = TRUE)
    idx <- order_indices(small$values, target)
    if (length(idx) < k) {
      stop("LOBPCG trial subspace rank dropped below requested k.", call. = FALSE)
    }
    idx <- idx[seq_len(k)]
    X_next <- Q %*% small$vectors[, idx, drop = FALSE]
    P <- X_next - X
    orth <- lobpcg_b_orthonormalize(X_next, Bop)
    record_orthogonalization(orth)
    X <- orth$Q
    values <- small$values[idx]
  }

  cert <- certify_eigen_operator(op, values, X, Bop = Bop, tol = tol)
  history <- data.frame(
    iteration = seq_len(iterations),
    max_relative_residual = history_max_relative_residual[seq_len(iterations)],
    nconv = history_nconv[seq_len(iterations)]
  )
  list(
    values = values,
    vectors = X,
    residuals = cert$residuals,
    backward_error = cert$backward_error,
    orthogonality = cert$orthogonality,
    certificate = cert,
    iterations = iterations,
    matvecs = matvecs,
    preconditioner_calls = preconditioner_calls,
    convergence_history = history,
    preconditioned = !is.null(preconditioner),
    preconditioner = preconditioner_info,
    generalized = !is.null(Bop),
    orthogonalization = list(
      native = orth_native,
      all_native = orth_all_native,
      methods = orth_methods,
      generalized = !is.null(Bop)
    )
  )
}

#' @keywords internal
native_lobpcg_hermitian <- function(op, k, target = smallest(), tol = 1e-8,
                                    maxit = 200L, preconditioner = NULL,
                                    seed = NULL) {
  op <- as_operator(op)
  if (op$dim[1L] != op$dim[2L]) {
    stop("Native LOBPCG requires a square operator.", call. = FALSE)
  }
  if (!identical(op$structure$kind, "hermitian")) {
    stop("Native LOBPCG requires a Hermitian operator.", call. = FALSE)
  }
  if (!native_lobpcg_target_supported(target)) {
    stop("Native LOBPCG does not support target ", target_label(target), ".", call. = FALSE)
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }

  n <- op$dim[1L]
  k <- as.integer(k)
  maxit <- as.integer(maxit)
  preconditioner_info <- eigencore_preconditioner_info(preconditioner, include_arrays = FALSE)
  preconditioner_args <- native_lobpcg_preconditioner_args(preconditioner)
  start <- matrix(stats::rnorm(n * k), nrow = n, ncol = k)
  storage <- op$metadata$storage %||% NULL
  source <- source_or_null(op)
  iter <- if (identical(storage, "dgCMatrix")) {
    A <- op$metadata$matrix
    .Call(
      "eigencore_lobpcg_csc",
      methods::slot(A, "i"),
      methods::slot(A, "p"),
      methods::slot(A, "x"),
      methods::slot(A, "Dim"),
      as.integer(k),
      as.integer(maxit),
      as.integer(lanczos_target_kind(target)),
      as.numeric(tol),
      start,
      preconditioner_args$lower,
      preconditioner_args$diag,
      preconditioner_args$upper,
      PACKAGE = "eigencore"
    )
  } else if (is.matrix(source) && is.double(source)) {
    .Call(
      "eigencore_lobpcg_dense",
      source,
      as.integer(k),
      as.integer(maxit),
      as.integer(lanczos_target_kind(target)),
      as.numeric(tol),
      start,
      preconditioner_args$lower,
      preconditioner_args$diag,
      preconditioner_args$upper,
      PACKAGE = "eigencore"
    )
  } else {
    stop("Native LOBPCG requires a built-in dense or dgCMatrix operator.", call. = FALSE)
  }

  cert <- certify_eigen_operator_residuals(op, iter$values, iter$vectors,
                                           iter$residuals, tol = tol)
  history <- data.frame(
    iteration = seq_len(iter$iterations),
    max_relative_residual = iter$history_max_relative_residual[seq_len(iter$iterations)],
    nconv = iter$history_nconv[seq_len(iter$iterations)]
  )
  list(
    values = iter$values,
    vectors = iter$vectors,
    residuals = cert$residuals,
    backward_error = cert$backward_error,
    orthogonality = cert$orthogonality,
    certificate = cert,
    iterations = iter$iterations,
    matvecs = iter$matvecs,
    preconditioner_calls = iter$preconditioner_calls,
    convergence_history = history,
    preconditioned = isTRUE(preconditioner_info$supplied),
    preconditioner = preconditioner_info,
    generalized = FALSE,
    orthogonalization = list(
      native = TRUE,
      all_native = TRUE,
      methods = "native_mgs2",
      generalized = FALSE
    ),
    native = TRUE,
    q_rank_final = iter$q_rank_final
  )
}

#' @keywords internal
native_lobpcg_tridiagonal_hermitian <- function(op, k, target = smallest(),
                                                tol = 1e-8, maxit = 80L,
                                                shift = 1e-3, seed = NULL) {
  op <- as_operator(op)
  if (op$dim[1L] != op$dim[2L]) {
    stop("Native tridiagonal LOBPCG requires a square operator.", call. = FALSE)
  }
  if (!identical(op$structure$kind, "hermitian")) {
    stop("Native tridiagonal LOBPCG requires a Hermitian operator.", call. = FALSE)
  }
  if (!native_lobpcg_target_supported(target)) {
    stop("Native tridiagonal LOBPCG does not support target ", target_label(target), ".", call. = FALSE)
  }
  storage <- op$metadata$storage %||% NULL
  if (!identical(storage, "dgCMatrix")) {
    stop("Native tridiagonal LOBPCG currently requires a dgCMatrix operator.", call. = FALSE)
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }

  n <- op$dim[1L]
  k <- as.integer(k)
  maxit <- as.integer(maxit)
  shift <- as.numeric(shift)
  if (length(shift) != 1L || is.na(shift) || shift < 0) {
    stop("shift must be a single non-negative number.", call. = FALSE)
  }
  A <- op$metadata$matrix
  start <- matrix(stats::rnorm(n * k), nrow = n, ncol = k)
  iter <- .Call(
    "eigencore_lobpcg_csc_shifted_tridiagonal",
    methods::slot(A, "i"),
    methods::slot(A, "p"),
    methods::slot(A, "x"),
    methods::slot(A, "Dim"),
    as.integer(k),
    as.integer(maxit),
    as.integer(lanczos_target_kind(target)),
    as.numeric(tol),
    start,
    as.numeric(shift),
    PACKAGE = "eigencore"
  )

  cert <- certify_eigen_operator_residuals(op, iter$values, iter$vectors,
                                           iter$residuals, tol = tol)
  history <- data.frame(
    iteration = seq_len(iter$iterations),
    max_relative_residual = iter$history_max_relative_residual[seq_len(iter$iterations)],
    nconv = iter$history_nconv[seq_len(iter$iterations)]
  )
  list(
    values = iter$values,
    vectors = iter$vectors,
    residuals = cert$residuals,
    backward_error = cert$backward_error,
    orthogonality = cert$orthogonality,
    certificate = cert,
    iterations = iter$iterations,
    matvecs = iter$matvecs,
    preconditioner_calls = iter$preconditioner_calls,
    convergence_history = history,
    preconditioned = TRUE,
    preconditioner = list(
      supplied = TRUE,
      typed = TRUE,
      kind = "shifted_tridiagonal",
      native = TRUE,
      n = n,
      shift = shift,
      factorization = "tridiagonal_thomas"
    ),
    generalized = FALSE,
    orthogonalization = list(
      native = TRUE,
      all_native = TRUE,
      methods = "native_mgs2",
      generalized = FALSE
    ),
    native = TRUE,
    q_rank_final = iter$q_rank_final
  )
}

#' @keywords internal
native_generalized_lobpcg_hermitian <- function(op, Bop, k, target = smallest(),
                                                tol = 1e-8, maxit = 200L,
                                                seed = NULL) {
  op <- as_operator(op)
  Bop <- as_operator(Bop)
  if (op$dim[1L] != op$dim[2L] || !identical(Bop$dim, op$dim)) {
    stop("Native generalized LOBPCG requires square A and conformable B.", call. = FALSE)
  }
  if (!identical(op$structure$kind, "hermitian") ||
      !identical(Bop$structure$kind, "hermitian")) {
    stop("Native generalized LOBPCG requires Hermitian A and B.", call. = FALSE)
  }
  if (!native_lobpcg_target_supported(target)) {
    stop("Native generalized LOBPCG does not support target ", target_label(target), ".", call. = FALSE)
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }

  n <- op$dim[1L]
  k <- as.integer(k)
  maxit <- as.integer(maxit)
  start <- matrix(stats::rnorm(n * k), nrow = n, ncol = k)
  source <- source_or_null(op)
  Bsource <- source_or_null(Bop)
  storage <- op$metadata$storage %||% NULL
  Bstorage <- Bop$metadata$storage %||% NULL
  iter <- if (is.matrix(source) && is.double(source) &&
    is.matrix(Bsource) && is.double(Bsource)) {
    .Call(
      "eigencore_lobpcg_dense_dense_b",
      source,
      Bsource,
      as.integer(k),
      as.integer(maxit),
      as.integer(lanczos_target_kind(target)),
      as.numeric(tol),
      start,
      PACKAGE = "eigencore"
    )
  } else if (is.matrix(source) && is.double(source) &&
    identical(Bstorage, "ddiMatrix")) {
    B <- Bop$metadata$matrix
    unit <- identical(methods::slot(B, "diag"), "U")
    diagonal <- if (unit) numeric(0) else methods::slot(B, "x")
    .Call(
      "eigencore_lobpcg_dense_diagonal_b",
      source,
      as.numeric(diagonal),
      isTRUE(unit),
      as.integer(k),
      as.integer(maxit),
      as.integer(lanczos_target_kind(target)),
      as.numeric(tol),
      start,
      PACKAGE = "eigencore"
    )
  } else if (is.matrix(source) && is.double(source) &&
    identical(Bstorage, "dgCMatrix")) {
    B <- Bop$metadata$matrix
    .Call(
      "eigencore_lobpcg_dense_csc_b",
      source,
      methods::slot(B, "i"),
      methods::slot(B, "p"),
      methods::slot(B, "x"),
      methods::slot(B, "Dim"),
      as.integer(k),
      as.integer(maxit),
      as.integer(lanczos_target_kind(target)),
      as.numeric(tol),
      start,
      PACKAGE = "eigencore"
    )
  } else if (identical(storage, "dgCMatrix") &&
    identical(Bstorage, "ddiMatrix")) {
    A <- op$metadata$matrix
    B <- Bop$metadata$matrix
    unit <- identical(methods::slot(B, "diag"), "U")
    diagonal <- if (unit) numeric(0) else methods::slot(B, "x")
    .Call(
      "eigencore_lobpcg_csc_diagonal_b",
      methods::slot(A, "i"),
      methods::slot(A, "p"),
      methods::slot(A, "x"),
      methods::slot(A, "Dim"),
      as.numeric(diagonal),
      isTRUE(unit),
      as.integer(k),
      as.integer(maxit),
      as.integer(lanczos_target_kind(target)),
      as.numeric(tol),
      start,
      PACKAGE = "eigencore"
    )
  } else if (identical(storage, "dgCMatrix") &&
    identical(Bstorage, "dgCMatrix")) {
    A <- op$metadata$matrix
    B <- Bop$metadata$matrix
    .Call(
      "eigencore_lobpcg_csc_csc_b",
      methods::slot(A, "i"),
      methods::slot(A, "p"),
      methods::slot(A, "x"),
      methods::slot(A, "Dim"),
      methods::slot(B, "i"),
      methods::slot(B, "p"),
      methods::slot(B, "x"),
      methods::slot(B, "Dim"),
      as.integer(k),
      as.integer(maxit),
      as.integer(lanczos_target_kind(target)),
      as.numeric(tol),
      start,
      PACKAGE = "eigencore"
    )
  } else if (identical(storage, "ddiMatrix") &&
    identical(Bstorage, "ddiMatrix")) {
    A <- op$metadata$matrix
    B <- Bop$metadata$matrix
    a_unit <- identical(methods::slot(A, "diag"), "U")
    b_unit <- identical(methods::slot(B, "diag"), "U")
    a_diagonal <- if (a_unit) numeric(0) else methods::slot(A, "x")
    b_diagonal <- if (b_unit) numeric(0) else methods::slot(B, "x")
    .Call(
      "eigencore_lobpcg_diagonal_diagonal_b",
      as.numeric(a_diagonal),
      isTRUE(a_unit),
      methods::slot(A, "Dim"),
      as.numeric(b_diagonal),
      isTRUE(b_unit),
      as.integer(k),
      as.integer(maxit),
      as.integer(lanczos_target_kind(target)),
      as.numeric(tol),
      start,
      PACKAGE = "eigencore"
    )
  } else {
    stop("Native generalized LOBPCG slice currently supports dense/CSC/diagonal A with dense, diagonal, or CSC B.", call. = FALSE)
  }

  cert <- certify_eigen_operator_residuals(op, iter$values, iter$vectors,
                                           iter$residuals, Bop = Bop, tol = tol)
  history <- data.frame(
    iteration = seq_len(iter$iterations),
    max_relative_residual = iter$history_max_relative_residual[seq_len(iter$iterations)],
    nconv = iter$history_nconv[seq_len(iter$iterations)]
  )
  list(
    values = iter$values,
    vectors = iter$vectors,
    residuals = cert$residuals,
    backward_error = cert$backward_error,
    orthogonality = cert$orthogonality,
    certificate = cert,
    iterations = iter$iterations,
    matvecs = iter$matvecs,
    preconditioner_calls = iter$preconditioner_calls,
    convergence_history = history,
    preconditioned = FALSE,
    preconditioner = eigencore_preconditioner_info(NULL, include_arrays = FALSE),
    generalized = TRUE,
    orthogonalization = list(
      native = TRUE,
      all_native = TRUE,
      methods = if (identical(Bop$metadata$storage %||% NULL, "ddiMatrix")) {
        "native_diagonal_b_mgs2"
      } else if (identical(Bop$metadata$storage %||% NULL, "dgCMatrix")) {
        "native_csc_b_mgs2"
      } else {
        "native_dense_b_mgs2"
      },
      generalized = TRUE
    ),
    native = TRUE,
    q_rank_final = iter$q_rank_final
  )
}

#' @keywords internal
native_lobpcg_supported <- function(op, target = smallest(), preconditioner = NULL,
                                    Bop = NULL) {
  op <- as_operator(op)
  if (!is.null(Bop) || !identical(op$structure$kind, "hermitian") ||
      !native_lobpcg_target_supported(target)) {
    return(FALSE)
  }
  storage <- op$metadata$storage %||% NULL
  source <- source_or_null(op)
  built_in <- identical(storage, "dgCMatrix") || (is.matrix(source) && is.double(source))
  built_in && native_lobpcg_preconditioner_supported(preconditioner)
}

#' @keywords internal
native_generalized_lobpcg_supported <- function(op, Bop, target = smallest(),
                                                preconditioner = NULL) {
  op <- as_operator(op)
  Bop <- as_operator(Bop)
  if (!is.null(preconditioner) ||
      !identical(op$structure$kind, "hermitian") ||
      !identical(Bop$structure$kind, "hermitian") ||
      !native_lobpcg_target_supported(target)) {
    return(FALSE)
  }
  source <- source_or_null(op)
  Bsource <- source_or_null(Bop)
  storage <- op$metadata$storage %||% NULL
  Bstorage <- Bop$metadata$storage %||% NULL
  dense_A <- is.matrix(source) && is.double(source)
  dense_B <- is.matrix(Bsource) && is.double(Bsource)
  csc_A <- identical(storage, "dgCMatrix")
  diagonal_A <- identical(storage, "ddiMatrix")
  diagonal_B <- identical(Bstorage, "ddiMatrix")
  csc_B <- identical(Bstorage, "dgCMatrix")
  (dense_A && (dense_B || diagonal_B || csc_B)) ||
    (csc_A && (diagonal_B || csc_B)) ||
    (diagonal_A && diagonal_B)
}

#' @keywords internal
native_lobpcg_target_supported <- function(target) {
  kind <- if (inherits(target, "eigencore_target")) target$kind else "largest"
  kind %in% c("largest", "smallest", "largest_magnitude", "smallest_magnitude")
}

#' @keywords internal
native_lobpcg_preconditioner_supported <- function(preconditioner) {
  if (is.null(preconditioner)) {
    return(TRUE)
  }
  info <- eigencore_preconditioner_info(preconditioner, include_arrays = FALSE)
  isTRUE(info$native) && identical(info$kind, "shifted_tridiagonal")
}

#' @keywords internal
native_lobpcg_preconditioner_args <- function(preconditioner) {
  if (is.null(preconditioner)) {
    return(list(lower = numeric(0), diag = numeric(0), upper = numeric(0)))
  }
  info <- attr(preconditioner, "eigencore_preconditioner", exact = TRUE)
  if (!native_lobpcg_preconditioner_supported(preconditioner)) {
    stop("Native LOBPCG only supports NULL or shifted_tridiagonal preconditioners.", call. = FALSE)
  }
  list(
    lower = info$lower,
    diag = info$diag,
    upper = info$upper
  )
}

#' @keywords internal
lobpcg_b_orthonormalize <- function(X, Bop = NULL, tol = sqrt(.Machine$double.eps)) {
  X <- as.matrix(X)
  if (is.null(Bop)) {
    Q <- tryCatch(native_cholqr2(X)$Q, error = function(e) NULL)
    if (is.null(Q)) {
      qrX <- qr(X, tol = tol)
      if (qrX$rank < 1L) {
        stop("LOBPCG trial subspace is numerically degenerate.", call. = FALSE)
      }
      Q <- qr.Q(qrX)[, seq_len(qrX$rank), drop = FALSE]
      return(list(Q = Q, BQ = Q, native = FALSE, method = "qr"))
    }
    return(list(Q = Q, BQ = Q, native = TRUE, method = "cholqr2"))
  }

  Bsrc <- source_or_null(Bop)
  if (is.matrix(Bsrc) && is.double(Bsrc)) {
    native <- tryCatch(native_b_cholqr2(X, Bsrc), error = function(e) NULL)
    if (!is.null(native)) {
      BQ <- Bsrc %*% native$Q
      return(list(Q = native$Q, BQ = BQ, native = TRUE, method = "b_cholqr2"))
    }
  }
  if (identical(Bop$metadata$storage %||% NULL, "ddiMatrix")) {
    B <- Bop$metadata$matrix
    unit <- identical(methods::slot(B, "diag"), "U")
    diagonal <- if (unit) numeric(0) else methods::slot(B, "x")
    native <- tryCatch(native_diagonal_b_cholqr2(X, diagonal, unit = unit),
                       error = function(e) NULL)
    if (!is.null(native)) {
      return(list(
        Q = native$Q,
        BQ = native$BQ,
        native = TRUE,
        method = "diagonal_b_cholqr2"
      ))
    }
  }

  qrX <- qr(X, tol = tol)
  if (qrX$rank < 1L) {
    stop("LOBPCG trial subspace is numerically B-degenerate.", call. = FALSE)
  }
  X <- qr.Q(qrX)[, seq_len(qrX$rank), drop = FALSE]
  BQ <- apply_operator(Bop, X)
  gram <- crossprod(X, BQ)
  gram <- (gram + t(gram)) / 2
  chol_try <- chol(gram)
  invR <- backsolve(chol_try, diag(ncol(chol_try)))
  Q <- X %*% invR
  BQ <- BQ %*% invR
  list(Q = Q, BQ = BQ, native = FALSE, method = "operator_b_qr_chol")
}
