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
  preconditioner_info <- eigencore_preconditioner_info(preconditioner)
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
