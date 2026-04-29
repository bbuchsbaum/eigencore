#' @keywords internal
shift_invert_plan_label <- function(problem, has_metric, is_hermitian,
                                    is_dense_source, is_native_csc) {
  user_solve <- problem$transform$solve
  if (!is_hermitian) {
    return("shift-invert requested (only Hermitian shift-invert is implemented)")
  }
  if (has_metric) {
    return("shift-invert requested (generalized SPD shift-invert is on the V1 roadmap)")
  }
  if (!is.null(user_solve)) {
    return("reference Hermitian Lanczos shift-invert (user solve)")
  }
  csc_available <- inherits(problem$A$metadata$matrix, "CsparseMatrix")
  if (is_native_csc || csc_available) {
    return("reference Hermitian Lanczos shift-invert (sparse LU)")
  }
  if (is_dense_source) {
    return("reference Hermitian Lanczos shift-invert (dense LU)")
  }
  "shift-invert requested (provide method$solve for matrix-free A)"
}

#' @keywords internal
shift_invert_apply_factory <- function(solve_fn) {
  function(X, alpha = 1, beta = 0, Y = NULL) {
    Z <- solve_fn(X)
    if (is.null(Y) || beta == 0) {
      if (alpha == 1) Z else alpha * Z
    } else {
      alpha * Z + beta * Y
    }
  }
}

#' @keywords internal
shift_invert_operator_fingerprint <- function(op) {
  op <- as_operator(op)
  source <- source_or_null(op)
  storage <- op$metadata$storage %||% NULL
  matrix <- op$metadata$matrix %||% NULL
  if (is.matrix(source)) {
    return(list(
      kind = "dense",
      dim = dim(source),
      storage_mode = storage.mode(source),
      values = as.vector(source)
    ))
  }
  if (inherits(matrix, "sparseMatrix")) {
    if (!inherits(matrix, "CsparseMatrix")) {
      matrix <- methods::as(matrix, "CsparseMatrix")
    }
    return(list(
      kind = storage %||% class(matrix)[[1L]],
      class = class(matrix),
      dim = methods::slot(matrix, "Dim"),
      i = methods::slot(matrix, "i"),
      p = methods::slot(matrix, "p"),
      x = if ("x" %in% methods::slotNames(matrix)) {
        methods::slot(matrix, "x")
      } else {
        numeric(0)
      }
    ))
  }
  list(
    kind = "operator",
    dim = op$dim,
    name = op$name %||% NA_character_,
    storage = storage,
    source_available = !is.null(source)
  )
}

#' @keywords internal
shift_invert_factorization_cache_key <- function(Aop, sigma, Bop = NULL) {
  sigma <- as.numeric(sigma)
  if (length(sigma) != 1L || !is.finite(sigma)) {
    stop("shift-invert cache key requires a single finite sigma.", call. = FALSE)
  }
  key <- list(
    transform = "shift_invert",
    sigma = sigma,
    A = shift_invert_operator_fingerprint(Aop),
    B = if (is.null(Bop)) NULL else shift_invert_operator_fingerprint(Bop)
  )
  class(key) <- "eigencore_shift_invert_cache_key"
  key
}

#' @keywords internal
shift_invert_factorization_cache_info <- function(Aop, sigma, Bop = NULL,
                                                  label_kind = NA_character_) {
  list(
    key = shift_invert_factorization_cache_key(Aop, sigma, Bop = Bop),
    label_kind = label_kind,
    native = FALSE,
    reusable_within_operator = TRUE,
    external_cache = FALSE
  )
}

#' @keywords internal
shift_invert_solver_dense <- function(A, sigma) {
  M <- A - sigma * diag(nrow(A))
  cond <- tryCatch(base::rcond(M), error = function(e) NA_real_)
  if (!is.finite(cond) || cond <= sqrt(.Machine$double.eps)) {
    stop(
      "shift_invert(sigma = ", sigma, ") produced a singular or near-singular ",
      "dense shifted operator; perturb sigma or supply a stable solve function.",
      call. = FALSE
    )
  }
  # Use base::solve at apply time so we don't materialize M^{-1}
  list(
    solve_fn = function(X) base::solve(M, X),
    label    = "dense_lu",
    M        = M
  )
}

#' @keywords internal
shift_invert_solver_csc <- function(A, sigma) {
  n <- nrow(A)
  M <- methods::as(A - sigma * Matrix::Diagonal(n), "CsparseMatrix")
  list(
    solve_fn = function(X) {
      Z <- Matrix::solve(M, X)
      if (inherits(Z, "Matrix")) as.matrix(Z) else Z
    },
    label = "sparse_lu",
    M     = M
  )
}

#' @keywords internal
prepare_shift_invert_operator <- function(problem, sigma, user_solve = NULL) {
  Aop <- problem$A
  Bop <- problem$metric
  n <- Aop$dim[1L]
  cache_info <- shift_invert_factorization_cache_info(Aop, sigma, Bop = Bop)

  if (!is.null(Bop)) {
    stop(
      "shift_invert() currently supports standard Hermitian eigenproblems (B = NULL). ",
      "Generalized SPD shift-invert is on the V1 roadmap (plan_v1.md milestone L).",
      call. = FALSE
    )
  }

  if (!identical(problem$structure$kind, "hermitian")) {
    stop("shift_invert() currently requires a Hermitian eigenproblem.", call. = FALSE)
  }

  # User-supplied solve operator for matrix-free A
  if (!is.null(user_solve)) {
    if (!is.function(user_solve)) {
      stop("shift_invert(solve = ...) must be a function.", call. = FALSE)
    }
    op <- linear_operator(
      dim = c(n, n),
      apply = shift_invert_apply_factory(user_solve),
      apply_adjoint = NULL,
      structure = hermitian(),
      name = "shift_invert_user_solve",
      metadata = list(
        native = FALSE,
        factorization_cache = modifyList(cache_info, list(label_kind = "user_solve"))
      )
    )
    return(list(
      operator = op,
      label_kind = "user_solve",
      factorization_cache = modifyList(cache_info, list(label_kind = "user_solve"))
    ))
  }

  # Build factorization from source matrix; CSC operators carry the matrix in
  # metadata$matrix rather than metadata$source.
  source_A <- source_or_null(Aop)
  csc_A    <- if (identical(Aop$metadata$storage %||% NULL, "dgCMatrix")) {
    Aop$metadata$matrix
  } else if (inherits(Aop$metadata$matrix, "CsparseMatrix")) {
    Aop$metadata$matrix
  } else {
    NULL
  }

  prep <- if (is.matrix(source_A) && is.double(source_A)) {
    shift_invert_solver_dense(source_A, sigma)
  } else if (!is.null(csc_A)) {
    shift_invert_solver_csc(methods::as(csc_A, "generalMatrix"), sigma)
  } else if (inherits(source_A, "CsparseMatrix") || inherits(source_A, "dgCMatrix")) {
    shift_invert_solver_csc(source_A, sigma)
  } else {
    stop(
      "shift_invert() supports dense double matrices and dgCMatrix/dsCMatrix sources, ",
      "or a user-supplied solve operator.",
      call. = FALSE
    )
  }

  op <- linear_operator(
    dim = c(n, n),
    apply = shift_invert_apply_factory(prep$solve_fn),
    apply_adjoint = NULL,
    structure = hermitian(),
    name = paste0("shift_invert_", prep$label),
    metadata = list(
      native = FALSE,
      factorization_cache = modifyList(cache_info, list(label_kind = prep$label))
    )
  )
  list(
    operator = op,
    label_kind = prep$label,
    factorization_cache = modifyList(cache_info, list(label_kind = prep$label))
  )
}

#' @keywords internal
solve_shift_invert_hermitian <- function(problem, k, method, tol, maxit,
                                          vectors, certify, plan) {
  sigma <- method$sigma
  if (length(sigma) != 1L || !is.finite(sigma)) {
    stop("shift_invert(sigma) requires a single finite numeric shift.", call. = FALSE)
  }

  prep <- prepare_shift_invert_operator(problem, sigma, user_solve = method$solve)
  M <- prep$operator
  n <- M$dim[1L]

  effective_maxit <- maxit %||% min(n, max(20L, 4L * as.integer(k) + 20L))

  iter <- reference_lanczos_hermitian(
    M,
    k = k,
    target = largest_magnitude(),
    tol = tol,
    maxit = effective_maxit,
    vectors = TRUE,
    reorthogonalize = TRUE
  )

  mu <- iter$values
  vec <- iter$vectors
  if (any(abs(mu) < .Machine$double.eps)) {
    stop(
      "shift_invert(sigma = ", sigma, ") produced a zero-magnitude eigenvalue ",
      "of the inverted operator; sigma is too close to a true eigenvalue. ",
      "Perturb sigma or use a tighter tolerance.",
      call. = FALSE
    )
  }
  lambda <- sigma + 1 / mu

  Aop <- problem$A
  cert_pre <- certify_eigen_operator(Aop, lambda, vec, tol = tol)

  ord <- order_indices(lambda, problem$target)
  if (length(ord) > k) ord <- ord[seq_len(k)]
  lambda <- lambda[ord]
  vec <- vec[, ord, drop = FALSE]

  cert <- if (isTRUE(certify) && !is.null(vec) && ncol(vec) > 0L) {
    certify_eigen_operator(Aop, lambda, vec, tol = tol)
  } else {
    empty_certificate(tol)
  }

  result <- list(
    values = lambda,
    vectors = if (isTRUE(vectors)) vec else NULL,
    residuals = cert$residuals,
    backward_error = cert$backward_error,
    orthogonality = cert$orthogonality,
    nconv = sum(cert$converged),
    requested = k,
    iterations = iter$iterations,
    matvecs = iter$matvecs,
    method = plan$method,
    target = target_label(problem$target),
    plan = plan,
    certificate = cert,
    sigma = sigma,
    transform = list(
      kind = "shift_invert",
      sigma = sigma,
      label_kind = prep$label_kind,
      factorization_cache = prep$factorization_cache
    ),
    warnings = paste0(
      "using reference Hermitian Lanczos shift-invert (",
      prep$label_kind, "); native shift-invert hot loop not yet implemented"
    )
  )
  class(result) <- "eigencore_eigen_result"
  result
}
