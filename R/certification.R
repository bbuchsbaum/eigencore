#' Extract a result certificate.
certificate <- function(x, ...) {
  x$certificate
}

#' Extract diagnostics.
diagnostics <- function(x, ...) {
  list(
    residuals = x$residuals,
    backward_error = x$backward_error,
    orthogonality = x$orthogonality,
    nconv = x$nconv,
    iterations = x$iterations,
    matvecs = x$matvecs,
    preconditioner_calls = x$preconditioner_calls,
    convergence_history = x$convergence_history,
    restart = x$restart,
    stage_seconds = x$stage_seconds %||% x$restart$stage_seconds %||% numeric(),
    preconditioner = x$preconditioner %||% x$restart$preconditioner %||% NULL,
    locked = x$locked,
    method = x$method,
    plan = x$plan,
    warnings = x$warnings
  )
}

#' Extract computed values.
values <- function(x, ...) {
  x$values
}

#' Extract eigenvectors.
vectors <- function(x, ...) {
  x$vectors
}

#' Extract left singular vectors.
left_vectors <- function(x, ...) {
  x$u
}

#' Extract right singular vectors.
right_vectors <- function(x, ...) {
  x$v
}

#' Extract residual diagnostics.
residuals <- function(x, ...) {
  x$residuals
}

#' Extract backward-error diagnostics.
backward_error <- function(x, ...) {
  x$backward_error
}

#' @keywords internal
certify_eigen <- function(A, values, vectors, B = NULL, tol = 1e-8,
                          require_orthogonality = TRUE) {
  diag <- native_dense_eigen_certificate(A, values, vectors, B = B, tol = tol)
  new_certificate(
    tol = tol,
    residuals = diag$residuals,
    backward_error = diag$backward_error,
    orthogonality = diag$orthogonality,
    converged = diag$converged,
    scale = diag$scale,
    norm_bound_type = "frobenius_exact",
    require_orthogonality = require_orthogonality
  )
}

#' @keywords internal
certify_eigen_operator <- function(Aop, values, vectors, Bop = NULL, tol = 1e-8) {
  native <- native_builtin_eigen_certificate(Aop, values, vectors, Bop = Bop, tol = tol)
  if (!is.null(native)) {
    diag <- native$diagnostics
    return(new_certificate(
      tol = tol,
      residuals = diag$residuals,
      backward_error = diag$backward_error,
      orthogonality = diag$orthogonality,
      converged = diag$converged,
      scale = diag$scale,
      norm_bound_type = native$norm_bound_type
    ))
  }

  k <- length(values)
  Av <- apply_operator(Aop, vectors)
  Bv <- if (is.null(Bop)) vectors else apply_operator(Bop, vectors)
  residual_matrix <- Av - sweep(Bv, 2L, values, `*`)
  residuals <- col_norms(residual_matrix)
  norm_A <- operator_norm_for_certificate_info(Aop)
  norm_B <- if (is.null(Bop)) {
    list(value = 1, norm_bound_type = "identity_exact", scale_is_estimate = FALSE)
  } else {
    operator_norm_for_certificate_info(Bop)
  }
  scale <- eigen_backward_scale(
    norm_A$value,
    norm_B$value,
    values,
    vectors
  )
  backward <- residuals / pmax(scale, .Machine$double.eps)
  gram <- if (is.null(Bop)) crossprod(vectors) else crossprod(vectors, Bv)
  orth <- max(abs(gram - diag(k)))
  new_certificate(
    tol = tol,
    residuals = residuals,
    backward_error = backward,
    orthogonality = orth,
    converged = backward <= tol,
    scale = scale,
    norm_bound_type = paste(c(norm_A$norm_bound_type, norm_B$norm_bound_type), collapse = "+"),
    scale_is_estimate = isTRUE(norm_A$scale_is_estimate) || isTRUE(norm_B$scale_is_estimate)
  )
}

#' @keywords internal
certify_eigen_operator_residuals <- function(Aop, values, vectors, residuals,
                                             Bop = NULL, tol = 1e-8) {
  norm_A <- operator_norm_for_certificate_info(Aop)
  norm_B <- if (is.null(Bop)) {
    list(value = 1, norm_bound_type = "identity_exact", scale_is_estimate = FALSE)
  } else {
    operator_norm_for_certificate_info(Bop)
  }
  scale <- eigen_backward_scale(norm_A$value, norm_B$value, values, vectors)
  backward <- residuals / pmax(scale, .Machine$double.eps)
  if (is.null(Bop)) {
    orth <- orthogonality_loss(vectors)
  } else {
    Bv <- apply_operator(Bop, vectors)
    orth <- max(abs(crossprod(vectors, Bv) - diag(length(values))))
  }
  new_certificate(
    tol = tol,
    residuals = residuals,
    backward_error = backward,
    orthogonality = orth,
    converged = backward <= tol,
    scale = scale,
    norm_bound_type = paste(c(norm_A$norm_bound_type, norm_B$norm_bound_type), collapse = "+"),
    scale_is_estimate = isTRUE(norm_A$scale_is_estimate) || isTRUE(norm_B$scale_is_estimate)
  )
}

#' @keywords internal
certify_svd <- function(A, d, u, v, tol = 1e-8) {
  diag <- native_dense_svd_certificate(A, d, u, v, tol = tol)
  new_certificate(
    tol = tol,
    residuals = list(left = diag$left, right = diag$right, combined = diag$combined),
    backward_error = diag$backward_error,
    orthogonality = diag$orthogonality,
    converged = diag$converged,
    scale = diag$scale,
    norm_bound_type = "frobenius_exact"
  )
}

#' @keywords internal
certify_svd_operator <- function(Aop, d, u, v, tol = 1e-8) {
  native <- native_builtin_svd_certificate(Aop, d, u, v, tol = tol)
  if (!is.null(native)) {
    diag <- native$diagnostics
    return(new_certificate(
      tol = tol,
      residuals = list(left = diag$left, right = diag$right, combined = diag$combined),
      backward_error = diag$backward_error,
      orthogonality = diag$orthogonality,
      converged = diag$converged,
      scale = diag$scale,
      norm_bound_type = native$norm_bound_type
    ))
  }

  left_residual_matrix <- apply_operator(Aop, v) - sweep(u, 2L, d, `*`)
  right_residual_matrix <- apply_adjoint_operator(Aop, u) - sweep(v, 2L, d, `*`)
  left <- col_norms(left_residual_matrix)
  right <- col_norms(right_residual_matrix)
  combined <- sqrt(left^2 + right^2)
  norm_A <- operator_norm_for_certificate_info(Aop)
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
certify_svd_operator_cached_av <- function(Aop, d, u, v, Av, tol = 1e-8,
                                           return_residual_vectors = FALSE) {
  if (is.null(Av)) {
    cert <- certify_svd_operator(Aop, d, u, v, tol = tol)
    if (isTRUE(return_residual_vectors)) {
      return(list(certificate = cert, right_residual_vectors = NULL))
    }
    return(cert)
  }
  if (!isTRUE(return_residual_vectors)) {
    native <- native_builtin_svd_certificate_cached_av(Aop, d, u, v, Av, tol = tol)
    if (!is.null(native)) {
      diag <- native$diagnostics
      return(new_certificate(
        tol = tol,
        residuals = list(left = diag$left, right = diag$right, combined = diag$combined),
        backward_error = diag$backward_error,
        orthogonality = diag$orthogonality,
        converged = diag$converged,
        scale = diag$scale,
        norm_bound_type = native$norm_bound_type
      ))
    }
    return(certify_svd_operator(Aop, d, u, v, tol = tol))
  }
  Av <- as.matrix(Av)
  if (nrow(Av) != Aop$dim[[1L]] || ncol(Av) != length(d)) {
    stop("Av must have nrow equal to nrow(Aop) and one column per singular value.",
         call. = FALSE)
  }
  left_residual_matrix <- Av - sweep(u, 2L, d, `*`)
  right_residual_matrix <- apply_adjoint_operator(Aop, u) - sweep(v, 2L, d, `*`)
  left <- col_norms(left_residual_matrix)
  right <- col_norms(right_residual_matrix)
  combined <- sqrt(left^2 + right^2)
  norm_A <- operator_norm_for_certificate_info(Aop)
  scale <- svd_backward_scale(norm_A$value, d)
  backward <- combined / scale
  orth_u <- max(abs(crossprod(u) - diag(length(d))))
  orth_v <- max(abs(crossprod(v) - diag(length(d))))
  cert <- new_certificate(
    tol = tol,
    residuals = list(left = left, right = right, combined = combined),
    backward_error = backward,
    orthogonality = c(U = orth_u, V = orth_v),
    converged = backward <= tol,
    scale = scale,
    norm_bound_type = norm_A$norm_bound_type,
    scale_is_estimate = isTRUE(norm_A$scale_is_estimate)
  )
  list(certificate = cert, right_residual_vectors = right_residual_matrix)
}

#' @keywords internal
certify_svd_operator_cached_sides <- function(Aop, d, u, v, Av, Atu,
                                             tol = 1e-8) {
  Av <- as.matrix(Av)
  Atu <- as.matrix(Atu)
  if (nrow(Av) != Aop$dim[[1L]] || ncol(Av) != length(d)) {
    stop("Av must have nrow equal to nrow(Aop) and one column per singular value.",
         call. = FALSE)
  }
  if (nrow(Atu) != Aop$dim[[2L]] || ncol(Atu) != length(d)) {
    stop("Atu must have nrow equal to ncol(Aop) and one column per singular value.",
         call. = FALSE)
  }
  left_residual_matrix <- Av - sweep(u, 2L, d, `*`)
  right_residual_matrix <- Atu - sweep(v, 2L, d, `*`)
  left <- col_norms(left_residual_matrix)
  right <- col_norms(right_residual_matrix)
  combined <- sqrt(left^2 + right^2)
  norm_A <- operator_norm_for_certificate_info(Aop)
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
new_certificate <- function(tol, residuals, backward_error, orthogonality,
                            converged, scale, notes = character(),
                            certificate_type = "residual_backward_error",
                            norm_bound_type = "unspecified",
                            scale_is_estimate = FALSE,
                            require_orthogonality = TRUE) {
  if (isTRUE(scale_is_estimate)) {
    notes <- c(notes, "certificate scale uses a stochastic norm estimate; passed is withheld")
  }
  orthogonality_tolerance <- max(tol, sqrt(.Machine$double.eps))
  max_orthogonality_loss <- if (length(orthogonality)) max(orthogonality) else NA_real_
  orthogonality_passed <- !isTRUE(require_orthogonality) ||
    is.na(max_orthogonality_loss) ||
    max_orthogonality_loss <= orthogonality_tolerance
  if (!isTRUE(orthogonality_passed)) {
    notes <- c(notes, "orthogonality loss exceeds certificate tolerance")
  }
  cert <- list(
    passed = all(converged) && orthogonality_passed && !isTRUE(scale_is_estimate),
    tolerance = tol,
    orthogonality_tolerance = orthogonality_tolerance,
    orthogonality_required = isTRUE(require_orthogonality),
    certificate_type = certificate_type,
    norm_bound_type = norm_bound_type,
    scale_is_estimate = isTRUE(scale_is_estimate),
    max_backward_error = if (length(backward_error)) max(backward_error) else NA_real_,
    max_residual = max_residual_value(residuals),
    max_orthogonality_loss = max_orthogonality_loss,
    orthogonality_passed = orthogonality_passed,
    failed_indices = which(!converged),
    scale = scale,
    notes = notes,
    residuals = residuals,
    backward_error = backward_error,
    orthogonality = orthogonality,
    converged = converged
  )
  class(cert) <- "eigencore_certificate"
  cert
}

#' @keywords internal
empty_certificate <- function(tol, note) {
  new_certificate(
    tol = tol,
    residuals = numeric(),
    backward_error = numeric(),
    orthogonality = numeric(),
    converged = FALSE,
    scale = NA_real_,
    notes = note,
    certificate_type = "uncomputed",
    norm_bound_type = "none"
  )
}

#' @export
print.eigencore_certificate <- function(x, ...) {
  cat("eigencore certificate\n")
  cat("  passed:", x$passed, "\n")
  cat("  tolerance:", format(x$tolerance), "\n")
  cat("  type:", x$certificate_type, "\n")
  cat("  norm bound:", x$norm_bound_type, "\n")
  cat("  scale estimated:", x$scale_is_estimate, "\n")
  cat("  max residual:", format(x$max_residual), "\n")
  cat("  max backward error:", format(x$max_backward_error), "\n")
  cat("  max orthogonality loss:", format(x$max_orthogonality_loss), "\n")
  cat("  orthogonality tolerance:", format(x$orthogonality_tolerance), "\n")
  cat("  orthogonality required:", x$orthogonality_required, "\n")
  if (length(x$failed_indices)) {
    cat("  failed indices:", paste(x$failed_indices, collapse = ", "), "\n")
  }
  if (length(x$notes)) {
    cat("  notes:", paste(x$notes, collapse = "; "), "\n")
  }
  invisible(x)
}

#' @keywords internal
col_norms <- function(x) {
  .Call("eigencore_col_norms", as.matrix(x), PACKAGE = "eigencore")
}

#' @keywords internal
dense_eigen_residuals <- function(A, values, vectors, B = NULL) {
  .Call(
    "eigencore_dense_eigen_residuals",
    as.matrix(A),
    as.numeric(values),
    as.matrix(vectors),
    if (is.null(B)) NULL else as.matrix(B),
    PACKAGE = "eigencore"
  )
}

#' @keywords internal
native_dense_eigen_certificate <- function(A, values, vectors, B = NULL, tol = 1e-8) {
  .Call(
    "eigencore_dense_eigen_certificate",
    as.matrix(A),
    as.numeric(values),
    as.matrix(vectors),
    if (is.null(B)) NULL else as.matrix(B),
    as.numeric(tol),
    PACKAGE = "eigencore"
  )
}

#' @keywords internal
dense_svd_residuals <- function(A, d, u, v) {
  .Call(
    "eigencore_dense_svd_residuals",
    as.matrix(A),
    as.numeric(d),
    as.matrix(u),
    as.matrix(v),
    PACKAGE = "eigencore"
  )
}

#' @keywords internal
native_dense_svd_certificate <- function(A, d, u, v, tol = 1e-8) {
  .Call(
    "eigencore_dense_svd_certificate",
    as.matrix(A),
    as.numeric(d),
    as.matrix(u),
    as.matrix(v),
    as.numeric(tol),
    PACKAGE = "eigencore"
  )
}

#' @keywords internal
native_dense_svd_certificate_cached_av <- function(A, d, u, v, Av, tol = 1e-8) {
  .Call(
    "eigencore_dense_svd_certificate_cached_av",
    as.matrix(A),
    as.numeric(d),
    as.matrix(u),
    as.matrix(v),
    as.matrix(Av),
    as.numeric(tol),
    PACKAGE = "eigencore"
  )
}

#' @keywords internal
native_builtin_eigen_certificate <- function(Aop, values, vectors, Bop = NULL, tol = 1e-8) {
  storage <- Aop$metadata$storage %||% NULL
  source <- source_or_null(Aop)
  if (!is.null(Bop)) {
    B_source <- source_or_null(Bop)
    if (is.matrix(source) && is.double(source) &&
        is.matrix(B_source) && is.double(B_source)) {
      return(list(
        diagnostics = native_dense_eigen_certificate(source, values, vectors,
                                                    B = B_source, tol = tol),
        norm_bound_type = "frobenius_exact+frobenius_exact"
      ))
    }
    return(NULL)
  }
  if (is.matrix(source) && is.double(source)) {
    return(list(
      diagnostics = native_dense_eigen_certificate(source, values, vectors, tol = tol),
      norm_bound_type = "frobenius_exact+identity_exact"
    ))
  }
  if (identical(storage, "dgCMatrix")) {
    A <- Aop$metadata$matrix
    return(list(
      diagnostics = .Call(
        "eigencore_csc_eigen_certificate",
        methods::slot(A, "i"),
        methods::slot(A, "p"),
        methods::slot(A, "x"),
        methods::slot(A, "Dim"),
        as.numeric(values),
        as.matrix(vectors),
        as.numeric(Aop$metadata$frobenius_norm),
        as.numeric(tol),
        PACKAGE = "eigencore"
      ),
      norm_bound_type = "frobenius_metadata+identity_exact"
    ))
  }
  if (identical(storage, "ddiMatrix")) {
    A <- Aop$metadata$matrix
    return(list(
      diagnostics = .Call(
        "eigencore_diagonal_eigen_certificate",
        methods::slot(A, "x"),
        methods::slot(A, "Dim"),
        identical(methods::slot(A, "diag"), "U"),
        as.numeric(values),
        as.matrix(vectors),
        as.numeric(Aop$metadata$frobenius_norm),
        as.numeric(tol),
        PACKAGE = "eigencore"
      ),
      norm_bound_type = "frobenius_metadata+identity_exact"
    ))
  }
  NULL
}

#' @keywords internal
native_builtin_svd_certificate <- function(Aop, d, u, v, tol = 1e-8) {
  storage <- Aop$metadata$storage %||% NULL
  source <- source_or_null(Aop)
  if (is.matrix(source) && is.double(source)) {
    return(list(
      diagnostics = native_dense_svd_certificate(source, d, u, v, tol = tol),
      norm_bound_type = "frobenius_exact"
    ))
  }
  if (identical(storage, "dgCMatrix")) {
    A <- Aop$metadata$matrix
    return(list(
      diagnostics = .Call(
        "eigencore_csc_svd_certificate",
        methods::slot(A, "i"),
        methods::slot(A, "p"),
        methods::slot(A, "x"),
        methods::slot(A, "Dim"),
        as.numeric(d),
        as.matrix(u),
        as.matrix(v),
        as.numeric(Aop$metadata$frobenius_norm),
        as.numeric(tol),
        PACKAGE = "eigencore"
      ),
      norm_bound_type = "frobenius_metadata"
    ))
  }
  if (identical(storage, "ddiMatrix")) {
    A <- Aop$metadata$matrix
    return(list(
      diagnostics = .Call(
        "eigencore_diagonal_svd_certificate",
        methods::slot(A, "x"),
        methods::slot(A, "Dim"),
        identical(methods::slot(A, "diag"), "U"),
        as.numeric(d),
        as.matrix(u),
        as.matrix(v),
        as.numeric(Aop$metadata$frobenius_norm),
        as.numeric(tol),
        PACKAGE = "eigencore"
      ),
      norm_bound_type = "frobenius_metadata"
    ))
  }
  NULL
}

#' @keywords internal
native_builtin_svd_certificate_cached_av <- function(Aop, d, u, v, Av, tol = 1e-8) {
  storage <- Aop$metadata$storage %||% NULL
  source <- source_or_null(Aop)
  if (is.matrix(source) && is.double(source)) {
    return(list(
      diagnostics = native_dense_svd_certificate_cached_av(source, d, u, v, Av, tol = tol),
      norm_bound_type = "frobenius_exact"
    ))
  }
  if (identical(storage, "dgCMatrix")) {
    A <- Aop$metadata$matrix
    return(list(
      diagnostics = .Call(
        "eigencore_csc_svd_certificate_cached_av",
        methods::slot(A, "i"),
        methods::slot(A, "p"),
        methods::slot(A, "x"),
        methods::slot(A, "Dim"),
        as.numeric(d),
        as.matrix(u),
        as.matrix(v),
        as.matrix(Av),
        as.numeric(Aop$metadata$frobenius_norm),
        as.numeric(tol),
        PACKAGE = "eigencore"
      ),
      norm_bound_type = "frobenius_metadata"
    ))
  }
  if (identical(storage, "ddiMatrix")) {
    A <- Aop$metadata$matrix
    return(list(
      diagnostics = .Call(
        "eigencore_diagonal_svd_certificate_cached_av",
        methods::slot(A, "x"),
        methods::slot(A, "Dim"),
        identical(methods::slot(A, "diag"), "U"),
        as.numeric(d),
        as.matrix(u),
        as.matrix(v),
        as.matrix(Av),
        as.numeric(Aop$metadata$frobenius_norm),
        as.numeric(tol),
        PACKAGE = "eigencore"
      ),
      norm_bound_type = "frobenius_metadata"
    ))
  }
  NULL
}

#' @keywords internal
matrix_norm <- function(x) {
  norm(x, type = "F")
}

#' @keywords internal
max_residual_value <- function(x) {
  if (is.list(x)) {
    vals <- unlist(x, use.names = FALSE)
  } else {
    vals <- x
  }
  if (!length(vals)) NA_real_ else max(vals)
}

#' @keywords internal
eigen_backward_scale <- function(norm_A, norm_B, values, vectors) {
  pmax((norm_A + abs(values) * norm_B) * pmax(col_norms(vectors), .Machine$double.eps),
       .Machine$double.eps)
}

#' @keywords internal
svd_backward_scale <- function(norm_A, d) {
  rep(max(norm_A, .Machine$double.eps), length(d))
}

#' @keywords internal
operator_norm_for_certificate <- function(op) {
  operator_norm_for_certificate_info(op)$value
}

#' @keywords internal
operator_norm_for_certificate_info <- function(op) {
  if (!is.null(op$metadata$frobenius_norm)) {
    return(list(
      value = op$metadata$frobenius_norm,
      norm_bound_type = "frobenius_metadata",
      scale_is_estimate = FALSE
    ))
  }
  src <- source_or_null(op)
  if (!is.null(src)) {
    return(list(
      value = matrix_norm(as.matrix(src)),
      norm_bound_type = "frobenius_exact",
      scale_is_estimate = FALSE
    ))
  }
  list(
    value = estimate_operator_frobenius_norm(op),
    norm_bound_type = "frobenius_hutchinson_estimate",
    scale_is_estimate = TRUE
  )
}

#' @keywords internal
estimate_operator_frobenius_norm <- function(op) {
  # Hutchinson-style Frobenius estimate used only when a matrix-free operator
  # has no explicit norm metadata. It keeps certificate scaling consistent in
  # form, but native solvers should pass exact or bounded operator norms.
  probes <- min(8L, max(2L, op$dim[2L]))
  norms <- numeric(probes)
  for (i in seq_len(probes)) {
    z <- matrix(sample(c(-1, 1), op$dim[2L], replace = TRUE), op$dim[2L], 1L)
    norms[[i]] <- sum(apply_operator(op, z)^2)
  }
  sqrt(mean(norms))
}
