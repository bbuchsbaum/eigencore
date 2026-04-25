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
certify_eigen <- function(A, values, vectors, B = NULL, tol = 1e-8) {
  k <- length(values)
  residuals <- dense_eigen_residuals(A, values, vectors, B)
  scale <- eigen_backward_scale(matrix_norm(A), if (is.null(B)) 1 else matrix_norm(B), values, vectors)
  backward <- residuals / pmax(scale, .Machine$double.eps)
  orth <- orthogonality_loss(vectors, B = B)
  new_certificate(
    tol = tol,
    residuals = residuals,
    backward_error = backward,
    orthogonality = orth,
    converged = backward <= tol,
    scale = scale
  )
}

#' @keywords internal
certify_eigen_operator <- function(Aop, values, vectors, Bop = NULL, tol = 1e-8) {
  k <- length(values)
  Av <- apply_operator(Aop, vectors)
  Bv <- if (is.null(Bop)) vectors else apply_operator(Bop, vectors)
  residual_matrix <- Av - sweep(Bv, 2L, values, `*`)
  residuals <- col_norms(residual_matrix)
  scale <- eigen_backward_scale(
    operator_norm_for_certificate(Aop),
    if (is.null(Bop)) 1 else operator_norm_for_certificate(Bop),
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
    scale = scale
  )
}

#' @keywords internal
certify_svd <- function(A, d, u, v, tol = 1e-8) {
  residuals <- dense_svd_residuals(A, d, u, v)
  scale <- svd_backward_scale(matrix_norm(A), d)
  backward <- residuals$combined / scale
  orth_u <- orthogonality_loss(u)
  orth_v <- orthogonality_loss(v)
  new_certificate(
    tol = tol,
    residuals = residuals,
    backward_error = backward,
    orthogonality = c(U = orth_u, V = orth_v),
    converged = backward <= tol,
    scale = scale
  )
}

#' @keywords internal
certify_svd_operator <- function(Aop, d, u, v, tol = 1e-8) {
  left_residual_matrix <- apply_operator(Aop, v) - sweep(u, 2L, d, `*`)
  right_residual_matrix <- apply_adjoint_operator(Aop, u) - sweep(v, 2L, d, `*`)
  left <- col_norms(left_residual_matrix)
  right <- col_norms(right_residual_matrix)
  combined <- sqrt(left^2 + right^2)
  scale <- svd_backward_scale(operator_norm_for_certificate(Aop), d)
  backward <- combined / scale
  orth_u <- max(abs(crossprod(u) - diag(length(d))))
  orth_v <- max(abs(crossprod(v) - diag(length(d))))
  new_certificate(
    tol = tol,
    residuals = list(left = left, right = right, combined = combined),
    backward_error = backward,
    orthogonality = c(U = orth_u, V = orth_v),
    converged = backward <= tol,
    scale = scale
  )
}

#' @keywords internal
new_certificate <- function(tol, residuals, backward_error, orthogonality,
                            converged, scale, notes = character()) {
  cert <- list(
    passed = all(converged),
    tolerance = tol,
    max_backward_error = if (length(backward_error)) max(backward_error) else NA_real_,
    max_residual = max_residual_value(residuals),
    max_orthogonality_loss = if (length(orthogonality)) max(orthogonality) else NA_real_,
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
    notes = note
  )
}

#' @export
print.eigencore_certificate <- function(x, ...) {
  cat("eigencore certificate\n")
  cat("  passed:", x$passed, "\n")
  cat("  tolerance:", format(x$tolerance), "\n")
  cat("  max residual:", format(x$max_residual), "\n")
  cat("  max backward error:", format(x$max_backward_error), "\n")
  cat("  max orthogonality loss:", format(x$max_orthogonality_loss), "\n")
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
  if (!is.null(op$metadata$frobenius_norm)) {
    return(op$metadata$frobenius_norm)
  }
  src <- source_or_null(op)
  if (!is.null(src)) {
    return(matrix_norm(as.matrix(src)))
  }
  estimate_operator_frobenius_norm(op)
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
