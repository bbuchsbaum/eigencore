#' @keywords internal
generalized_pencil_values <- function(alpha, beta,
                                      tol = sqrt(.Machine$double.eps)) {
  alpha <- as.vector(alpha)
  beta <- as.vector(beta)
  if (length(alpha) != length(beta)) {
    stop("alpha and beta must have the same length.", call. = FALSE)
  }
  if (!length(alpha)) {
    stop("alpha and beta must contain at least one pair.", call. = FALSE)
  }
  if (length(tol) != 1L || is.na(tol) || tol < 0) {
    stop("tol must be a single non-negative number.", call. = FALSE)
  }
  if (any(!is.finite(Mod(alpha))) || any(!is.finite(Mod(beta)))) {
    stop("alpha and beta must contain finite homogeneous coordinates.", call. = FALSE)
  }

  scale <- pmax(1, Mod(alpha), Mod(beta))
  alpha_zero <- Mod(alpha) <= tol * scale
  beta_zero <- Mod(beta) <= tol * scale
  finite <- !beta_zero
  infinite <- beta_zero & !alpha_zero
  undefined <- beta_zero & alpha_zero
  classification <- ifelse(finite, "finite",
                           ifelse(infinite, "infinite", "undefined"))

  complex_values <- is.complex(alpha) || is.complex(beta)
  values <- if (complex_values) {
    rep(NA_complex_, length(alpha))
  } else {
    rep(NA_real_, length(alpha))
  }
  values[finite] <- alpha[finite] / beta[finite]
  values[infinite] <- if (complex_values) {
    complex(real = Inf, imaginary = 0)
  } else {
    Inf
  }

  structure(
    list(
      values = values,
      alpha = alpha,
      beta = beta,
      classification = classification,
      finite = finite,
      infinite = infinite,
      undefined = undefined,
      alpha_zero = alpha_zero,
      beta_zero = beta_zero,
      tolerance = tol
    ),
    class = "eigencore_generalized_pencil_values"
  )
}

#' @keywords internal
certify_dense_generalized_pencil <- function(A, B, alpha, beta, vectors,
                                             tol = 1e-8,
                                             beta_tol = sqrt(.Machine$double.eps)) {
  A <- as.matrix(A)
  B <- as.matrix(B)
  vectors <- as.matrix(vectors)
  pencil <- generalized_pencil_values(alpha, beta, tol = beta_tol)
  generalized_pencil_validate_dimensions(A, B, vectors, length(pencil$values))

  residuals <- generalized_pencil_dense_residuals(A, B, pencil, vectors)
  generalized_pencil_certificate_from_residuals(
    pencil = pencil,
    residuals = residuals,
    norm_A = matrix_norm(A),
    norm_B = matrix_norm(B),
    norm_bound_type = "frobenius_exact+frobenius_exact",
    scale_is_estimate = FALSE,
    vectors = vectors,
    tol = tol
  )
}

#' @keywords internal
certify_generalized_pencil_operator <- function(Aop, Bop, alpha, beta, vectors,
                                                tol = 1e-8,
                                                beta_tol = sqrt(.Machine$double.eps)) {
  Aop <- as_operator(Aop)
  Bop <- as_operator(Bop)
  vectors <- as.matrix(vectors)
  pencil <- generalized_pencil_values(alpha, beta, tol = beta_tol)
  if (Aop$dim[1L] != Bop$dim[1L] || Aop$dim[2L] != Bop$dim[2L] ||
      Aop$dim[1L] != Aop$dim[2L] || nrow(vectors) != Aop$dim[2L] ||
      ncol(vectors) != length(pencil$values)) {
    stop("A, B, vectors, alpha, and beta must have compatible dimensions.", call. = FALSE)
  }

  residuals <- generalized_pencil_operator_residuals(Aop, Bop, pencil, vectors)
  norm_A <- operator_norm_for_certificate_info(Aop)
  norm_B <- operator_norm_for_certificate_info(Bop)
  generalized_pencil_certificate_from_residuals(
    pencil = pencil,
    residuals = residuals,
    norm_A = norm_A$value,
    norm_B = norm_B$value,
    norm_bound_type = paste(c(norm_A$norm_bound_type, norm_B$norm_bound_type),
                            collapse = "+"),
    scale_is_estimate = isTRUE(norm_A$scale_is_estimate) ||
      isTRUE(norm_B$scale_is_estimate),
    vectors = vectors,
    tol = tol
  )
}

#' @keywords internal
generalized_pencil_validate_dimensions <- function(A, B, vectors, k) {
  if (nrow(A) != ncol(A) || nrow(B) != ncol(B) ||
      nrow(A) != nrow(B) || nrow(vectors) != ncol(A) ||
      ncol(vectors) != k) {
    stop("A, B, vectors, alpha, and beta must have compatible dimensions.", call. = FALSE)
  }
  invisible(TRUE)
}

#' @keywords internal
generalized_pencil_dense_residuals <- function(A, B, pencil, vectors) {
  residuals <- rep(Inf, length(pencil$values))
  if (any(pencil$finite)) {
    idx <- which(pencil$finite)
    V <- vectors[, idx, drop = FALSE]
    Bv <- B %*% V
    residual_matrix <- A %*% V - sweep(Bv, 2L, pencil$values[idx], `*`)
    residuals[idx] <- col_norms(residual_matrix)
  }
  residuals
}

#' @keywords internal
generalized_pencil_operator_residuals <- function(Aop, Bop, pencil, vectors) {
  residuals <- rep(Inf, length(pencil$values))
  if (any(pencil$finite)) {
    idx <- which(pencil$finite)
    V <- vectors[, idx, drop = FALSE]
    Av <- generalized_pencil_apply_for_residual(Aop, V)
    Bv <- generalized_pencil_apply_for_residual(Bop, V)
    residual_matrix <- Av - sweep(Bv, 2L, pencil$values[idx], `*`)
    residuals[idx] <- col_norms(residual_matrix)
  }
  residuals
}

#' @keywords internal
generalized_pencil_apply_for_residual <- function(op, vectors) {
  op <- as_operator(op)
  source <- source_or_null(op) %||% op$metadata$matrix %||% NULL
  if (is.complex(vectors) && !is.null(source)) {
    return(as.matrix(as.matrix(source) %*% vectors))
  }
  apply_operator(op, vectors)
}

#' @keywords internal
generalized_pencil_certificate_from_residuals <- function(pencil, residuals,
                                                          norm_A, norm_B,
                                                          norm_bound_type,
                                                          scale_is_estimate,
                                                          vectors, tol) {
  scale <- rep(Inf, length(pencil$values))
  backward <- rep(Inf, length(pencil$values))
  if (any(pencil$finite)) {
    idx <- which(pencil$finite)
    scale[idx] <- eigen_backward_scale(
      norm_A,
      norm_B,
      pencil$values[idx],
      vectors[, idx, drop = FALSE]
    )
    backward[idx] <- residuals[idx] / pmax(scale[idx], .Machine$double.eps)
  }
  converged <- pencil$finite & backward <= tol
  notes <- generalized_pencil_certificate_notes(pencil)
  new_certificate(
    tol = tol,
    residuals = residuals,
    backward_error = backward,
    orthogonality = numeric(),
    converged = converged,
    scale = scale,
    notes = notes,
    certificate_type = "generalized_pencil_right_residual_backward_error",
    norm_bound_type = norm_bound_type,
    scale_is_estimate = isTRUE(scale_is_estimate),
    require_orthogonality = FALSE
  )
}

#' @keywords internal
generalized_pencil_certificate_notes <- function(pencil) {
  notes <- character()
  if (any(pencil$infinite)) {
    notes <- c(
      notes,
      "infinite generalized eigenvalues have beta equal to zero; finite residual certificate is unsupported"
    )
  }
  if (any(pencil$undefined)) {
    notes <- c(
      notes,
      "undefined generalized eigenvalues have alpha and beta equal to zero; residual certificate is unsupported"
    )
  }
  notes
}
