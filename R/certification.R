#' Extract a result certificate.
#'
#' @param x An eigencore result object.
#' @param ... Reserved for future methods.
#' @return The `eigencore_certificate` object stored on `x`, or `NULL` if the
#'   result does not carry a certificate field.
#' @details
#' Backward errors use the standard normwise 2-norm definition:
#' `||A x - lambda B x|| / ((||A||_2 + |lambda| ||B||_2) ||x||)` for eigenpairs
#' and `sqrt(||A v - sigma u||^2 + ||A^H u - sigma v||^2) / ||A||_2` for
#' singular triplets. The norms in the denominator are exact or LOWER bounds
#' (never estimates), so the reported backward error is never smaller than the
#' true one and `passed` is sound. Fields describing the scale:
#' \describe{
#'   \item{`norm_bound_type`}{`"two_norm_exact"` or `"two_norm_lower_bound"`
#'     (eigen certificates report `A+B` parts; `B = I` is
#'     `"identity_exact"`).}
#'   \item{`norm_source`}{Where each value came from: `"diagonal"`,
#'     `"full_spectrum"`, `"metadata"` (an operator's `metadata$two_norm`),
#'     `"column_norms"`, `"applied_vectors"` (`||A x|| / ||x||` for the
#'     certified vectors), `"ritz"` (residual-corrected Ritz values),
#'     `"frobenius_rank_bound"`, or `"lanczos"` (a short deterministic
#'     Krylov estimate run only when it could change the verdict).}
#'   \item{`norm_values`}{The values used, named `A` (and `B`).}
#'   \item{`scale_is_estimate`}{Always `FALSE` for built-in certificates.}
#' }
#' @examples
#' fit <- eig_partial(diag(c(3, 2, 1)), k = 1, target = largest())
#' cert <- certificate(fit)
#' cert$passed
#' cert$max_residual
certificate <- function(x, ...) {
  x$certificate
}

#' Extract diagnostics.
#'
#' @param x An eigencore result object.
#' @param ... Reserved for future methods.
#' @return A named list of diagnostic fields, including residuals, backward
#'   errors, orthogonality diagnostics, iteration counts, method/plan metadata,
#'   warnings, and any available left-eigenvector diagnostics.
#' @examples
#' fit <- eig_partial(diag(c(3, 2, 1)), k = 1, target = largest())
#' d <- diagnostics(fit)
#' d$nconv
#' d$method
diagnostics <- function(x, ...) {
  restart <- if (is.list(x$restart)) x$restart else NULL
  out <- list(
    residuals = x$residuals,
    backward_error = x$backward_error,
    orthogonality = x$orthogonality,
    nconv = x$nconv,
    iterations = x$iterations,
    matvecs = x$matvecs,
    work = work(x),
    preconditioner_calls = x$preconditioner_calls,
    convergence_history = x$convergence_history,
    restart = restart,
    stage_seconds = x$stage_seconds %||% restart$stage_seconds %||% numeric(),
    preconditioner = x$preconditioner %||% restart$preconditioner %||% NULL,
    locked = x$locked,
    method = x$method,
    plan = x$plan,
    warnings = x$warnings
  )
  # Warm-start provenance is specific to the Hermitian Lanczos eigen path;
  # append it only when present so SVD and other result schemas are unchanged.
  if (!is.null(x$start_source)) {
    out$start_source <- x$start_source
    out$initial_subspace <- x$initial_subspace
    out$operator_block_calls <- x$operator_block_calls
    out$operator_columns <- x$operator_columns
    out$certification_operator_columns <- x$certification_operator_columns
  }
  if (!is.null(x$left_vectors) && !is.null(x$left_certificate)) {
    out$left_vectors <- x$left_vectors
    out$left_certificate <- x$left_certificate
    out$biorthogonality <- x$biorthogonality
  }
  out
}

#' Extract computed values.
#'
#' @param x An eigencore result object.
#' @param ... Reserved for future methods.
#' @return A numeric or complex vector of computed eigenvalues, singular values,
#'   or generalized singular values.
#' @examples
#' fit <- eig_partial(diag(c(3, 2, 1)), k = 2, target = largest())
#' values(fit)
values <- function(x, ...) {
  x$values
}

#' Extract homogeneous generalized coordinates.
#'
#' Generalized dense-pencil and generalized Schur results store eigenvalues in
#' homogeneous form as `alpha / beta`; generalized SVD results use the same
#' alpha/beta convention for generalized singular values. `alpha_beta()`
#' exposes those coordinates along with the finite/infinite/undefined
#' classification computed by eigencore.
#'
#' @param x An eigencore result with homogeneous `alpha` and `beta` fields.
#' @param ... Reserved for future methods.
#' @return A list containing `alpha`, `beta`, and any available `values`,
#'   `classification`, `finite`, `infinite`, and `undefined` fields. Results
#'   that record how the finite/infinite/undefined labels were decided also
#'   include a `classification_policy` list with the policy name, the
#'   tolerance, the per-coordinate zero thresholds, and any pencil norms used
#'   for norm-scaled classification. A `reason` explains exact structural
#'   policies used by GSVD and transformed sparse-pencil results.
#' @examples
#' A <- diag(c(2, 3, 0))
#' B <- diag(c(1, 0, 0))
#' fit <- eig_full(A, B = B, structure = general())
#' alpha_beta(fit)$classification
#' @export
alpha_beta <- function(x, ...) {
  if (is.null(x$alpha) || is.null(x$beta)) {
    stop(
      "alpha/beta coordinates are not available for this result.",
      call. = FALSE
    )
  }
  out <- list(alpha = x$alpha, beta = x$beta)
  for (nm in c("values", "classification", "finite", "infinite", "undefined",
               "classification_policy")) {
    if (!is.null(x[[nm]])) {
      out[[nm]] <- x[[nm]]
    }
  }
  out
}

#' Extract eigenvectors.
#'
#' @param x An eigencore eigen result object.
#' @param ... Reserved for future methods.
#' @return A matrix whose columns are computed right eigenvectors, or `NULL`
#'   when vectors were not requested or are unavailable.
#' @examples
#' fit <- eig_partial(diag(c(3, 2, 1)), k = 2, target = largest())
#' dim(vectors(fit))
vectors <- function(x, ...) {
  x$vectors
}

#' Extract left singular vectors or left eigenvectors.
#'
#' For SVD results this returns the left singular vectors `U`. For
#' nonsymmetric and dense general-pencil eigen results this returns left
#' eigenvectors when the solver computed them (for example, the dense
#' general-pencil `eig_full()` path, which computes left generalized
#' eigenvectors satisfying `w^H A = lambda w^H B`).
#'
#' @param x An eigencore SVD or eigen result object.
#' @param ... Reserved for future methods.
#' @return A matrix of left singular vectors or left eigenvectors, or `NULL`
#'   when the result does not contain a left-vector field.
left_vectors <- function(x, ...) {
  result_field(x, c("left_vectors", "u", "U"))
}

#' Extract right singular vectors.
#'
#' @param x An eigencore SVD or nonsymmetric eigen result object.
#' @param ... Reserved for future methods.
#' @return A matrix of right singular vectors or right eigenvectors, or `NULL`
#'   when the result does not contain a right-vector field.
right_vectors <- function(x, ...) {
  result_field(x, c("right_vectors", "v", "V", "vectors"))
}

result_field <- function(x, names) {
  for (nm in names) {
    value <- x[[nm, exact = TRUE]]
    if (!is.null(value)) {
      return(value)
    }
  }
  NULL
}

#' Extract residual diagnostics.
#'
#' Methods for the [stats::residuals()] generic: return the per-pair (or
#' per-triplet) residual norms stored in an eigencore result or certificate.
#'
#' @param object An eigencore result or certificate object.
#' @param ... Reserved for future methods.
#' @name residuals
#' @rdname residuals
residuals.eigencore_eigen_result <- function(object, ...) {
  object$residuals
}

#' @rdname residuals
residuals.eigencore_svd_result <- function(object, ...) {
  object$residuals
}

#' @rdname residuals
residuals.eigencore_certificate <- function(object, ...) {
  object$residuals
}

#' Extract backward-error diagnostics.
#'
#' @param x An eigencore result object.
#' @param ... Reserved for future methods.
#' @return A numeric vector of per-pair or per-triplet backward-error
#'   estimates.
backward_error <- function(x, ...) {
  x$backward_error
}

#' @keywords internal
certify_eigen <- function(A, values, vectors, B = NULL, tol = 1e-8,
                          require_orthogonality = TRUE,
                          full_spectrum = NULL) {
  if (is.complex(A) || is.complex(values) || is.complex(vectors) ||
      (!is.null(B) && is.complex(B))) {
    return(certify_dense_eigen_r_residual(
      A,
      values,
      vectors,
      B = B,
      tol = tol,
      require_orthogonality = require_orthogonality,
      full_spectrum = full_spectrum
    ))
  }
  diag <- native_dense_eigen_certificate(A, values, vectors, B = B, tol = tol)
  norm <- eigen_two_norm_backward(
    A, values, diag$residuals, diag$vector_norms, tol,
    B = B,
    free_A = list(
      applied_bound(diag$norm_A_applied_bound),
      full_spectrum_bound(full_spectrum, B)
    ),
    free_B = list(applied_bound(diag$norm_B_applied_bound)),
    structural_A = native_column_bound(diag$norm_A_column_bound),
    structural_B = if (!is.null(B)) native_column_bound(diag$norm_B_column_bound)
  )
  eigen_certificate_from_norms(
    norm, tol, diag$residuals, diag$orthogonality,
    require_orthogonality = require_orthogonality
  )
}

# All eigenvalues of a Hermitian A (B = NULL) give ||A||_2 = max |lambda|
# exactly (up to rounding).
#' @keywords internal
full_spectrum_bound <- function(full_spectrum, B = NULL) {
  if (is.null(full_spectrum) || !is.null(B) || !length(full_spectrum) ||
      !all(is.finite(Mod(full_spectrum)))) {
    return(NULL)
  }
  cert_norm_bound(max(Mod(full_spectrum)), exact = TRUE,
                  source = "full_spectrum")
}

#' @keywords internal
eigen_certificate_from_norms <- function(norm, tol, residuals, orthogonality,
                                         converged = NULL,
                                         notes = character(),
                                         certificate_type = "residual_backward_error",
                                         require_orthogonality = TRUE) {
  if (is.null(converged)) {
    converged <- is.finite(norm$backward) & norm$backward <= tol
  }
  new_certificate(
    tol = tol,
    residuals = residuals,
    backward_error = norm$backward,
    orthogonality = orthogonality,
    converged = converged,
    scale = norm$scale,
    notes = notes,
    certificate_type = certificate_type,
    norm_bound_type = norm$norm_bound_type,
    norm_source = norm$norm_source,
    norm_values = norm$norms,
    frobenius_norm = norm$frobenius_norm,
    require_orthogonality = require_orthogonality
  )
}

#' @keywords internal
certify_dense_eigen_r_residual <- function(A, values, vectors, B = NULL,
                                           tol = 1e-8,
                                           require_orthogonality = TRUE,
                                           full_spectrum = NULL) {
  A <- as.matrix(A)
  vectors <- as.matrix(vectors)
  values <- as.vector(values)
  k <- length(values)
  if (ncol(vectors) != k) {
    stop("values and vectors must have compatible dimensions.", call. = FALSE)
  }
  Av <- A %*% vectors
  Bv <- if (is.null(B)) vectors else as.matrix(B) %*% vectors
  residual_matrix <- Av - sweep(Bv, 2L, values, `*`)
  residuals <- col_norms(residual_matrix)
  vec_norms <- col_norms(vectors)
  norm <- eigen_two_norm_backward(
    A, values, residuals, vec_norms, tol,
    B = B,
    free_A = list(
      bound_from_ratios(col_norms(Av), vec_norms, "applied_vectors"),
      full_spectrum_bound(full_spectrum, B)
    ),
    free_B = if (!is.null(B)) {
      list(bound_from_ratios(col_norms(Bv), vec_norms, "applied_vectors"))
    }
  )
  gram <- if (is.null(B)) certificate_gram(vectors) else certificate_gram(vectors, Bv)
  orth <- max(abs(gram - diag(k)))
  eigen_certificate_from_norms(
    norm, tol, residuals, orth,
    require_orthogonality = require_orthogonality
  )
}

#' @keywords internal
certify_eigen_operator <- function(Aop, values, vectors, Bop = NULL, tol = 1e-8) {
  work_phase <- work_phase_enter("certification")
  on.exit(work_phase_exit(work_phase), add = TRUE)
  native <- native_builtin_eigen_certificate(Aop, values, vectors, Bop = Bop, tol = tol)
  if (!is.null(native)) {
    diag <- native$diagnostics
    norm <- eigen_two_norm_backward(
      Aop, values, diag$residuals,
      diag$vector_norms %||% col_norms(vectors), tol,
      B = Bop,
      free_A = list(applied_bound(diag$norm_A_applied_bound)),
      free_B = list(applied_bound(diag$norm_B_applied_bound)),
      structural_A = native_column_bound(diag$norm_A_column_bound),
      structural_B = native_column_bound(diag$norm_B_column_bound)
    )
    return(eigen_certificate_from_norms(
      norm, tol, diag$residuals, diag$orthogonality
    ))
  }

  k <- length(values)
  Av <- apply_operator(Aop, vectors)
  Bv <- if (is.null(Bop)) vectors else apply_operator(Bop, vectors)
  residual_matrix <- Av - sweep(Bv, 2L, values, `*`)
  residuals <- col_norms(residual_matrix)
  vec_norms <- col_norms(vectors)
  norm <- eigen_two_norm_backward(
    Aop, values, residuals, vec_norms, tol,
    B = Bop,
    free_A = list(bound_from_ratios(col_norms(Av), vec_norms, "applied_vectors")),
    free_B = if (!is.null(Bop)) {
      list(bound_from_ratios(col_norms(Bv), vec_norms, "applied_vectors"))
    }
  )
  gram <- if (is.null(Bop)) certificate_gram(vectors) else certificate_gram(vectors, Bv)
  orth <- max(abs(gram - diag(k)))
  eigen_certificate_from_norms(norm, tol, residuals, orth)
}

#' @keywords internal
certify_dense_general_eigen <- function(A, values, vectors, tol = 1e-8) {
  A <- as.matrix(A)
  vectors <- as.matrix(vectors)
  values <- as.vector(values)
  k <- length(values)
  if (ncol(vectors) != k) {
    stop("values and vectors must have compatible dimensions.", call. = FALSE)
  }
  Av <- A %*% vectors
  residual_matrix <- Av - sweep(vectors, 2L, values, `*`)
  residuals <- col_norms(residual_matrix)
  vec_norms <- col_norms(vectors)
  norm <- eigen_two_norm_backward(
    A, values, residuals, vec_norms, tol,
    free_A = list(bound_from_ratios(col_norms(Av), vec_norms, "applied_vectors"))
  )
  gram <- certificate_gram(vectors)
  orth <- max(abs(gram - diag(k)))
  eigen_certificate_from_norms(
    norm, tol, residuals, orth,
    notes = "right residual certificate for dense general eigenpairs; eigenvector orthogonality is not required",
    certificate_type = "right_residual_backward_error",
    require_orthogonality = FALSE
  )
}

# Apply a real operator (or its adjoint) to a possibly complex block as two
# real applies, so sparse or matrix-free sources are never densified for
# complex eigenvectors (C39). Complex-typed operators get the block as is.
#' @keywords internal
apply_operator_split_complex <- function(op, X, adjoint = FALSE) {
  f <- if (isTRUE(adjoint)) apply_adjoint_operator else apply_operator
  if (!is.complex(X) || !identical(op$dtype %||% "double", "double")) {
    return(as.matrix(f(op, X)))
  }
  X <- as.matrix(X)
  p <- ncol(X)
  out <- as.matrix(f(op, cbind(Re(X), Im(X))))
  if (is.complex(out)) {
    return(out[, seq_len(p), drop = FALSE] +
             1i * out[, p + seq_len(p), drop = FALSE])
  }
  matrix(
    complex(
      real = out[, seq_len(p), drop = FALSE],
      imaginary = out[, p + seq_len(p), drop = FALSE]
    ),
    nrow(out),
    p
  )
}

#' @keywords internal
certify_general_eigen_operator <- function(Aop, values, vectors, tol = 1e-8) {
  work_phase <- work_phase_enter("certification")
  on.exit(work_phase_exit(work_phase), add = TRUE)
  values <- as.vector(values)
  vectors <- as.matrix(vectors)
  Av <- apply_operator_split_complex(Aop, vectors)
  residual_matrix <- Av - sweep(vectors, 2L, values, `*`)
  residuals <- col_norms(residual_matrix)
  vec_norms <- col_norms(vectors)
  norm <- eigen_two_norm_backward(
    Aop, values, residuals, vec_norms, tol,
    free_A = list(bound_from_ratios(col_norms(Av), vec_norms, "applied_vectors"))
  )
  gram <- certificate_gram(vectors)
  orth <- max(abs(gram - diag(length(values))))
  eigen_certificate_from_norms(
    norm, tol, residuals, orth,
    notes = "right residual certificate for general eigenpairs; eigenvector orthogonality is not required",
    certificate_type = "right_residual_backward_error",
    require_orthogonality = FALSE
  )
}

#' @keywords internal
certify_left_eigen_operator <- function(Aop, values, left_vectors,
                                        right_vectors = NULL, tol = 1e-8) {
  work_phase <- work_phase_enter("certification")
  on.exit(work_phase_exit(work_phase), add = TRUE)
  values <- as.vector(values)
  left_vectors <- as.matrix(left_vectors)
  if (ncol(left_vectors) != length(values)) {
    stop("values and left_vectors must have compatible dimensions.", call. = FALSE)
  }

  Astar_w <- apply_operator_split_complex(Aop, left_vectors, adjoint = TRUE)
  residual_matrix <- Astar_w - sweep(left_vectors, 2L, values, `*`)
  left_residuals <- col_norms(residual_matrix)
  vec_norms <- col_norms(left_vectors)
  # ||A^H w|| / ||w|| <= ||A^H||_2 = ||A||_2.
  norm <- eigen_two_norm_backward(
    Aop, values, left_residuals, vec_norms, tol,
    free_A = list(bound_from_ratios(col_norms(Astar_w), vec_norms, "applied_vectors"))
  )

  biorthogonality <- numeric()
  if (!is.null(right_vectors)) {
    right_vectors <- as.matrix(right_vectors)
    cross <- crossprod(left_vectors, right_vectors)
    biorthogonality <- max(abs(cross - diag(length(values))))
  }

  eigen_certificate_from_norms(
    norm, tol, list(left = left_residuals), biorthogonality,
    notes = "left residual and biorthogonality certificate for nonsymmetric eigenpairs",
    certificate_type = "left_residual_biorthogonal_backward_error",
    require_orthogonality = !is.null(right_vectors)
  )
}

#' @keywords internal
certify_eigen_operator_residuals <- function(Aop, values, vectors, residuals,
                                             Bop = NULL, tol = 1e-8) {
  work_phase <- work_phase_enter("certification")
  on.exit(work_phase_exit(work_phase), add = TRUE)
  vec_norms <- col_norms(vectors)
  Bv <- NULL
  if (is.null(Bop)) {
    orth <- orthogonality_loss(vectors)
  } else {
    Bv <- apply_operator(Bop, vectors)
    orth <- max(abs(crossprod(vectors, Bv) - diag(length(values))))
  }
  bv_norms <- if (is.null(Bv)) vec_norms else col_norms(Bv)
  norm <- eigen_two_norm_backward(
    Aop, values, residuals, vec_norms, tol,
    B = Bop,
    free_A = list(eigen_ritz_bound(values, residuals, vec_norms, bv_norms)),
    free_B = if (!is.null(Bop)) {
      list(bound_from_ratios(bv_norms, vec_norms, "applied_vectors"))
    }
  )
  eigen_certificate_from_norms(norm, tol, residuals, orth)
}

#' @keywords internal
certify_svd <- function(A, d, u, v, tol = 1e-8) {
  if (is.complex(A) || is.complex(u) || is.complex(v)) {
    A <- as.matrix(A)
    d <- as.vector(d)
    u <- as.matrix(u)
    v <- as.matrix(v)
    Av <- A %*% v
    Ahu <- Conj(t(A)) %*% u
    left <- col_norms(Av - sweep(u, 2L, d, `*`))
    right <- col_norms(Ahu - sweep(v, 2L, d, `*`))
    orth_u <- max(abs(certificate_gram(u) - diag(length(d))))
    orth_v <- max(abs(certificate_gram(v) - diag(length(d))))
    applied <- max(0, c(
      bound_from_ratios(col_norms(Av), col_norms(v), "applied_vectors")$value,
      bound_from_ratios(col_norms(Ahu), col_norms(u), "applied_vectors")$value
    ))
    return(svd_certificate_from_residuals(
      A, d, left, right, c(U = orth_u, V = orth_v), tol,
      u = u, v = v, applied = applied
    ))
  }
  diag <- native_dense_svd_certificate(A, d, u, v, tol = tol)
  svd_certificate_from_native_diagnostics(A, d, diag, tol, u = u, v = v)
}

#' @keywords internal
certify_svd_operator <- function(Aop, d, u, v, tol = 1e-8) {
  work_phase <- work_phase_enter("certification")
  on.exit(work_phase_exit(work_phase), add = TRUE)
  native <- native_builtin_svd_certificate(Aop, d, u, v, tol = tol)
  if (!is.null(native)) {
    return(svd_certificate_from_native_diagnostics(
      Aop, d, native$diagnostics, tol, u = u, v = v
    ))
  }

  Av <- apply_operator(Aop, v)
  Atu <- apply_adjoint_operator(Aop, u)
  svd_certificate_from_applied(Aop, d, u, v, Av, Atu, tol)
}

#' @keywords internal
svd_certificate_from_applied <- function(Aop, d, u, v, Av, Atu, tol) {
  left <- col_norms(Av - sweep(u, 2L, d, `*`))
  right <- col_norms(Atu - sweep(v, 2L, d, `*`))
  applied <- max(0, c(
    bound_from_ratios(col_norms(Av), col_norms(v), "applied_vectors")$value,
    bound_from_ratios(col_norms(Atu), col_norms(u), "applied_vectors")$value
  ))
  orth_u <- max(abs(certificate_gram(u) - diag(length(d))))
  orth_v <- max(abs(certificate_gram(v) - diag(length(d))))
  svd_certificate_from_residuals(
    Aop, d, left, right, c(U = orth_u, V = orth_v), tol,
    u = u, v = v, applied = applied
  )
}

#' @keywords internal
certify_svd_operator_cached_av <- function(Aop, d, u, v, Av, tol = 1e-8,
                                           return_residual_vectors = FALSE) {
  work_phase <- work_phase_enter("certification")
  on.exit(work_phase_exit(work_phase), add = TRUE)
  if (is.null(Av)) {
    cert <- certify_svd_operator(Aop, d, u, v, tol = tol)
    if (isTRUE(return_residual_vectors)) {
      return(list(certificate = cert, right_residual_vectors = NULL))
    }
    return(cert)
  }
  # The cache trusts the caller's Av to be apply_operator(Aop, v). A stale or
  # non-finite Av would silently produce wrong residuals and a passing
  # certificate. Guard the obvious failure modes (non-finite, mismatched shape)
  # at the entry; cache invalidation against operator/v fingerprints belongs in
  # a follow-up since it requires plumbing fingerprints through call sites.
  if (anyNA(Av) || any(!is.finite(Av))) {
    stop(
      "certify_svd_operator_cached_av: cached Av contains non-finite values; ",
      "the cache is stale or corrupted. Recompute apply_operator(Aop, v).",
      call. = FALSE
    )
  }
  if (nrow(Av) != Aop$dim[[1L]] || ncol(Av) != length(d)) {
    stop(
      "certify_svd_operator_cached_av: Av must have nrow == nrow(Aop) (",
      Aop$dim[[1L]], ") and one column per singular value (", length(d),
      "); got ", nrow(Av), " x ", ncol(Av), ".",
      call. = FALSE
    )
  }
  if (!isTRUE(return_residual_vectors)) {
    native <- native_builtin_svd_certificate_cached_av(Aop, d, u, v, Av, tol = tol)
    if (!is.null(native)) {
      return(svd_certificate_from_native_diagnostics(
        Aop, d, native$diagnostics, tol, u = u, v = v
      ))
    }
    return(certify_svd_operator(Aop, d, u, v, tol = tol))
  }
  Av <- as.matrix(Av)
  Atu <- apply_adjoint_operator(Aop, u)
  right_residual_matrix <- Atu - sweep(v, 2L, d, `*`)
  cert <- svd_certificate_from_applied(Aop, d, u, v, Av, Atu, tol)
  list(certificate = cert, right_residual_vectors = right_residual_matrix)
}

#' @keywords internal
certify_svd_operator_cached_sides <- function(Aop, d, u, v, Av, Atu,
                                             tol = 1e-8) {
  work_phase <- work_phase_enter("certification")
  on.exit(work_phase_exit(work_phase), add = TRUE)
  Av <- as.matrix(Av)
  Atu <- as.matrix(Atu)
  # Same trust contract as certify_svd_operator_cached_av: stale or non-finite
  # cached sides silently produce wrong residuals. Guard the obvious failure
  # modes at entry.
  if (anyNA(Av) || any(!is.finite(Av))) {
    stop(
      "certify_svd_operator_cached_sides: cached Av contains non-finite values; ",
      "the cache is stale or corrupted. Recompute apply_operator(Aop, v).",
      call. = FALSE
    )
  }
  if (anyNA(Atu) || any(!is.finite(Atu))) {
    stop(
      "certify_svd_operator_cached_sides: cached Atu contains non-finite values; ",
      "the cache is stale or corrupted. Recompute apply_adjoint_operator(Aop, u).",
      call. = FALSE
    )
  }
  if (nrow(Av) != Aop$dim[[1L]] || ncol(Av) != length(d)) {
    stop("Av must have nrow equal to nrow(Aop) and one column per singular value.",
         call. = FALSE)
  }
  if (nrow(Atu) != Aop$dim[[2L]] || ncol(Atu) != length(d)) {
    stop("Atu must have nrow equal to ncol(Aop) and one column per singular value.",
         call. = FALSE)
  }
  svd_certificate_from_applied(Aop, d, u, v, Av, Atu, tol)
}

#' @keywords internal
new_certificate <- function(tol, residuals, backward_error, orthogonality,
                            converged, scale, notes = character(),
                            certificate_type = "residual_backward_error",
                            norm_bound_type = "unspecified",
                            scale_is_estimate = FALSE,
                            require_orthogonality = TRUE,
                            norm_source = NA_character_,
                            norm_values = numeric(),
                            frobenius_norm = NA_real_) {
  # A stochastic (unbounded) scale would make the backward error neither an
  # upper nor a lower bound, so such certificates never pass. Every built-in
  # certificate now scales by a two-norm value that is exact or a lower bound
  # (C12), which keeps `passed` sound; the flag remains for callers that
  # construct certificates from their own estimates.
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
    norm_source = norm_source,
    norm_values = norm_values,
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
  if (length(frobenius_norm) == 1L && is.finite(frobenius_norm)) {
    cert$frobenius_norm <- frobenius_norm
  }
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
  if (!is.null(x$norm_source) && !is.na(x$norm_source[[1L]])) {
    cat("  norm source:", x$norm_source, "\n")
  }
  cat("  scale estimated:", x$scale_is_estimate, "\n")
  cat("  max residual:", format(x$max_residual), "\n")
  cat("  max backward error:", format(x$max_backward_error), "\n")
  cat("  max orthogonality loss:", format(x$max_orthogonality_loss), "\n")
  cat("  orthogonality tolerance:", format(x$orthogonality_tolerance), "\n")
  cat("  orthogonality required:", x$orthogonality_required, "\n")
  if (!is.null(x$target_completeness)) {
    cat("  target completeness:", x$target_completeness, "\n")
  }
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
  if (is.complex(x)) {
    return(sqrt(colSums(Mod(as.matrix(x))^2)))
  }
  .Call("eigencore_col_norms", as.matrix(x), PACKAGE = "eigencore")
}

# Certificate Gram V^H W. For W = V the symmetric product uses crossprod(),
# which R evaluates with the BLAS symmetric rank-k update (dsyrk) (C32).
#' @keywords internal
certificate_gram <- function(x, y = NULL) {
  x <- as.matrix(x)
  if (is.null(y)) {
    if (is.complex(x)) {
      return(Conj(t(x)) %*% x)
    }
    return(crossprod(x))
  }
  y <- as.matrix(y)
  if (is.complex(x) || is.complex(y)) {
    return(Conj(t(x)) %*% y)
  }
  crossprod(x, y)
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

# Native residual kernels for built-in storages. Each returns the residual
# diagnostics plus the applied-vector norm bound computed in the same call;
# the backward error is (re)assembled in R with the two-norm bounds.
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
                                                    B = B_source, tol = tol)
      ))
    }
    return(NULL)
  }
  if (is.matrix(source) && is.double(source)) {
    return(list(
      diagnostics = native_dense_eigen_certificate(source, values, vectors, tol = tol)
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
        as.numeric(two_norm_structural_bound(Aop)$value),
        as.numeric(tol),
        PACKAGE = "eigencore"
      )
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
        as.numeric(two_norm_structural_bound(Aop)$value),
        as.numeric(tol),
        PACKAGE = "eigencore"
      )
    ))
  }
  NULL
}

#' @keywords internal
native_tridiagonal_eigen_certificate <- function(Aop, parts, values, vectors, tol = 1e-8) {
  diag <- .Call(
    "eigencore_tridiagonal_eigen_certificate",
    as.numeric(parts$diag),
    as.numeric(parts$upper),
    as.numeric(values),
    as.matrix(vectors),
    as.numeric(two_norm_structural_bound(Aop)$value),
    as.numeric(tol),
    PACKAGE = "eigencore"
  )
  norm <- eigen_two_norm_backward(
    Aop, values, diag$residuals, col_norms(vectors), tol,
    free_A = list(
      cert_norm_bound(diag$norm_A_column_bound, source = "column_norms"),
      applied_bound(diag$norm_A_applied_bound)
    )
  )
  eigen_certificate_from_norms(norm, tol, diag$residuals, diag$orthogonality)
}

#' @keywords internal
native_builtin_svd_certificate <- function(Aop, d, u, v, tol = 1e-8) {
  storage <- Aop$metadata$storage %||% NULL
  source <- source_or_null(Aop)
  if (is.matrix(source) && is.double(source)) {
    return(list(
      diagnostics = native_dense_svd_certificate(source, d, u, v, tol = tol)
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
        as.numeric(two_norm_structural_bound(Aop)$value),
        as.numeric(tol),
        PACKAGE = "eigencore"
      )
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
        as.numeric(two_norm_structural_bound(Aop)$value),
        as.numeric(tol),
        PACKAGE = "eigencore"
      )
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
      diagnostics = native_dense_svd_certificate_cached_av(source, d, u, v, Av, tol = tol)
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
        as.numeric(two_norm_structural_bound(Aop)$value),
        as.numeric(tol),
        PACKAGE = "eigencore"
      )
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
        as.numeric(two_norm_structural_bound(Aop)$value),
        as.numeric(tol),
        PACKAGE = "eigencore"
      )
    ))
  }
  NULL
}

#' @keywords internal
matrix_norm <- function(x) {
  # base::norm() drops imaginary parts, so complex inputs use Mod().
  if (inherits(x, "sparseMatrix")) {
    x <- methods::as(x, "generalMatrix")
    if (methods::.hasSlot(x, "x")) {
      return(sqrt(sum(Mod(x@x)^2)))
    }
    return(sqrt(sum(Mod(as.matrix(x))^2)))
  }
  if (is.complex(x)) {
    return(sqrt(sum(Mod(x)^2)))
  }
  norm(as.matrix(x), type = "F")
}

#' @keywords internal
matrix_norm_one <- function(x) {
  # base::norm() drops imaginary parts of complex matrices.
  x <- as.matrix(x)
  if (!length(x)) {
    return(0)
  }
  max(colSums(Mod(x)))
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

# Pre-solve spectral-norm value for an operator: exact where cheap (identity,
# diagonal, asserted metadata) and otherwise a structural LOWER bound (largest
# column norm, or ||A||_F / sqrt(min(m, n))). Native kernels use it as the
# floor of their relative convergence and breakdown scales; because it never
# exceeds ||A||_2 those decisions are at least as strict as the certificate.
# No operator applies and no randomness (the former Hutchinson Frobenius
# estimate is gone, C12).
#' @keywords internal
operator_norm_for_certificate_info <- function(op) {
  info <- two_norm_structural_bound(op)
  list(
    value = info$value,
    norm_bound_type = if (isTRUE(info$exact)) "two_norm_exact" else "two_norm_lower_bound",
    norm_source = info$source,
    scale_is_estimate = FALSE
  )
}

# Session-level memo for per-operator values such as the structural and
# Krylov two-norm bounds. Entries are keyed by the
# construction token linear_operator() assigns, and reused only while the
# operator still carries its construction-time metadata and apply closure, so
# an operator whose fields were replaced after construction recomputes. The
# memo lives outside the operator so solving never mutates a plan's
# serialised state. `value` is a promise, evaluated only on a miss.
.eigencore_operator_memo <- new.env(parent = emptyenv())

#' @keywords internal
operator_memoised_value <- function(op, key, value) {
  cache <- if (is.list(op)) op$cache else NULL
  valid <- is.environment(cache) &&
    is.character(cache$memo_token) &&
    identical(op$metadata, cache$metadata) &&
    identical(op$apply, cache$wrapped_apply) &&
    identical(as.integer(op$dim), as.integer(cache$dim))
  if (!valid) {
    return(value)
  }
  memo_key <- paste(cache$memo_token, key, sep = "|")
  hit <- .eigencore_operator_memo[[memo_key]]
  if (!is.null(hit)) {
    return(hit)
  }
  if (length(.eigencore_operator_memo) >= 512L) {
    rm(list = ls(.eigencore_operator_memo, all.names = TRUE),
       envir = .eigencore_operator_memo)
  }
  assign(memo_key, value, envir = .eigencore_operator_memo)
  value
}
