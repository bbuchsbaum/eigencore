#' @keywords internal
#'
#' Classify homogeneous generalized eigenvalue pairs.
#'
#' Two tolerance policies are supported. The default `per_pair_magnitude`
#' policy compares each `|alpha|`/`|beta|` against `tol * max(1, |alpha|,
#' |beta|)`. When `norm_A` and `norm_B` are supplied (LAPACK balanced
#' one-norms from DGGEVX, or plain one-norms of the input pencil), the
#' `pencil_norm_scaled` policy compares `|alpha|` against `tol * norm_A` and
#' `|beta|` against `tol * norm_B` instead. LAPACK guarantees `alpha` is
#' bounded by (and usually comparable with) `norm(A)` and `beta` by
#' `norm(B)`, so the norm-scaled policy is invariant under joint rescaling
#' `(A, B) -> (c A, c B)` and does not misclassify well-defined eigenvalues
#' of uniformly small pencils as infinite.
generalized_pencil_values <- function(alpha, beta,
                                      tol = sqrt(.Machine$double.eps),
                                      norm_A = NULL, norm_B = NULL) {
  alpha <- as.vector(alpha)
  beta <- as.vector(beta)
  if (length(alpha) != length(beta)) {
    stop("alpha and beta must have the same length.", call. = FALSE)
  }
  if (length(tol) != 1L || is.na(tol) || tol < 0) {
    stop("tol must be a single non-negative number.", call. = FALSE)
  }
  if (any(!is.finite(Mod(alpha))) || any(!is.finite(Mod(beta)))) {
    stop("alpha and beta must contain finite homogeneous coordinates.", call. = FALSE)
  }

  norm_scaled <- !is.null(norm_A) && !is.null(norm_B) &&
    is.finite(norm_A) && is.finite(norm_B) && norm_A >= 0 && norm_B >= 0
  if (norm_scaled) {
    policy <- "pencil_norm_scaled"
    alpha_threshold <- tol * norm_A
    beta_threshold <- tol * norm_B
    alpha_zero <- Mod(alpha) <= alpha_threshold
    beta_zero <- Mod(beta) <= beta_threshold
  } else {
    policy <- "per_pair_magnitude"
    scale <- pmax(1, Mod(alpha), Mod(beta))
    alpha_threshold <- tol * scale
    beta_threshold <- alpha_threshold
    alpha_zero <- Mod(alpha) <= alpha_threshold
    beta_zero <- Mod(beta) <= beta_threshold
  }
  finite <- !beta_zero
  infinite <- beta_zero & !alpha_zero
  undefined <- beta_zero & alpha_zero
  classification <- rep("finite", length(alpha))
  classification[infinite] <- "infinite"
  classification[undefined] <- "undefined"

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
      tolerance = tol,
      tolerance_policy = policy,
      alpha_threshold = alpha_threshold,
      beta_threshold = beta_threshold,
      norm_A = if (norm_scaled) norm_A else NULL,
      norm_B = if (norm_scaled) norm_B else NULL
    ),
    class = "eigencore_generalized_pencil_values"
  )
}

#' @keywords internal
certify_dense_generalized_pencil <- function(A, B, alpha, beta, vectors,
                                             tol = 1e-8,
                                             beta_tol = sqrt(.Machine$double.eps),
                                             norm_A = NULL, norm_B = NULL) {
  A <- as.matrix(A)
  B <- as.matrix(B)
  vectors <- as.matrix(vectors)
  pencil <- generalized_pencil_values(alpha, beta, tol = beta_tol,
                                      norm_A = norm_A, norm_B = norm_B)
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
certify_dense_generalized_pencil_left <- function(A, B, pencil,
                                                   left_vectors,
                                                   right_vectors,
                                                   tol = 1e-8) {
  A <- as.matrix(A)
  B <- as.matrix(B)
  left_vectors <- as.matrix(left_vectors)
  right_vectors <- as.matrix(right_vectors)
  k <- length(pencil$values)
  generalized_pencil_validate_dimensions(A, B, left_vectors, k)
  generalized_pencil_validate_dimensions(A, B, right_vectors, k)

  # LAPACK normalizes left and right vectors independently. Biorthonormalize
  # each finite left eigenspace against its paired B V block. Treating a
  # repeated eigenvalue as a block matters because LAPACK may choose different
  # bases for its left and right eigenspaces. The transformation preserves the
  # eigenspace and makes the generalized biorthogonality contract inspectable.
  Bv <- B %*% right_vectors
  finite_idx <- which(pencil$finite)
  stable <- rep(FALSE, k)
  if (length(finite_idx)) {
    clusters <- generalized_pencil_value_clusters(
      pencil$values,
      finite_idx
    )
    for (cluster in clusters) {
      W <- left_vectors[, cluster, drop = FALSE]
      Bv_cluster <- Bv[, cluster, drop = FALSE]
      gram <- certificate_gram(W, Bv_cluster)
      singular_values <- tryCatch(
        svd(gram, nu = 0L, nv = 0L)$d,
        error = function(e) numeric()
      )
      pairing_scale <- matrix_norm(W) * matrix_norm(Bv_cluster)
      reciprocal_condition <- tryCatch(
        rcond(gram),
        error = function(e) 0
      )
      pairing_is_resolvable <- length(singular_values) &&
        all(is.finite(singular_values)) &&
        min(singular_values) > 100 * .Machine$double.eps *
          pmax(pairing_scale, .Machine$double.xmin)
      transform <- if (pairing_is_resolvable &&
                       is.finite(reciprocal_condition) &&
                       reciprocal_condition > 100 * .Machine$double.eps) {
        tryCatch(Conj(t(solve(gram))), error = function(e) NULL)
      } else {
        NULL
      }
      if (!is.null(transform) && all(is.finite(Mod(transform)))) {
        left_vectors[, cluster] <- W %*% transform
        stable[cluster] <- TRUE
      }
    }
  }

  left_residuals <- rep(Inf, k)
  scale <- rep(Inf, k)
  backward <- rep(Inf, k)
  if (length(finite_idx)) {
    W <- left_vectors[, finite_idx, drop = FALSE]
    residual_matrix <- Conj(t(A)) %*% W - sweep(
      Conj(t(B)) %*% W,
      2L,
      Conj(pencil$values[finite_idx]),
      `*`
    )
    left_residuals[finite_idx] <- col_norms(residual_matrix)
    scale[finite_idx] <- eigen_backward_scale(
      matrix_norm(A),
      matrix_norm(B),
      pencil$values[finite_idx],
      W
    )
    backward[finite_idx] <- left_residuals[finite_idx] /
      pmax(scale[finite_idx], .Machine$double.eps)
  }

  biorthogonality <- certificate_gram(left_vectors, Bv)
  biorthogonality_loss <- numeric()
  if (length(finite_idx)) {
    finite_cross <- biorthogonality[finite_idx, finite_idx, drop = FALSE]
    biorthogonality_loss <- max(
      abs(finite_cross - diag(length(finite_idx)))
    )
  }
  notes <- c(
    generalized_pencil_certificate_notes(pencil),
    paste(
      "left residual uses A^H w - conj(lambda) B^H w in original",
      "coordinates; biorthogonality uses W^H B V"
    )
  )
  if (any(pencil$finite & !stable)) {
    notes <- c(
      notes,
      paste(
        "one or more finite left/right eigenspaces could not be stably",
        "biorthonormalized"
      )
    )
  }
  cert <- new_certificate(
    tol = tol,
    residuals = list(left = left_residuals),
    backward_error = backward,
    orthogonality = biorthogonality_loss,
    converged = pencil$finite & stable & backward <= tol,
    scale = scale,
    notes = notes,
    certificate_type = "generalized_pencil_left_residual_biorthogonal_backward_error",
    norm_bound_type = "frobenius_exact+frobenius_exact",
    scale_is_estimate = FALSE,
    require_orthogonality = TRUE
  )

  list(
    vectors = left_vectors,
    certificate = cert,
    biorthogonality = biorthogonality
  )
}

#' @keywords internal
generalized_pencil_value_clusters <- function(values, finite_idx) {
  clusters <- list()
  for (idx in finite_idx) {
    matched <- FALSE
    for (cluster_id in seq_along(clusters)) {
      representative <- clusters[[cluster_id]][[1L]]
      scale <- max(1, Mod(values[[idx]]), Mod(values[[representative]]))
      if (Mod(values[[idx]] - values[[representative]]) <=
          100 * .Machine$double.eps * scale) {
        clusters[[cluster_id]] <- c(clusters[[cluster_id]], idx)
        matched <- TRUE
        break
      }
    }
    if (!matched) {
      clusters[[length(clusters) + 1L]] <- idx
    }
  }
  clusters
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
    if (inherits(source, "sparseMatrix") && !is.complex(source)) {
      # Supported Matrix releases differ in whether a real sparse matrix can
      # multiply a complex dense block directly. Split the block so the source
      # remains sparse without requiring a zgeMatrix representation.
      real_part <- as.matrix(source %*% Re(vectors))
      imaginary_part <- as.matrix(source %*% Im(vectors))
      return(real_part + 1i * imaginary_part)
    }
    return(as.matrix(source %*% vectors))
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
