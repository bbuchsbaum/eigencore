#' Compute a dense generalized singular value decomposition
#'
#' `generalized_svd()` is eigencore's dense GSVD compatibility surface for
#' matrix pairs with the same number of columns. The current native path uses
#' the real LAPACK GSVD routine available through R. Sparse/operator inputs are
#' not silently densified, and complex GSVD remains explicit future scope until
#' a bundled or platform `ZGGSVD3`-equivalent native path is available.
#'
#' @param A Base dense real matrix with `m` rows and `n` columns.
#' @param B Base dense real matrix with `p` rows and `n` columns.
#' @param tol Reconstruction and orthogonality certification tolerance.
#' @param ... Reserved for future options.
#' @return An `eigencore_gsvd_result` with fields `alpha`, `beta`, `values`,
#'   `classification`, `U`, `V`, `Q`, `D1`, `D2`, `R`, `zero_R`, `A_factor`,
#'   `B_factor`, `k`, `l`, `rank`, `method`, `plan`, and `certificate`.
#' @examples
#' A <- matrix(c(1, 2, 3, 3, 2, 1), nrow = 2, byrow = TRUE)
#' B <- matrix(1:9, nrow = 3)
#' fit <- generalized_svd(A, B)
#' alpha_beta(fit)$values
#' certificate(fit)$passed
#' @export
generalized_svd <- function(A, B, tol = 1e-8, ...) {
  dots <- list(...)
  if (length(dots)) {
    stop(
      "unused generalized_svd() options: ",
      paste(names(dots), collapse = ", "),
      call. = FALSE
    )
  }
  A <- generalized_svd_dense_input(A, "A")
  B <- generalized_svd_dense_input(B, "B")
  generalized_svd_validate_dimensions(A, B)
  generalized_svd_validate_tol(tol)

  if (is.complex(A) || is.complex(B)) {
    stop(
      "native complex GSVD requires ZGGSVD3 or equivalent; this R LAPACK ",
      "does not export a complex GSVD routine, so complex generalized_svd() ",
      "is not promoted yet",
      call. = FALSE
    )
  }

  raw <- native_dense_generalized_svd(A, B)
  generalized_svd_result(raw, A = A, B = B, tol = tol)
}

#' @keywords internal
generalized_svd_result <- function(raw, A, B, tol) {
  k <- as.integer(raw$k)
  l <- as.integer(raw$l)
  rank <- k + l
  pencil <- generalized_pencil_values(raw$alpha, raw$beta)
  factors <- generalized_svd_factors(raw$A_factor, raw$B_factor, raw$alpha,
                                     raw$beta, k = k, l = l)
  method <- native_dense_generalized_svd_label()
  plan <- list(
    problem_type = "generalized_svd",
    requested = ncol(A),
    method = method,
    target = "all",
    reasons = c("full dense real generalized SVD requested"),
    fallback = "none",
    controls = list(
      full = TRUE,
      dense = TRUE,
      gsvd = TRUE,
      real = TRUE,
      sparse_densified = FALSE
    )
  )
  class(plan) <- "eigencore_plan"
  certificate <- certify_generalized_svd(
    A = A,
    B = B,
    U = raw$U,
    V = raw$V,
    Q = raw$Q,
    D1 = factors$D1,
    D2 = factors$D2,
    zero_R = factors$zero_R,
    tol = tol
  )
  out <- list(
    alpha = raw$alpha,
    beta = raw$beta,
    values = pencil$values,
    classification = pencil$classification,
    finite = pencil$finite,
    infinite = pencil$infinite,
    undefined = pencil$undefined,
    U = raw$U,
    V = raw$V,
    Q = raw$Q,
    D1 = factors$D1,
    D2 = factors$D2,
    R = factors$R,
    zero_R = factors$zero_R,
    A_factor = raw$A_factor,
    B_factor = raw$B_factor,
    k = k,
    l = l,
    rank = rank,
    dimensions = c(m = nrow(A), n = ncol(A), p = nrow(B)),
    residuals = certificate$residuals,
    backward_error = certificate$backward_error,
    orthogonality = certificate$orthogonality,
    nconv = rank,
    requested = ncol(A),
    iterations = 1L,
    matvecs = 0L,
    method = method,
    target = "all",
    plan = plan,
    certificate = certificate,
    warnings = "using native dense real LAPACK GSVD full decomposition"
  )
  class(out) <- "eigencore_gsvd_result"
  out
}

#' @keywords internal
generalized_svd_factors <- function(A_factor, B_factor, alpha, beta, k, l) {
  m <- nrow(A_factor)
  n <- ncol(A_factor)
  p <- nrow(B_factor)
  rank <- k + l
  R <- generalized_svd_R(A_factor, B_factor, k = k, l = l)
  zero_R <- if (n > rank) {
    cbind(matrix(0, nrow = rank, ncol = n - rank), R)
  } else {
    R
  }
  list(
    R = R,
    zero_R = zero_R,
    D1 = generalized_svd_D1(m, k = k, l = l, alpha = alpha),
    D2 = generalized_svd_D2(p, m, k = k, l = l, beta = beta)
  )
}

#' @keywords internal
generalized_svd_R <- function(A_factor, B_factor, k, l) {
  m <- nrow(A_factor)
  n <- ncol(A_factor)
  rank <- k + l
  if (rank == 0L) {
    return(matrix(0, nrow = 0L, ncol = 0L))
  }
  if (m - rank < 0L) {
    lower_rows <- seq.int(m - k + 1L, l)
    combined <- rbind(A_factor, B_factor[lower_rows, , drop = FALSE])
    R <- combined[seq_len(rank), , drop = FALSE]
    if (n > rank) {
      R <- R[, seq.int(n - rank + 1L, n), drop = FALSE]
    }
    return(R)
  }
  A_factor[seq_len(rank), seq.int(n - rank + 1L, n), drop = FALSE]
}

#' @keywords internal
generalized_svd_D1 <- function(m, k, l, alpha) {
  rank <- k + l
  D1 <- matrix(0, nrow = m, ncol = rank)
  if (k > 0L) {
    D1[seq_len(k), seq_len(k)] <- diag(1, k)
  }
  if (m - rank >= 0L) {
    if (l > 0L) {
      rows <- seq.int(k + 1L, rank)
      D1[rows, rows] <- diag(alpha[rows], nrow = l, ncol = l)
    }
  } else if (m > k) {
    rows <- seq.int(k + 1L, m)
    width <- length(rows)
    D1[rows, rows] <- diag(alpha[rows], nrow = width, ncol = width)
  }
  D1
}

#' @keywords internal
generalized_svd_D2 <- function(p, m, k, l, beta) {
  rank <- k + l
  D2 <- matrix(0, nrow = p, ncol = rank)
  if (rank == 0L) {
    return(D2)
  }
  if (m - rank >= 0L) {
    if (l > 0L) {
      rows <- seq_len(l)
      cols <- seq.int(k + 1L, rank)
      D2[rows, cols] <- diag(beta[cols], nrow = l, ncol = l)
    }
  } else {
    if (m > k) {
      rows <- seq_len(m - k)
      cols <- seq.int(k + 1L, m)
      width <- length(cols)
      D2[rows, cols] <- diag(beta[cols], nrow = width, ncol = width)
    }
    if (rank > m) {
      rows <- seq.int(m - k + 1L, l)
      cols <- seq.int(m + 1L, rank)
      width <- length(cols)
      D2[rows, cols] <- diag(1, width)
    }
  }
  D2
}

#' @keywords internal
certify_generalized_svd <- function(A, B, U, V, Q, D1, D2, zero_R, tol) {
  Q_adj <- t(Q)
  A_hat <- U %*% D1 %*% zero_R %*% Q_adj
  B_hat <- V %*% D2 %*% zero_R %*% Q_adj
  residuals <- c(
    A = matrix_norm(A - A_hat),
    B = matrix_norm(B - B_hat)
  )
  scales <- c(A = max(1, matrix_norm(A)), B = max(1, matrix_norm(B)))
  backward_error <- residuals / scales
  orthogonality <- c(
    U = matrix_norm(crossprod(U) - diag(ncol(U))),
    V = matrix_norm(crossprod(V) - diag(ncol(V))),
    Q = matrix_norm(crossprod(Q) - diag(ncol(Q)))
  )
  new_certificate(
    tol = tol,
    residuals = residuals,
    backward_error = backward_error,
    orthogonality = orthogonality,
    converged = c(backward_error <= tol, orthogonality <= max(tol, sqrt(.Machine$double.eps))),
    scale = scales,
    notes = "GSVD reconstruction and orthogonality certificate",
    certificate_type = "gsvd_reconstruction",
    norm_bound_type = "frobenius_exact+frobenius_exact",
    require_orthogonality = TRUE
  )
}

#' @keywords internal
generalized_svd_dense_input <- function(x, name) {
  if (!is.matrix(x)) {
    stop(
      name, " must be a base dense matrix for generalized_svd(); ",
      "sparse/operator GSVD inputs are not silently densified",
      call. = FALSE
    )
  }
  if (any(!is.finite(x))) {
    stop(name, " must contain only finite values.", call. = FALSE)
  }
  if (any(dim(x) <= 0L)) {
    stop(name, " must have positive dimensions.", call. = FALSE)
  }
  if (is.complex(x)) {
    storage.mode(x) <- "complex"
  } else {
    storage.mode(x) <- "double"
  }
  x
}

#' @keywords internal
generalized_svd_validate_dimensions <- function(A, B) {
  if (ncol(A) != ncol(B)) {
    stop("A and B must have the same number of columns.", call. = FALSE)
  }
  invisible(TRUE)
}

#' @keywords internal
generalized_svd_validate_tol <- function(tol) {
  if (!is.numeric(tol) || length(tol) != 1L || is.na(tol) || tol < 0) {
    stop("tol must be a single non-negative number.", call. = FALSE)
  }
  invisible(TRUE)
}

#' @keywords internal
native_dense_generalized_svd <- function(A, B) {
  .Call(
    "eigencore_dense_generalized_svd",
    as.matrix(A),
    as.matrix(B),
    PACKAGE = "eigencore"
  )
}

#' @keywords internal
native_dense_generalized_svd_label <- function() {
  "native dense real LAPACK GSVD full"
}
