# Native route for both_ends(k_low, k_high) on standard real Hermitian
# operators (dense, dgCMatrix, matrix-free callbacks).
#
# The native thick-restart Lanczos kernels select one end of the spectrum,
# so both_ends() used to fall back to the unrestarted reference Lanczos, which
# does not converge on spectra whose ends are not well separated (O10). This
# route runs two native thick-restart solves (k_low smallest, k_high largest),
# merges them by a Rayleigh-Ritz step on the span of both blocks (which
# restores cross-orthogonality between the halves) and certifies the merged
# set from scratch. Missing copies of a repeated eigenvalue in either half are
# caught by the target-completeness check (two parts, R/completeness_hermitian.R).
# A warm start (initial_subspace, e.g. eigs_sym(opts = list(initvec =))) is
# the start block of both solves.

#' @keywords internal
native_both_ends_lanczos_label <- function() {
  "native Hermitian Lanczos both ends (two thick-restart solves)"
}

# Can the both-ends route serve this problem?
#' @keywords internal
native_both_ends_lanczos_supported <- function(problem) {
  op <- problem$A
  target <- problem$target
  if (!inherits(target, "eigencore_target") || !identical(target$kind, "both_ends") ||
      !is.null(problem$metric) || !inherits(op, "eigencore_operator") ||
      !identical(op$structure$kind, "hermitian") ||
      !identical(op$dtype %||% "double", "double")) {
    return(FALSE)
  }
  source <- source_or_null(op)
  (is.matrix(source) && is.double(source)) ||
    identical(op$metadata$storage %||% NULL, "dgCMatrix") ||
    native_matrix_free_block_lanczos_available(op)
}

#' @keywords internal
native_both_ends_lanczos_hermitian <- function(op, k, target, tol = 1e-8,
                                               maxit = NULL, block = 1L,
                                               max_restarts = 100L,
                                               vectors = TRUE,
                                               check_stride = 0L,
                                               start = NULL) {
  op <- as_operator(op)
  n <- as.integer(op$dim[[1L]])
  kl <- as.integer(target$value$k_low)
  kh <- as.integer(target$value$k_high)
  block <- max(1L, as.integer(block %||% 1L))
  max_restarts <- as.integer(max_restarts %||% 100L)
  run <- function(kk, end_target) {
    if (kk < 1L) {
      return(NULL)
    }
    m <- if (is.null(maxit)) {
      if (block > 1L) {
        min(n, default_block_lanczos_max_subspace(kk, block))
      } else {
        default_lanczos_max_subspace(kk, n, op = op)
      }
    } else {
      min(n, max(as.integer(maxit), kk + block))
    }
    if (block > 1L) {
      native_block_lanczos_hermitian(
        op, k = kk, target = end_target, tol = tol, maxit = m, block = block,
        max_restarts = max_restarts, vectors = TRUE,
        # A warm start must be iterated, not bypassed by the dense
        # full-subspace shortcut (as on the single-end route).
        full_subspace = is.null(start), start = start,
        check_stride = check_stride
      )
    } else {
      native_lanczos_hermitian(
        op, k = kk, target = end_target, tol = tol, maxit = m,
        max_restarts = max_restarts, vectors = TRUE, start = start,
        check_stride = check_stride
      )
    }
  }
  low <- run(kl, smallest())
  high <- run(kh, largest())
  halves <- Filter(Negate(is.null), list(low, high))
  count <- function(field) {
    sum(vapply(halves, function(h) {
      as.numeric(h[[field]] %||% h$matvecs %||% 0)
    }, numeric(1L)))
  }
  V <- do.call(cbind, lapply(halves, function(h) as.matrix(h$vectors)))
  # Rayleigh-Ritz on span(V_low, V_high): restores orthogonality between the
  # halves (their eigenvalues differ unless the ends overlap) and keeps the
  # converged values.
  Q <- qr.Q(qr(V))
  AQ <- as.matrix(apply_operator(op, Q))
  H <- crossprod(Q, AQ)
  H <- (H + t(H)) / 2
  eig <- eigen(H, symmetric = TRUE)
  ord <- order_indices(eig$values, target)
  values <- eig$values[ord]
  X <- Q %*% eig$vectors[, ord, drop = FALSE]
  cert <- certify_eigen_operator(op, values, X, tol = tol)
  rr_columns <- ncol(Q)
  list(
    values = values,
    vectors = if (isTRUE(vectors)) X else NULL,
    residuals = cert$residuals,
    backward_error = cert$backward_error,
    orthogonality = cert$orthogonality,
    certificate = cert,
    iterations = as.integer(count("iterations")),
    matvecs = as.integer(count("matvecs")),
    restarts = as.integer(count("restarts")),
    operator_block_calls = as.integer(count("operator_block_calls") + 2L),
    operator_columns = as.integer(count("operator_columns") + rr_columns + k),
    certification_operator_columns = as.integer(
      count("certification_operator_columns") + k
    ),
    locked = which(cert$converged),
    block = block,
    restart = list(
      kind = "both_ends_two_thick_restart_solves",
      implemented = TRUE,
      native = TRUE,
      block = block,
      max_restarts = max_restarts,
      max_subspace = max(vapply(halves, function(h) {
        as.numeric(h$restart$max_subspace %||% NA_real_)
      }, numeric(1L))),
      restarts_used = as.integer(count("restarts")),
      ends = c(k_low = kl, k_high = kh),
      merge = "rayleigh_ritz_on_both_halves"
    )
  )
}
