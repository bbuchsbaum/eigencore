#' RSpectra-compatible eigen shim.
#'
#' @param A Matrix or eigencore operator.
#' @param k Number of eigenpairs to compute.
#' @param which RSpectra-style target selector.
#' @param opts Compatibility options list; currently accepted for API
#'   compatibility and not interpreted directly.
#' @param ... Additional arguments passed to [eig_partial()].
#' @return A list compatible with `RSpectra::eigs()`, including `values`,
#'   `vectors`, convergence counts, operation counts, certificate diagnostics,
#'   and left/right vector fields when available.
#' @examples
#' A <- diag(c(5, 4, 3, 2, 1))
#' A[1, 2] <- 0.5
#' res <- eigs(A, k = 2, which = "LM")
#' res$values
eigs <- function(A, k, which = "LM", opts = list(), ...) {
  target <- target_from_which(which, k = k)
  fit <- eig_partial(A, k = k, target = target, ...)
  list(
    values = fit$values,
    vectors = fit$vectors,
    left_vectors = left_vectors(fit),
    right_vectors = right_vectors(fit),
    nconv = fit$nconv,
    niter = fit$iterations,
    nops = fit$matvecs,
    left_certificate = fit$left_certificate,
    biorthogonality = fit$biorthogonality,
    certificate = fit$certificate,
    diagnostics = diagnostics(fit)
  )
}

#' RSpectra-compatible symmetric eigen shim.
#'
#' @param A Matrix or eigencore operator.
#' @param k Number of eigenpairs to compute.
#' @param which RSpectra-style target selector.
#' @param opts Compatibility options list; currently accepted for API
#'   compatibility and not interpreted directly.
#' @param ... Additional arguments passed to [solve.eigencore_eigen_problem()].
#' @return A list compatible with `RSpectra::eigs_sym()`, including `values`,
#'   `vectors`, convergence counts, operation counts, certificate diagnostics,
#'   and eigencore diagnostics.
#' @examples
#' A <- diag(c(5, 4, 3, 2, 1))
#' res <- eigs_sym(A, k = 2, which = "LA")
#' res$values
eigs_sym <- function(A, k, which = "LA", opts = list(), ...) {
  target <- target_from_which(which, k = k)
  P <- eigen_problem(A, structure = hermitian(), target = target)
  fit <- solve(P, k = k, ...)
  list(
    values = fit$values,
    vectors = fit$vectors,
    nconv = fit$nconv,
    niter = fit$iterations,
    nops = fit$matvecs,
    certificate = fit$certificate,
    diagnostics = diagnostics(fit)
  )
}

#' RSpectra-compatible SVD shim.
#'
#' @param A Matrix or eigencore operator.
#' @param k Number of singular values to compute.
#' @param nu Number of left singular vectors requested.
#' @param nv Number of right singular vectors requested.
#' @param opts Compatibility options list; currently accepted for API
#'   compatibility and not interpreted directly.
#' @param ... Additional arguments passed to [svd_partial()].
#' @return A list compatible with `RSpectra::svds()`, including `d`, optional
#'   `u` and `v`, convergence counts, operation counts, certificate
#'   diagnostics, and eigencore diagnostics.
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(60), 10, 6)
#' res <- svds(X, k = 2)
#' res$d
svds <- function(A, k, nu = k, nv = k, opts = list(), ...) {
  vector_mode <- if (nu > 0 && nv > 0) {
    "both"
  } else if (nu > 0) {
    "left"
  } else if (nv > 0) {
    "right"
  } else {
    "none"
  }
  fit <- svd_partial(A, rank = k, vectors = vector_mode, ...)
  list(
    d = fit$d,
    u = fit$u,
    v = fit$v,
    nconv = fit$nconv,
    niter = fit$iterations,
    nops = fit$matvecs,
    certificate = fit$certificate,
    diagnostics = diagnostics(fit)
  )
}

#' @keywords internal
target_from_which <- function(which, k = NULL) {
  both_ends_from_k <- function(k) {
    k <- as.integer(k)
    if (length(k) != 1L || is.na(k) || k < 1L) {
      stop("ARPACK which = 'BE' requires a positive k.", call. = FALSE)
    }
    k_low <- k %/% 2L
    k_high <- k - k_low
    both_ends(k_low, k_high)
  }
  switch(
    toupper(which),
    LM = largest_magnitude(),
    SM = smallest_magnitude(),
    LA = largest(),
    SA = smallest(),
    LR = largest_real(),
    SR = smallest_real(),
    LI = largest_imaginary(),
    SI = smallest_imaginary(),
    BE = both_ends_from_k(k),
    largest()
  )
}
