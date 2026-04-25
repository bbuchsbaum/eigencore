#' RSpectra-compatible eigen shim.
eigs <- function(A, k, which = "LM", opts = list(), ...) {
  target <- target_from_which(which, k = k)
  fit <- eig_partial(A, k = k, target = target, ...)
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

#' RSpectra-compatible symmetric eigen shim.
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
