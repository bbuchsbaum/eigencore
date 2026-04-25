#' Compute a partial eigendecomposition.
eig_partial <- function(A, k, target = largest(), B = NULL, method = auto(),
                        tol = 1e-8, maxit = NULL, vectors = TRUE, seed = NULL,
                        certify = TRUE) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  P <- eigen_problem(A, metric = B, target = target, transform = if (inherits(method, "eigencore_method") && method$kind == "shift_invert") method else NULL)
  solve(P, k = k, method = method, tol = tol, maxit = maxit, vectors = vectors, certify = certify)
}

#' Compute a partial singular-value decomposition.
svd_partial <- function(A, rank, target = largest(), method = auto(), tol = 1e-8,
                        vectors = c("both", "left", "right", "none"),
                        seed = NULL, certify = TRUE) {
  vectors <- match.arg(vectors)
  if (!is.null(seed)) {
    set.seed(seed)
  }
  P <- svd_problem(A, target = target)
  solve(P, rank = rank, method = method, tol = tol, vectors = vectors, certify = certify)
}

#' @export
solve.eigencore_eigen_problem <- function(a, b, k, method = auto(), tol = 1e-8,
                                          maxit = NULL, vectors = TRUE,
                                          certify = TRUE, ...) {
  plan <- plan_solver(a, k = k, method = method)
  if (inherits(a$transform, "eigencore_method") && identical(a$transform$kind, "shift_invert")) {
    stop("shift_invert() is declared in the API but not implemented in the solver core yet.", call. = FALSE)
  }
  use_lanczos <- should_use_lanczos(a, method)
  if (use_lanczos) {
    method_maxit <- if (inherits(method, "eigencore_method") && identical(method$kind, "lanczos")) method$max_subspace else NULL
    method_reorth <- if (inherits(method, "eigencore_method") && identical(method$kind, "lanczos")) method$reorthogonalize else TRUE
    iter <- reference_lanczos_hermitian(
      a$A,
      k = k,
      target = a$target,
      tol = tol,
      maxit = maxit %||% method_maxit,
      vectors = vectors,
      reorthogonalize = method_reorth
    )

    result <- list(
      values = iter$values,
      vectors = iter$vectors,
      residuals = iter$residuals,
      backward_error = iter$backward_error,
      orthogonality = iter$orthogonality,
      nconv = sum(iter$certificate$converged),
      requested = k,
      iterations = iter$iterations,
      matvecs = iter$matvecs,
      method = plan$method,
      target = target_label(a$target),
      plan = plan,
      certificate = iter$certificate,
      warnings = if (identical(plan$method, "native CSC-backed prototype Hermitian Lanczos")) {
        "using R-level prototype Lanczos over native CSC block apply; native hot loop not yet implemented"
      } else {
        "using R-level prototype Lanczos; native hot loop not yet implemented"
      }
    )
    class(result) <- "eigencore_eigen_result"
    return(result)
  }

  if (should_use_native_dense_hermitian(a, method)) {
    A <- operator_source_matrix(a$A)
    eig <- native_dense_symmetric_eigen(A)
    idx <- order_indices(eig$values, a$target)
    idx <- idx[seq_len(min(k, length(idx)))]
    vals <- eig$values[idx]
    vecs <- if (vectors) eig$vectors[, idx, drop = FALSE] else NULL
    cert <- if (certify && !is.null(vecs)) {
      certify_eigen(A, vals, vecs, tol = tol)
    } else {
      empty_certificate(tol, note = "vectors not returned; residual certificate not computed")
    }

    result <- list(
      values = vals,
      vectors = vecs,
      residuals = cert$residuals,
      backward_error = cert$backward_error,
      orthogonality = cert$orthogonality,
      nconv = sum(cert$converged),
      requested = k,
      iterations = 1L,
      matvecs = 0L,
      method = plan$method,
      target = target_label(a$target),
      plan = plan,
      certificate = cert,
      warnings = "using native dense Hermitian LAPACK fallback; iterative engine not yet implemented"
    )
    class(result) <- "eigencore_eigen_result"
    return(result)
  }

  A <- operator_source_matrix(a$A)
  B <- if (is.null(a$metric)) NULL else operator_source_matrix(a$metric)

  if (!is.null(B)) {
    eig <- dense_generalized_spd_eigen(A, B, vectors = vectors)
  } else {
    eig <- eigen(A, symmetric = identical(a$structure$kind, "hermitian"))
    if (!vectors) {
      eig$vectors <- NULL
    }
  }

  idx <- order_indices(eig$values, a$target)
  idx <- idx[seq_len(min(k, length(idx)))]
  vals <- eig$values[idx]
  vecs <- if (vectors && !is.null(eig$vectors)) eig$vectors[, idx, drop = FALSE] else NULL

  cert <- if (certify && !is.null(vecs)) {
    certify_eigen(A, vals, vecs, B = B, tol = tol)
  } else {
    empty_certificate(tol, note = "vectors not returned; residual certificate not computed")
  }

  result <- list(
    values = vals,
    vectors = vecs,
    residuals = cert$residuals,
    backward_error = cert$backward_error,
    orthogonality = cert$orthogonality,
    nconv = sum(cert$converged),
    requested = k,
    iterations = 1L,
    matvecs = 0L,
    method = plan$method,
    target = target_label(a$target),
    plan = plan,
    certificate = cert,
    warnings = if (identical(plan$method, "native dense generalized SPD LAPACK fallback")) {
      "using native dense generalized SPD LAPACK fallback; iterative engine not yet implemented"
    } else if (identical(plan$fallback, "dense oracle prototype")) {
      "using dense oracle prototype solver"
    } else {
      character()
    }
  )
  class(result) <- "eigencore_eigen_result"
  result
}

#' @export
solve.eigencore_svd_problem <- function(a, b, rank, method = auto(), tol = 1e-8,
                                        vectors = c("both", "left", "right", "none"),
                                        certify = TRUE, ...) {
  vectors <- match.arg(vectors)
  plan <- plan_solver(a, rank = rank, method = method)
  use_gk <- should_use_golub_kahan(a, method)
  if (use_gk) {
    method_maxit <- if (inherits(method, "eigencore_method") && identical(method$kind, "golub_kahan")) method$max_subspace else NULL
    method_reorth <- if (inherits(method, "eigencore_method") && identical(method$kind, "golub_kahan")) method$reorthogonalize else TRUE
    iter <- reference_golub_kahan_svd(
      a$A,
      rank = rank,
      target = a$target,
      tol = tol,
      maxit = method_maxit,
      vectors = vectors,
      reorthogonalize = method_reorth
    )

    cert <- if (isTRUE(certify) && !is.null(iter$u) && !is.null(iter$v)) {
      iter$certificate
    } else {
      empty_certificate(tol, note = "both left and right vectors are required for full SVD certification")
    }
    result <- list(
      d = iter$d,
      u = iter$u,
      v = iter$v,
      values = iter$d,
      residuals = cert$residuals,
      backward_error = cert$backward_error,
      orthogonality = cert$orthogonality,
      nconv = sum(cert$converged),
      requested = rank,
      iterations = iter$iterations,
      matvecs = iter$matvecs,
      method = plan$method,
      target = target_label(a$target),
      plan = plan,
      certificate = cert,
      warnings = if (identical(plan$method, "native CSC-backed prototype Golub-Kahan")) {
        "using R-level prototype Golub-Kahan over native CSC block apply; native hot loop not yet implemented"
      } else {
        "using R-level prototype Golub-Kahan; native hot loop not yet implemented"
      }
    )
    class(result) <- "eigencore_svd_result"
    return(result)
  }

  A <- operator_source_matrix(a$A)
  decomp <- if (should_use_native_dense_svd(a, method)) {
    native_dense_svd(A)
  } else {
    nu <- if (vectors %in% c("both", "left")) min(rank, nrow(A)) else 0L
    nv <- if (vectors %in% c("both", "right")) min(rank, ncol(A)) else 0L
    svd(A, nu = nu, nv = nv)
  }
  idx <- order_indices(decomp$d, a$target)
  idx <- idx[seq_len(min(rank, length(idx)))]

  d <- decomp$d[idx]
  u <- if (vectors %in% c("both", "left")) decomp$u[, idx, drop = FALSE] else NULL
  v <- if (vectors %in% c("both", "right")) decomp$v[, idx, drop = FALSE] else NULL

  cert <- if (certify && !is.null(u) && !is.null(v)) {
    certify_svd(A, d, u, v, tol = tol)
  } else {
    empty_certificate(tol, note = "both left and right vectors are required for full SVD certification")
  }

  result <- list(
    d = d,
    u = u,
    v = v,
    values = d,
    residuals = cert$residuals,
    backward_error = cert$backward_error,
    orthogonality = cert$orthogonality,
    nconv = sum(cert$converged),
    requested = rank,
    iterations = 1L,
    matvecs = 0L,
    method = plan$method,
    target = target_label(a$target),
    plan = plan,
    certificate = cert,
    warnings = if (should_use_native_dense_svd(a, method)) {
      "using native dense LAPACK SVD fallback; iterative engine not yet implemented"
    } else {
      "using dense oracle prototype solver"
    }
  )
  class(result) <- "eigencore_svd_result"
  result
}

#' @keywords internal
dense_generalized_spd_eigen <- function(A, B, vectors = TRUE) {
  eig <- native_dense_generalized_spd_eigen(A, B)
  if (!vectors) eig$vectors <- NULL
  eig
}

#' @keywords internal
native_dense_generalized_spd_eigen <- function(A, B) {
  .Call(
    "eigencore_dense_generalized_spd_eigen",
    as.matrix(A),
    as.matrix(B),
    PACKAGE = "eigencore"
  )
}

#' @keywords internal
order_indices <- function(x, target) {
  kind <- if (inherits(target, "eigencore_target")) target$kind else "largest"
  switch(
    kind,
    largest = order(Re(x), decreasing = TRUE),
    smallest = order(Re(x), decreasing = FALSE),
    largest_magnitude = order(Mod(x), decreasing = TRUE),
    smallest_magnitude = order(Mod(x), decreasing = FALSE),
    largest_real = order(Re(x), decreasing = TRUE),
    smallest_real = order(Re(x), decreasing = FALSE),
    largest_imaginary = order(Im(x), decreasing = TRUE),
    smallest_imaginary = order(Im(x), decreasing = FALSE),
    nearest = order(abs(x - target$value), decreasing = FALSE),
    order(Re(x), decreasing = TRUE)
  )
}

#' @keywords internal
should_use_lanczos <- function(problem, method) {
  if (!is.null(problem$metric)) {
    return(FALSE)
  }
  if (!identical(problem$structure$kind, "hermitian")) {
    return(FALSE)
  }
  if (inherits(method, "eigencore_method") && identical(method$kind, "lanczos")) {
    return(TRUE)
  }
  inherits(method, "eigencore_method") &&
    identical(method$kind, "auto") &&
    is.null(source_or_null(problem$A))
}

#' @keywords internal
should_use_native_dense_hermitian <- function(problem, method) {
  if (!is.null(problem$metric)) {
    return(FALSE)
  }
  if (!identical(problem$structure$kind, "hermitian")) {
    return(FALSE)
  }
  if (inherits(method, "eigencore_method") && !identical(method$kind, "auto")) {
    return(FALSE)
  }
  src <- source_or_null(problem$A)
  is.matrix(src) && is.double(src)
}

#' @keywords internal
native_dense_symmetric_eigen <- function(A) {
  .Call("eigencore_dense_symmetric_eigen", as.matrix(A), PACKAGE = "eigencore")
}

#' @keywords internal
should_use_native_dense_svd <- function(problem, method) {
  if (inherits(method, "eigencore_method") && !identical(method$kind, "auto")) {
    return(FALSE)
  }
  src <- source_or_null(problem$A)
  is.matrix(src) && is.double(src)
}

#' @keywords internal
native_dense_svd <- function(A) {
  .Call("eigencore_dense_svd", as.matrix(A), PACKAGE = "eigencore")
}

#' @keywords internal
should_use_golub_kahan <- function(problem, method) {
  if (is.null(problem$A$apply_adjoint)) {
    return(FALSE)
  }
  if (inherits(method, "eigencore_method") && identical(method$kind, "golub_kahan")) {
    return(TRUE)
  }
  inherits(method, "eigencore_method") &&
    identical(method$kind, "auto") &&
    is.null(source_or_null(problem$A))
}
