# Per-method runtime helpers for solve.eigencore_eigen_problem.
# Each helper takes the eigen problem `a`, the requested k, and the resolved
# `plan` produced by plan_solver(); it returns a classed eigencore_eigen_result.
# The dispatcher in R/solve.R selects between these based on plan$method labels.

#' @keywords internal
solve_eigen_lobpcg <- function(a, k, method, tol, maxit, vectors, certify, plan) {
  if (!identical(a$structure$kind, "hermitian")) {
    stop("lobpcg() prototype currently requires a Hermitian eigenproblem.", call. = FALSE)
  }
  controls <- plan$controls %||% list()
  method_is_lobpcg <- inherits(method, "eigencore_method") && identical(method$kind, "lobpcg")
  method_preconditioner <- if (method_is_lobpcg) method$preconditioner else NULL
  method_constraints <- if (method_is_lobpcg) method$constraints else NULL
  method_maxit <- maxit %||% controls$maxit %||%
    if (method_is_lobpcg) method$maxit else getOption("eigencore.lobpcg_maxit", 200L)
  iter <- if (identical(plan$method, native_standard_lobpcg_label())) {
    native_lobpcg_hermitian(
      a$A,
      k = k,
      target = a$target,
      tol = tol,
      maxit = method_maxit,
      preconditioner = method_preconditioner,
      constraints = method_constraints
    )
  } else if (identical(plan$method, native_generalized_lobpcg_label())) {
    native_generalized_lobpcg_hermitian(
      a$A,
      a$metric,
      k = k,
      target = a$target,
      tol = tol,
      maxit = method_maxit,
      preconditioner = method_preconditioner,
      constraints = method_constraints
    )
  } else {
    reference_lobpcg_hermitian(
      a$A,
      k = k,
      target = a$target,
      tol = tol,
      maxit = method_maxit,
      preconditioner = method_preconditioner,
      Bop = a$metric,
      constraints = method_constraints
    )
  }
  warning_msg <- if (!isTRUE(iter$certificate$passed)) {
    paste0(plan$method, " exhausted ", method_maxit,
           " iterations before all ", k, " requested pairs converged")
  } else if (isTRUE(iter$native)) {
    if (isTRUE(iter$generalized)) {
      character()
    } else {
      "using native standard LOBPCG prototype; locking/generalized production path not yet implemented"
    }
  } else {
    "using R-level reference LOBPCG prototype; native hot loop not yet implemented"
  }
  make_eigen_result(
    values = iter$values,
    vectors = if (isTRUE(vectors)) iter$vectors else NULL,
    certificate = iter$certificate,
    iter = iter,
    requested = k,
    method_label = plan$method,
    target_label_value = target_label(a$target),
    plan = plan,
    warnings = warning_msg,
    extras = list(
      preconditioner_calls = iter$preconditioner_calls,
      convergence_history = iter$convergence_history,
      preconditioner = iter$preconditioner,
      restart = list(
        kind = "lobpcg",
        implemented = TRUE,
        native = isTRUE(iter$native),
        native_kernels = any(c(
          isTRUE(iter$native),
          isTRUE(iter$orthogonalization$native),
          isTRUE(iter$preconditioner$native)
        )),
        orthogonalization = iter$orthogonalization,
        orthogonalization_native = isTRUE(iter$orthogonalization$native),
        orthogonalization_methods = iter$orthogonalization$methods,
        preconditioned = isTRUE(iter$preconditioned),
        preconditioner_kind = iter$preconditioner$kind,
        preconditioner_native = isTRUE(iter$preconditioner$native),
        preconditioner_calls = iter$preconditioner_calls,
        preconditioner = iter$preconditioner,
        constrained = isTRUE(iter$constrained),
        constraints_rank = iter$constraints_rank %||% 0L,
        generalized = isTRUE(iter$generalized),
        maxit = method_maxit,
        q_rank_final = iter$q_rank_final %||% NA_integer_
      ),
      locked = which(iter$certificate$converged)
    )
  )
}

#' @keywords internal
solve_eigen_lanczos <- function(a, k, method, tol, maxit, vectors, certify, plan) {
  controls <- plan$controls %||% list()
  method_maxit <- controls$max_subspace %||%
    if (inherits(method, "eigencore_method") && identical(method$kind, "lanczos")) method$max_subspace else NULL
  method_reorth <- controls$reorthogonalize %||%
    if (inherits(method, "eigencore_method") && identical(method$kind, "lanczos")) method$reorthogonalize else TRUE
  method_max_restarts <- controls$max_restarts %||%
    if (inherits(method, "eigencore_method") && identical(method$kind, "lanczos")) method$max_restarts else NULL
  method_block <- controls$block %||%
    if (inherits(method, "eigencore_method") && identical(method$kind, "lanczos")) method$block else 1L
  iter <- if (plan_dispatches_native_lanczos(plan)) {
    if (method_block > 1L) {
      native_block_lanczos_hermitian(
        a$A,
        k = k,
        target = a$target,
        tol = tol,
        maxit = maxit %||% method_maxit,
        block = method_block,
        max_restarts = method_max_restarts,
        vectors = vectors
      )
    } else {
      native_lanczos_hermitian(
        a$A,
        k = k,
        target = a$target,
        tol = tol,
        maxit = maxit %||% method_maxit,
        max_restarts = method_max_restarts,
        vectors = vectors
      )
    }
  } else {
    reference_lanczos_hermitian(
      a$A,
      k = k,
      target = a$target,
      tol = tol,
      maxit = maxit %||% method_maxit,
      vectors = vectors,
      reorthogonalize = method_reorth
    )
  }

  warning_msg <- if (plan_dispatches_native_lanczos(plan)) {
    if (!isTRUE(iter$certificate$passed)) {
      paste0(
        plan$method, " exhausted its current budget and did not converge all ",
        k,
        " requested pairs; restart budget/subspace: ",
        iter$restart$max_restarts %||% NA_integer_,
        "/",
        iter$restart$max_subspace %||% NA_integer_
      )
    } else if (isTRUE(iter$restart$fallback_used)) {
      paste0(
        "native block Lanczos failed certification; used ",
        iter$restart$kind,
        " for a certified result"
      )
    } else {
      character()
    }
  } else {
    "using R-level prototype Lanczos; native hot loop not yet implemented"
  }

  make_eigen_result(
    values = iter$values,
    vectors = iter$vectors,
    certificate = iter$certificate,
    iter = iter,
    requested = k,
    method_label = plan$method,
    target_label_value = target_label(a$target),
    plan = plan,
    warnings = warning_msg,
    extras = list(
      residuals = iter$residuals %||% iter$certificate$residuals,
      backward_error = iter$backward_error %||% iter$certificate$backward_error,
      orthogonality = iter$orthogonality %||% iter$certificate$orthogonality,
      restarts = iter$restarts %||% iter$restart$restarts_used %||% NA_integer_,
      ortho_passes = iter$ortho_passes %||% iter$restart$ortho_passes %||% NA_integer_,
      locking_events = iter$locking_events %||% iter$restart$locking_events %||% NA_integer_,
      block = iter$block %||% iter$restart$block %||% NA_integer_,
      operator_allocations = iter$operator_allocations %||%
        iter$restart$operator_allocations %||%
        if (identical(iter$restart$kind, "block_full_subspace_dense_lapack")) 0 else NA_real_,
      operator_bytes_allocated = iter$operator_bytes_allocated %||%
        iter$restart$operator_bytes_allocated %||%
        if (identical(iter$restart$kind, "block_full_subspace_dense_lapack")) 0 else NA_real_,
      stage_seconds = iter$stage_seconds %||% iter$restart$stage_seconds %||% numeric(),
      convergence_history = iter$convergence_history %||% NULL,
      restart = iter$restart %||% NULL,
      locked = iter$locked %||% which(iter$certificate$converged)
    )
  )
}

#' @keywords internal
solve_eigen_native_dense_hermitian <- function(a, k, tol, vectors, certify,
                                                allow_dense_fallback, plan) {
  A <- materialize_dense_fallbacks(list(A = a$A), allow = allow_dense_fallback)$A
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
  make_eigen_result(
    values = vals,
    vectors = vecs,
    certificate = cert,
    iter = list(iterations = 1L, matvecs = 0L),
    requested = k,
    method_label = plan$method,
    target_label_value = target_label(a$target),
    plan = plan,
    warnings = "using native dense Hermitian LAPACK fallback; iterative engine not yet implemented"
  )
}

#' @keywords internal
solve_eigen_dense_oracle <- function(a, k, tol, vectors, certify,
                                     allow_dense_fallback, plan) {
  dense_inputs <- materialize_dense_fallbacks(
    if (is.null(a$metric)) list(A = a$A) else list(A = a$A, B = a$metric),
    allow = allow_dense_fallback
  )
  A <- dense_inputs$A
  B <- dense_inputs$B

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

  # The SPD certification gate only certifies real eigenpairs. Non-Hermitian
  # general eigenproblems can return complex eigenpairs from base::eigen(),
  # which certify_eigen has no documented contract for. Skip certification
  # with an explicit note rather than silently passing complex inputs through.
  complex_eigenpairs <- is.complex(vals) ||
    (is.complex(eig$values) && any(abs(Im(vals)) > 0)) ||
    (!is.null(vecs) && is.complex(vecs))
  cert <- if (certify && !is.null(vecs) && !complex_eigenpairs) {
    certify_eigen(
      A,
      vals,
      vecs,
      B = B,
      tol = tol,
      require_orthogonality = identical(a$structure$kind, "hermitian") || !is.null(B)
    )
  } else if (complex_eigenpairs) {
    empty_certificate(
      tol,
      note = paste0(
        "general dense eigen oracle returned complex eigenpairs; ",
        "the SPD certification gate has no complex contract. ",
        "Inspect $values/$vectors directly."
      )
    )
  } else {
    empty_certificate(tol, note = "vectors not returned; residual certificate not computed")
  }

  warnings_msg <- if (identical(plan$method, "native dense generalized SPD LAPACK fallback")) {
    "using native dense generalized SPD LAPACK fallback; iterative engine not yet implemented"
  } else if (identical(plan$method, "dense LAPACK general eigen oracle (prototype fallback)")) {
    if (complex_eigenpairs) {
      "using dense general eigen oracle prototype; LAPACK returned complex eigenpairs which are not yet certified"
    } else {
      "using dense general eigen oracle prototype (non-Hermitian)"
    }
  } else if (identical(plan$method, "dense LAPACK Hermitian eigen oracle (prototype fallback)")) {
    "using dense Hermitian eigen oracle prototype"
  } else if (identical(plan$fallback, "dense oracle prototype")) {
    "using dense oracle prototype solver"
  } else {
    character()
  }

  make_eigen_result(
    values = vals,
    vectors = vecs,
    certificate = cert,
    iter = list(iterations = 1L, matvecs = 0L),
    requested = k,
    method_label = plan$method,
    target_label_value = target_label(a$target),
    plan = plan,
    warnings = warnings_msg
  )
}
