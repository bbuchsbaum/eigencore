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
  generalized_path <- identical(plan$method, generalized_lanczos_label())
  iter <- if (generalized_path) {
    reference_generalized_lanczos_hermitian(
      a$A,
      a$metric,
      k = k,
      target = a$target,
      tol = tol,
      maxit = maxit %||% method_maxit,
      vectors = vectors,
      reorthogonalize = method_reorth
    )
  } else if (plan_dispatches_native_lanczos(plan)) {
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

  warning_msg <- if (generalized_path) {
    if (!isTRUE(iter$certificate$passed)) {
      paste0(
        plan$method, " exhausted its current subspace before all ",
        k, " requested generalized pairs converged"
      )
    } else {
      "using reference generalized SPD B-orthogonal Lanczos refinement; native generalized Lanczos hot loop not yet implemented"
    }
  } else if (plan_dispatches_native_lanczos(plan)) {
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
      generalized = isTRUE(iter$generalized),
      metric_solve = iter$metric_solve %||% NULL,
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
solve_eigen_arnoldi <- function(a, k, method, tol, maxit, vectors, certify, plan) {
  controls <- plan$controls %||% list()
  native_path <- identical(plan$method, native_arnoldi_label())
  default_maxit <- if (native_path) {
    native_arnoldi_default_max_subspace(a$A$dim[[1L]], k)
  } else {
    max(k + 8L, 2L * k + 4L)
  }
  method_maxit <- maxit %||% controls$max_subspace %||% default_maxit
  method_max_restarts <- controls$max_restarts %||% 0L
  arnoldi_solver <- if (native_path) native_arnoldi_general else reference_arnoldi_general
  iter <- arnoldi_solver(
    a$A,
    k = k,
    target = a$target,
    tol = tol,
    maxit = method_maxit,
    max_restarts = method_max_restarts,
    vectors = vectors
  )
  warning_msg <- if (native_path && isTRUE(iter$certificate$passed)) {
    "using native Arnoldi cycle with native Ritz extraction; right residuals certified"
  } else if (native_path) {
    "using native Arnoldi cycle with native Ritz extraction; result did not pass certificate"
  } else if (isTRUE(iter$certificate$passed)) {
    "using R-level reference Arnoldi prototype; native nonsymmetric Arnoldi not yet implemented"
  } else {
    "using R-level reference Arnoldi prototype; result did not pass certificate"
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
      residuals = iter$certificate$residuals,
      backward_error = iter$certificate$backward_error,
      orthogonality = iter$certificate$orthogonality,
      restarts = iter$restart$restart_count,
      restart = iter$restart,
      locked = which(iter$certificate$converged)
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
solve_eigen_native_tridiagonal_hermitian <- function(a, k, tol, vectors, certify, plan) {
  A <- a$A$metadata$matrix %||% source_or_null(a$A)
  parts <- shift_invert_tridiagonal_parts(A, shift = 0)
  if (is.null(parts)) {
    stop("native tridiagonal Hermitian eigensolver requires a symmetric tridiagonal source.", call. = FALSE)
  }
  eig <- native_tridiagonal_eigen_selected(parts$diag, parts$upper, k, a$target)
  vals <- eig$values
  vecs_for_cert <- eig$vectors
  cert <- if (isTRUE(certify) && ncol(vecs_for_cert) > 0L) {
    native_tridiagonal_eigen_certificate(a$A, parts, vals, vecs_for_cert, tol = tol)
  } else {
    empty_certificate(
      tol,
      note = if (!isTRUE(certify)) {
        "native tridiagonal eigensolver: certification disabled by caller"
      } else {
        "native tridiagonal eigensolver: no eigenpairs returned; residual certificate not computed"
      }
    )
  }
  restart <- list(
    kind = "tridiagonal_lapack_selected",
    native = TRUE,
    implemented = TRUE,
    selected = length(vals)
  )
  make_eigen_result(
    values = vals,
    vectors = if (vectors) vecs_for_cert else NULL,
    certificate = cert,
    iter = list(iterations = 1L, matvecs = 0L, restart = restart),
    requested = k,
    method_label = plan$method,
    target_label_value = target_label(a$target),
    plan = plan,
    warnings = character(),
    extras = list(
      residuals = cert$residuals,
      backward_error = cert$backward_error,
      orthogonality = cert$orthogonality,
      restarts = 0L,
      restart = restart,
      locked = which(cert$converged)
    )
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

  general_dense <- is.null(B) && !identical(a$structure$kind, "hermitian")
  cert <- if (certify && !is.null(vecs) && isTRUE(general_dense)) {
    certify_dense_general_eigen(A, vals, vecs, tol = tol)
  } else if (certify && !is.null(vecs)) {
    certify_eigen(
      A,
      vals,
      vecs,
      B = B,
      tol = tol,
      require_orthogonality = identical(a$structure$kind, "hermitian") || !is.null(B)
    )
  } else {
    empty_certificate(tol, note = "vectors not returned; residual certificate not computed")
  }

  warnings_msg <- if (identical(plan$method, "native dense generalized SPD LAPACK fallback")) {
    "using native dense generalized SPD LAPACK fallback; iterative engine not yet implemented"
  } else if (identical(plan$method, "dense LAPACK general eigen oracle (prototype fallback)")) {
    if (isTRUE(certify) && !is.null(vecs)) {
      "using dense general eigen oracle prototype (non-Hermitian); right residuals certified"
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
