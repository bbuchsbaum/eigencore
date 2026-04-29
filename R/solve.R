#' Compute a partial eigendecomposition.
eig_partial <- function(A, k, target = largest(), B = NULL, method = auto(),
                        tol = 1e-8, maxit = NULL, vectors = TRUE, seed = NULL,
                        certify = TRUE,
                        allow_dense_fallback = c("auto", "never", "always")) {
  allow_dense_fallback <- match.arg(allow_dense_fallback)
  if (!is.null(seed)) {
    set.seed(seed)
  }
  P <- eigen_problem(A, metric = B, target = target,
                     transform = if (is_transform_method(method)) method else NULL)
  solve(P, k = k, method = method, tol = tol, maxit = maxit, vectors = vectors,
        certify = certify, allow_dense_fallback = allow_dense_fallback)
}

#' Compute a partial singular-value decomposition.
svd_partial <- function(A, rank, target = largest(), method = auto(), tol = 1e-8,
                        vectors = c("both", "left", "right", "none"),
                        seed = NULL, certify = TRUE,
                        allow_dense_fallback = c("auto", "never", "always")) {
  vectors <- match.arg(vectors)
  allow_dense_fallback <- match.arg(allow_dense_fallback)
  if (!is.null(seed)) {
    set.seed(seed)
  }
  P <- svd_problem(A, target = target)
  solve(P, rank = rank, method = method, tol = tol, vectors = vectors,
        certify = certify, allow_dense_fallback = allow_dense_fallback)
}

#' @export
solve.eigencore_eigen_problem <- function(a, b, k, method = auto(), tol = 1e-8,
                                          maxit = NULL, vectors = TRUE,
                                          certify = TRUE,
                                          allow_dense_fallback = c("auto", "never", "always"), ...) {
  allow_dense_fallback <- match.arg(allow_dense_fallback)
  plan <- plan_solver(a, k = k, method = method)
  if (is_transform_method(a$transform)) {
    if (identical(a$transform$kind, "shift_invert")) {
      return(solve_shift_invert_hermitian(
        a, k = k, method = a$transform, tol = tol,
        maxit = maxit, vectors = vectors, certify = certify, plan = plan
      ))
    }
    stop(
      "transform method '", a$transform$kind,
      "' is registered in transform_method_kinds() but has no solver dispatch wired in solve.eigencore_eigen_problem.",
      call. = FALSE
    )
  }
  use_lobpcg_method <- inherits(method, "eigencore_method") && identical(method$kind, "lobpcg")
  use_auto_generalized_lobpcg <- inherits(method, "eigencore_method") &&
    identical(method$kind, "auto") &&
    identical(plan$method, native_generalized_lobpcg_label())
  use_auto_reference_generalized_lobpcg <- inherits(method, "eigencore_method") &&
    identical(method$kind, "auto") &&
    identical(plan$method, reference_generalized_lobpcg_label())
  if (use_lobpcg_method || use_auto_generalized_lobpcg || use_auto_reference_generalized_lobpcg) {
    if (!identical(a$structure$kind, "hermitian")) {
      stop("lobpcg() prototype currently requires a Hermitian eigenproblem.", call. = FALSE)
    }
    controls <- plan$controls %||% list()
    method_preconditioner <- if (use_lobpcg_method) method$preconditioner else NULL
    method_constraints <- if (use_lobpcg_method) method$constraints else NULL
    method_maxit <- maxit %||% controls$maxit %||%
      if (use_lobpcg_method) method$maxit else getOption("eigencore.lobpcg_maxit", 200L)
    use_native_lobpcg <- native_lobpcg_supported(
      a$A,
      target = a$target,
      preconditioner = method_preconditioner,
      Bop = a$metric,
      constraints = method_constraints
    )
    use_native_generalized_lobpcg <- !use_native_lobpcg && !is.null(a$metric) &&
      native_generalized_lobpcg_supported(
        a$A,
        a$metric,
        target = a$target,
        preconditioner = method_preconditioner,
        constraints = method_constraints
      )
    iter <- if (use_native_lobpcg) {
      native_lobpcg_hermitian(
        a$A,
        k = k,
        target = a$target,
        tol = tol,
        maxit = method_maxit,
        preconditioner = method_preconditioner,
        constraints = method_constraints
      )
    } else if (use_native_generalized_lobpcg) {
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
    result <- list(
      values = iter$values,
      vectors = if (isTRUE(vectors)) iter$vectors else NULL,
      residuals = iter$residuals,
      backward_error = iter$backward_error,
      orthogonality = iter$orthogonality,
      nconv = sum(iter$certificate$converged),
      requested = k,
      iterations = iter$iterations,
      matvecs = iter$matvecs,
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
      locked = which(iter$certificate$converged),
      method = plan$method,
      target = target_label(a$target),
      plan = plan,
      certificate = iter$certificate,
      warnings = warning_msg
    )
    class(result) <- "eigencore_eigen_result"
    return(result)
  }
  use_lanczos <- should_use_lanczos(a, method, k = k)
  if (use_lanczos) {
    controls <- plan$controls %||% list()
    method_maxit <- controls$max_subspace %||%
      if (inherits(method, "eigencore_method") && identical(method$kind, "lanczos")) method$max_subspace else NULL
    method_reorth <- controls$reorthogonalize %||%
      if (inherits(method, "eigencore_method") && identical(method$kind, "lanczos")) method$reorthogonalize else TRUE
    method_max_restarts <- controls$max_restarts %||%
      if (inherits(method, "eigencore_method") && identical(method$kind, "lanczos")) method$max_restarts else NULL
    method_block <- controls$block %||%
      if (inherits(method, "eigencore_method") && identical(method$kind, "lanczos")) method$block else 1L
    iter <- if (should_use_native_lanczos(a, method, k = k)) {
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

    native_path <- plan$method %in% c(
      "native scalar thick-restart Hermitian Lanczos",
      "native block Hermitian Lanczos thick-restart candidate",
      "native block Hermitian Lanczos (thick restart, locking)"
    )
    warning_msg <- if (native_path) {
      if (!isTRUE(iter$certificate$passed)) {
        paste0(
          plan$method, " exhausted its current budget and did not converge all ",
          k,
          " requested pairs; restart budget/subspace: ",
          iter$restart$max_restarts %||% NA_integer_,
          "/",
          iter$restart$max_subspace %||% NA_integer_
        )
      } else {
        character()
      }
    } else {
      "using R-level prototype Lanczos; native hot loop not yet implemented"
    }

    result <- list(
      values = iter$values,
      vectors = iter$vectors,
      residuals = iter$residuals %||% iter$certificate$residuals,
      backward_error = iter$backward_error %||% iter$certificate$backward_error,
      orthogonality = iter$orthogonality %||% iter$certificate$orthogonality,
      nconv = sum(iter$certificate$converged),
      requested = k,
      iterations = iter$iterations,
      matvecs = iter$matvecs,
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
      locked = iter$locked %||% which(iter$certificate$converged),
      method = plan$method,
      target = target_label(a$target),
      plan = plan,
      certificate = iter$certificate,
      warnings = warning_msg
    )
    class(result) <- "eigencore_eigen_result"
    return(result)
  }

  if (should_use_native_dense_hermitian(a, method, k = k)) {
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
  )
  class(result) <- "eigencore_eigen_result"
  result
}

#' @export
solve.eigencore_svd_problem <- function(a, b, rank, method = auto(), tol = 1e-8,
                                        vectors = c("both", "left", "right", "none"),
                                        certify = TRUE,
                                        allow_dense_fallback = c("auto", "never", "always"), ...) {
  vectors <- match.arg(vectors)
  allow_dense_fallback <- match.arg(allow_dense_fallback)
  plan <- plan_solver(a, rank = rank, method = method)
  if (inherits(method, "eigencore_method") && identical(method$kind, "randomized")) {
    controls <- plan$controls %||% list()
    iter <- reference_randomized_svd(
      a$A,
      rank = rank,
      target = a$target,
      tol = tol,
      oversample = controls$oversample %||% method$oversample,
      n_iter = controls$n_iter %||% method$n_iter,
      vectors = vectors,
      refine = controls$refine %||% method$refine,
      normalizer = controls$normalizer %||% method$normalizer
    )
    cert <- if (isTRUE(certify) && !is.null(iter[["u"]]) && !is.null(iter[["v"]])) {
      iter$certificate
    } else {
      empty_certificate(tol, note = "both left and right vectors are required for full SVD certification")
    }
    result <- list(
      d = iter$d,
      u = iter[["u"]],
      v = iter[["v"]],
      values = iter$d,
      residuals = cert$residuals,
      backward_error = cert$backward_error,
      orthogonality = cert$orthogonality,
      nconv = sum(cert$converged),
      requested = rank,
      iterations = iter$iterations,
      matvecs = iter$matvecs,
      stage_seconds = iter$stage_seconds %||% iter$restart$stage_seconds %||% numeric(),
      method = plan$method,
      target = target_label(a$target),
      plan = plan,
      certificate = cert,
      restart = iter$restart,
      warnings = if (isTRUE(iter$restart$refinement_passed)) {
        "using reference randomized SVD prototype with certified native refinement"
      } else if (isTRUE(cert$passed)) {
        "using reference randomized SVD prototype with residual certification"
      } else {
        "using reference randomized SVD prototype; residual certificate did not meet tolerance"
      }
    )
    class(result) <- "eigencore_svd_result"
    return(result)
  }
  if (should_use_native_gram_svd(a, method, rank = rank)) {
    iter <- native_gram_svd(
      a$A,
      rank = rank,
      target = a$target,
      tol = tol,
      vectors = vectors
    )
    cert <- if (isTRUE(certify) && !is.null(iter[["u"]]) && !is.null(iter[["v"]])) {
      iter$certificate
    } else {
      empty_certificate(tol, note = "both left and right vectors are required for full SVD certification")
    }
    gram_cert <- cert
    gram_restart <- iter$restart
    fallback_attempted <- FALSE
    fallback_used <- FALSE
    fallback_iter <- NULL
    fallback_cert <- NULL
    if (isTRUE(certify) && !isTRUE(cert$passed) && vectors == "both") {
      fallback_attempted <- TRUE
      fallback_iter <- tryCatch(
        native_golub_kahan_svd(
          a$A,
          rank = rank,
          target = a$target,
          tol = tol,
          maxit = NULL,
          vectors = vectors
        ),
        error = function(e) {
          structure(list(error = conditionMessage(e)), class = "eigencore_fallback_error")
        }
      )
      if (!inherits(fallback_iter, "eigencore_fallback_error")) {
        fallback_cert <- fallback_iter$certificate
        fallback_used <- isTRUE(fallback_cert$passed) ||
          (is.finite(fallback_cert$max_backward_error) &&
            (!is.finite(cert$max_backward_error) ||
              fallback_cert$max_backward_error < cert$max_backward_error))
        if (isTRUE(fallback_used)) {
          iter <- fallback_iter
          cert <- fallback_cert
        }
      }
    }
    restart <- iter$restart
    restart$fallback_attempted <- fallback_attempted
    restart$fallback_used <- fallback_used
    restart$fallback_method <- if (fallback_attempted) "native prototype Golub-Kahan" else NA_character_
    restart$fallback_error <- if (inherits(fallback_iter, "eigencore_fallback_error")) {
      fallback_iter$error
    } else {
      NA_character_
    }
    restart$gram_certificate_passed <- isTRUE(gram_cert$passed)
    restart$gram_max_backward_error <- gram_cert$max_backward_error
    restart$gram_restart <- gram_restart
    if (fallback_attempted && !is.null(fallback_cert)) {
      restart$fallback_max_backward_error <- fallback_cert$max_backward_error
    }
    result <- list(
      d = iter$d,
      u = iter[["u"]],
      v = iter[["v"]],
      values = iter$d,
      residuals = cert$residuals,
      backward_error = cert$backward_error,
      orthogonality = cert$orthogonality,
      nconv = sum(cert$converged),
      requested = rank,
      iterations = iter$iterations,
      matvecs = iter$matvecs,
      stage_seconds = iter$stage_seconds %||% iter$restart$stage_seconds %||% numeric(),
      method = if (isTRUE(fallback_used)) {
        "native prototype Golub-Kahan fallback from Gram SVD"
      } else {
        plan$method
      },
      target = target_label(a$target),
      plan = plan,
      certificate = cert,
      restart = restart,
      warnings = if (isTRUE(fallback_used)) {
        "native certified Gram SVD special case failed certification; using native Golub-Kahan fallback"
      } else if (fallback_attempted) {
        "native certified Gram SVD special case failed certification; native Golub-Kahan fallback was not better"
      } else {
        "using native certified Gram SVD special case; residuals certified in original coordinates"
      }
    )
    class(result) <- "eigencore_svd_result"
    return(result)
  }
  if (should_use_native_retained_golub_kahan(a, method, rank = rank)) {
    iter <- native_block_golub_kahan_retained_cycle_svd(
      a$A,
      rank = rank,
      target = a$target,
      tol = tol,
      vectors = vectors
    )
    cert <- if (isTRUE(certify) && !is.null(iter[["u"]]) && !is.null(iter[["v"]])) {
      iter$certificate
    } else {
      empty_certificate(tol, note = "both left and right vectors are required for full SVD certification")
    }
    result <- list(
      d = iter$d,
      u = iter[["u"]],
      v = iter[["v"]],
      values = iter$d,
      residuals = cert$residuals,
      backward_error = cert$backward_error,
      orthogonality = cert$orthogonality,
      nconv = sum(cert$converged),
      requested = rank,
      iterations = iter$iterations,
      matvecs = iter$matvecs,
      stage_seconds = iter$stage_seconds %||% iter$restart$stage_seconds %||% numeric(),
      method = plan$method,
      target = target_label(a$target),
      plan = plan,
      certificate = cert,
      restart = iter$restart,
      warnings = if (isTRUE(cert$passed)) {
        "using native retained Golub-Kahan SVD with thick restart"
      } else {
        "using native retained Golub-Kahan SVD; adaptive subspace budget exhausted before full certification"
      }
    )
    class(result) <- "eigencore_svd_result"
    return(result)
  }
  use_gk <- should_use_golub_kahan(a, method)
  if (use_gk) {
    controls <- plan$controls %||% list()
    method_maxit <- controls$max_subspace %||%
      if (inherits(method, "eigencore_method") && identical(method$kind, "golub_kahan")) method$max_subspace else NULL
    method_reorth <- controls$reorthogonalize %||%
      if (inherits(method, "eigencore_method") && identical(method$kind, "golub_kahan")) method$reorthogonalize else TRUE
    iter <- if (should_use_native_golub_kahan(a, method)) {
      native_golub_kahan_svd(
        a$A,
        rank = rank,
        target = a$target,
        tol = tol,
        maxit = method_maxit,
        vectors = vectors
      )
    } else {
      reference_golub_kahan_svd(
        a$A,
        rank = rank,
        target = a$target,
        tol = tol,
        maxit = method_maxit,
        vectors = vectors,
        reorthogonalize = method_reorth
      )
    }

    cert <- if (isTRUE(certify) && !is.null(iter[["u"]]) && !is.null(iter[["v"]])) {
      iter$certificate
    } else {
      empty_certificate(tol, note = "both left and right vectors are required for full SVD certification")
    }
    result <- list(
      d = iter$d,
      u = iter[["u"]],
      v = iter[["v"]],
      values = iter$d,
      residuals = cert$residuals,
      backward_error = cert$backward_error,
      orthogonality = cert$orthogonality,
      nconv = sum(cert$converged),
      requested = rank,
      iterations = iter$iterations,
      matvecs = iter$matvecs,
      stage_seconds = iter$stage_seconds %||% iter$restart$stage_seconds %||% numeric(),
      method = plan$method,
      target = target_label(a$target),
      plan = plan,
      certificate = cert,
      restart = iter$restart,
      warnings = if (identical(plan$method, "native prototype Golub-Kahan")) {
        if (isTRUE(iter$restart$converged)) {
          "using native prototype Golub-Kahan iteration with adaptive subspace growth"
        } else {
          "using native prototype Golub-Kahan iteration; adaptive subspace budget exhausted before full certification"
        }
      } else {
        "using R-level prototype Golub-Kahan; native hot loop not yet implemented"
      }
    )
    class(result) <- "eigencore_svd_result"
    return(result)
  }

  A <- materialize_dense_fallbacks(list(A = a$A), allow = allow_dense_fallback)$A
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
    both_ends = {
      low <- order(Re(x), decreasing = FALSE)
      high <- order(Re(x), decreasing = TRUE)
      unique(c(utils::head(low, target$value$k_low), utils::head(high, target$value$k_high)))
    },
    order(Re(x), decreasing = TRUE)
  )
}

#' @keywords internal
should_use_lanczos <- function(problem, method, k = NULL) {
  if (!is.null(problem$metric)) {
    return(FALSE)
  }
  if (!identical(problem$structure$kind, "hermitian")) {
    return(FALSE)
  }
  if (inherits(method, "eigencore_method") && identical(method$kind, "lanczos")) {
    return(TRUE)
  }
  if (!(inherits(method, "eigencore_method") && identical(method$kind, "auto"))) {
    return(FALSE)
  }
  # NN-3: any sparse Hermitian operator must route to the Lanczos branch
  # (native or reference) before the dense fallback. A sparse operator
  # carrying a sparseMatrix in metadata$matrix or metadata$source must
  # NEVER silently densify under allow_dense_fallback = "always" when an
  # honest reference_lanczos_hermitian path is already available.
  src <- source_or_null(problem$A)
  matrix_meta <- problem$A$metadata$matrix
  has_sparse_payload <- inherits(matrix_meta, "sparseMatrix") ||
    inherits(src, "sparseMatrix")
  if (isTRUE(has_sparse_payload)) {
    return(TRUE)
  }
  is.null(src) || auto_dense_partial_lanczos(problem, k)
}

#' @keywords internal
should_use_native_lanczos <- function(problem, method, k = NULL) {
  if (!is.null(problem$metric)) {
    return(FALSE)
  }
  if (!identical(problem$structure$kind, "hermitian")) {
    return(FALSE)
  }
  if (!has_native_kernel(problem$A)) {
    return(FALSE)
  }
  if (!native_lanczos_target_supported(problem$target)) {
    return(FALSE)
  }
  if (inherits(method, "eigencore_method") && identical(method$kind, "lanczos")) {
    return(TRUE)
  }
  inherits(method, "eigencore_method") &&
    identical(method$kind, "auto") &&
    (identical(native_kernel_kind(problem$A), "csc") ||
       auto_dense_partial_lanczos(problem, k))
}

#' @keywords internal
should_use_native_dense_hermitian <- function(problem, method, k = NULL) {
  if (!is.null(problem$metric)) {
    return(FALSE)
  }
  if (!identical(problem$structure$kind, "hermitian")) {
    return(FALSE)
  }
  if (inherits(method, "eigencore_method") && !identical(method$kind, "auto")) {
    return(FALSE)
  }
  if (auto_dense_partial_lanczos(problem, k)) {
    return(FALSE)
  }
  src <- source_or_null(problem$A)
  is.matrix(src) && is.double(src)
}

#' @keywords internal
auto_dense_partial_lanczos <- function(problem, k = NULL) {
  if (is.null(k)) {
    k <- problem$requested %||% NA_integer_
  }
  if (is.na(k) || length(k) != 1L) {
    return(FALSE)
  }
  k <- as.integer(k)
  if (!is.null(problem$metric) || !identical(problem$structure$kind, "hermitian")) {
    return(FALSE)
  }
  if (!native_lanczos_target_supported(problem$target)) {
    return(FALSE)
  }
  source <- source_or_null(problem$A)
  if (!(is.matrix(source) && is.double(source))) {
    return(FALSE)
  }
  n <- as.integer(problem$A$dim[1L])
  min_n <- as.integer(getOption("eigencore.dense_partial_lanczos_min_n", 128L))
  max_fraction <- as.numeric(getOption("eigencore.dense_partial_lanczos_max_fraction", 0.25))
  if (length(min_n) != 1L || is.na(min_n) || min_n < 1L) {
    min_n <- 128L
  }
  if (length(max_fraction) != 1L || is.na(max_fraction) || max_fraction <= 0 || max_fraction > 1) {
    max_fraction <- 0.25
  }
  n >= min_n && k >= 1L && k < n && (k / n) <= max_fraction
}

#' @keywords internal
native_dense_symmetric_eigen <- function(A) {
  .Call("eigencore_dense_symmetric_eigen", as.matrix(A), PACKAGE = "eigencore")
}

#' @keywords internal
native_dense_symmetric_eigen_selected <- function(A, k, target) {
  .Call(
    "eigencore_dense_symmetric_eigen_selected",
    as.matrix(A),
    as.integer(k),
    as.integer(lanczos_target_kind(target)),
    PACKAGE = "eigencore"
  )
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
dense_fallback_budget_bytes <- function() {
  mb <- getOption("eigencore.dense_fallback_mb", NULL)
  if (is.null(mb)) {
    legacy <- getOption("eigencore.max_dense_fallback_bytes", NULL)
    if (!is.null(legacy)) {
      if (!is.numeric(legacy) || length(legacy) != 1L || is.na(legacy) || legacy < 0) {
        stop("Option eigencore.max_dense_fallback_bytes must be one non-negative numeric value.", call. = FALSE)
      }
      return(legacy)
    }
    mb <- 256
  }
  if (!is.numeric(mb) || length(mb) != 1L || is.na(mb) || mb < 0) {
    stop("Option eigencore.dense_fallback_mb must be one non-negative numeric value.", call. = FALSE)
  }
  mb * 1e6
}

#' @keywords internal
materialize_dense_fallbacks <- function(ops, budget = dense_fallback_budget_bytes(),
                                        allow = c("auto", "never", "always")) {
  allow <- match.arg(allow)
  if (identical(allow, "never")) {
    stop("Dense fallback disabled by allow_dense_fallback = 'never'.", call. = FALSE)
  }
  stopifnot(is.list(ops), length(ops) > 0L)
  roles <- names(ops)
  if (is.null(roles) || any(!nzchar(roles))) {
    roles <- paste0("operator_", seq_along(ops))
  }

  allow_sparse <- identical(allow, "always")
  if (!identical(allow, "always")) {
    sizes <- vapply(
      seq_along(ops),
      function(i) dense_fallback_bytes(ops[[i]], roles[[i]], allow_sparse = allow_sparse),
      numeric(1)
    )
    total <- sum(sizes)
    if (is.finite(budget) && total > budget) {
      stop(
        "Dense fallback would materialize ", format_bytes(total),
        " across ", paste(roles, collapse = ", "),
        ", exceeding eigencore.dense_fallback_mb = ", format_bytes(budget), ".",
        call. = FALSE
      )
    }
  }

  out <- Map(function(op, role) {
    materialize_dense_fallback(op, role, allow_sparse = allow_sparse)
  }, ops, roles)
  names(out) <- roles
  out
}

#' @keywords internal
dense_fallback_bytes <- function(op, role = "operator", allow_sparse = FALSE) {
  src <- dense_fallback_source(op, role, allow_sparse = allow_sparse)
  dims <- dim(src)
  if (length(dims) != 2L) {
    stop("Dense fallback source for ", role, " is not matrix-like.", call. = FALSE)
  }
  prod(as.numeric(dims)) * 8
}

#' @keywords internal
materialize_dense_fallback <- function(op, role = "operator", allow_sparse = FALSE) {
  src <- dense_fallback_source(op, role, allow_sparse = allow_sparse)
  as.matrix(src)
}

#' @keywords internal
dense_fallback_source <- function(op, role = "operator", allow_sparse = FALSE) {
  src <- if (inherits(op, "eigencore_operator")) {
    op$metadata$source %||% op$metadata$matrix
  } else {
    op
  }
  if (is.null(src)) {
    stop("Dense fallback for ", role, " needs an explicit matrix-backed operator.", call. = FALSE)
  }
  if (inherits(src, "sparseMatrix") && !isTRUE(allow_sparse)) {
    stop(
      "Refusing to densify sparse ", role,
      " in a solver path. Use a native sparse iterative path, pass as.matrix(",
      role, ") explicitly, or set allow_dense_fallback = 'always' to opt into dense fallback.",
      call. = FALSE
    )
  }
  src
}

#' @keywords internal
format_bytes <- function(bytes) {
  units <- c("B", "KB", "MB", "GB", "TB")
  value <- as.numeric(bytes)
  unit <- 1L
  while (is.finite(value) && abs(value) >= 1024 && unit < length(units)) {
    value <- value / 1024
    unit <- unit + 1L
  }
  paste0(format(round(value, 2), trim = TRUE), " ", units[[unit]])
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

#' @keywords internal
should_use_native_golub_kahan <- function(problem, method) {
  if (is.null(problem$A$apply_adjoint)) {
    return(FALSE)
  }
  if (!has_native_kernel(problem$A)) {
    return(FALSE)
  }
  if (inherits(method, "eigencore_method") && identical(method$kind, "golub_kahan")) {
    return(TRUE)
  }
  inherits(method, "eigencore_method") &&
    identical(method$kind, "auto") &&
    identical(native_kernel_kind(problem$A), "csc")
}

#' @keywords internal
should_use_native_retained_golub_kahan <- function(problem, method, rank = NULL) {
  if (is.null(problem$A$apply_adjoint)) {
    return(FALSE)
  }
  if (!inherits(method, "eigencore_method") || !identical(method$kind, "auto")) {
    return(FALSE)
  }
  if (should_use_native_gram_svd(problem, method, rank = rank)) {
    return(FALSE)
  }
  storage <- problem$A$metadata$storage %||% NULL
  identical(storage, "dgCMatrix")
}

#' @keywords internal
should_use_native_gram_svd <- function(problem, method, rank = NULL) {
  if (!inherits(method, "eigencore_method") || !identical(method$kind, "auto")) {
    return(FALSE)
  }
  if (!native_gram_svd_target_supported(problem$target)) {
    return(FALSE)
  }
  if (!has_native_kernel(problem$A)) {
    return(FALSE)
  }
  dims <- problem$A$dim
  reduced <- min(dims)
  full <- max(dims)
  rank <- as.integer(rank %||% problem$rank %||% 1L)
  reduced <= getOption("eigencore.gram_svd_max_dimension", 512L) &&
    full >= 2L * reduced &&
    rank <= reduced / 2
}
