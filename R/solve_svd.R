# Per-method runtime helpers for solve.eigencore_svd_problem.
# Each helper takes the svd problem `a`, the requested rank, and the resolved
# `plan` produced by plan_solver(); it returns a classed eigencore_svd_result.
# The dispatcher in R/solve.R selects between these based on plan$method labels.

#' @keywords internal
solve_svd_randomized <- function(a, rank, method, tol, vectors, certify, plan) {
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
  warnings_msg <- if (isTRUE(iter$restart$refinement_passed)) {
    "using reference randomized SVD prototype with certified native refinement"
  } else if (isTRUE(cert$passed)) {
    "using reference randomized SVD prototype with residual certification"
  } else {
    "using reference randomized SVD prototype; residual certificate did not meet tolerance"
  }
  make_svd_result(
    d = iter$d,
    u = iter[["u"]],
    v = iter[["v"]],
    certificate = cert,
    iter = iter,
    requested = rank,
    method_label = plan$method,
    target_label_value = target_label(a$target),
    plan = plan,
    warnings = warnings_msg
  )
}

#' @keywords internal
try_svd_partial_native_gram_fastpath <- function(A, rank, target, method, tol,
                                                 vectors, certify) {
  if (!inherits(method, "eigencore_method") || !identical(method$kind, "auto")) {
    return(NULL)
  }
  if (!native_gram_svd_target_supported(target)) {
    return(NULL)
  }
  if (!inherits(A, "dgCMatrix")) {
    return(NULL)
  }
  dims <- dim(A)
  if (length(dims) != 2L || any(!is.finite(dims)) || dims[1L] >= dims[2L]) {
    return(NULL)
  }
  rank <- as.integer(rank)
  if (length(rank) != 1L || is.na(rank) || rank < 1L) {
    return(NULL)
  }
  reduced <- min(dims)
  full <- max(dims)
  if (rank > reduced / 2) {
    return(NULL)
  }
  gram_max <- as.integer(getOption("eigencore.gram_svd_max_dimension", 512L))
  if (reduced > gram_max || full < 2L * reduced) {
    return(NULL)
  }

  target_kind <- if (inherits(target, "eigencore_target")) target$kind else "largest"
  if (identical(vectors, "both") && isTRUE(certify) && identical(target_kind, "largest")) {
    fast_native <- .Call(
      "eigencore_csc_left_gram_svd_fast_result",
      methods::slot(A, "i"),
      methods::slot(A, "p"),
      methods::slot(A, "x"),
      methods::slot(A, "Dim"),
      as.integer(rank),
      as.numeric(tol),
      PACKAGE = "eigencore"
    )
    if (!is.null(fast_native)) {
      return(fast_native)
    }
  }

  native <- .Call(
    "eigencore_csc_left_gram_svd",
    methods::slot(A, "i"),
    methods::slot(A, "p"),
    methods::slot(A, "x"),
    methods::slot(A, "Dim"),
    as.integer(rank),
    as.numeric(tol),
    PACKAGE = "eigencore"
  )
  zero_tol <- gram_svd_zero_tolerance(native$d, tol)
  if (any(native$d <= zero_tol)) {
    op <- as_operator(A)
    plan <- native_gram_svd_fast_plan(op, rank, target)
    return(solve_svd_gram(
      list(A = op, target = target),
      rank = rank,
      tol = tol,
      vectors = vectors,
      certify = certify,
      plan = plan
    ))
  }

  plan <- native_gram_svd_fast_plan_from_dims(dim(A), rank, target)
  cert <- if (isTRUE(certify) && identical(vectors, "both")) {
    new_certificate(
      tol = tol,
      residuals = list(
        left = native$diagnostics$left,
        right = native$diagnostics$right,
        combined = native$diagnostics$combined
      ),
      backward_error = native$diagnostics$backward_error,
      orthogonality = native$diagnostics$orthogonality,
      converged = native$diagnostics$converged,
      scale = native$diagnostics$scale,
      norm_bound_type = "frobenius_exact"
    )
  } else {
    empty_certificate(tol, note = "both left and right vectors are required for full SVD certification")
  }
  u <- if (vectors %in% c("both", "left")) native$u else NULL
  v <- if (vectors %in% c("both", "right")) native$v else NULL
  restart <- list(
    kind = "gram_svd_special_case",
    implemented = TRUE,
    native = TRUE,
    gram_side = "left",
    gram_dimension = min(dim(A)),
    native_gram_kernel = "csc_left_gram",
    native_gram_eigensolver = native$eigensolver %||% "lapack_dsyevr",
    native_gram_subspace_max_backward_error =
      native$subspace_max_backward_error %||% NA_real_,
    native_implicit_normal_lanczos_max_backward_error =
      native$implicit_lanczos_max_backward_error %||% NA_real_,
    native_implicit_normal_lanczos_iterations =
      native$implicit_lanczos_iterations %||% 0L,
    native_gram_krylov_iterations =
      native$gram_krylov_iterations %||% 0L,
    normal_operator_implicit =
      identical(native$eigensolver %||% "", "implicit_normal_lanczos"),
    materialized_gram =
      !identical(native$eigensolver %||% "", "implicit_normal_lanczos"),
    stage_seconds = native$stage_seconds,
    zero_singular_completion = FALSE,
    zero_singular_threshold = zero_tol,
    certificate_reuses_gram_sides = TRUE,
    certified_in_original_coordinates = TRUE,
    fallback_attempted = FALSE,
    fallback_used = FALSE,
    fallback_method = NA_character_,
    fallback_error = NA_character_,
    gram_certificate_passed = isTRUE(cert$passed),
    gram_max_backward_error = cert$max_backward_error
  )
  out <- list(
    d = native$d,
    u = u,
    v = v,
    values = native$d,
    residuals = cert$residuals,
    backward_error = cert$backward_error,
    orthogonality = cert$orthogonality,
    nconv = sum(cert$converged),
    requested = rank,
    iterations = 1L,
    matvecs = 1L,
    stage_seconds = native$stage_seconds,
    method = plan$method,
    target = target_label(target),
    plan = plan,
    certificate = cert,
    restart = restart,
    warnings = "using native certified Gram SVD special case; residuals certified in original coordinates"
  )
  class(out) <- "eigencore_svd_result"
  out
}

#' @keywords internal
native_gram_svd_fast_plan <- function(op, rank, target) {
  dims <- op$dim
  full <- max(dims)
  problem <- list(type = "svd", A = op, target = target)
  chosen <- "native certified Gram SVD special case"
  controls <- svd_plan_controls(problem, rank = rank, method = auto(), chosen = chosen)
  controls$svd_partial_fastpath <- TRUE
  controls$full_dimension <- as.integer(full)
  new_plan(
    problem,
    k = as.integer(rank),
    method = chosen,
    reasons = c(
      paste0("target: ", target_label(target)),
      "rectangular SVD problem",
      "adjoint is available",
      "small rectangular sparse problem: materializes the smaller Gram matrix as an explicit certified special case",
      "built-in sparse CSC operator has native block apply",
      "direct svd_partial() fast path avoids S3 dispatch overhead"
    ),
    fallback = "native Golub-Kahan if Gram special case is disabled or uncertified",
    controls = controls
  )
}

#' @keywords internal
native_gram_svd_fast_plan_from_dims <- function(dims, rank, target) {
  dims <- as.integer(dims)
  full <- max(dims)
  reduced <- min(dims)
  problem <- list(type = "svd", target = target)
  controls <- list(
    gram_side = if (dims[2L] <= dims[1L]) "right" else "left",
    gram_dimension = reduced,
    gram_max_dimension = as.integer(getOption("eigencore.gram_svd_max_dimension", 512L)),
    rank_fraction_limit = 0.5,
    certified_in_original_coordinates = TRUE,
    materializes = "smaller Gram matrix only",
    fallback_policy = "certification-gated",
    runtime_fallback = "native Golub-Kahan if original-coordinate certificate is weaker",
    fallback_requires_vectors = "both",
    svd_partial_fastpath = TRUE,
    full_dimension = as.integer(full)
  )
  new_plan(
    problem,
    k = as.integer(rank),
    method = "native certified Gram SVD special case",
    reasons = c(
      paste0("target: ", target_label(target)),
      "rectangular SVD problem",
      "adjoint is available",
      "small rectangular sparse problem: materializes the smaller Gram matrix as an explicit certified special case",
      "built-in sparse CSC operator has native block apply",
      "direct svd_partial() fast path avoids S3 dispatch overhead"
    ),
    fallback = "native Golub-Kahan if Gram special case is disabled or uncertified",
    controls = controls
  )
}

#' @keywords internal
solve_svd_gram <- function(a, rank, tol, vectors, certify, plan) {
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
  warnings_msg <- if (isTRUE(fallback_used)) {
    "native certified Gram SVD special case failed certification; using native Golub-Kahan fallback"
  } else if (fallback_attempted) {
    "native certified Gram SVD special case failed certification; native Golub-Kahan fallback was not better"
  } else {
    "using native certified Gram SVD special case; residuals certified in original coordinates"
  }
  method_label_used <- if (isTRUE(fallback_used)) {
    "native prototype Golub-Kahan fallback from Gram SVD"
  } else {
    plan$method
  }
  iter$restart <- restart
  make_svd_result(
    d = iter$d,
    u = iter[["u"]],
    v = iter[["v"]],
    certificate = cert,
    iter = iter,
    requested = rank,
    method_label = method_label_used,
    target_label_value = target_label(a$target),
    plan = plan,
    warnings = warnings_msg
  )
}

#' @keywords internal
solve_svd_retained_golub_kahan <- function(a, rank, tol, vectors, certify, plan) {
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
  warnings_msg <- if (isTRUE(cert$passed)) {
    "using native retained Golub-Kahan SVD with thick restart"
  } else {
    "using native retained Golub-Kahan SVD; adaptive subspace budget exhausted before full certification"
  }
  make_svd_result(
    d = iter$d,
    u = iter[["u"]],
    v = iter[["v"]],
    certificate = cert,
    iter = iter,
    requested = rank,
    method_label = plan$method,
    target_label_value = target_label(a$target),
    plan = plan,
    warnings = warnings_msg
  )
}

#' @keywords internal
solve_svd_golub_kahan <- function(a, rank, method, tol, vectors, certify, plan) {
  controls <- plan$controls %||% list()
  method_maxit <- controls$max_subspace %||%
    if (inherits(method, "eigencore_method") && identical(method$kind, "golub_kahan")) method$max_subspace else NULL
  method_reorth <- controls$reorthogonalize %||%
    if (inherits(method, "eigencore_method") && identical(method$kind, "golub_kahan")) method$reorthogonalize else TRUE
  iter <- if (identical(plan$method, "native prototype Golub-Kahan")) {
    native_golub_kahan_svd(
      a$A,
      rank = rank,
      target = a$target,
      tol = tol,
      maxit = method_maxit,
      vectors = vectors,
      reorthogonalize = method_reorth
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
  warnings_msg <- if (identical(plan$method, "native prototype Golub-Kahan")) {
    if (isTRUE(iter$restart$converged)) {
      "using native prototype Golub-Kahan iteration with adaptive subspace growth"
    } else {
      "using native prototype Golub-Kahan iteration; adaptive subspace budget exhausted before full certification"
    }
  } else {
    "using R-level prototype Golub-Kahan; native hot loop not yet implemented"
  }
  make_svd_result(
    d = iter$d,
    u = iter[["u"]],
    v = iter[["v"]],
    certificate = cert,
    iter = iter,
    requested = rank,
    method_label = plan$method,
    target_label_value = target_label(a$target),
    plan = plan,
    warnings = warnings_msg
  )
}

#' @keywords internal
solve_svd_dense <- function(a, rank, tol, vectors, certify, allow_dense_fallback, plan) {
  A <- materialize_dense_fallbacks(list(A = a$A), allow = allow_dense_fallback)$A
  decomp <- if (identical(plan$method, "native dense LAPACK SVD fallback")) {
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

  warnings_msg <- if (identical(plan$method, "native dense LAPACK SVD fallback")) {
    "using native dense LAPACK SVD fallback; iterative engine not yet implemented"
  } else {
    "using dense oracle prototype solver"
  }
  # The dense path does not surface a restart record; build a minimal iter that
  # the shared make_svd_result builder consumes.
  make_svd_result(
    d = d,
    u = u,
    v = v,
    certificate = cert,
    iter = list(iterations = 1L, matvecs = 0L, restart = NULL,
                stage_seconds = numeric()),
    requested = rank,
    method_label = plan$method,
    target_label_value = target_label(a$target),
    plan = plan,
    warnings = warnings_msg,
    extras = list(restart = NULL)
  )
}
