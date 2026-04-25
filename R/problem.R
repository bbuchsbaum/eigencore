#' Define an eigenproblem.
eigen_problem <- function(A, metric = NULL, structure = NULL, target = largest(),
                          transform = NULL) {
  Aop <- as_operator(A)
  if (is.null(structure)) {
    structure <- Aop$structure
  }
  Bop <- if (is.null(metric)) NULL else as_operator(metric)
  problem <- list(
    type = "eigen",
    A = Aop,
    metric = Bop,
    structure = structure,
    target = target,
    transform = transform
  )
  class(problem) <- "eigencore_eigen_problem"
  problem
}

#' Define an SVD problem.
svd_problem <- function(A, domain = NULL, codomain = NULL, target = largest()) {
  Aop <- as_operator(A)
  if (is.null(domain)) {
    domain <- euclidean(Aop$dim[2L])
  }
  if (is.null(codomain)) {
    codomain <- euclidean(Aop$dim[1L])
  }
  problem <- list(
    type = "svd",
    A = Aop,
    domain = domain,
    codomain = codomain,
    target = target
  )
  class(problem) <- "eigencore_svd_problem"
  problem
}

#' Plan a solver for a problem.
plan_solver <- function(problem, ...) {
  UseMethod("plan_solver")
}

#' @export
plan_solver.eigencore_eigen_problem <- function(problem, k, method = auto(), ...) {
  has_metric <- !is.null(problem$metric)
  has_shift <- inherits(problem$transform, "eigencore_method") &&
    identical(problem$transform$kind, "shift_invert")
  is_hermitian <- identical(problem$structure$kind, "hermitian")
  source_matrix <- source_or_null(problem$A)
  metric_matrix <- if (has_metric) source_or_null(problem$metric) else NULL
  is_dense_source <- is.matrix(source_matrix) && is.double(source_matrix)
  is_dense_metric <- is.matrix(metric_matrix) && is.double(metric_matrix)
  is_native_csc <- identical(problem$A$metadata$storage, "dgCMatrix")
  preconditioner_reason <- if (inherits(method, "eigencore_method") &&
    identical(method$kind, "lobpcg")) {
    preconditioner_plan_reason(method$preconditioner)
  } else {
    NULL
  }
  native_lanczos_supported <- is_hermitian &&
    native_lanczos_target_supported(problem$target) &&
    (is_native_csc || is_dense_source)
  lanczos_block <- if (inherits(method, "eigencore_method") &&
    identical(method$kind, "lanczos")) method$block %||% 1L else 1L

  chosen <- if (inherits(method, "eigencore_method") && method$kind != "auto") {
    if (identical(method$kind, "lanczos")) {
      if (!is_hermitian) {
        "dense LAPACK eigen oracle (Lanczos requires Hermitian structure)"
      } else if (native_lanczos_supported && lanczos_block > 1L) {
        "native block Hermitian Lanczos prototype"
      } else if (native_lanczos_supported) {
        "native scalar thick-restart Hermitian Lanczos"
      } else {
        "reference Hermitian Lanczos (prototype/oracle fallback)"
      }
    } else if (identical(method$kind, "lobpcg")) {
      if (is_hermitian && is.null(problem$metric)) {
        "reference LOBPCG prototype"
      } else if (is_hermitian && !is.null(problem$metric)) {
        "reference generalized SPD LOBPCG prototype"
      } else {
        "dense LAPACK eigen oracle (LOBPCG prototype requires Hermitian structure)"
      }
    } else {
      method_label(method)
    }
  } else if (has_shift) {
    "shift-invert requested (not implemented)"
  } else if (has_metric && is_hermitian && is_dense_source && is_dense_metric) {
    "native dense generalized SPD LAPACK fallback"
  } else if (has_metric && is_hermitian) {
    "dense generalized SPD LAPACK oracle (prototype fallback)"
  } else if (is_hermitian && is_native_csc && native_lanczos_target_supported(problem$target)) {
    "native scalar thick-restart Hermitian Lanczos"
  } else if (is_hermitian && is_native_csc) {
    "reference Hermitian Lanczos (target unsupported by native path)"
  } else if (is_hermitian && is.null(source_or_null(problem$A))) {
    "reference Hermitian Lanczos (prototype/oracle fallback)"
  } else if (auto_dense_partial_lanczos(problem, k)) {
    "native scalar thick-restart Hermitian Lanczos"
  } else if (is_hermitian && is_dense_source) {
    "native dense Hermitian LAPACK fallback"
  } else if (is_hermitian) {
    "dense LAPACK eigen oracle (prototype fallback)"
  } else {
    "dense LAPACK eigen oracle (prototype fallback)"
  }

  reasons <- c(
    paste0("structure: ", problem$structure$kind),
    paste0("target: ", target_label(problem$target)),
    if (has_metric) "metric/operator B supplied" else "standard eigenproblem",
    if (has_shift) "shift-invert transform requested" else NULL,
    preconditioner_reason,
    operator_kernel_reason(problem$A)
  )

  fallback <- if (grepl("Hermitian Lanczos", chosen, fixed = TRUE) ||
    grepl("LOBPCG", chosen, fixed = TRUE)) {
    "dense oracle prototype if unsupported"
  } else {
    "dense oracle prototype"
  }
  new_plan(problem, k = k, method = chosen, reasons = reasons, fallback = fallback)
}

#' @export
plan_solver.eigencore_svd_problem <- function(problem, rank, method = auto(), ...) {
  source_matrix <- source_or_null(problem$A)
  is_dense_source <- is.matrix(source_matrix) && is.double(source_matrix)
  is_native_csc <- identical(problem$A$metadata$storage, "dgCMatrix")
  chosen <- if (inherits(method, "eigencore_method") && method$kind != "auto") {
    if (identical(method$kind, "golub_kahan")) {
      if (is_native_csc || is_dense_source) "native prototype Golub-Kahan" else "prototype Golub-Kahan"
    } else if (identical(method$kind, "randomized")) {
      "reference randomized SVD prototype"
    } else {
      method_label(method)
    }
  } else if (is_native_csc && !is.null(problem$A$apply_adjoint)) {
    "native prototype Golub-Kahan"
  } else if (is.null(source_or_null(problem$A)) && !is.null(problem$A$apply_adjoint)) {
    "prototype Golub-Kahan"
  } else if (is_dense_source) {
    "native dense LAPACK SVD fallback"
  } else {
    "dense LAPACK SVD oracle (prototype fallback)"
  }
  reasons <- c(
    paste0("target: ", target_label(problem$target)),
    "rectangular SVD problem",
    if (!is.null(problem$A$apply_adjoint)) "adjoint is available" else "adjoint is missing",
    "default avoids normal equations",
    operator_kernel_reason(problem$A)
  )
  fallback <- if (grepl("prototype Golub-Kahan", chosen, fixed = TRUE)) {
    "dense oracle prototype if unsupported"
  } else {
    "dense oracle prototype"
  }
  new_plan(problem, k = rank, method = chosen, reasons = reasons, fallback = fallback)
}

#' @keywords internal
operator_kernel_reason <- function(op) {
  storage <- op$metadata$storage %||% NULL
  if (identical(storage, "dgCMatrix")) {
    "built-in sparse CSC operator has native block apply"
  } else if (identical(storage, "ddiMatrix")) {
    "built-in diagonal operator has native block apply"
  } else if (isTRUE(op$metadata$native)) {
    "built-in dense operator has native block apply"
  } else {
    "operator uses R-level apply path in current prototype"
  }
}

#' @keywords internal
new_plan <- function(problem, k, method, reasons, fallback = "dense oracle prototype") {
  plan <- list(
    problem_type = problem$type,
    requested = k,
    method = method,
    target = target_label(problem$target),
    reasons = reasons,
    fallback = fallback
  )
  class(plan) <- "eigencore_plan"
  plan
}

#' @export
print.eigencore_plan <- function(x, ...) {
  cat("eigencore solver plan\n")
  cat("  problem:", x$problem_type, "\n")
  cat("  requested:", x$requested, "\n")
  cat("  target:", x$target, "\n")
  cat("  method:", x$method, "\n")
  cat("  reasons:\n")
  for (reason in x$reasons) {
    cat("   -", reason, "\n")
  }
  cat("  fallback:", x$fallback, "\n")
  invisible(x)
}
