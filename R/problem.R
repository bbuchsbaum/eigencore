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
  constraints_reason <- if (inherits(method, "eigencore_method") &&
    identical(method$kind, "lobpcg")) {
    lobpcg_constraints_plan_reason(method$constraints)
  } else {
    NULL
  }
  native_lanczos_supported <- is_hermitian &&
    native_lanczos_target_supported(problem$target) &&
    (is_native_csc || is_dense_source)
  lanczos_block <- if (inherits(method, "eigencore_method") &&
    identical(method$kind, "lanczos")) method$block %||% 1L else 1L

  chosen <- if (inherits(method, "eigencore_method") &&
                identical(method$kind, "shift_invert")) {
    shift_invert_plan_label(problem, has_metric, is_hermitian,
                             is_dense_source, is_native_csc)
  } else if (inherits(method, "eigencore_method") && method$kind != "auto") {
    if (identical(method$kind, "lanczos")) {
      if (!is_hermitian) {
        "dense LAPACK eigen oracle (Lanczos requires Hermitian structure)"
      } else if (native_lanczos_supported && lanczos_block > 1L) {
        "native block Hermitian Lanczos (thick restart, locking)"
      } else if (native_lanczos_supported) {
        "native scalar thick-restart Hermitian Lanczos"
      } else {
        "reference Hermitian Lanczos (prototype/oracle fallback)"
      }
    } else if (identical(method$kind, "lobpcg")) {
      if (native_lobpcg_supported(
          problem$A,
          target = problem$target,
          preconditioner = method$preconditioner,
          Bop = problem$metric,
          constraints = method$constraints
      )) {
        "native standard Hermitian LOBPCG prototype"
      } else if (is_hermitian && !is.null(problem$metric) &&
        native_generalized_lobpcg_supported(
          problem$A,
          problem$metric,
          target = problem$target,
          preconditioner = method$preconditioner,
          constraints = method$constraints
        )) {
        native_generalized_lobpcg_label()
      } else if (is_hermitian && is.null(problem$metric)) {
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
    shift_invert_plan_label(problem, has_metric, is_hermitian,
                             is_dense_source, is_native_csc)
  } else if (has_metric && should_auto_native_generalized_lobpcg(problem, k)) {
    native_generalized_lobpcg_label()
  } else if (has_metric && should_auto_reference_generalized_lobpcg(problem, k)) {
    reference_generalized_lobpcg_label()
  } else if (has_metric && is_hermitian && is_dense_source && is_dense_metric) {
    "native dense generalized SPD LAPACK fallback"
  } else if (has_metric && is_hermitian) {
    "dense generalized SPD LAPACK oracle (prototype fallback)"
  } else if (!is.null(promoted_block_lanczos_controls(problem, k))) {
    "native block Hermitian Lanczos (thick restart, locking)"
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
    constraints_reason,
    operator_kernel_reason(problem$A)
  )

  fallback <- if (grepl("Hermitian Lanczos", chosen, fixed = TRUE) ||
    grepl("LOBPCG", chosen, fixed = TRUE)) {
    "dense oracle prototype if unsupported"
  } else {
    "dense oracle prototype"
  }
  controls <- lanczos_plan_controls(problem, k = k, method = method, chosen = chosen)
  if (grepl("LOBPCG", chosen, fixed = TRUE)) {
    controls <- lobpcg_plan_controls(method)
  }
  new_plan(
    problem,
    k = k,
    method = chosen,
    reasons = reasons,
    fallback = fallback,
    controls = controls
  )
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
  } else if (should_use_native_gram_svd(problem, method, rank = rank)) {
    "native certified Gram SVD special case"
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
    if (identical(chosen, "native certified Gram SVD special case")) {
      "small rectangular sparse problem: materializes the smaller Gram matrix as an explicit certified special case"
    } else {
      "default avoids normal equations"
    },
    operator_kernel_reason(problem$A)
  )
  fallback <- if (identical(chosen, "native certified Gram SVD special case")) {
    "native Golub-Kahan if Gram special case is disabled or uncertified"
  } else if (grepl("prototype Golub-Kahan", chosen, fixed = TRUE)) {
    "dense oracle prototype if unsupported"
  } else {
    "dense oracle prototype"
  }
  controls <- svd_plan_controls(problem, rank = rank, method = method, chosen = chosen)
  new_plan(
    problem,
    k = rank,
    method = chosen,
    reasons = reasons,
    fallback = fallback,
    controls = controls
  )
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
new_plan <- function(problem, k, method, reasons, fallback = "dense oracle prototype",
                     controls = list()) {
  plan <- list(
    problem_type = problem$type,
    requested = k,
    method = method,
    target = target_label(problem$target),
    reasons = reasons,
    fallback = fallback,
    controls = controls
  )
  class(plan) <- "eigencore_plan"
  plan
}

#' @keywords internal
lanczos_plan_controls <- function(problem, k, method, chosen) {
  if (!grepl("Lanczos", chosen, fixed = TRUE)) {
    return(list())
  }
  is_lanczos_method <- inherits(method, "eigencore_method") &&
    identical(method$kind, "lanczos")
  promoted <- if (is_lanczos_method) NULL else promoted_block_lanczos_controls(problem, k)
  block <- if (is_lanczos_method) {
    method$block %||% 1L
  } else {
    promoted$block %||% 1L
  }
  block <- as.integer(block)
  n <- problem$A$dim[1L]
  max_restarts <- if (is_lanczos_method) method$max_restarts else NULL
  max_restarts <- as.integer(max_restarts %||% 100L)
  max_subspace <- if (is_lanczos_method) method$max_subspace else NULL
  max_subspace <- if (is.null(max_subspace)) {
    if (!is.null(promoted$max_subspace)) {
      promoted$max_subspace
    } else if (block > 1L) {
      source <- source_or_null(problem$A)
      dense_full_n <- as.integer(getOption("eigencore.block_dense_full_subspace_max_n", 256L))
      if (length(dense_full_n) != 1L || is.na(dense_full_n) || dense_full_n < 1L) {
        dense_full_n <- 256L
      }
      if (is.matrix(source) && is.double(source) && n <= dense_full_n) {
        n
      } else {
        default_block_lanczos_max_subspace(k, block)
      }
    } else {
      max(as.integer(k) + 1L, 3L * as.integer(k) + 20L)
    }
  } else {
    as.integer(max_subspace)
  }
  max_subspace <- min(n, max_subspace)
  reorthogonalize <- if (is_lanczos_method) {
    isTRUE(method$reorthogonalize)
  } else {
    TRUE
  }
  controls <- list(
    block = block,
    max_subspace = max_subspace,
    max_restarts = max_restarts,
    reorthogonalize = reorthogonalize
  )
  if (inherits(problem$transform, "eigencore_method") &&
      identical(problem$transform$kind, "shift_invert")) {
    controls <- c(controls, list(
      transform = "shift_invert",
      transformed_operator_target = "largest_magnitude",
      eigenvalue_recovery = "lambda = sigma + 1 / mu",
      certified_in_original_coordinates = TRUE,
      certification_policy = "residual certificate on original eigenproblem"
    ))
  }
  controls
}

#' @keywords internal
lobpcg_plan_controls <- function(method) {
  is_lobpcg_method <- inherits(method, "eigencore_method") &&
    identical(method$kind, "lobpcg")
  maxit <- if (is_lobpcg_method) {
    method$maxit
  } else {
    getOption("eigencore.lobpcg_maxit", 200L)
  }
  maxit <- as.integer(maxit)
  if (length(maxit) != 1L || is.na(maxit) || maxit < 1L) {
    maxit <- 200L
  }
  list(maxit = maxit)
}

#' @keywords internal
promoted_block_lanczos_controls <- function(problem, k) {
  if (is.null(k) || is.na(k) || length(k) != 1L) {
    return(NULL)
  }
  k <- as.integer(k)
  if (k < 16L || !is.null(problem$metric) ||
      !identical(problem$structure$kind, "hermitian") ||
      !native_lanczos_target_supported(problem$target)) {
    return(NULL)
  }

  n <- as.integer(problem$A$dim[1L])
  source <- source_or_null(problem$A)
  storage <- problem$A$metadata$storage %||% NULL
  dense_full_n <- as.integer(getOption("eigencore.block_dense_full_subspace_max_n", 256L))
  if (length(dense_full_n) != 1L || is.na(dense_full_n) || dense_full_n < 1L) {
    dense_full_n <- 256L
  }

  dense_max_fraction <- as.numeric(getOption("eigencore.dense_partial_lanczos_max_fraction", 0.25))
  if (length(dense_max_fraction) != 1L || is.na(dense_max_fraction) ||
      dense_max_fraction <= 0 || dense_max_fraction > 1) {
    dense_max_fraction <- 0.25
  }

  if (is.matrix(source) && is.double(source) && n <= dense_full_n &&
      (k / n) <= dense_max_fraction) {
    return(list(block = 2L, max_subspace = n))
  }
  if (identical(storage, "dgCMatrix") && n >= 5000L) {
    block <- 4L
    return(list(
      block = block,
      max_subspace = min(n, max(default_block_lanczos_max_subspace(k, block), 16L * k))
    ))
  }
  if (identical(storage, "dgCMatrix") && n >= 500L) {
    block <- 2L
    return(list(
      block = block,
      max_subspace = min(n, default_block_lanczos_max_subspace(k, block))
    ))
  }
  NULL
}

#' @keywords internal
svd_plan_controls <- function(problem, rank, method, chosen) {
  dims <- as.integer(problem$A$dim)
  rank <- min(as.integer(rank), min(dims))
  if (inherits(method, "eigencore_method") && identical(method$kind, "randomized")) {
    oversample <- max(0L, as.integer(method$oversample %||% 10L))
    n_iter <- max(0L, as.integer(method$n_iter %||% 2L))
    normalizer <- method$normalizer %||% "qr"
    if (!normalizer %in% c("qr", "lu", "none")) {
      normalizer <- "qr"
    }
    sample_dimension <- min(min(dims), rank + oversample)
    return(list(
      oversample = oversample,
      n_iter = n_iter,
      sample_dimension = sample_dimension,
      normalizer = normalizer,
      approximate = TRUE,
      auto_selected = FALSE,
      refine = isTRUE(method$refine),
      certification_policy = if (isTRUE(method$refine)) {
        "residual certificate, with deterministic refinement when needed"
      } else {
        "residual certificate only; stochastic sketch is not sufficient to pass"
      },
      certification_refinement = "native Gram SVD when available and randomized certificate fails"
    ))
  }

  if (identical(chosen, "native certified Gram SVD special case")) {
    return(list(
      gram_side = if (dims[2L] <= dims[1L]) "right" else "left",
      gram_dimension = min(dims),
      gram_max_dimension = as.integer(getOption("eigencore.gram_svd_max_dimension", 512L)),
      rank_fraction_limit = 0.5,
      certified_in_original_coordinates = TRUE,
      materializes = "smaller Gram matrix only",
      fallback_policy = "certification-gated",
      runtime_fallback = "native Golub-Kahan if original-coordinate certificate is weaker",
      fallback_requires_vectors = "both"
    ))
  }

  if (grepl("Golub-Kahan", chosen, fixed = TRUE)) {
    is_gk <- inherits(method, "eigencore_method") && identical(method$kind, "golub_kahan")
    requested_max_subspace <- if (is_gk) method$max_subspace else NULL
    default_initial_subspace <- default_golub_kahan_initial_subspace(dims, rank)
    initial_max_subspace <- if (is.null(requested_max_subspace)) {
      default_initial_subspace
    } else {
      min(min(dims), as.integer(requested_max_subspace))
    }
    controls <- list(
      adaptive_subspace = is.null(requested_max_subspace),
      initial_max_subspace = initial_max_subspace,
      reorthogonalize = if (is_gk) isTRUE(method$reorthogonalize) else TRUE,
      requires_adjoint = TRUE,
      default_normal_equations = FALSE
    )
    if (!is.null(requested_max_subspace)) {
      controls$max_subspace <- min(min(dims), as.integer(requested_max_subspace))
    }
    return(controls)
  }

  if (grepl("dense", chosen, fixed = TRUE)) {
    return(list(
      dense_fallback_mb = getOption("eigencore.dense_fallback_mb", 256),
      certification = "original dense coordinates"
    ))
  }

  list()
}

#' @keywords internal
default_golub_kahan_initial_subspace <- function(dims, rank) {
  dims <- as.integer(dims)
  rank <- as.integer(rank)
  limit <- min(dims)
  base <- max(rank + 1L, 4L * rank + 20L)
  if (dims[[2L]] > dims[[1L]]) {
    base <- max(base, 8L * rank + 20L)
  }
  min(limit, as.integer(base))
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
  if (length(x$controls)) {
    cat("  controls:\n")
    for (nm in names(x$controls)) {
      cat("   -", nm, ":", x$controls[[nm]], "\n")
    }
  }
  cat("  fallback:", x$fallback, "\n")
  invisible(x)
}
