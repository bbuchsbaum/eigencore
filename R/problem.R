#' Define an eigenproblem.
#'
#' @param A Matrix or operator defining the linear map.
#' @param metric Optional metric operator for generalized eigenproblems.
#' @param structure Optional structure descriptor; defaults to the operator
#'   structure.
#' @param target Eigencore target descriptor.
#' @param transform Optional transform method such as [shift_invert()].
#' @return An `eigencore_eigen_problem` object containing the operator, optional
#'   metric, structure, target, and transform metadata consumed by
#'   [plan_solver()] and [solve()].
#' @examples
#' A <- diag(c(4, 3, 2, 1))
#' P <- eigen_problem(A, target = largest())
#' fit <- solve(P, k = 2)
#' values(fit)
eigen_problem <- function(A, metric = NULL, structure = NULL, target = largest(),
                          transform = NULL) {
  Aop <- as_operator(A)
  if (is.null(structure)) {
    structure <- Aop$structure
  }
  Bop <- if (is.null(metric)) NULL else as_operator(metric)
  if (!is.null(metric)) {
    validate_metric_symmetric(metric, Bop, structure = structure)
  }
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

#' @keywords internal
#' Enforce the Hermitian metric contract.
#'
#' Hermitian generalized eigenproblems interpret `metric =` as an
#' SPD/Hermitian metric. General-structure problems may use `metric =` as the
#' right-hand pencil only when a later planner branch explicitly accepts that
#' sparse partial boundary; unsupported cases fail before any SPD fallback.
validate_metric_symmetric <- function(metric, Bop, structure) {
  if (!identical(structure$kind, "hermitian")) {
    return(invisible(TRUE))
  }
  src <- if (is.matrix(metric) || inherits(metric, "Matrix")) {
    metric
  } else {
    source_or_null(Bop) %||% Bop$metadata$matrix %||% NULL
  }
  if (is.null(src) || !(is.matrix(src) || inherits(src, "Matrix"))) {
    return(invisible(TRUE))
  }
  if (!is_square_symmetric(src)) {
    stop(
      "metric B must be symmetric (real) or Hermitian (complex) for a ",
      "generalized eigenproblem; nonsymmetric or general right-hand pencils ",
      "are not accepted through Hermitian metric=. Use structure = general() ",
      "for supported sparse partial pencil boundaries or eig_full(A, B = ...) ",
      "for general dense pencils.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Define an SVD problem.
#'
#' @param A Matrix or operator defining the rectangular linear map.
#' @param domain Optional domain-space descriptor.
#' @param codomain Optional codomain-space descriptor.
#' @param target Eigencore singular-value target descriptor.
#' @return An `eigencore_svd_problem` object containing the operator, domain and
#'   codomain descriptors, and singular-value target consumed by [plan_solver()]
#'   and [solve()].
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(40), 8, 5)
#' S <- svd_problem(X, target = largest())
#' fit <- solve(S, rank = 2)
#' values(fit)
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
#'
#' @param problem Eigencore eigen or SVD problem object.
#' @param ... Additional planning arguments passed to methods.
#' @return An `eigencore_plan` list describing the requested problem, chosen
#'   method label, target, planner reasons, fallback label, and control
#'   metadata used by solver dispatch.
#' @examples
#' A <- diag(c(4, 3, 2, 1))
#' plan <- plan_solver(eigen_problem(A), k = 2)
#' plan$method
#' plan$reasons
plan_solver <- function(problem, ...) {
  UseMethod("plan_solver")
}

#' @export
plan_solver.eigencore_eigen_problem <- function(problem, k, method = auto(), ...) {
  auto_shift <- auto_shift_invert_route(problem, method)
  problem <- auto_shift$problem
  method <- auto_shift$method

  has_metric <- !is.null(problem$metric)
  has_shift <- inherits(problem$transform, "eigencore_method") &&
    identical(problem$transform$kind, "shift_invert")
  is_hermitian <- identical(problem$structure$kind, "hermitian")
  source_matrix <- source_or_null(problem$A)
  metric_matrix <- if (has_metric) source_or_null(problem$metric) else NULL
  is_dense_source <- is.matrix(source_matrix) && is.double(source_matrix)
  is_complex_dense_source <- is.matrix(source_matrix) && is.complex(source_matrix)
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
      if (has_metric && is_hermitian &&
          native_generalized_lanczos_supported(problem$A, problem$metric, target = problem$target)) {
        native_generalized_lanczos_label()
      } else if (has_metric && is_hermitian &&
          generalized_lanczos_supported(problem$A, problem$metric, target = problem$target)) {
        generalized_lanczos_label()
      } else if (has_metric && is_hermitian) {
        "reference generalized SPD LOBPCG prototype"
      } else if (!is_hermitian) {
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
  } else if (has_metric && !is_hermitian &&
      sparse_general_pencil_diagonal_arnoldi_supported(problem)) {
    sparse_general_pencil_diagonal_arnoldi_label()
  } else if (has_metric && !is_hermitian) {
    sparse_general_pencil_unsupported_label()
  } else if (has_metric && should_auto_native_generalized_lobpcg(problem, k)) {
    native_generalized_lobpcg_label()
  } else if (has_metric && should_auto_reference_generalized_lobpcg(problem, k)) {
    reference_generalized_lobpcg_label()
  } else if (has_metric && is_hermitian && is_dense_source && is_dense_metric) {
    "native dense generalized SPD LAPACK fallback"
  } else if (has_metric && is_hermitian) {
    "dense generalized SPD LAPACK oracle (prototype fallback)"
  } else if (is_hermitian && is_complex_dense_source) {
    native_dense_complex_hermitian_label()
  } else if (structured_grid_laplacian_2d_supported(problem, method, k)) {
    structured_grid_laplacian_2d_label()
  } else if (should_use_native_tridiagonal_hermitian(problem, k)) {
    native_tridiagonal_hermitian_label()
  } else if (!is.null(promoted_block_lanczos_controls(problem, k))) {
    "native block Hermitian Lanczos (thick restart, locking)"
  } else if (is_hermitian && is_native_csc && native_lanczos_target_supported(problem$target)) {
    "native scalar thick-restart Hermitian Lanczos"
  } else if (!is_hermitian && !has_metric &&
      reference_arnoldi_target_supported(problem$target) &&
      native_arnoldi_available(problem$A) &&
      (is_native_csc || is_dense_source)) {
    native_refined_arnoldi_label()
  } else if (!is_hermitian && !has_metric &&
      reference_arnoldi_target_supported(problem$target) &&
      native_matrix_free_arnoldi_available(problem$A)) {
    native_matrix_free_arnoldi_label()
  } else if (!is_hermitian && !has_metric &&
      reference_arnoldi_target_supported(problem$target) &&
      (is_native_csc || is.null(source_or_null(problem$A)))) {
    reference_arnoldi_label()
  } else if (is_hermitian && is_native_csc) {
    "reference Hermitian Lanczos (target unsupported by native path)"
  } else if (is_hermitian && is.null(source_or_null(problem$A))) {
    "reference Hermitian Lanczos (prototype/oracle fallback)"
  } else if (auto_dense_partial_lanczos(problem, k)) {
    "native scalar thick-restart Hermitian Lanczos"
  } else if (is_hermitian && is_dense_source) {
    "native dense Hermitian LAPACK fallback"
  } else if (!is_hermitian && !has_metric && is_complex_dense_source) {
    native_dense_complex_general_label()
  } else if (is_hermitian) {
    "dense LAPACK Hermitian eigen oracle (prototype fallback)"
  } else {
    "dense LAPACK general eigen oracle (prototype fallback)"
  }

  reasons <- c(
    paste0("structure: ", problem$structure$kind),
    paste0("target: ", target_label(problem$target)),
    if (has_metric) "metric/operator B supplied" else "standard eigenproblem",
    if (auto_shift$nearest_implicit) "nearest target routed through shift_invert(sigma)" else NULL,
    if (auto_shift$tridiagonal_edge_implicit) {
      paste0(
        "tridiagonal ", target_label(problem$target),
        " target auto-routed through factorized shift_invert(sigma = ",
        format(auto_shift$sigma, digits = 8), ")"
      )
    } else {
      NULL
    },
    if (has_shift) "shift-invert transform requested" else NULL,
    preconditioner_reason,
    constraints_reason,
    operator_kernel_reason(problem$A),
    if (identical(chosen, structured_grid_laplacian_2d_label())) {
      "explicit 2D path-grid Laplacian metadata enables separable diagnostic prototype; no arbitrary sparse recognition"
    } else {
      NULL
    }
  )

  fallback <- if (identical(chosen, sparse_general_pencil_diagonal_arnoldi_label())) {
    "none; sparse general-pencil boundary is explicit"
  } else if (identical(chosen, sparse_general_pencil_unsupported_label())) {
    "none; unsupported sparse general-pencil boundary fails before dense fallback"
  } else if (grepl("Hermitian Lanczos", chosen, fixed = TRUE) ||
    grepl("LOBPCG", chosen, fixed = TRUE)) {
    "dense oracle prototype if unsupported"
  } else {
    "dense oracle prototype"
  }
  controls <- lanczos_plan_controls(problem, k = k, method = method, chosen = chosen)
  if (identical(chosen, native_arnoldi_label()) ||
      identical(chosen, native_refined_arnoldi_label()) ||
      identical(chosen, native_matrix_free_arnoldi_label()) ||
      identical(chosen, reference_arnoldi_label())) {
    controls <- arnoldi_plan_controls(problem, k = k, chosen = chosen)
  }
  if (identical(chosen, generalized_lanczos_label()) ||
      identical(chosen, native_generalized_lanczos_label())) {
    controls <- generalized_lanczos_plan_controls(problem, k = k, method = method)
  }
  if (identical(chosen, sparse_general_pencil_diagonal_arnoldi_label())) {
    controls <- sparse_general_pencil_arnoldi_plan_controls(problem, k = k)
  }
  if (identical(chosen, structured_grid_laplacian_2d_label())) {
    controls <- structured_grid_laplacian_2d_controls(problem, k = k)
  }
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

#' @keywords internal
sparse_general_pencil_diagonal_arnoldi_label <- function() {
  "native transformed sparse general-pencil Arnoldi (diagonal B)"
}

#' @keywords internal
sparse_general_pencil_unsupported_label <- function() {
  "unsupported sparse general-pencil partial solver"
}

#' @keywords internal
sparse_general_pencil_diagonal_values <- function(Bop, require_nonsingular = TRUE) {
  Bop <- as_operator(Bop)
  storage <- Bop$metadata$storage %||% NULL
  values <- NULL
  if (identical(storage, "ddiMatrix")) {
    B <- Bop$metadata$matrix
    values <- if (identical(methods::slot(B, "diag"), "U")) {
      rep(1, Bop$dim[1L])
    } else {
      methods::slot(B, "x")
    }
  } else {
    source <- source_or_null(Bop)
    if (is.matrix(source) && nrow(source) == ncol(source) &&
        isTRUE(all(source[row(source) != col(source)] == 0))) {
      values <- diag(source)
    }
  }
  if (is.null(values)) {
    return(NULL)
  }
  values <- as.numeric(values)
  if (length(values) != Bop$dim[1L] || any(!is.finite(values))) {
    return(NULL)
  }
  if (isTRUE(require_nonsingular)) {
    scale <- max(1, abs(values))
    if (any(abs(values) <= sqrt(.Machine$double.eps) * scale)) {
      return(NULL)
    }
  }
  values
}

#' @keywords internal
sparse_general_pencil_diagonal_arnoldi_supported <- function(problem) {
  if (is.null(problem$metric) || identical(problem$structure$kind, "hermitian")) {
    return(FALSE)
  }
  if (!reference_arnoldi_target_supported(problem$target)) {
    return(FALSE)
  }
  B_values <- sparse_general_pencil_diagonal_values(problem$metric)
  if (is.null(B_values)) {
    return(FALSE)
  }
  A_matrix <- problem$A$metadata$matrix %||% NULL
  inherits(A_matrix, "CsparseMatrix") &&
    identical(problem$A$metadata$storage %||% NULL, "dgCMatrix")
}

#' @keywords internal
sparse_general_pencil_arnoldi_plan_controls <- function(problem, k) {
  controls <- arnoldi_plan_controls(
    problem,
    k = k,
    chosen = native_refined_arnoldi_label()
  )
  controls$transformed_operator <- "B^{-1} A"
  controls$right_hand_pencil <- TRUE
  controls$metric_solve <- "nonsingular diagonal B row scaling"
  controls$certification_policy <- "generalized right residual A * x - lambda * B * x in original coordinates"
  controls$unsupported_sparse_qz <- TRUE
  controls$target_family <- "sparse_general_pencil_partial"
  controls$promotion_gate <- "sparse_general_pencil_partial:diagonal_B"
  controls$dense_fallback_policy <- "none; sparse general-pencil boundary is explicit"
  controls$alpha_beta_semantics <- "finite alpha/beta with beta = 1 from transformed Arnoldi"
  n <- as.integer(problem$A$dim[1L])
  controls$max_subspace <- sparse_general_pencil_default_max_subspace(n, k)
  controls
}

#' @keywords internal
is_auto_method <- function(method) {
  inherits(method, "eigencore_method") && identical(method$kind, "auto")
}

#' @keywords internal
is_nearest_target <- function(target) {
  inherits(target, "eigencore_target") && identical(target$kind, "nearest")
}

#' @keywords internal
auto_nearest_shift_invert <- function(problem, method) {
  if (!is_auto_method(method) ||
      is_transform_method(problem$transform) ||
      !is_nearest_target(problem$target)) {
    return(list(problem = problem, method = method, implicit = FALSE))
  }

  sigma <- problem$target$value
  if (length(sigma) != 1L || !is.finite(sigma)) {
    stop("nearest(sigma) auto-routing requires a single finite numeric sigma.", call. = FALSE)
  }

  transform <- shift_invert(sigma)
  problem$transform <- transform
  list(problem = problem, method = transform, implicit = TRUE)
}

#' @keywords internal
auto_shift_invert_route <- function(problem, method) {
  nearest_shift <- auto_nearest_shift_invert(problem, method)
  if (isTRUE(nearest_shift$implicit)) {
    nearest_shift$nearest_implicit <- TRUE
    nearest_shift$tridiagonal_edge_implicit <- FALSE
    nearest_shift$sigma <- nearest_shift$method$sigma
    return(nearest_shift)
  }

  edge_shift <- auto_tridiagonal_edge_shift_invert(problem, method)
  edge_shift$nearest_implicit <- FALSE
  edge_shift$tridiagonal_edge_implicit <- isTRUE(edge_shift$implicit)
  edge_shift
}

#' @keywords internal
auto_tridiagonal_edge_shift_invert <- function(problem, method) {
  if (!is_auto_method(method) ||
      is_transform_method(problem$transform) ||
      !is.null(problem$metric) ||
      !identical(problem$structure$kind, "hermitian")) {
    return(list(problem = problem, method = method, implicit = FALSE, sigma = NULL))
  }

  target_kind <- if (inherits(problem$target, "eigencore_target")) {
    problem$target$kind
  } else {
    "largest"
  }
  if (!target_kind %in% c("largest", "smallest")) {
    return(list(problem = problem, method = method, implicit = FALSE, sigma = NULL))
  }

  A <- problem$A$metadata$matrix %||% source_or_null(problem$A)
  if (!(inherits(A, "CsparseMatrix") || inherits(A, "diagonalMatrix"))) {
    return(list(problem = problem, method = method, implicit = FALSE, sigma = NULL))
  }
  parts <- shift_invert_tridiagonal_parts(A, shift = 0)
  if (is.null(parts)) {
    return(list(problem = problem, method = method, implicit = FALSE, sigma = NULL))
  }

  sigma <- tridiagonal_edge_shift_sigma(parts, target_kind)
  if (is.null(sigma) || length(sigma) != 1L || !is.finite(sigma)) {
    return(list(problem = problem, method = method, implicit = FALSE, sigma = NULL))
  }

  transform <- shift_invert(sigma)
  problem$transform <- transform
  list(problem = problem, method = transform, implicit = TRUE, sigma = sigma)
}

#' @keywords internal
tridiagonal_gershgorin_bounds <- function(parts) {
  diag <- as.numeric(parts$diag)
  n <- length(diag)
  radius <- numeric(n)
  if (n > 1L) {
    offdiag <- abs(as.numeric(parts$upper))
    radius[seq_len(n - 1L)] <- radius[seq_len(n - 1L)] + offdiag
    radius[seq.int(2L, n)] <- radius[seq.int(2L, n)] + offdiag
  }
  lower <- min(diag - radius)
  upper <- max(diag + radius)
  list(lower = lower, upper = upper)
}

#' @keywords internal
tridiagonal_shift_pivot_ratio <- function(parts, sigma) {
  diag <- as.numeric(parts$diag) - as.numeric(sigma)
  lower <- as.numeric(parts$lower)
  upper <- as.numeric(parts$upper)
  n <- length(diag)
  denom <- numeric(n)
  cprime <- numeric(max(n - 1L, 0L))
  if (n < 1L || length(sigma) != 1L || !is.finite(sigma)) {
    return(NA_real_)
  }
  if (abs(diag[[1L]]) <= .Machine$double.eps) {
    return(0)
  }
  denom[[1L]] <- diag[[1L]]
  if (n > 1L) {
    cprime[[1L]] <- upper[[1L]] / denom[[1L]]
    for (i in seq.int(2L, n)) {
      denom[[i]] <- diag[[i]] - lower[[i - 1L]] * cprime[[i - 1L]]
      if (abs(denom[[i]]) <= .Machine$double.eps) {
        return(0)
      }
      if (i < n) {
        cprime[[i]] <- upper[[i]] / denom[[i]]
      }
    }
  }
  abs_denom <- abs(denom)
  max_abs <- max(abs_denom)
  if (!is.finite(max_abs) || max_abs <= 0) {
    return(NA_real_)
  }
  min(abs_denom) / max_abs
}

#' @keywords internal
tridiagonal_edge_shift_sigma <- function(parts, target_kind) {
  bounds <- tridiagonal_gershgorin_bounds(parts)
  diag <- as.numeric(parts$diag)
  offdiag <- abs(as.numeric(parts$upper))
  scale <- max(1, abs(bounds$lower), abs(bounds$upper), abs(diag), offdiag)
  span <- max(bounds$upper - bounds$lower, scale)
  margin <- max(1e-8 * span, 10 * sqrt(.Machine$double.eps) * scale)
  multipliers <- c(1, 10, 100, 1000, 10000)
  stable_enough <- function(sigma) {
    ratio <- tridiagonal_shift_pivot_ratio(parts, sigma)
    is.finite(ratio) && ratio > sqrt(.Machine$double.eps)
  }

  candidates <- if (identical(target_kind, "smallest")) {
    c(
      if (bounds$lower >= 0 && abs(bounds$lower) <= margin) 0 else numeric(),
      bounds$lower - margin * multipliers
    )
  } else if (identical(target_kind, "largest")) {
    c(
      if (bounds$upper <= 0 && abs(bounds$upper) <= margin) 0 else numeric(),
      bounds$upper + margin * multipliers
    )
  } else {
    numeric()
  }

  candidates <- unique(as.numeric(candidates[is.finite(candidates)]))
  for (sigma in candidates) {
    if (stable_enough(sigma)) {
      return(sigma)
    }
  }
  NULL
}

#' @export
plan_solver.eigencore_svd_problem <- function(problem, rank, method = auto(), ...) {
  source_matrix <- source_or_null(problem$A)
  is_dense_source <- is.matrix(source_matrix) && is.double(source_matrix)
  is_complex_dense_source <- is.matrix(source_matrix) && is.complex(source_matrix)
  is_native_csc <- identical(problem$A$metadata$storage, "dgCMatrix")
  is_native_matrix_free_gk <- native_matrix_free_golub_kahan_available(problem$A)
  is_smallest_svd_target <- svd_target_is_smallest(problem$target)
  is_interior_svd_target <- svd_target_is_interior(problem$target)
  has_exact_operator_scale <- operator_has_nonestimated_norm_provenance(problem$A)
  chosen <- if (inherits(method, "eigencore_method") && method$kind != "auto") {
    if (identical(method$kind, "golub_kahan")) {
      if (is_smallest_svd_target && is_native_csc) {
        native_smallest_golub_kahan_label()
      } else if (is_smallest_svd_target && is_native_matrix_free_gk && has_exact_operator_scale) {
        native_matrix_free_smallest_golub_kahan_label()
      } else if (is_interior_svd_target && is_native_csc) {
        native_interior_golub_kahan_label()
      } else if (is_interior_svd_target && is_native_matrix_free_gk && has_exact_operator_scale) {
        native_matrix_free_interior_golub_kahan_label()
      } else if (is_native_csc || is_dense_source) {
        "native prototype Golub-Kahan"
      } else if (is_native_matrix_free_gk) {
        native_matrix_free_golub_kahan_label()
      } else {
        "prototype Golub-Kahan"
      }
    } else if (identical(method$kind, "randomized")) {
      if (native_dense_randomized_svd_supported(problem, method)) {
        native_dense_randomized_svd_label()
      } else if (native_csc_randomized_svd_supported(problem, method)) {
        native_csc_randomized_svd_label()
      } else {
        "reference randomized SVD prototype"
      }
    } else {
      method_label(method)
    }
  } else if (should_use_native_gram_svd(problem, method, rank = rank)) {
    "native certified Gram SVD special case"
  } else if (should_use_native_retained_golub_kahan(problem, method, rank = rank)) {
    native_retained_golub_kahan_diagnostic_label()
  } else if (is_smallest_svd_target && is_native_csc && !is.null(problem$A$apply_adjoint)) {
    native_smallest_golub_kahan_label()
  } else if (is_smallest_svd_target && is_native_matrix_free_gk && has_exact_operator_scale) {
    native_matrix_free_smallest_golub_kahan_label()
  } else if (is_interior_svd_target && is_native_csc && !is.null(problem$A$apply_adjoint)) {
    native_interior_golub_kahan_label()
  } else if (is_interior_svd_target && is_native_matrix_free_gk && has_exact_operator_scale) {
    native_matrix_free_interior_golub_kahan_label()
  } else if (is_native_csc && !is.null(problem$A$apply_adjoint)) {
    "native prototype Golub-Kahan"
  } else if (is_native_matrix_free_gk) {
    native_matrix_free_golub_kahan_label()
  } else if (is.null(source_or_null(problem$A)) && !is.null(problem$A$apply_adjoint)) {
    "prototype Golub-Kahan"
  } else if (is_complex_dense_source) {
    native_dense_complex_svd_label()
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
      "bounded smaller-side normal problem with exact original-coordinate certification; explicit Gram is the production default"
    } else if (identical(chosen, native_retained_golub_kahan_diagnostic_label())) {
      "sparse explicit operator uses diagnostic native retained block Golub-Kahan with thick restart; not production-promoted"
    } else if (identical(chosen, native_matrix_free_golub_kahan_label()) ||
        identical(chosen, native_matrix_free_smallest_golub_kahan_label()) ||
        identical(chosen, native_matrix_free_interior_golub_kahan_label())) {
      "matrix-free operator uses a native Golub-Kahan callback cycle with native Ritz extraction; sparse/matrix-free performance promotion remains gated"
    } else if (identical(chosen, native_smallest_golub_kahan_label())) {
      "smallest-SVD target uses a native Golub-Kahan cycle with exact two-sided certification"
    } else if (identical(chosen, native_interior_golub_kahan_label())) {
      "nearest-SVD target uses a native full-subspace Golub-Kahan cycle without densifying the original operator"
    } else if (identical(chosen, native_dense_randomized_svd_label())) {
      "dense randomized request uses the native QR range/subspace/projected-core controller with exact residual certification"
    } else {
      "default avoids normal equations"
    },
    operator_kernel_reason(problem$A)
  )
  fallback <- if (identical(chosen, "native certified Gram SVD special case")) {
    paste(
      "opt-in implicit candidate retries explicit Gram;",
      "native Golub-Kahan if Gram remains uncertified"
    )
  } else if (identical(chosen, native_matrix_free_golub_kahan_label()) ||
      identical(chosen, native_matrix_free_smallest_golub_kahan_label()) ||
      identical(chosen, native_matrix_free_interior_golub_kahan_label())) {
    "dense oracle prototype if native callback boundary is unsupported"
  } else if (grepl("prototype Golub-Kahan", chosen, fixed = TRUE)) {
    "dense oracle prototype if unsupported"
  } else if (identical(chosen, native_retained_golub_kahan_diagnostic_label())) {
    "native prototype Golub-Kahan if retained restart is unsupported"
  } else if (identical(chosen, native_dense_randomized_svd_label())) {
    "reference randomized SVD prototype if dense native controller is unsupported"
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
arnoldi_plan_controls <- function(problem, k, chosen) {
  n <- as.integer(problem$A$dim[1L])
  k <- as.integer(k)
  native_path <- identical(chosen, native_arnoldi_label()) ||
    identical(chosen, native_refined_arnoldi_label()) ||
    identical(chosen, native_matrix_free_arnoldi_label())
  refined_native_path <- identical(chosen, native_refined_arnoldi_label())
  matrix_free_native_path <- identical(chosen, native_matrix_free_arnoldi_label())
  source_matrix <- source_or_null(problem$A)
  dense_native_path <- native_path && is.matrix(source_matrix) && is.double(source_matrix)
  default_restarts <- if (native_path) 5L else 0L
  max_restarts <- getOption("eigencore.arnoldi_max_restarts", default_restarts)
  max_restarts <- as.integer(max_restarts)
  if (length(max_restarts) != 1L || is.na(max_restarts) || max_restarts < 0L) {
    max_restarts <- default_restarts
  }
  max_subspace <- if (dense_native_path) {
    n
  } else if (native_path) {
    native_arnoldi_default_max_subspace(n, k)
  } else {
    min(n, max(k + 8L, 2L * k + 4L))
  }
  list(
    max_subspace = max_subspace,
    max_restarts = max_restarts,
    restart = if (matrix_free_native_path) {
      "native matrix-free Arnoldi callback restart budget"
    } else if (native_path) {
      "native Arnoldi cycle restart budget"
    } else {
      "reference Arnoldi restart budget"
    },
    ritz_extraction_native = native_path,
    arnoldi_extraction = if (refined_native_path) "refined_ritz" else "projected_ritz",
    refined_extraction_native = refined_native_path,
    krylov_schur = FALSE,
    krylov_schur_status = if (refined_native_path) {
      "not implemented; V2 tranche promotes native refined Ritz extraction only"
    } else {
      "not requested"
    },
    v2_issue = if (refined_native_path) "bd-01KTF6H41S9XDN286TR3V184P4" else NULL,
    certification_policy = "right residual certificate on original nonsymmetric eigenproblem"
  )
}

#' @keywords internal
operator_kernel_reason <- function(op) {
  storage <- op$metadata$storage %||% NULL
  if (identical(storage, "dgCMatrix")) {
    "built-in sparse CSC operator has native block apply"
  } else if (identical(storage, "ddiMatrix")) {
    "built-in diagonal operator has native block apply"
  } else if (identical(storage, "complex_dense_matrix")) {
    "base complex dense source has native dense LAPACK decomposition kernels and native zgemm block apply"
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
  # Sparse block auto-promotion is diagnostic-only until non-quick installed G1
  # gates are green. Explicit lanczos(block > 1) remains available.
  promote_sparse <- isTRUE(getOption("eigencore.promote_sparse_block_lanczos", FALSE))
  if (!promote_sparse) {
    return(NULL)
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
native_retained_golub_kahan_diagnostic_label <- function() {
  "diagnostic native retained Golub-Kahan SVD (thick restart; not production-promoted)"
}

#' @keywords internal
native_smallest_golub_kahan_label <- function() {
  "native certified smallest Golub-Kahan SVD"
}

#' @keywords internal
native_matrix_free_smallest_golub_kahan_label <- function() {
  "native certified smallest matrix-free Golub-Kahan callback SVD (exact-norm boundary)"
}

#' @keywords internal
native_interior_golub_kahan_label <- function() {
  "native full-subspace interior Golub-Kahan SVD"
}

#' @keywords internal
native_matrix_free_interior_golub_kahan_label <- function() {
  "native full-subspace interior matrix-free Golub-Kahan SVD (exact-norm boundary)"
}

#' @keywords internal
native_dense_randomized_svd_label <- function() {
  "native dense randomized SVD controller (QR, exact-certificate boundary)"
}

#' @keywords internal
native_csc_randomized_svd_label <- function() {
  "native sparse CSC randomized SVD controller (QR, exact-certificate boundary)"
}

#' @keywords internal
native_dense_randomized_svd_supported <- function(problem, method) {
  source_matrix <- source_or_null(problem$A)
  is_dense_source <- is.matrix(source_matrix) && is.double(source_matrix)
  kind <- svd_target_kind(problem$target)
  inherits(method, "eigencore_method") &&
    identical(method$kind, "randomized") &&
    is_dense_source &&
    identical(method$normalizer %||% "qr", "qr") &&
    kind %in% c("largest", "largest_magnitude")
}

#' @keywords internal
native_csc_randomized_svd_supported <- function(problem, method) {
  source_matrix <- source_or_null(problem$A) %||% problem$A$metadata$matrix
  is_csc_source <- inherits(source_matrix, "dgCMatrix")
  kind <- svd_target_kind(problem$target)
  inherits(method, "eigencore_method") &&
    identical(method$kind, "randomized") &&
    is_csc_source &&
    identical(method$normalizer %||% "qr", "qr") &&
    kind %in% c("largest", "largest_magnitude")
}

#' @keywords internal
svd_target_kind <- function(target) {
  if (inherits(target, "eigencore_target")) target$kind else "largest"
}

#' @keywords internal
svd_target_is_smallest <- function(target) {
  svd_target_kind(target) %in% c("smallest", "smallest_magnitude")
}

#' @keywords internal
svd_target_is_interior <- function(target) {
  identical(svd_target_kind(target), "nearest")
}

#' @keywords internal
svd_native_iterative_plan <- function(method_label) {
  method_label %in% c(
    "native prototype Golub-Kahan",
    native_smallest_golub_kahan_label(),
    native_interior_golub_kahan_label(),
    native_retained_golub_kahan_diagnostic_label(),
    native_matrix_free_golub_kahan_label(),
    native_matrix_free_smallest_golub_kahan_label(),
    native_matrix_free_interior_golub_kahan_label()
  )
}

#' @keywords internal
operator_has_nonestimated_norm_provenance <- function(op) {
  !is.null(op$metadata$frobenius_norm) || !is.null(source_or_null(op))
}

#' @keywords internal
svd_target_plan_controls <- function(problem, chosen) {
  target <- problem$target %||% largest()
  kind <- svd_target_kind(target)
  family <- if (kind %in% c("largest", "largest_magnitude")) {
    "largest"
  } else if (kind %in% c("smallest", "smallest_magnitude")) {
    "smallest"
  } else if (identical(kind, "nearest")) {
    "interior"
  } else {
    "unsupported"
  }

  boundary <- if (identical(family, "interior") && grepl("dense", chosen, fixed = TRUE)) {
    "dense exact fallback"
  } else if (identical(chosen, "native certified Gram SVD special case") &&
      identical(family, "smallest")) {
    "native certified smallest SVD production boundary"
  } else if (identical(chosen, native_interior_golub_kahan_label()) ||
      identical(chosen, native_matrix_free_interior_golub_kahan_label())) {
    "native full-subspace interior SVD boundary"
  } else if (identical(family, "interior") && svd_native_iterative_plan(chosen)) {
    "unsupported native iterative interior SVD"
  } else if (identical(family, "interior")) {
    "reference/prototype interior selection"
  } else if (identical(family, "smallest") && grepl("dense", chosen, fixed = TRUE)) {
    "dense exact fallback"
  } else if (identical(chosen, native_smallest_golub_kahan_label()) ||
      identical(chosen, native_matrix_free_smallest_golub_kahan_label())) {
    "native certified smallest SVD production boundary"
  } else if (identical(family, "smallest")) {
    "exact-certificate iterative smallest SVD boundary"
  } else {
    "largest singular-value surface"
  }

  list(
    svd_target_kind = kind,
    svd_target_family = family,
    svd_target_boundary = boundary,
    svd_target_issue = if (family %in% c("smallest", "interior")) {
      "bd-01KTE8G6RYE4RD5F6CN7SNKKC6"
    } else {
      NULL
    },
    svd_target_closed_decision_issue = if (family %in% c("smallest", "interior")) {
      "bd-01KTEH6862GB19JJWX2M3FQP6T"
    } else {
      NULL
    },
    svd_target_certificate_policy = "exact two-sided residual certificate in original coordinates"
  )
}

#' @keywords internal
svd_plan_controls <- function(problem, rank, method, chosen) {
  dims <- as.integer(problem$A$dim)
  rank <- min(as.integer(rank), min(dims))
  target_controls <- svd_target_plan_controls(problem, chosen)
  if (inherits(method, "eigencore_method") && identical(method$kind, "randomized")) {
    oversample <- max(0L, as.integer(method$oversample %||% 10L))
    n_iter <- max(0L, as.integer(method$n_iter %||% 2L))
    normalizer <- method$normalizer %||% "qr"
    if (!normalizer %in% c("qr", "lu", "none")) {
      normalizer <- "qr"
    }
    sample_dimension <- min(min(dims), rank + oversample)
    return(c(target_controls, list(
      oversample = oversample,
      n_iter = n_iter,
      sample_dimension = sample_dimension,
      normalizer = normalizer,
      approximate = TRUE,
      auto_selected = FALSE,
      refine = isTRUE(method$refine),
      randomized_controller = if (identical(chosen, native_dense_randomized_svd_label())) {
        "native_dense_qr"
      } else if (identical(chosen, native_csc_randomized_svd_label())) {
        "native_csc_qr"
      } else {
        "reference_control"
      },
      randomized_controller_native = chosen %in% c(
        native_dense_randomized_svd_label(),
        native_csc_randomized_svd_label()
      ),
      randomized_controller_issue = "bd-01KTEPZ7TP4Q3J1WA83XZH6A05",
      randomized_controller_boundary = if (identical(chosen, native_dense_randomized_svd_label())) {
        "dense QR randomized native controller with R result construction"
      } else if (identical(chosen, native_csc_randomized_svd_label())) {
        "sparse CSC QR randomized native controller with R result construction"
      } else {
        "reference-control randomized SVD prototype"
      },
      certification_policy = if (isTRUE(method$refine)) {
        "residual certificate, with deterministic refinement when needed"
      } else {
        "residual certificate only; stochastic sketch is not sufficient to pass"
      },
      certification_refinement = "native Gram SVD when available and randomized certificate fails"
    )))
  }

  if (identical(chosen, "native certified Gram SVD special case")) {
    controls <- c(
      target_controls,
      gram_svd_plan_controls(dims, rank, problem$target)
    )
    if (identical(target_controls$svd_target_family, "smallest")) {
      controls$promotion_status <- "production_smallest_gram_certified"
      controls$promotion_gate <- "post_v1_svd_smallest_surface"
      controls$promotion_gate_issue <- "bd-01KTE8G6RYE4RD5F6CN7SNKKC6"
      controls$closed_decision_issue <- "bd-01KTEH6862GB19JJWX2M3FQP6T"
      controls$promotion_requires <- c(
        "small-side Gram dimension passes explicit materialization gate",
        "exact two-sided certificate passes",
        "original sparse operator is not densified"
      )
    }
    return(controls)
  }

  if (grepl("Golub-Kahan", chosen, fixed = TRUE)) {
    is_gk <- inherits(method, "eigencore_method") && identical(method$kind, "golub_kahan")
    is_auto <- inherits(method, "eigencore_method") && identical(method$kind, "auto")
    is_native_csc <- identical(problem$A$metadata$storage %||% NULL, "dgCMatrix")
    is_native_matrix_free <- identical(chosen, native_matrix_free_golub_kahan_label()) ||
      identical(chosen, native_matrix_free_smallest_golub_kahan_label()) ||
      identical(chosen, native_matrix_free_interior_golub_kahan_label())
    requested_max_subspace <- if (is_gk) method$max_subspace else NULL
    method_reorthogonalize <- if (is_gk) {
      isTRUE(method$reorthogonalize)
    } else if (is_auto && is_native_csc && identical(chosen, "native prototype Golub-Kahan")) {
      FALSE
    } else {
      TRUE
    }
    default_initial_subspace <- default_golub_kahan_initial_subspace(
      dims,
      rank,
      reorthogonalize = method_reorthogonalize
    )
    initial_max_subspace <- if (is.null(requested_max_subspace)) {
      default_initial_subspace
    } else {
      min(min(dims), as.integer(requested_max_subspace))
    }
    controls <- c(target_controls, list(
      adaptive_subspace = is.null(requested_max_subspace),
      initial_max_subspace = initial_max_subspace,
      reorthogonalize = method_reorthogonalize,
      requires_adjoint = TRUE,
      default_normal_equations = FALSE
    ))
    if (identical(chosen, native_retained_golub_kahan_diagnostic_label())) {
      controls$retained_restart <- TRUE
      controls$thick_restart <- TRUE
      controls$restart_policy <- "native Ritz-plus-random retained restart"
      controls$cached_av_retention <- FALSE
      controls$promotion_status <- "diagnostic_only"
      controls$promotion_gate <- "post_v1_svd_hard_surface"
      controls$promotion_gate_issue <- "bd-01KTE8G6RYE4RD5F6CN7SNKKC6"
      controls$closed_decision_issue <- "bd-01KTE8J9SF16Y1832D8HQQ9KEC"
      controls$promotion_requires <- c(
        "general sparse and matrix-free rows pass strict installed gates",
        "time-to-certified-answer beats certified RSpectra/irlba baselines",
        "memory is no worse than the best certified reference",
        "exact two-sided SVD certificates pass",
        "planner label remains diagnostic until gates pass"
      )
    }
    if (is_native_matrix_free) {
      controls$matrix_free_native <- TRUE
      controls$callback_boundary <- TRUE
      controls$promotion_status <- "callback_boundary_native"
      controls$promotion_gate <- "post_v1_operator_sidecars"
      controls$promotion_gate_issue <- "bd-01KTE8G6RYE4RD5F6CN7SNKKC6"
      controls$closed_decision_issue <- "bd-01KTE8J9SF16Y1832D8HQQ9KEC"
      controls$promotion_requires <- c(
        "operator-sidecar gate passes with native label and exact two-sided certificate",
        "general sparse and matrix-free rows still require strict installed time and memory gates",
        "do not reuse this callback-boundary label as a production sparse SVD promotion"
      )
    }
    if (identical(chosen, native_smallest_golub_kahan_label())) {
      controls$promotion_status <- "production_smallest_certified"
      controls$promotion_gate <- "post_v1_svd_smallest_surface"
      controls$promotion_gate_issue <- "bd-01KTE8G6RYE4RD5F6CN7SNKKC6"
      controls$closed_decision_issue <- "bd-01KTEH6862GB19JJWX2M3FQP6T"
      controls$promotion_requires <- c(
        "exact two-sided certificate passes",
        "no sparse densification is required",
        "nearest/interior SVD remains a separate future-scope boundary"
      )
    }
    if (identical(chosen, native_interior_golub_kahan_label())) {
      controls$promotion_status <- "production_interior_full_subspace"
      controls$promotion_gate <- "post_v1_svd_interior_surface"
      controls$promotion_gate_issue <- "bd-01KTE8G6RYE4RD5F6CN7SNKKC6"
      controls$closed_decision_issue <- "bd-01KTEH6862GB19JJWX2M3FQP6T"
      controls$full_subspace_interior <- TRUE
      controls$promotion_requires <- c(
        "exact two-sided certificate passes",
        "native Golub-Kahan reaches the full smaller subspace",
        "no sparse densification is required"
      )
    }
    if (identical(chosen, native_matrix_free_smallest_golub_kahan_label())) {
      controls$promotion_status <- "production_smallest_exact_norm_callback"
      controls$promotion_gate <- "post_v1_svd_smallest_surface"
      controls$promotion_gate_issue <- "bd-01KTE8G6RYE4RD5F6CN7SNKKC6"
      controls$closed_decision_issue <- "bd-01KTEH6862GB19JJWX2M3FQP6T"
      controls$requires_nonestimated_norm_scale <- TRUE
      controls$promotion_requires <- c(
        "operator supplies exact Frobenius norm metadata",
        "exact two-sided certificate passes with scale_is_estimate = FALSE",
        "nearest/interior SVD remains a separate future-scope boundary"
      )
    }
    if (identical(chosen, native_matrix_free_interior_golub_kahan_label())) {
      controls$promotion_status <- "production_interior_exact_norm_callback_full_subspace"
      controls$promotion_gate <- "post_v1_svd_interior_surface"
      controls$promotion_gate_issue <- "bd-01KTE8G6RYE4RD5F6CN7SNKKC6"
      controls$closed_decision_issue <- "bd-01KTEH6862GB19JJWX2M3FQP6T"
      controls$requires_nonestimated_norm_scale <- TRUE
      controls$full_subspace_interior <- TRUE
      controls$promotion_requires <- c(
        "operator supplies exact Frobenius norm metadata",
        "native Golub-Kahan reaches the full smaller subspace",
        "exact two-sided certificate passes with scale_is_estimate = FALSE"
      )
    }
    if (!is.null(requested_max_subspace)) {
      controls$max_subspace <- min(min(dims), as.integer(requested_max_subspace))
    }
    return(controls)
  }

  if (grepl("dense", chosen, fixed = TRUE)) {
    return(c(target_controls, list(
      dense_fallback_mb = getOption("eigencore.dense_fallback_mb", 256),
      certification = "original dense coordinates",
      svd_dense_fallback_exact = TRUE
    )))
  }

  target_controls
}

#' @keywords internal
default_golub_kahan_initial_subspace <- function(dims, rank, reorthogonalize = TRUE) {
  dims <- as.integer(dims)
  rank <- as.integer(rank)
  limit <- min(dims)
  base <- max(rank + 1L, 4L * rank + 20L)
  if (!isTRUE(reorthogonalize)) {
    base <- max(base, 9L * rank)
  } else if (dims[[2L]] > dims[[1L]]) {
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
