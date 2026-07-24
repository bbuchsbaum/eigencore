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
solve_eigen_lanczos <- function(a, k, method, tol, maxit, vectors, certify, plan,
                                initial_subspace = NULL) {
  controls <- plan$controls %||% list()
  method_maxit <- controls$max_subspace %||%
    if (inherits(method, "eigencore_method") && identical(method$kind, "lanczos")) method$max_subspace else NULL
  method_reorth <- controls$reorthogonalize %||%
    if (inherits(method, "eigencore_method") && identical(method$kind, "lanczos")) method$reorthogonalize else TRUE
  method_max_restarts <- controls$max_restarts %||%
    if (inherits(method, "eigencore_method") && identical(method$kind, "lanczos")) method$max_restarts else NULL
  method_block <- controls$block %||%
    if (inherits(method, "eigencore_method") && identical(method$kind, "lanczos")) method$block else 1L
  method_check_stride <- controls$check_stride %||%
    if (inherits(method, "eigencore_method") && identical(method$kind, "lanczos")) method$check_stride %||% 0L else 0L
  method_check_stride <- as.integer(method_check_stride %||% 0L)

  # Warm start: fit a user-supplied initial_subspace to the method's start
  # block at the solver boundary. The plan-support guard in
  # solve.eigencore_eigen_problem() has already rejected unsupported paths, so
  # a non-NULL start here is guaranteed to reach a standard Hermitian Lanczos
  # dispatch (native dense/CSC, native block matrix-free callback, or scalar
  # matrix-free reference); the defensive stop keeps that invariant local.
  warm <- !is.null(initial_subspace)
  start_block <- NULL
  start_provenance <- warm_start_cold_provenance()
  if (warm) {
    rng_state_before_prepare <- if (exists(".Random.seed", envir = .GlobalEnv,
                                           inherits = FALSE)) {
      get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    } else {
      NULL
    }
    native_warm_path <- plan_dispatches_native_lanczos(plan)
    if (!native_warm_path && !plan_dispatches_reference_hermitian_lanczos(plan)) {
      stop("Internal error: initial_subspace reached an unsupported Lanczos dispatch.",
           call. = FALSE)
    }
    start_width <- if (native_warm_path) as.integer(method_block) else 1L
    prep <- prepare_initial_subspace(initial_subspace, n = a$A$dim[1L],
                                     width = start_width)
    start_block <- prep$start
    start_provenance <- prep[c("start_source", "supplied", "accepted",
                               "rejected", "augmented", "rank", "compressed",
                               "invariant_guard_used",
                               "invariant_relative_residual",
                               "guard_operator_block_calls",
                               "guard_operator_columns")]
    if (prep$augmented == 0L && prep$rank > 0L) {
      guard <- warm_start_invariant_guard(a$A, prep$accepted_basis, tol = tol)
      start_provenance$invariant_guard_used <- TRUE
      start_provenance$invariant_relative_residual <- guard$relative_residual
      start_provenance$guard_operator_block_calls <- guard$operator_block_calls
      start_provenance$guard_operator_columns <- guard$operator_columns
      if (isTRUE(guard$discard)) {
        start_block <- NULL
        if (!is.null(rng_state_before_prepare)) {
          assign(".Random.seed", rng_state_before_prepare, envir = .GlobalEnv)
        }
        start_provenance$start_source <-
          "user_supplied_discarded_invariant_guard"
      }
    }
  }

  native_generalized_path <- identical(plan$method, native_generalized_lanczos_label())
  reference_generalized_path <- identical(plan$method, generalized_lanczos_label())
  iter <- if (native_generalized_path) {
    native_generalized_lanczos_hermitian(
      a$A,
      a$metric,
      k = k,
      target = a$target,
      tol = tol,
      maxit = maxit %||% method_maxit,
      vectors = vectors,
      block = method_block,
      max_restarts = method_max_restarts
    )
  } else if (reference_generalized_path) {
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
        vectors = vectors,
        # A warm start must actually be iterated: disable the exact
        # full-subspace dsyev shortcut so the supplied block is consumed
        # rather than silently ignored.
        full_subspace = !warm,
        start = start_block,
        check_stride = method_check_stride
      )
    } else {
      native_lanczos_hermitian(
        a$A,
        k = k,
        target = a$target,
        tol = tol,
        maxit = maxit %||% method_maxit,
        max_restarts = method_max_restarts,
        vectors = vectors,
        start = start_block,
        check_stride = method_check_stride
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
      reorthogonalize = method_reorth,
      start = start_block
    )
  }

  # The certificate-rescue fallback recomputes from a fresh start, discarding
  # the supplied subspace; report that honestly rather than claiming it was used.
  if (warm && isTRUE(iter$restart$fallback_used)) {
    start_provenance$start_source <- "user_supplied_discarded_on_fallback"
  }

  guard_block_calls <- start_provenance$guard_operator_block_calls %||% 0L
  guard_columns <- start_provenance$guard_operator_columns %||% 0L
  iter$operator_block_calls <- as.integer(
    (iter$operator_block_calls %||% iter$matvecs %||% 0L) +
      guard_block_calls
  )
  iter$operator_columns <- as.integer(
    (iter$operator_columns %||% iter$matvecs %||% 0L) + guard_columns
  )
  iter$certification_operator_columns <- as.integer(
    iter$certification_operator_columns %||% 0L
  )
  plan$initial_subspace <- start_provenance
  plan$controls$initial_subspace_supported <- TRUE

  warning_msg <- if (native_generalized_path) {
    if (!isTRUE(iter$certificate$passed)) {
      paste0(
        plan$method, " exhausted its current transformed subspace before all ",
        k, " requested generalized pairs converged"
      )
    } else {
      "using native transformed generalized SPD B-orthogonal Lanczos; original generalized residuals certified"
    }
  } else if (reference_generalized_path) {
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
      operator_block_calls = iter$operator_block_calls,
      operator_columns = iter$operator_columns,
      certification_operator_columns = iter$certification_operator_columns,
      convergence_history = iter$convergence_history %||% NULL,
      restart = iter$restart %||% NULL,
      locked = iter$locked %||% which(iter$certificate$converged),
      start_source = start_provenance$start_source,
      initial_subspace = start_provenance
    )
  )
}

#' @keywords internal
solve_eigen_arnoldi <- function(a, k, method, tol, maxit, vectors, certify, plan) {
  controls <- plan$controls %||% list()
  native_path <- identical(plan$method, native_arnoldi_label()) ||
    identical(plan$method, native_refined_arnoldi_label()) ||
    identical(plan$method, native_matrix_free_arnoldi_label())
  refined_native_path <- identical(plan$method, native_refined_arnoldi_label())
  matrix_free_native_path <- identical(plan$method, native_matrix_free_arnoldi_label())
  default_maxit <- if (native_path) {
    native_arnoldi_default_max_subspace(a$A$dim[[1L]], k)
  } else {
    max(k + 8L, 2L * k + 4L)
  }
  method_maxit <- maxit %||% controls$max_subspace %||% default_maxit
  method_max_restarts <- controls$max_restarts %||% 0L
  method_extraction <- controls$arnoldi_extraction %||%
    if (refined_native_path) "refined_ritz" else "projected_ritz"
  arnoldi_solver <- if (native_path) native_arnoldi_general else reference_arnoldi_general
  iter <- arnoldi_solver(
    a$A,
    k = k,
    target = a$target,
    tol = tol,
    maxit = method_maxit,
    max_restarts = method_max_restarts,
    vectors = vectors,
    extraction = method_extraction
  )
  left_contract <- arnoldi_left_eigen_contract(
    a$A,
    iter$values,
    iter$vectors,
    target = a$target,
    tol = tol,
    maxit = method_maxit,
    max_restarts = method_max_restarts,
    extraction = method_extraction
  )
  warning_msg <- if (matrix_free_native_path && isTRUE(iter$certificate$passed)) {
    "using native matrix-free Arnoldi callback cycle with native Ritz extraction; right residuals certified"
  } else if (matrix_free_native_path) {
    "using native matrix-free Arnoldi callback cycle with native Ritz extraction; result did not pass certificate"
  } else if (refined_native_path && isTRUE(iter$certificate$passed)) {
    "using native Arnoldi cycle with native refined Ritz extraction; right residuals certified"
  } else if (refined_native_path) {
    "using native Arnoldi cycle with native refined Ritz extraction; result did not pass certificate"
  } else if (native_path && isTRUE(iter$certificate$passed)) {
    "using native Arnoldi cycle with native Ritz extraction; right residuals certified"
  } else if (native_path) {
    "using native Arnoldi cycle with native Ritz extraction; result did not pass certificate"
  } else if (isTRUE(iter$certificate$passed)) {
    "using R-level reference Arnoldi prototype; native nonsymmetric Arnoldi not yet implemented"
  } else {
    "using R-level reference Arnoldi prototype; result did not pass certificate"
  }
  if (isTRUE(left_contract$supported) && isTRUE(left_contract$certificate$passed)) {
    warning_msg <- paste(
      warning_msg,
      "left eigenvectors computed from adjoint Arnoldi; left residuals and biorthogonality certified",
      sep = "; "
    )
  } else if (isTRUE(left_contract$supported)) {
    warning_msg <- paste(
      warning_msg,
      "left eigenvectors computed from adjoint Arnoldi; left residuals or biorthogonality did not pass certificate",
      sep = "; "
    )
  } else {
    warning_msg <- paste(
      warning_msg,
      paste0("left eigenvectors unavailable: ", left_contract$reason),
      sep = "; "
    )
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
      locked = which(iter$certificate$converged),
      right_vectors = iter$vectors,
      left_vectors = left_contract$vectors,
      left_certificate = left_contract$certificate,
      biorthogonality = left_contract$biorthogonality,
      left_eigenvectors = left_contract
    )
  )
}

#' @keywords internal
sparse_general_pencil_unsupported_message <- function(problem) {
  paste(
    "sparse general-pencil partial support currently requires general structure,",
    "a sparse dgCMatrix A, and a nonsingular diagonal B so eigencore can run",
    "Arnoldi on B^{-1} A and certify A * x - lambda * B * x in original coordinates.",
    "Use eig_full(A, B = ..., structure = general()) or generalized_schur(A, B)",
    "for dense full general pencils; sparse QZ and arbitrary sparse/factorized",
    "general pencils are not implemented."
  )
}

#' @keywords internal
sparse_general_pencil_transformed_operator <- function(Aop, Bop) {
  A_matrix <- Aop$metadata$matrix %||% NULL
  if (!inherits(A_matrix, "CsparseMatrix")) {
    stop("sparse general-pencil Arnoldi requires a sparse CSC A.", call. = FALSE)
  }
  B_values <- sparse_general_pencil_diagonal_values(Bop)
  if (is.null(B_values)) {
    stop(
      "sparse general-pencil Arnoldi requires a nonsingular finite diagonal B.",
      call. = FALSE
    )
  }
  if (length(B_values) != Aop$dim[1L]) {
    stop("A and B must have compatible dimensions.", call. = FALSE)
  }
  transformed <- Matrix::Diagonal(x = 1 / B_values) %*%
    methods::as(A_matrix, "generalMatrix")
  transformed <- methods::as(transformed, "CsparseMatrix")
  Cop <- as_operator(transformed)
  Cop$name <- "B_inverse_A_sparse_general_pencil"
  Cop$metadata$general_pencil_transform <- list(
    kind = "diagonal_left_scaling",
    transformed_operator = "B^{-1} A",
    materialized_sparse_operator = TRUE,
    materialized_dense_operator = FALSE,
    metric_solve = "nonsingular diagonal B row scaling",
    min_abs_B_diagonal = min(abs(B_values)),
    max_abs_B_diagonal = max(abs(B_values))
  )
  Cop
}

#' @keywords internal
solve_eigen_sparse_general_pencil_arnoldi <- function(a, k, method, tol, maxit,
                                                      vectors, certify, plan) {
  controls <- plan$controls %||% list()
  Cop <- sparse_general_pencil_transformed_operator(a$A, a$metric)
  method_maxit <- maxit %||% controls$max_subspace %||%
    sparse_general_pencil_default_max_subspace(a$A$dim[[1L]], k)
  method_max_restarts <- controls$max_restarts %||% 5L
  method_extraction <- controls$arnoldi_extraction %||% "refined_ritz"

  iter <- native_arnoldi_general(
    Cop,
    k = k,
    target = a$target,
    tol = tol,
    maxit = method_maxit,
    max_restarts = method_max_restarts,
    vectors = TRUE,
    extraction = method_extraction
  )
  vals <- iter$values
  vecs_for_cert <- iter$vectors
  alpha <- vals
  beta <- rep(1, length(vals))
  pencil <- generalized_pencil_values(alpha, beta)
  cert <- if (isTRUE(certify) && !is.null(vecs_for_cert) &&
      ncol(vecs_for_cert) > 0L) {
    certify_generalized_pencil_operator(
      a$A,
      a$metric,
      alpha,
      beta,
      vecs_for_cert,
      tol = tol
    )
  } else {
    empty_certificate(
      tol,
      note = if (!isTRUE(certify)) {
        "sparse general-pencil Arnoldi: certification disabled by caller"
      } else {
        "sparse general-pencil Arnoldi: no eigenpairs returned; residual certificate not computed"
      }
    )
  }

  restart <- iter$restart
  restart$kind <- "native_transformed_sparse_general_pencil_arnoldi"
  restart$generalized <- TRUE
  restart$right_hand_pencil <- TRUE
  restart$native <- TRUE
  restart$transformed_operator <- "B^{-1} A"
  restart$transformed_operator_storage <- Cop$metadata$storage %||% NA_character_
  restart$metric_solve <- "nonsingular diagonal B row scaling"
  restart$certificate_problem <- "original_generalized_pencil"
  restart$materialized_dense_operator <- FALSE
  restart$materialized_sparse_operator <- TRUE
  restart$general_pencil_transform <- Cop$metadata$general_pencil_transform

  warning_msg <- if (isTRUE(cert$passed)) {
    "using native transformed sparse general-pencil Arnoldi with diagonal B; original generalized residuals certified"
  } else {
    "using native transformed sparse general-pencil Arnoldi with diagonal B; result did not pass original generalized residual certificate"
  }

  make_eigen_result(
    values = vals,
    vectors = if (isTRUE(vectors)) vecs_for_cert else NULL,
    certificate = cert,
    iter = iter,
    requested = k,
    method_label = plan$method,
    target_label_value = target_label(a$target),
    plan = plan,
    warnings = warning_msg,
    extras = list(
      residuals = cert$residuals,
      backward_error = cert$backward_error,
      orthogonality = cert$orthogonality,
      restarts = restart$restart_count,
      restart = restart,
      locked = which(cert$converged),
      generalized = TRUE,
      right_hand_pencil = TRUE,
      alpha = pencil$alpha,
      beta = pencil$beta,
      classification = pencil$classification,
      finite = pencil$finite,
      infinite = pencil$infinite,
      undefined = pencil$undefined,
      transform = list(
        kind = "sparse_general_pencil_diagonal_B",
        transformed_operator = "B^{-1} A",
        metric_solve = "nonsingular diagonal B row scaling",
        certification = list(
          problem = "original",
          residual_formula = "A * x - lambda * B * x",
          transformed_residuals_used = FALSE
        )
      )
    )
  )
}

#' @keywords internal
solve_eigen_grid_laplacian_2d <- function(a, k, tol, vectors, certify, plan) {
  meta <- structured_grid_laplacian_2d_metadata(a$A)
  if (is.null(meta)) {
    stop("structured 2D-grid Laplacian solver requires explicit grid metadata.", call. = FALSE)
  }
  eig <- structured_grid_laplacian_2d_eigen(meta$nx, meta$ny, k)
  vals <- eig$values
  vecs_for_cert <- eig$vectors
  residual_vectors <- apply_operator(a$A, vecs_for_cert) -
    sweep(vecs_for_cert, 2L, vals, `*`)
  residuals <- col_norms(residual_vectors)
  cert <- if (isTRUE(certify)) {
    certify_eigen_operator_residuals(a$A, vals, vecs_for_cert, residuals, tol = tol)
  } else {
    empty_certificate(tol, note = "certification disabled by caller")
  }
  vecs <- if (isTRUE(vectors)) vecs_for_cert else NULL
  restart <- list(
    kind = "separable_2d_grid_laplacian",
    implemented = TRUE,
    native = FALSE,
    prototype = TRUE,
    grid_nx = meta$nx,
    grid_ny = meta$ny,
    mode_pairs = eig$mode_pairs,
    materialized_dense_operator = FALSE,
    certificate_in_original_coordinates = TRUE
  )
  make_eigen_result(
    values = vals,
    vectors = vecs,
    certificate = cert,
    iter = list(iterations = 1L, matvecs = 0L),
    requested = k,
    method_label = plan$method,
    target_label_value = target_label(a$target),
    plan = plan,
    warnings = "using diagnostic separable 2D-grid Laplacian prototype; explicit metadata only",
    extras = list(
      restart = restart,
      locked = which(cert$converged)
    )
  )
}

#' @keywords internal
solve_eigen_native_dense_hermitian <- function(a, k, tol, vectors, certify,
                                                allow_dense_fallback, plan) {
  A <- materialize_dense_fallbacks(list(A = a$A), allow = allow_dense_fallback)$A
  target_kind <- if (inherits(a$target, "eigencore_target")) {
    a$target$kind
  } else {
    "largest"
  }
  selected_range <- !is.complex(A) && k < nrow(A) &&
    target_kind %in% c("largest", "smallest")
  eig <- if (is.complex(A)) {
    native_dense_complex_hermitian_eigen(A)
  } else if (selected_range) {
    native_dense_symmetric_eigen_selected(A, k, a$target)
  } else {
    native_dense_symmetric_eigen(A)
  }
  if (selected_range) {
    vals <- eig$values
    vecs <- if (vectors) eig$vectors else NULL
  } else {
    idx <- order_indices(eig$values, a$target)
    idx <- idx[seq_len(min(k, length(idx)))]
    vals <- eig$values[idx]
    vecs <- if (vectors) eig$vectors[, idx, drop = FALSE] else NULL
  }
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
    warnings = if (identical(plan$method, native_dense_complex_hermitian_label())) {
      "using native dense complex Hermitian LAPACK fallback; iterative engine not yet implemented"
    } else if (selected_range) {
      "using native dense Hermitian LAPACK fallback (selected-range dsyevr); iterative engine not yet implemented"
    } else {
      "using native dense Hermitian LAPACK fallback; iterative engine not yet implemented"
    },
    extras = list(
      restart = list(
        kind = "dense_hermitian_lapack",
        implemented = TRUE,
        native = TRUE,
        eigensolver = if (is.complex(A)) {
          "native_dense_complex_hermitian"
        } else if (selected_range) {
          "lapack_dsyevr_selected"
        } else {
          "lapack_dsyev_full"
        },
        selected_range = selected_range,
        selected_count = if (selected_range) length(vals) else NA_integer_,
        full_dimension = nrow(A),
        materialized_dense_operator = TRUE,
        certified_in_original_coordinates = isTRUE(certify) && !is.null(vecs)
      )
    )
  )
}

#' @keywords internal
solve_eigen_native_dense_general <- function(a, k, tol, vectors, certify,
                                             allow_dense_fallback, plan) {
  A <- materialize_dense_fallbacks(list(A = a$A), allow = allow_dense_fallback)$A
  if (!is.complex(A)) {
    stop("native dense complex general eigensolver requires a complex dense matrix.", call. = FALSE)
  }
  eig <- native_dense_complex_general_eigen(A)
  idx <- order_indices(eig$values, a$target)
  idx <- idx[seq_len(min(k, length(idx)))]
  vals <- eig$values[idx]
  vecs <- if (vectors) eig$vectors[, idx, drop = FALSE] else NULL
  cert <- if (certify && !is.null(vecs)) {
    certify_dense_general_eigen(A, vals, vecs, tol = tol)
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
    warnings = if (isTRUE(certify) && !is.null(vecs)) {
      "using native dense complex general LAPACK fallback; right residuals certified"
    } else {
      "using native dense complex general LAPACK fallback"
    }
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
