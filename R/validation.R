#' Validate eigencore eigen results against a dense oracle.
#'
#' @param A Matrix or eigencore operator to validate.
#' @param k Number of eigenpairs to validate.
#' @param target Eigencore eigenvalue target descriptor.
#' @param B Optional metric matrix or operator for generalized problems.
#' @param fit Optional precomputed eigencore eigen result.
#' @param tol Validation tolerance.
#' @keywords internal
validate_eigen_accuracy <- function(A, k, target = largest(), B = NULL,
                                    fit = NULL, tol = 1e-8) {
  if (is.null(fit)) {
    fit <- eig_partial(A, k = k, target = target, B = B, tol = tol)
  }

  A_dense <- operator_source_matrix(as_operator(A))
  B_dense <- if (is.null(B)) NULL else operator_source_matrix(as_operator(B))
  oracle <- if (is.null(B_dense)) {
    eigen(A_dense, symmetric = is_square_symmetric(A_dense))
  } else {
    dense_generalized_spd_eigen(A_dense, B_dense, vectors = TRUE)
  }
  idx <- order_indices(oracle$values, target)
  idx <- idx[seq_len(min(k, length(idx)))]
  oracle_values <- oracle$values[idx]

  value_error <- value_error_summary(values(fit), oracle_values)
  cert <- if (!is.null(fit$vectors)) {
    certify_eigen(A_dense, values(fit), fit$vectors, B = B_dense, tol = tol)
  } else {
    empty_certificate(tol, "vectors not available for independent recomputation")
  }

  out <- list(
    kind = "eigen",
    requested = k,
    target = target_label(target),
    values = values(fit),
    oracle_values = oracle_values,
    value_abs_error = value_error$abs,
    value_rel_error = value_error$rel,
    max_value_abs_error = max(value_error$abs),
    max_value_rel_error = max(value_error$rel),
    recomputed_certificate = cert,
    certificate_agrees = certificates_agree(certificate(fit), cert),
    passed = max(value_error$rel) <= max(tol, 100 * .Machine$double.eps) &&
      isTRUE(cert$passed)
  )
  class(out) <- "eigencore_validation"
  out
}

#' Validate eigencore SVD results against base::svd().
#'
#' @param A Matrix or eigencore operator to validate.
#' @param rank Number of singular values to validate.
#' @param target Eigencore singular-value target descriptor.
#' @param fit Optional precomputed eigencore SVD result.
#' @param tol Validation tolerance.
#' @keywords internal
validate_svd_accuracy <- function(A, rank, target = largest(), fit = NULL,
                                  tol = 1e-8) {
  if (is.null(fit)) {
    fit <- svd_partial(A, rank = rank, target = target, tol = tol)
  }
  A_dense <- operator_source_matrix(as_operator(A))
  oracle <- svd(A_dense, nu = rank, nv = rank)
  idx <- order_indices(oracle$d, target)
  idx <- idx[seq_len(min(rank, length(idx)))]
  oracle_values <- oracle$d[idx]

  value_error <- value_error_summary(values(fit), oracle_values)
  cert <- if (!is.null(fit$u) && !is.null(fit$v)) {
    certify_svd(A_dense, fit$d, fit$u, fit$v, tol = tol)
  } else {
    empty_certificate(tol, "both singular-vector sides are required for independent recomputation")
  }

  out <- list(
    kind = "svd",
    requested = rank,
    target = target_label(target),
    values = values(fit),
    oracle_values = oracle_values,
    value_abs_error = value_error$abs,
    value_rel_error = value_error$rel,
    max_value_abs_error = max(value_error$abs),
    max_value_rel_error = max(value_error$rel),
    recomputed_certificate = cert,
    certificate_agrees = certificates_agree(certificate(fit), cert),
    passed = max(value_error$rel) <= max(tol, 100 * .Machine$double.eps) &&
      isTRUE(cert$passed)
  )
  class(out) <- "eigencore_validation"
  out
}

#' Benchmark eigen methods against base and optional references.
#'
#' @param A Matrix or eigencore operator to benchmark.
#' @param k Number of eigenpairs to compute.
#' @param target Eigencore eigenvalue target descriptor.
#' @param repeats Number of timing repetitions.
#' @param include Character vector of method labels to include when available.
#' @param tol Solver tolerance.
#' @keywords internal
benchmark_eigen_methods <- function(A, k, target = largest(), repeats = 3L,
                                    include = c("eigencore", "base", "RSpectra"),
                                    tol = 1e-8) {
  include <- intersect(include, available_eigen_methods())
  rows <- lapply(include, function(method) {
    timed <- time_repeated(repeats, {
      run_eigen_method(method, A, k = k, target = target, tol = tol)
    })
    fit <- timed$value
    list(
      method = method,
      median_seconds = stats::median(timed$times),
      min_seconds = min(timed$times),
      max_seconds = max(timed$times),
      mem_alloc = timed$mem_alloc,
      values = method_values(fit, kind = "eigen"),
      certificate_passed = method_certificate_passed(fit),
      certificate_type = method_certificate_field(fit, "certificate_type"),
      norm_bound_type = method_certificate_field(fit, "norm_bound_type"),
      scale_is_estimate = method_certificate_field(fit, "scale_is_estimate"),
      preconditioner_kind = method_preconditioner_field(fit, "kind"),
      preconditioner_native = method_preconditioner_field(fit, "native"),
      preconditioner_calls = method_preconditioner_calls(fit)
    )
  })
  class(rows) <- c("eigencore_benchmark", "list")
  rows
}

#' Benchmark SVD methods against base and optional references.
#'
#' @param A Matrix or eigencore operator to benchmark.
#' @param rank Number of singular values to compute.
#' @param repeats Number of timing repetitions.
#' @param include Character vector of method labels to include when available.
#' @param tol Solver tolerance.
#' @keywords internal
benchmark_svd_methods <- function(A, rank, repeats = 3L,
                                  include = c("eigencore", "base", "RSpectra", "irlba", "rsvd"),
                                  tol = 1e-8) {
  include <- intersect(include, available_svd_methods())
  rows <- lapply(include, function(method) {
    timed <- time_repeated(repeats, {
      run_svd_method(method, A, rank = rank, tol = tol)
    })
    fit <- timed$value
    list(
      method = method,
      median_seconds = stats::median(timed$times),
      min_seconds = min(timed$times),
      max_seconds = max(timed$times),
      mem_alloc = timed$mem_alloc,
      values = method_values(fit, kind = "svd"),
      certificate_passed = method_certificate_passed(fit),
      certificate_type = method_certificate_field(fit, "certificate_type"),
      norm_bound_type = method_certificate_field(fit, "norm_bound_type"),
      scale_is_estimate = method_certificate_field(fit, "scale_is_estimate"),
      fallback_attempted = method_restart_field(fit, "fallback_attempted"),
      fallback_used = method_restart_field(fit, "fallback_used"),
      fallback_method = method_restart_field(fit, "fallback_method"),
      gram_max_backward_error = method_restart_field(fit, "gram_max_backward_error"),
      fallback_max_backward_error = method_restart_field(fit, "fallback_max_backward_error"),
      preconditioner_kind = method_preconditioner_field(fit, "kind"),
      preconditioner_native = method_preconditioner_field(fit, "native"),
      preconditioner_calls = method_preconditioner_calls(fit)
    )
  })
  class(rows) <- c("eigencore_benchmark", "list")
  rows
}

#' @export
print.eigencore_validation <- function(x, ...) {
  cat("eigencore", x$kind, "validation\n")
  cat("  requested:", x$requested, "\n")
  cat("  target:", x$target, "\n")
  cat("  max value abs error:", format(x$max_value_abs_error), "\n")
  cat("  max value rel error:", format(x$max_value_rel_error), "\n")
  cat("  certificate agrees:", x$certificate_agrees, "\n")
  cat("  passed:", x$passed, "\n")
  invisible(x)
}

#' @export
print.eigencore_benchmark <- function(x, ...) {
  cat("eigencore benchmark\n")
  for (row in x) {
    cat("  ", row$method, ": median ", format(row$median_seconds),
        "s, certificate ", row$certificate_passed, "\n", sep = "")
  }
  invisible(x)
}

#' @keywords internal
available_eigen_methods <- function() {
  c(
    "eigencore",
    "eigencore_block_candidate",
    "eigencore_block",
    "eigencore_lobpcg",
    "eigencore_lobpcg_preconditioned",
    "eigencore_lobpcg_tridiagonal",
    "base",
    if (requireNamespace("RSpectra", quietly = TRUE)) "RSpectra"
  )
}

#' @keywords internal
available_svd_methods <- function() {
  c(
    "eigencore",
    "eigencore_smallest",
    "eigencore_interior",
    "eigencore_golub_kahan",
    "eigencore_golub_kahan_smallest",
    "eigencore_golub_kahan_interior",
    "eigencore_golub_kahan_one_sided",
    "eigencore_irlba_lbd_one_sided",
    "eigencore_irlba_lbd_retained_native",
    "eigencore_irlba_lbd_retained_bpro",
    "eigencore_irlba_lbd_retained_bpro_one_sided_guarded",
    "eigencore_irlba_lbd_retained_bpro_block_guarded",
    "eigencore_irlba_lbd_normal_scout",
    "eigencore_golub_kahan_projected",
    "eigencore_implicit_normal_lanczos",
    "eigencore_gram_dsyevx",
    "eigencore_block_golub_kahan_cycle",
    "eigencore_block_golub_kahan_cycle_cached",
    "eigencore_block_golub_kahan_cycle_cached_random",
    "eigencore_block_golub_kahan_cycle_residual",
    "eigencore_block_golub_kahan_cycle_lean",
    "eigencore_block_golub_kahan_retained",
    "eigencore_block_golub_kahan_retained_cached",
    "eigencore_block_golub_kahan_retained_deflated",
    "eigencore_randomized",
    "base",
    "base_smallest",
    "base_interior",
    if (requireNamespace("RSpectra", quietly = TRUE)) "RSpectra",
    if (requireNamespace("irlba", quietly = TRUE)) "irlba",
    if (requireNamespace("rsvd", quietly = TRUE)) "rsvd"
  )
}

#' @keywords internal
benchmark_block_candidate_lanczos_method <- function(A, k) {
  n <- nrow(A)
  k <- as.integer(k)
  if (k >= 16L && n >= 5000L) {
    lanczos(
      block = 4L,
      max_subspace = max(default_block_lanczos_max_subspace(k, 4L), 16L * k),
      max_restarts = 100L
    )
  } else {
    lanczos(block = 2L, max_restarts = 100L)
  }
}

#' @keywords internal
benchmark_rspectra_eigen_opts <- function(A, k, target) {
  n <- nrow(A)
  k <- as.integer(k)
  target_kind <- if (inherits(target, "eigencore_target")) target$kind else "largest"
  storage <- tryCatch(as_operator(A)$metadata$storage %||% NULL, error = function(e) NULL)
  if (identical(storage, "dgCMatrix") && identical(target_kind, "smallest") &&
      k >= 16L && n >= 5000L) {
    return(list(ncv = min(n, max(2L * k + 1L, 6L * k)), maxitr = 20000L))
  }
  list()
}

#' @keywords internal
run_eigen_method <- function(method, A, k, target, tol) {
  switch(
    method,
    eigencore = eig_partial(A, k = k, target = target, tol = tol),
    eigencore_block_candidate = eig_partial(
      A,
      k = k,
      target = target,
      method = benchmark_block_candidate_lanczos_method(A, k),
      tol = tol
    ),
    eigencore_block = eig_partial(
      A,
      k = k,
      target = target,
      method = benchmark_block_candidate_lanczos_method(A, k),
      tol = tol
    ),
    eigencore_lobpcg = eig_partial(
      A,
      k = k,
      target = target,
      method = lobpcg(maxit = 200L),
      tol = tol
    ),
    eigencore_lobpcg_preconditioned = {
      preconditioner <- shifted_cholesky_preconditioner(A, shift = 1e-3)
      eig_partial(
        A,
        k = k,
        target = target,
        method = lobpcg(maxit = 80L, preconditioner = preconditioner),
        tol = tol
      )
    },
    eigencore_lobpcg_tridiagonal = {
      native_lobpcg_tridiagonal_hermitian(
        as_operator(A),
        k = k,
        target = target,
        tol = tol,
        maxit = 80L,
        shift = 1e-3
      )
    },
    base = {
      eig <- eigen(as.matrix(A), symmetric = is_square_symmetric(as.matrix(A)))
      idx <- order_indices(eig$values, target)
      list(values = eig$values[idx[seq_len(k)]], vectors = eig$vectors[, idx[seq_len(k)], drop = FALSE])
    },
    RSpectra = {
      which <- target_to_rspectra_which(target, symmetric = TRUE)
      RSpectra::eigs_sym(
        A,
        k = k,
        which = which,
        opts = benchmark_rspectra_eigen_opts(A, k, target)
      )
    },
    stop("Unsupported eigen benchmark method: ", method, call. = FALSE)
  )
}

#' @keywords internal
irlba_lbd_one_sided_warm_start <- function(A, small) {
  dims <- as_operator(A)$dim
  if (length(dims) != 2L || is.null(small$u) || is.null(small$v)) {
    return(NULL)
  }
  start <- if (dims[[1L]] < dims[[2L]]) {
    small$u[, 1L]
  } else {
    small$v[, 1L]
  }
  start <- as.numeric(start)
  start_norm <- sqrt(sum(start^2))
  if (!is.finite(start_norm) || start_norm <= 100 * .Machine$double.eps) {
    return(NULL)
  }
  start / start_norm
}

#' @keywords internal
run_irlba_lbd_one_sided_method <- function(A, rank, tol, seed = NULL) {
  initial_work <- max(rank + 7L, 2L * rank + 1L)
  small <- svd_partial(
    A,
    rank = rank,
    method = golub_kahan(
      max_subspace = initial_work,
      reorthogonalize = FALSE
    ),
    tol = tol,
    seed = seed
  )
  small$restart$irlba_lbd_policy <- "small work one-sided LBD with certified adaptive fallback"
  small$restart$irlba_lbd_small_work_attempted <- TRUE
  small$restart$irlba_lbd_small_work_max_subspace <- initial_work
  small$restart$irlba_lbd_fallback_attempted <- FALSE
  small$restart$irlba_lbd_fallback_used <- FALSE
  small$restart$fallback_attempted <- FALSE
  small$restart$fallback_used <- FALSE
  if (isTRUE(small$certificate$passed)) {
    return(small)
  }

  warm_start <- irlba_lbd_one_sided_warm_start(A, small)
  fallback <- tryCatch(
    native_golub_kahan_svd(
      A,
      rank = rank,
      target = largest(),
      tol = tol,
      vectors = "both",
      reorthogonalize = FALSE,
      internal_start = warm_start
    ),
    error = function(e) {
      svd_partial(
        A,
        rank = rank,
        method = golub_kahan(reorthogonalize = FALSE),
        tol = tol,
        seed = seed
      )
    }
  )
  small_row <- data.frame(
    attempt = 1L,
    max_subspace = initial_work,
    iterations = small$iterations,
    matvecs = small$matvecs,
    accounted_seconds = sum(small$stage_seconds %||% NA_real_, na.rm = TRUE),
    warm_started = FALSE,
    certificate_passed = isTRUE(small$certificate$passed),
    max_backward_error = small$certificate$max_backward_error,
    max_residual = small$certificate$max_residual,
    stringsAsFactors = FALSE
  )
  fallback_row <- data.frame(
    attempt = 2L,
    max_subspace = fallback$restart$final_max_subspace %||%
      fallback$restart$max_subspace %||% fallback$iterations,
    iterations = fallback$iterations,
    matvecs = fallback$matvecs,
    accounted_seconds = sum(fallback$stage_seconds %||% NA_real_, na.rm = TRUE),
    warm_started = !is.null(warm_start),
    certificate_passed = isTRUE(fallback$certificate$passed),
    max_backward_error = fallback$certificate$max_backward_error,
    max_residual = fallback$certificate$max_residual,
    stringsAsFactors = FALSE
  )
  fallback$restart$irlba_lbd_policy <- "small work one-sided LBD with certified adaptive fallback"
  fallback$restart$irlba_lbd_small_work_attempted <- TRUE
  fallback$restart$irlba_lbd_small_work_max_subspace <- initial_work
  fallback$restart$irlba_lbd_small_work_certificate_passed <- FALSE
  fallback$restart$irlba_lbd_small_work_max_backward_error <-
    small$certificate$max_backward_error
  fallback$restart$irlba_lbd_small_work_matvecs <- small$matvecs
  fallback$restart$irlba_lbd_small_work_iterations <- small$iterations
  fallback$restart$irlba_lbd_small_work_accounted_seconds <-
    small_row$accounted_seconds[[1L]]
  fallback$restart$irlba_lbd_fallback_attempted <- TRUE
  fallback$restart$irlba_lbd_fallback_used <- TRUE
  fallback$restart$irlba_lbd_fallback_method <- "adaptive one-sided Golub-Kahan"
  fallback$restart$irlba_lbd_fallback_warm_started <- !is.null(warm_start)
  fallback$restart$irlba_lbd_fallback_matvecs <- fallback$matvecs
  fallback$restart$irlba_lbd_fallback_iterations <- fallback$iterations
  fallback$restart$irlba_lbd_fallback_accounted_seconds <-
    fallback_row$accounted_seconds[[1L]]
  fallback$restart$irlba_lbd_scout_matvec_overhead_fraction <-
    small$matvecs / max(1L, small$matvecs + fallback$matvecs)
  fallback$restart$fallback_attempted <- TRUE
  fallback$restart$fallback_used <- TRUE
  fallback$restart$fallback_method <- "adaptive one-sided Golub-Kahan"
  fallback$restart$irlba_lbd_attempt_history <- rbind(small_row, fallback_row)
  fallback$restart$attempt_history <- fallback$restart$irlba_lbd_attempt_history
  fallback$restart$attempted_subspaces <- fallback$restart$attempt_history$max_subspace
  fallback$restart$certified_attempt <- if (isTRUE(fallback$certificate$passed)) 2L else NA_integer_
  fallback
}

#' @keywords internal
irlba_lbd_normal_scout_operator <- function(A) {
  op <- as_operator(A)
  dims <- op$dim
  if (is.null(op$apply_adjoint)) {
    stop("normal-scout IRLBA/LBD requires an adjoint operator.", call. = FALSE)
  }
  if (dims[[1L]] < dims[[2L]]) {
    apply_normal <- function(X, alpha = 1, beta = 0, Y = NULL) {
      AX <- apply_adjoint_operator(op, X)
      apply_operator(op, AX, alpha = alpha, beta = beta, Y = Y)
    }
    side <- "left"
    dim <- dims[[1L]]
  } else {
    apply_normal <- function(X, alpha = 1, beta = 0, Y = NULL) {
      AX <- apply_operator(op, X)
      apply_adjoint_operator(op, AX, alpha = alpha, beta = beta, Y = Y)
    }
    side <- "right"
    dim <- dims[[2L]]
  }
  list(
    op = linear_operator(
      dim = c(dim, dim),
      apply = apply_normal,
      apply_adjoint = apply_normal,
      structure = hermitian(),
      name = paste0("normal_scout_", side, "(", op$name, ")"),
      metadata = list(
        parent = op,
        normal_scout = TRUE,
        normal_scout_side = side,
        materialized_normal = FALSE
      )
    ),
    side = side,
    internal_transposed = dims[[1L]] < dims[[2L]]
  )
}

#' @keywords internal
irlba_lbd_normal_scout_start <- function(A, rank, tol, steps, seed = NULL) {
  normal <- irlba_lbd_normal_scout_operator(A)
  if (!is.null(seed)) {
    set.seed(seed)
  }
  steps <- as.integer(steps)
  if (length(steps) != 1L || is.na(steps) || steps < rank + 1L) {
    stop("normal-scout steps must be at least rank + 1.", call. = FALSE)
  }
  fit <- eig_partial(
    normal$op,
    k = rank,
    target = largest(),
    method = lanczos(max_subspace = steps, max_restarts = 0L),
    tol = max(tol, 1e-6),
    vectors = TRUE,
    seed = seed,
    certify = FALSE,
    allow_dense_fallback = "never"
  )
  if (is.null(fit$vectors) || ncol(fit$vectors) < 1L) {
    stop("normal-scout Lanczos did not return a start vector.", call. = FALSE)
  }
  start <- as.numeric(fit$vectors[, 1L])
  start_norm <- sqrt(sum(start^2))
  if (!is.finite(start_norm) || start_norm <= 100 * .Machine$double.eps) {
    stop("normal-scout start vector has zero norm.", call. = FALSE)
  }
  list(
    start = start / start_norm,
    fit = fit,
    side = normal$side,
    internal_transposed = normal$internal_transposed,
    steps = steps
  )
}

#' @keywords internal
run_irlba_lbd_normal_scout_method <- function(A, rank, tol, seed = NULL,
                                              scout_steps = c(8L, 12L, 16L, 20L)) {
  limit <- min(as_operator(A)$dim)
  rank <- as.integer(rank)
  scout_steps <- sort(unique(as.integer(scout_steps)))
  scout_steps <- scout_steps[is.finite(scout_steps)]
  scout_steps <- scout_steps[scout_steps >= rank + 1L & scout_steps <= limit]
  if (!length(scout_steps)) {
    scout_steps <- min(limit, rank + 1L)
  }

  scouts <- lapply(seq_along(scout_steps), function(i) {
    steps <- scout_steps[[i]]
    tryCatch(
      irlba_lbd_normal_scout_start(
        A,
        rank = rank,
        tol = tol,
        steps = steps,
        seed = (seed %||% 1L) + i - 1L
      ),
      error = function(e) {
        structure(
          list(error = conditionMessage(e), steps = steps),
          class = "eigencore_irlba_lbd_normal_scout_error"
        )
      }
    )
  })
  ok <- which(!vapply(scouts, inherits, logical(1), "eigencore_irlba_lbd_normal_scout_error"))
  if (!length(ok)) {
    stop("normal-scout IRLBA/LBD could not produce a start vector.", call. = FALSE)
  }
  chosen <- tail(ok, 1L)[[1L]]
  scout <- scouts[[chosen]]

  if (!is.null(seed)) {
    set.seed(seed + 101L)
  }
  polished <- native_golub_kahan_svd(
    A,
    rank = rank,
    target = largest(),
    tol = tol,
    vectors = "both",
    reorthogonalize = FALSE,
    internal_start = scout$start
  )

  scout_rows <- lapply(seq_along(scouts), function(i) {
    item <- scouts[[i]]
    if (inherits(item, "eigencore_irlba_lbd_normal_scout_error")) {
      return(data.frame(
        attempt = i,
        max_subspace = item$steps,
        iterations = NA_integer_,
        matvecs = NA_integer_,
        accounted_seconds = NA_real_,
        warm_started = FALSE,
        certificate_passed = FALSE,
        max_backward_error = NA_real_,
        max_residual = NA_real_,
        scout_error = item$error,
        stringsAsFactors = FALSE
      ))
    }
    data.frame(
      attempt = i,
      max_subspace = item$steps,
      iterations = item$fit$iterations,
      matvecs = item$fit$matvecs,
      accounted_seconds = sum(item$fit$stage_seconds %||% NA_real_, na.rm = TRUE),
      warm_started = FALSE,
      certificate_passed = FALSE,
      max_backward_error = item$fit$certificate$max_backward_error %||% NA_real_,
      max_residual = item$fit$certificate$max_residual %||% NA_real_,
      scout_error = NA_character_,
      stringsAsFactors = FALSE
    )
  })
  scout_history <- do.call(rbind, scout_rows)
  polish_row <- data.frame(
    attempt = nrow(scout_history) + 1L,
    max_subspace = polished$restart$final_max_subspace %||%
      polished$restart$max_subspace %||% polished$iterations,
    iterations = polished$iterations,
    matvecs = polished$matvecs,
    accounted_seconds = sum(polished$stage_seconds %||% NA_real_, na.rm = TRUE),
    warm_started = TRUE,
    certificate_passed = isTRUE(polished$certificate$passed),
    max_backward_error = polished$certificate$max_backward_error,
    max_residual = polished$certificate$max_residual,
    scout_error = NA_character_,
    stringsAsFactors = FALSE
  )
  history <- rbind(scout_history, polish_row)
  scout_matvecs <- sum(scout_history$matvecs, na.rm = TRUE)
  polish_matvecs <- polished$matvecs %||% 0L

  polished$restart$irlba_lbd_policy <-
    "implicit normal scout warm-start with certified one-sided LBD polish"
  polished$restart$irlba_lbd_normal_scout_attempted <- TRUE
  polished$restart$irlba_lbd_normal_scout_steps <- paste(scout_steps, collapse = ",")
  polished$restart$irlba_lbd_normal_scout_chosen_steps <- scout$steps
  polished$restart$irlba_lbd_normal_scout_count <- nrow(scout_history)
  polished$restart$irlba_lbd_normal_scout_side <- scout$side
  polished$restart$irlba_lbd_normal_scout_materialized <- FALSE
  polished$restart$irlba_lbd_normal_scout_certificate_trusted <- FALSE
  polished$restart$irlba_lbd_normal_scout_matvecs <- scout_matvecs
  polished$restart$irlba_lbd_normal_scout_operator_matvecs <- 2L * scout_matvecs
  polished$restart$irlba_lbd_normal_scout_iterations <-
    sum(scout_history$iterations, na.rm = TRUE)
  polished$restart$irlba_lbd_normal_scout_accounted_seconds <-
    sum(scout_history$accounted_seconds, na.rm = TRUE)
  polished$restart$irlba_lbd_normal_scout_polish_matvecs <- polish_matvecs
  polished$restart$irlba_lbd_normal_scout_polish_iterations <- polished$iterations
  polished$restart$irlba_lbd_normal_scout_polish_accounted_seconds <-
    polish_row$accounted_seconds[[1L]]
  polished$restart$irlba_lbd_fallback_attempted <- TRUE
  polished$restart$irlba_lbd_fallback_used <- TRUE
  polished$restart$irlba_lbd_fallback_method <- "adaptive one-sided Golub-Kahan polish"
  polished$restart$irlba_lbd_fallback_warm_started <- TRUE
  polished$restart$irlba_lbd_fallback_matvecs <- polish_matvecs
  polished$restart$irlba_lbd_fallback_iterations <- polished$iterations
  polished$restart$irlba_lbd_fallback_accounted_seconds <-
    polish_row$accounted_seconds[[1L]]
  polished$restart$irlba_lbd_scout_matvec_overhead_fraction <-
    scout_matvecs / max(1L, scout_matvecs + polish_matvecs)
  polished$restart$fallback_attempted <- TRUE
  polished$restart$fallback_used <- FALSE
  polished$restart$fallback_method <- NA_character_
  polished$restart$attempt_history <- history
  polished$restart$attempted_subspaces <- history$max_subspace
  polished$restart$certified_attempt <- if (isTRUE(polished$certificate$passed)) nrow(history) else NA_integer_
  polished
}

#' @keywords internal
run_rspectra_svd_method <- function(A, rank, tol, solver = RSpectra::svds) {
  solver(
    A,
    k = rank,
    nu = rank,
    nv = rank,
    opts = list(tol = tol, maxitr = 1000L)
  )
}

#' @keywords internal
run_irlba_svd_method <- function(A, rank, tol, solver = irlba::irlba) {
  solver(A, nv = rank, nu = rank, tol = tol)
}

#' @keywords internal
run_svd_method <- function(method, A, rank, tol, seed = NULL) {
  switch(
    method,
    eigencore = svd_partial(A, rank = rank, tol = tol, seed = seed),
    eigencore_smallest = svd_partial(
      A,
      rank = rank,
      target = smallest(),
      tol = tol,
      seed = seed,
      allow_dense_fallback = "never"
    ),
    eigencore_interior = svd_partial(
      A,
      rank = rank,
      target = nearest(1),
      tol = tol,
      seed = seed,
      allow_dense_fallback = "never"
    ),
    eigencore_golub_kahan = svd_partial(
      A,
      rank = rank,
      method = golub_kahan(),
      tol = tol,
      seed = seed
    ),
    eigencore_golub_kahan_smallest = svd_partial(
      A,
      rank = rank,
      target = smallest(),
      method = golub_kahan(),
      tol = tol,
      seed = seed,
      allow_dense_fallback = "never"
    ),
    eigencore_golub_kahan_interior = svd_partial(
      A,
      rank = rank,
      target = nearest(1),
      method = golub_kahan(),
      tol = tol,
      seed = seed,
      allow_dense_fallback = "never"
    ),
    eigencore_golub_kahan_one_sided = svd_partial(
      A,
      rank = rank,
      method = golub_kahan(reorthogonalize = FALSE),
      tol = tol,
      seed = seed
    ),
    eigencore_irlba_lbd_one_sided = run_irlba_lbd_one_sided_method(
      A,
      rank = rank,
      tol = tol,
      seed = seed
    ),
    eigencore_irlba_lbd_retained_native = {
      set.seed(seed %||% 1L)
      native_irlba_lbd_retained_svd(
        A,
        rank = rank,
        target = largest(),
        work = max(rank + 7L, 2L * rank + 1L),
        retained = min(max(rank, rank + 2L), max(rank + 6L, 2L * rank)),
        max_restarts = 7L,
        tol = tol,
        vectors = "both"
      )
    },
    eigencore_irlba_lbd_retained_bpro = {
      set.seed(seed %||% 1L)
      native_irlba_lbd_retained_svd(
        A,
        rank = rank,
        target = largest(),
        work = max(rank + 7L, 2L * rank + 1L),
        retained = min(max(rank, rank + 2L), max(rank + 6L, 2L * rank)),
        max_restarts = 7L,
        tol = tol,
        vectors = "both",
        reorth_policy = "bpro_two_sided"
      )
    },
    eigencore_irlba_lbd_retained_bpro_one_sided_guarded = {
      set.seed(seed %||% 1L)
      native_irlba_lbd_retained_svd(
        A,
        rank = rank,
        target = largest(),
        work = max(rank + 7L, 2L * rank + 1L),
        retained = min(max(rank, rank + 2L), max(rank + 6L, 2L * rank)),
        max_restarts = 7L,
        tol = tol,
        vectors = "both",
        reorth_policy = "bpro_one_sided_guarded"
      )
    },
    eigencore_irlba_lbd_retained_bpro_block_guarded = {
      set.seed(seed %||% 1L)
      native_irlba_lbd_retained_svd(
        A,
        rank = rank,
        target = largest(),
        work = max(rank + 7L, 2L * rank + 1L),
        retained = min(max(rank, rank + 2L), max(rank + 6L, 2L * rank)),
        max_restarts = 7L,
        tol = tol,
        vectors = "both",
        reorth_policy = "bpro_block_guarded"
      )
    },
    eigencore_irlba_lbd_normal_scout = run_irlba_lbd_normal_scout_method(
      A,
      rank = rank,
      tol = tol,
      seed = seed
    ),
    eigencore_golub_kahan_projected = {
      old_options <- options(eigencore.golub_kahan_projected_stop = TRUE)
      on.exit(options(old_options), add = TRUE)
      svd_partial(
        A,
        rank = rank,
        method = golub_kahan(),
        tol = tol,
        seed = seed
      )
    },
    eigencore_implicit_normal_lanczos = {
      old_options <- options(
        eigencore.csc_left_normal_lanczos_attempt = TRUE,
        eigencore.csc_right_normal_lanczos_attempt = TRUE
      )
      on.exit(options(old_options), add = TRUE)
      svd_partial(A, rank = rank, tol = tol, seed = seed)
    },
    eigencore_gram_dsyevx = {
      old_options <- options(eigencore.csc_left_gram_dsyevx_attempt = TRUE)
      on.exit(options(old_options), add = TRUE)
      svd_partial(A, rank = rank, tol = tol, seed = seed)
    },
    eigencore_block_golub_kahan_cycle = {
      set.seed(seed %||% 1L)
      native_block_golub_kahan_cycle_svd(
        A,
        rank = rank,
        target = largest(),
        tol = tol,
        vectors = "both"
      )
    },
    eigencore_block_golub_kahan_cycle_cached = {
      set.seed(seed %||% 1L)
      native_block_golub_kahan_cycle_svd(
        A,
        rank = rank,
        target = largest(),
        tol = tol,
        adaptive_start = "ritz_cached",
        vectors = "both"
      )
    },
    eigencore_block_golub_kahan_cycle_cached_random = {
      set.seed(seed %||% 1L)
      native_block_golub_kahan_cycle_svd(
        A,
        rank = rank,
        target = largest(),
        tol = tol,
        adaptive_start = "ritz_cached_random",
        vectors = "both"
      )
    },
    eigencore_block_golub_kahan_cycle_residual = {
      set.seed(seed %||% 1L)
      native_block_golub_kahan_cycle_svd(
        A,
        rank = rank,
        target = largest(),
        tol = tol,
        adaptive_start = "ritz_residual",
        vectors = "both"
      )
    },
    eigencore_block_golub_kahan_cycle_lean = {
      set.seed(seed %||% 1L)
      native_block_golub_kahan_cycle_svd(
        A,
        rank = rank,
        target = largest(),
        tol = tol,
        adaptive_start = "ritz_lean",
        vectors = "both"
      )
    },
    eigencore_block_golub_kahan_retained = {
      set.seed(seed %||% 1L)
      native_block_golub_kahan_retained_cycle_svd(
        A,
        rank = rank,
        target = largest(),
        tol = tol,
        retained_av_cache = FALSE,
        vectors = "both"
      )
    },
    eigencore_block_golub_kahan_retained_cached = {
      set.seed(seed %||% 1L)
      native_block_golub_kahan_retained_cycle_svd(
        A,
        rank = rank,
        target = largest(),
        tol = tol,
        retained_av_cache = TRUE,
        vectors = "both"
      )
    },
    eigencore_block_golub_kahan_retained_deflated = {
      set.seed(seed %||% 1L)
      native_block_golub_kahan_retained_cycle_svd(
        A,
        rank = rank,
        target = largest(),
        tol = tol,
        retained_av_cache = TRUE,
        retained_deflation = TRUE,
        vectors = "both"
      )
    },
    eigencore_randomized = svd_partial(
      A,
      rank = rank,
      method = randomized(oversample = 10L, n_iter = 2L),
      tol = tol,
      seed = seed
    ),
    base = {
      decomp <- svd(as.matrix(A), nu = rank, nv = rank)
      idx <- seq_len(min(rank, length(decomp$d)))
      list(
        d = decomp$d[idx],
        u = decomp$u[, idx, drop = FALSE],
        v = decomp$v[, idx, drop = FALSE],
        values = decomp$d[idx]
      )
    },
    base_smallest = {
      decomp <- svd(as.matrix(A))
      idx <- order_indices(decomp$d, smallest())
      idx <- idx[seq_len(min(rank, length(idx)))]
      list(
        d = decomp$d[idx],
        u = decomp$u[, idx, drop = FALSE],
        v = decomp$v[, idx, drop = FALSE],
        values = decomp$d[idx]
      )
    },
    base_interior = {
      decomp <- svd(as.matrix(A))
      idx <- order_indices(decomp$d, nearest(1))
      idx <- idx[seq_len(min(rank, length(idx)))]
      list(
        d = decomp$d[idx],
        u = decomp$u[, idx, drop = FALSE],
        v = decomp$v[, idx, drop = FALSE],
        values = decomp$d[idx]
      )
    },
    RSpectra = run_rspectra_svd_method(A, rank = rank, tol = tol),
    irlba = run_irlba_svd_method(A, rank = rank, tol = tol),
    rsvd = rsvd::rsvd(A, k = rank, nu = rank, nv = rank),
    stop("Unsupported SVD benchmark method: ", method, call. = FALSE)
  )
}

#' @keywords internal
time_repeated <- function(repeats, expr) {
  expr <- substitute(expr)
  env <- parent.frame()
  if (requireNamespace("bench", quietly = TRUE)) {
    value <- NULL
    mark <- bench::mark(
      value <- eval(expr, envir = env),
      iterations = repeats,
      check = FALSE,
      time_unit = "s",
      memory = TRUE,
      filter_gc = FALSE
    )
    return(list(
      times = as.numeric(mark$time[[1L]]),
      value = value,
      mem_alloc = as.numeric(mark$mem_alloc[[1L]]),
      bench = mark
    ))
  }

  times <- numeric(repeats)
  value <- NULL
  for (i in seq_len(repeats)) {
    gc(FALSE)
    elapsed <- system.time(value <- eval(expr, envir = env))[["elapsed"]]
    times[[i]] <- elapsed
  }
  list(times = times, value = value, mem_alloc = NA_real_, bench = NULL)
}

#' @keywords internal
method_values <- function(x, kind) {
  if (kind == "svd") {
    return(x$d %||% x$values)
  }
  x$values
}

#' @keywords internal
method_certificate_passed <- function(x) {
  cert <- x$certificate
  if (is.null(cert$passed)) NA else cert$passed
}

#' @keywords internal
method_certificate_field <- function(x, field) {
  cert <- x$certificate
  if (is.null(cert[[field]])) NA else cert[[field]]
}

#' @keywords internal
method_restart_field <- function(x, field) {
  restart <- x$restart %||% NULL
  if (is.null(restart) || is.null(restart[[field]])) {
    return(NA)
  }
  restart[[field]]
}

#' @keywords internal
method_preconditioner_field <- function(x, field) {
  info <- x$preconditioner %||% x$restart$preconditioner %||% NULL
  if (is.null(info) || is.null(info[[field]])) {
    return(NA)
  }
  info[[field]]
}

#' @keywords internal
method_preconditioner_calls <- function(x) {
  x$preconditioner_calls %||% x$restart$preconditioner_calls %||% NA_integer_
}

#' @keywords internal
value_error_summary <- function(values, oracle_values) {
  abs_err <- abs(values - oracle_values)
  rel_err <- abs_err / pmax(abs(oracle_values), .Machine$double.eps)
  list(abs = abs_err, rel = rel_err)
}

#' @keywords internal
certificates_agree <- function(a, b, tol = 100 * .Machine$double.eps) {
  if (is.null(a) || is.null(b)) {
    return(FALSE)
  }
  both <- c(a$max_residual, b$max_residual, a$max_backward_error, b$max_backward_error)
  if (anyNA(both)) {
    return(FALSE)
  }
  abs(a$max_residual - b$max_residual) <= tol * max(1, abs(b$max_residual)) &&
    abs(a$max_backward_error - b$max_backward_error) <= tol * max(1, abs(b$max_backward_error))
}

#' @keywords internal
target_to_rspectra_which <- function(target, symmetric = TRUE) {
  kind <- if (inherits(target, "eigencore_target")) target$kind else "largest"
  switch(
    kind,
    largest = if (symmetric) "LA" else "LR",
    smallest = if (symmetric) "SA" else "SR",
    largest_magnitude = "LM",
    smallest_magnitude = "SM",
    largest_real = "LR",
    smallest_real = "SR",
    largest_imaginary = "LI",
    smallest_imaginary = "SI",
    both_ends = "BE",
    "LM"
  )
}
