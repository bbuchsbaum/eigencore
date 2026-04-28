#' Validate eigencore eigen results against a dense oracle.
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
    if (requireNamespace("RSpectra", quietly = TRUE)) "RSpectra",
    if (requireNamespace("PRIMME", quietly = TRUE)) "PRIMME"
  )
}

#' @keywords internal
available_svd_methods <- function() {
  c(
    "eigencore",
    "eigencore_golub_kahan",
    "eigencore_golub_kahan_projected",
    "eigencore_block_golub_kahan_cycle",
    "eigencore_block_golub_kahan_cycle_lean",
    "eigencore_randomized",
    "base",
    if (requireNamespace("RSpectra", quietly = TRUE)) "RSpectra",
    if (requireNamespace("PRIMME", quietly = TRUE)) "PRIMME",
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
    PRIMME = {
      which <- target_to_rspectra_which(target, symmetric = TRUE)
      PRIMME::eigs_sym(A, NEig = k, which = which, tol = tol)
    },
    stop("Unsupported eigen benchmark method: ", method, call. = FALSE)
  )
}

#' @keywords internal
run_svd_method <- function(method, A, rank, tol, seed = NULL) {
  switch(
    method,
    eigencore = svd_partial(A, rank = rank, tol = tol, seed = seed),
    eigencore_golub_kahan = svd_partial(
      A,
      rank = rank,
      method = golub_kahan(),
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
    RSpectra = RSpectra::svds(A, k = rank),
    PRIMME = PRIMME::svds(A, NSvals = rank, which = "L", tol = tol),
    irlba = irlba::irlba(A, nv = rank, nu = rank),
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
