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
      values = method_values(fit, kind = "eigen"),
      certificate_passed = method_certificate_passed(fit)
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
      values = method_values(fit, kind = "svd"),
      certificate_passed = method_certificate_passed(fit)
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
  c("eigencore", "base", if (requireNamespace("RSpectra", quietly = TRUE)) "RSpectra")
}

#' @keywords internal
available_svd_methods <- function() {
  c(
    "eigencore",
    "base",
    if (requireNamespace("RSpectra", quietly = TRUE)) "RSpectra",
    if (requireNamespace("irlba", quietly = TRUE)) "irlba",
    if (requireNamespace("rsvd", quietly = TRUE)) "rsvd"
  )
}

#' @keywords internal
run_eigen_method <- function(method, A, k, target, tol) {
  switch(
    method,
    eigencore = eig_partial(A, k = k, target = target, tol = tol),
    base = {
      eig <- eigen(as.matrix(A), symmetric = is_square_symmetric(as.matrix(A)))
      idx <- order_indices(eig$values, target)
      list(values = eig$values[idx[seq_len(k)]], vectors = eig$vectors[, idx[seq_len(k)], drop = FALSE])
    },
    RSpectra = {
      which <- target_to_rspectra_which(target, symmetric = TRUE)
      RSpectra::eigs_sym(A, k = k, which = which)
    },
    stop("Unsupported eigen benchmark method: ", method, call. = FALSE)
  )
}

#' @keywords internal
run_svd_method <- function(method, A, rank, tol) {
  switch(
    method,
    eigencore = svd_partial(A, rank = rank, tol = tol),
    base = svd(as.matrix(A), nu = rank, nv = rank),
    RSpectra = RSpectra::svds(A, k = rank),
    irlba = irlba::irlba(A, nv = rank, nu = rank),
    rsvd = rsvd::rsvd(A, k = rank, nu = rank, nv = rank),
    stop("Unsupported SVD benchmark method: ", method, call. = FALSE)
  )
}

#' @keywords internal
time_repeated <- function(repeats, expr) {
  expr <- substitute(expr)
  env <- parent.frame()
  times <- numeric(repeats)
  value <- NULL
  for (i in seq_len(repeats)) {
    gc(FALSE)
    elapsed <- system.time(value <- eval(expr, envir = env))[["elapsed"]]
    times[[i]] <- elapsed
  }
  list(times = times, value = value)
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
    "LM"
  )
}
