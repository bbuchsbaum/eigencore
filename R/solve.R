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
  fast <- try_svd_partial_native_gram_fastpath(
    A = A,
    rank = rank,
    target = target,
    method = method,
    tol = tol,
    vectors = vectors,
    certify = certify
  )
  if (!is.null(fast)) {
    return(fast)
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
  if (plan_dispatches_lobpcg(plan)) {
    return(solve_eigen_lobpcg(a, k, method, tol, maxit, vectors, certify, plan))
  }
  if (plan_dispatches_lanczos(plan)) {
    return(solve_eigen_lanczos(a, k, method, tol, maxit, vectors, certify, plan))
  }
  if (identical(plan$method, "native dense Hermitian LAPACK fallback")) {
    return(solve_eigen_native_dense_hermitian(a, k, tol, vectors, certify,
                                              allow_dense_fallback, plan))
  }
  solve_eigen_dense_oracle(a, k, tol, vectors, certify, allow_dense_fallback, plan)
}

#' @export
solve.eigencore_svd_problem <- function(a, b, rank, method = auto(), tol = 1e-8,
                                        vectors = c("both", "left", "right", "none"),
                                        certify = TRUE,
                                        allow_dense_fallback = c("auto", "never", "always"), ...) {
  vectors <- match.arg(vectors)
  allow_dense_fallback <- match.arg(allow_dense_fallback)
  plan <- plan_solver(a, rank = rank, method = method)
  if (identical(plan$method, "reference randomized SVD prototype")) {
    return(solve_svd_randomized(a, rank, method, tol, vectors, certify, plan))
  }
  if (identical(plan$method, "native certified Gram SVD special case")) {
    return(solve_svd_gram(a, rank, tol, vectors, certify, plan))
  }
  if (identical(plan$method, "native retained Golub-Kahan SVD (thick restart)")) {
    return(solve_svd_retained_golub_kahan(a, rank, tol, vectors, certify, plan))
  }
  if (plan_dispatches_golub_kahan(plan)) {
    return(solve_svd_golub_kahan(a, rank, method, tol, vectors, certify, plan))
  }
  solve_svd_dense(a, rank, tol, vectors, certify, allow_dense_fallback, plan)
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
native_standard_lobpcg_label <- function() {
  "native standard Hermitian LOBPCG prototype"
}

#' @keywords internal
plan_dispatches_lobpcg <- function(plan) {
  plan$method %in% c(
    native_standard_lobpcg_label(),
    native_generalized_lobpcg_label(),
    "reference LOBPCG prototype",
    "reference generalized SPD LOBPCG prototype",
    reference_generalized_lobpcg_label()
  )
}

#' @keywords internal
plan_dispatches_lanczos <- function(plan) {
  plan$method %in% c(
    "native scalar thick-restart Hermitian Lanczos",
    "native block Hermitian Lanczos thick-restart candidate",
    "native block Hermitian Lanczos (thick restart, locking)",
    "reference Hermitian Lanczos (target unsupported by native path)",
    "reference Hermitian Lanczos (prototype/oracle fallback)"
  )
}

#' @keywords internal
plan_dispatches_native_lanczos <- function(plan) {
  plan$method %in% c(
    "native scalar thick-restart Hermitian Lanczos",
    "native block Hermitian Lanczos thick-restart candidate",
    "native block Hermitian Lanczos (thick restart, locking)"
  )
}

#' @keywords internal
plan_dispatches_golub_kahan <- function(plan) {
  plan$method %in% c("native prototype Golub-Kahan", "prototype Golub-Kahan")
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
native_dense_symmetric_eigen_dsyevd <- function(A) {
  .Call("eigencore_dense_symmetric_eigen_dsyevd", as.matrix(A), PACKAGE = "eigencore")
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
  if (!isTRUE(getOption("eigencore.promote_retained_golub_kahan", FALSE))) {
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
