#' Compute a partial eigendecomposition.
#'
#' @param A Matrix or eigencore operator.
#' @param k Number of eigenpairs to compute.
#' @param target Eigencore eigenvalue target descriptor.
#' @param B Optional metric matrix or operator for generalized problems.
#' @param method Solver method descriptor.
#' @param tol Convergence and certification tolerance.
#' @param maxit Optional iteration limit.
#' @param vectors Whether to compute vectors.
#' @param seed Optional random seed for stochastic solver components.
#' @param certify Whether to compute certification diagnostics.
#' @param allow_dense_fallback Dense fallback policy.
#' @examples
#' A <- diag(c(5, 4, 3, 2, 1))
#' A[1, 2] <- A[2, 1] <- 0.1
#' fit <- eig_partial(A, k = 2, target = largest())
#' values(fit)
#' certificate(fit)$passed
#'
#' # Generalized SPD problem A x = lambda B x
#' B <- diag(c(2, 1, 1, 1, 1))
#' gfit <- eig_partial(A, B = B, k = 2, target = smallest())
#' values(gfit)
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
#'
#' @param A Matrix or eigencore operator.
#' @param rank Number of singular values to compute.
#' @param target Eigencore singular-value target descriptor.
#' @param method Solver method descriptor.
#' @param tol Convergence and certification tolerance.
#' @param vectors Which singular-vector sides to compute.
#' @param seed Optional random seed for stochastic solver components.
#' @param certify Whether to compute certification diagnostics.
#' @param allow_dense_fallback Dense fallback policy.
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(60), 10, 6)
#' fit <- svd_partial(X, rank = 3)
#' values(fit)
#' certificate(fit)$passed
svd_partial <- function(A, rank, target = largest(), method = auto(), tol = 1e-8,
                        vectors = c("both", "left", "right", "none"),
                        seed = NULL, certify = TRUE,
                        allow_dense_fallback = c("auto", "never", "always")) {
  vectors <- match.arg(vectors)
  allow_dense_fallback <- match.arg(allow_dense_fallback)
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
  if (!is.null(seed)) {
    set.seed(seed)
  }
  P <- svd_problem(A, target = target)
  solve(P, rank = rank, method = method, tol = tol, vectors = vectors,
        certify = certify, allow_dense_fallback = allow_dense_fallback)
}

#' Solve a planned eigenproblem.
#'
#' S3 method that runs the planned solver for an eigenproblem built by
#' [eigen_problem()]. Most users call [eig_partial()], which constructs the
#' problem and dispatches here; call `solve()` directly when you want to build
#' a problem once and reuse or inspect it. Returns a certified partial
#' eigendecomposition.
#' @param a Eigencore eigen problem object.
#' @param b Unused second argument reserved by the base [solve()] generic.
#' @param k Number of eigenpairs to compute.
#' @param method Solver method descriptor.
#' @param tol Convergence and certification tolerance.
#' @param maxit Optional iteration limit.
#' @param vectors Whether to compute vectors.
#' @param certify Whether to compute certification diagnostics.
#' @param allow_dense_fallback Dense fallback policy.
#' @param ... Reserved for future solver options.
#' @export
solve.eigencore_eigen_problem <- function(a, b, k, method = auto(), tol = 1e-8,
                                          maxit = NULL, vectors = TRUE,
                                          certify = TRUE,
                                          allow_dense_fallback = c("auto", "never", "always"), ...) {
  allow_dense_fallback <- match.arg(allow_dense_fallback)
  auto_shift <- auto_shift_invert_route(a, method)
  a <- auto_shift$problem
  method <- auto_shift$method
  plan <- plan_solver(a, k = k, method = method)
  plan <- validate_complex_eigen_plan(a, plan)
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
  if (plan_dispatches_structured_grid_laplacian_2d(plan)) {
    return(solve_eigen_grid_laplacian_2d(a, k, tol, vectors, certify, plan))
  }
  if (plan_dispatches_native_tridiagonal_eigen(plan)) {
    return(solve_eigen_native_tridiagonal_hermitian(a, k, tol, vectors, certify, plan))
  }
  if (plan_dispatches_lanczos(plan)) {
    return(solve_eigen_lanczos(a, k, method, tol, maxit, vectors, certify, plan))
  }
  if (plan_dispatches_arnoldi(plan)) {
    return(solve_eigen_arnoldi(a, k, method, tol, maxit, vectors, certify, plan))
  }
  if (plan$method %in% c("native dense Hermitian LAPACK fallback",
                         native_dense_complex_hermitian_label())) {
    return(solve_eigen_native_dense_hermitian(a, k, tol, vectors, certify,
                                              allow_dense_fallback, plan))
  }
  if (identical(plan$method, native_dense_complex_general_label())) {
    return(solve_eigen_native_dense_general(a, k, tol, vectors, certify,
                                            allow_dense_fallback, plan))
  }
  solve_eigen_dense_oracle(a, k, tol, vectors, certify, allow_dense_fallback, plan)
}

#' Solve a planned SVD problem.
#'
#' S3 method that runs the planned solver for an SVD problem built by
#' [svd_problem()]. Most users call [svd_partial()], which constructs the
#' problem and dispatches here; call `solve()` directly when you want to build
#' a problem once and reuse or inspect it. Returns a certified partial
#' singular-value decomposition.
#' @param a Eigencore SVD problem object.
#' @param b Unused second argument reserved by the base [solve()] generic.
#' @param rank Number of singular values to compute.
#' @param method Solver method descriptor.
#' @param tol Convergence and certification tolerance.
#' @param vectors Which singular-vector sides to compute.
#' @param certify Whether to compute certification diagnostics.
#' @param allow_dense_fallback Dense fallback policy.
#' @param ... Reserved for future solver options.
#' @export
solve.eigencore_svd_problem <- function(a, b, rank, method = auto(), tol = 1e-8,
                                        vectors = c("both", "left", "right", "none"),
                                        certify = TRUE,
                                        allow_dense_fallback = c("auto", "never", "always"), ...) {
  vectors <- match.arg(vectors)
  allow_dense_fallback <- match.arg(allow_dense_fallback)
  plan <- plan_solver(a, rank = rank, method = method)
  plan <- validate_complex_svd_plan(a, plan)
  plan <- validate_or_route_svd_target_plan(
    a,
    plan = plan,
    rank = rank,
    method = method,
    allow_dense_fallback = allow_dense_fallback
  )
  if (identical(plan$method, "reference randomized SVD prototype") ||
      identical(plan$method, native_dense_randomized_svd_label()) ||
      identical(plan$method, native_csc_randomized_svd_label())) {
    return(solve_svd_randomized(a, rank, method, tol, vectors, certify, plan))
  }
  if (identical(plan$method, "native certified Gram SVD special case")) {
    return(solve_svd_gram(a, rank, tol, vectors, certify, plan))
  }
  if (identical(plan$method, native_retained_golub_kahan_diagnostic_label())) {
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
native_dense_generalized_pencil_eigen <- function(A, B) {
  .Call(
    "eigencore_dense_generalized_pencil_eigen",
    as.matrix(A),
    as.matrix(B),
    PACKAGE = "eigencore"
  )
}

#' @keywords internal
native_dense_complex_generalized_hpd_eigen <- function(A, B) {
  .Call(
    "eigencore_dense_complex_generalized_hpd_eigen",
    as.matrix(A),
    as.matrix(B),
    PACKAGE = "eigencore"
  )
}

#' @keywords internal
native_dense_complex_generalized_pencil_eigen <- function(A, B) {
  .Call(
    "eigencore_dense_complex_generalized_pencil_eigen",
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
native_tridiagonal_hermitian_label <- function() {
  "native tridiagonal Hermitian LAPACK selected eigensolver"
}

#' @keywords internal
structured_grid_laplacian_2d_label <- function() {
  "diagnostic separable 2D-grid Laplacian eigensolver"
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
plan_dispatches_native_tridiagonal_eigen <- function(plan) {
  identical(plan$method, native_tridiagonal_hermitian_label())
}

#' @keywords internal
plan_dispatches_structured_grid_laplacian_2d <- function(plan) {
  identical(plan$method, structured_grid_laplacian_2d_label())
}

#' @keywords internal
plan_dispatches_lanczos <- function(plan) {
  plan$method %in% c(
    "native scalar thick-restart Hermitian Lanczos",
    "native block Hermitian Lanczos thick-restart candidate",
    "native block Hermitian Lanczos (thick restart, locking)",
    native_generalized_lanczos_label(),
    generalized_lanczos_label(),
    "reference Hermitian Lanczos (target unsupported by native path)",
    "reference Hermitian Lanczos (prototype/oracle fallback)"
  )
}

#' @keywords internal
plan_dispatches_arnoldi <- function(plan) {
  identical(plan$method, reference_arnoldi_label()) ||
    identical(plan$method, native_arnoldi_label()) ||
    identical(plan$method, native_refined_arnoldi_label()) ||
    identical(plan$method, native_matrix_free_arnoldi_label())
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
native_dense_complex_hermitian_label <- function() {
  "native dense complex Hermitian LAPACK fallback"
}

#' @keywords internal
native_dense_complex_general_label <- function() {
  "native dense complex general LAPACK fallback"
}

#' @keywords internal
native_dense_generalized_spd_full_label <- function() {
  "native dense generalized SPD/Hermitian LAPACK full"
}

#' @keywords internal
native_dense_generalized_pencil_full_label <- function() {
  "native dense general pencil LAPACK full"
}

#' @keywords internal
native_dense_complex_svd_label <- function() {
  "native dense complex LAPACK SVD fallback"
}

#' @keywords internal
plan_dispatches_golub_kahan <- function(plan) {
  plan$method %in% c(
    "native prototype Golub-Kahan",
    native_smallest_golub_kahan_label(),
    native_interior_golub_kahan_label(),
    native_matrix_free_golub_kahan_label(),
    native_matrix_free_smallest_golub_kahan_label(),
    native_matrix_free_interior_golub_kahan_label(),
    "prototype Golub-Kahan"
  )
}

#' @keywords internal
validate_complex_eigen_plan <- function(problem, plan) {
  if (!identical(problem$A$dtype, "complex")) {
    return(plan)
  }
  if (!is.null(source_or_null(problem$A))) {
    return(plan)
  }
  stop(
    "Complex matrix-free eigen operators are future scope in eigencore. ",
    "Use a base complex dense matrix for the native dense complex LAPACK path; ",
    "native complex callback/sparse operator support is not promoted yet.",
    call. = FALSE
  )
}

#' @keywords internal
validate_complex_svd_plan <- function(problem, plan) {
  if (!identical(problem$A$dtype, "complex")) {
    return(plan)
  }
  if (!is.null(source_or_null(problem$A))) {
    return(plan)
  }
  stop(
    "Complex matrix-free SVD operators are future scope in eigencore. ",
    "Use a base complex dense matrix for the native dense complex LAPACK path; ",
    "native complex callback/sparse SVD support is not promoted yet.",
    call. = FALSE
  )
}

#' @keywords internal
svd_plan_uses_unsupported_native_interior <- function(plan) {
  identical(plan$controls$svd_target_family %||% NULL, "interior") &&
    svd_native_iterative_plan(plan$method) &&
    !plan$method %in% c(
      native_interior_golub_kahan_label(),
      native_matrix_free_interior_golub_kahan_label()
    )
}

#' @keywords internal
svd_route_explicit_dense_interior_fallback <- function(problem, plan, rank,
                                                       method) {
  old_method <- plan$method
  source <- source_or_null(problem$A)
  plan$method <- if (is.matrix(source) && is.double(source)) {
    "native dense LAPACK SVD fallback"
  } else if (is.matrix(source) && is.complex(source)) {
    native_dense_complex_svd_label()
  } else {
    "dense LAPACK SVD oracle (prototype fallback)"
  }
  plan$reasons <- c(
    plan$reasons,
    paste0(
      "interior SVD target ", target_label(problem$target),
      " is unsupported by ", old_method,
      "; allow_dense_fallback = 'always' selected an explicit dense fallback"
    )
  )
  plan$fallback <- "explicit dense fallback for unsupported native interior SVD target"
  plan$controls <- svd_plan_controls(problem, rank = rank, method = method, chosen = plan$method)
  plan$controls$interior_svd_dense_fallback <- TRUE
  plan$controls$interior_svd_previous_method <- old_method
  plan
}

#' @keywords internal
validate_or_route_svd_target_plan <- function(problem, plan, rank, method,
                                              allow_dense_fallback) {
  if (!svd_plan_uses_unsupported_native_interior(plan)) {
    return(plan)
  }
  if (identical(allow_dense_fallback, "always")) {
    return(svd_route_explicit_dense_interior_fallback(problem, plan, rank, method))
  }
  stop(
    "Interior SVD target ", target_label(problem$target),
    " is not supported by ", plan$method,
    ". Use a dense matrix or set allow_dense_fallback = 'always' for an ",
    "explicit dense fallback; sparse and matrix-free shift-invert/refined ",
    "interior SVD is not production-promoted yet.",
    call. = FALSE
  )
}

#' @keywords internal
should_use_lanczos <- function(problem, method, k = NULL) {
  if (!identical(problem$structure$kind, "hermitian")) {
    return(FALSE)
  }
  if (!is.null(problem$metric)) {
    return(
      inherits(method, "eigencore_method") &&
        identical(method$kind, "lanczos") &&
        generalized_lanczos_supported(problem$A, problem$metric, target = problem$target)
    )
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
should_use_native_tridiagonal_hermitian <- function(problem, k = NULL) {
  if (!is.null(problem$metric) ||
      !identical(problem$structure$kind, "hermitian") ||
      !native_lanczos_target_supported(problem$target)) {
    return(FALSE)
  }
  target <- target_label(problem$target)
  if (!target %in% c("largest", "smallest")) {
    return(FALSE)
  }
  if (is.null(k) || length(k) != 1L || is.na(k) || as.integer(k) < 1L) {
    return(FALSE)
  }
  A <- problem$A$metadata$matrix %||% source_or_null(problem$A)
  if (!(inherits(A, "CsparseMatrix") || inherits(A, "diagonalMatrix"))) {
    return(FALSE)
  }
  !is.null(shift_invert_tridiagonal_parts(A, shift = 0))
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
native_dense_complex_hermitian_eigen <- function(A) {
  .Call("eigencore_dense_complex_hermitian_eigen", as.matrix(A), PACKAGE = "eigencore")
}

#' @keywords internal
native_dense_complex_general_eigen <- function(A) {
  .Call("eigencore_dense_complex_general_eigen", as.matrix(A), PACKAGE = "eigencore")
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
native_dense_complex_svd <- function(A) {
  .Call("eigencore_dense_complex_svd", as.matrix(A), PACKAGE = "eigencore")
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
gram_svd_numeric_option <- function(name, default, min = 0, allow_infinite = TRUE) {
  value <- getOption(name, default)
  if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
      value < min || (!allow_infinite && !is.finite(value))) {
    stop(
      "Option ", name, " must be one numeric value >= ", min,
      if (allow_infinite) " or Inf." else ".",
      call. = FALSE
    )
  }
  as.numeric(value)
}

#' @keywords internal
gram_svd_max_dimension <- function() {
  value <- gram_svd_numeric_option(
    "eigencore.gram_svd_max_dimension",
    512,
    min = 1
  )
  if (is.finite(value)) as.integer(value) else Inf
}

#' @keywords internal
gram_svd_max_dimension_wide <- function() {
  value <- gram_svd_numeric_option(
    "eigencore.gram_svd_max_dimension_wide",
    1024,
    min = 1
  )
  if (is.finite(value)) as.integer(value) else Inf
}

#' @keywords internal
gram_svd_memory_budget_bytes <- function() {
  gram_svd_numeric_option("eigencore.gram_svd_memory_mb", 64, min = 0) * 1e6
}

#' @keywords internal
gram_svd_rank_fraction_limit <- function() {
  value <- gram_svd_numeric_option(
    "eigencore.gram_svd_rank_fraction_limit",
    0.5,
    min = .Machine$double.eps,
    allow_infinite = FALSE
  )
  if (value > 1) {
    stop("Option eigencore.gram_svd_rank_fraction_limit must be <= 1.", call. = FALSE)
  }
  value
}

#' @keywords internal
gram_svd_min_aspect_ratio <- function() {
  gram_svd_numeric_option(
    "eigencore.gram_svd_min_aspect_ratio",
    2,
    min = 1,
    allow_infinite = FALSE
  )
}

#' @keywords internal
gram_svd_work_budget_units <- function() {
  gram_svd_numeric_option("eigencore.gram_svd_work_budget", Inf, min = 0)
}

#' @keywords internal
gram_svd_policy <- function(dims, rank, target = largest()) {
  dims <- as.numeric(dims)
  if (length(dims) != 2L || any(!is.finite(dims)) || any(dims < 1)) {
    return(list(eligible = FALSE, rejection_reasons = "invalid_dimensions"))
  }
  rank <- as.integer(rank %||% 1L)
  if (length(rank) != 1L || is.na(rank) || rank < 1L) {
    return(list(eligible = FALSE, rejection_reasons = "invalid_rank"))
  }

  reduced <- min(dims)
  full <- max(dims)
  gram_max <- gram_svd_max_dimension()
  # Wide (left-Gram) largest-target rows pass the installed cutoff speed gate
  # up to side 1024 (inst/benchmarks/results/20260611-svd-gram-cutoff-*);
  # tall (right-Gram) rows do not, so the wider default is scoped to that
  # promoted regime. max() keeps an explicitly raised global option in force.
  target_kind <- if (inherits(target, "eigencore_target")) target$kind else "largest"
  if (dims[2L] > dims[1L] && target_kind %in% c("largest", "largest_magnitude")) {
    gram_max <- max(gram_max, gram_svd_max_dimension_wide())
  }
  memory_budget <- gram_svd_memory_budget_bytes()
  rank_limit <- gram_svd_rank_fraction_limit()
  min_aspect <- gram_svd_min_aspect_ratio()
  work_budget <- gram_svd_work_budget_units()

  gram_bytes <- reduced * reduced * 8
  result_vector_bytes <- rank * sum(dims) * 8
  total_materialization_bytes <- gram_bytes + result_vector_bytes
  eigensolve_work <- reduced^3
  rank_fraction <- rank / reduced
  aspect_ratio <- full / reduced

  rejection_reasons <- character()
  if (reduced > gram_max) {
    rejection_reasons <- c(rejection_reasons, "gram_dimension_exceeds_max")
  }
  if (is.finite(memory_budget) && total_materialization_bytes > memory_budget) {
    rejection_reasons <- c(rejection_reasons, "gram_memory_budget_exceeded")
  }
  if (is.finite(work_budget) && eigensolve_work > work_budget) {
    rejection_reasons <- c(rejection_reasons, "gram_work_budget_exceeded")
  }
  if (aspect_ratio < min_aspect) {
    rejection_reasons <- c(rejection_reasons, "aspect_ratio_below_minimum")
  }
  if (rank_fraction > rank_limit) {
    rejection_reasons <- c(rejection_reasons, "rank_fraction_exceeds_limit")
  }

  list(
    eligible = !length(rejection_reasons),
    rejection_reasons = rejection_reasons,
    gram_dimension = as.integer(reduced),
    full_dimension = as.integer(full),
    gram_max_dimension = gram_max,
    gram_memory_budget_bytes = memory_budget,
    gram_memory_budget_mb = memory_budget / 1e6,
    estimated_gram_bytes = gram_bytes,
    estimated_result_vector_bytes = result_vector_bytes,
    estimated_total_materialization_bytes = total_materialization_bytes,
    estimated_gram_eigensolve_work_units = eigensolve_work,
    gram_work_budget_units = work_budget,
    rank_fraction = rank_fraction,
    rank_fraction_limit = rank_limit,
    aspect_ratio = aspect_ratio,
    min_aspect_ratio = min_aspect,
    materialization_policy = "budgeted smaller Gram materialization"
  )
}

#' @keywords internal
gram_svd_plan_controls <- function(dims, rank, target = largest(),
                                   svd_partial_fastpath = FALSE) {
  policy <- gram_svd_policy(dims, rank, target)
  c(list(
    gram_side = if (dims[2L] <= dims[1L]) "right" else "left",
    gram_dimension = policy$gram_dimension,
    full_dimension = policy$full_dimension,
    gram_max_dimension = policy$gram_max_dimension,
    rank_fraction = policy$rank_fraction,
    rank_fraction_limit = policy$rank_fraction_limit,
    aspect_ratio = policy$aspect_ratio,
    min_aspect_ratio = policy$min_aspect_ratio,
    gram_memory_budget_mb = policy$gram_memory_budget_mb,
    gram_memory_budget_bytes = policy$gram_memory_budget_bytes,
    estimated_gram_bytes = policy$estimated_gram_bytes,
    estimated_result_vector_bytes = policy$estimated_result_vector_bytes,
    estimated_total_materialization_bytes = policy$estimated_total_materialization_bytes,
    estimated_gram_eigensolve_work_units =
      policy$estimated_gram_eigensolve_work_units,
    gram_work_budget_units = policy$gram_work_budget_units,
    materialization_policy = policy$materialization_policy,
    gram_policy_passed = isTRUE(policy$eligible),
    gram_policy_rejection = paste(policy$rejection_reasons, collapse = ";"),
    certified_in_original_coordinates = TRUE,
    materializes = "smaller Gram matrix only",
    fallback_policy = "certification-gated",
    runtime_fallback = "native Golub-Kahan if original-coordinate certificate is weaker",
    fallback_requires_vectors = "both"
  ), if (isTRUE(svd_partial_fastpath)) {
    list(svd_partial_fastpath = TRUE)
  } else {
    list()
  })
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
  kind <- if (inherits(problem$target, "eigencore_target")) problem$target$kind else "largest"
  if (kind %in% c("smallest", "smallest_magnitude") &&
      !identical(native_kernel_kind(problem$A), "csc")) {
    return(FALSE)
  }
  dims <- problem$A$dim
  rank <- as.integer(rank %||% problem$rank %||% 1L)
  isTRUE(gram_svd_policy(dims, rank, problem$target)$eligible)
}
