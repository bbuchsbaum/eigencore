#' Compute a partial eigendecomposition.
#'
#' @param A Matrix or eigencore operator.
#' @param k Number of eigenpairs to compute. Required except for an
#'   [interval()] target, whose count comes from an inertia factorisation;
#'   there a supplied `k` is an upper bound (an error if the interval holds
#'   more eigenvalues).
#' @param target Eigencore eigenvalue target descriptor. With
#'   `method = shift_invert(sigma)` it defaults to `nearest(sigma)`; other
#'   targets (except `smallest_magnitude()` with `sigma = 0`) are an error.
#' @param B Optional metric matrix or operator for generalized problems.
#' @param method Solver method descriptor.
#' @param tol Convergence and certification tolerance.
#' @param maxit Optional iteration limit (`NULL` uses each route's default).
#'   It bounds outer iterations, never the Krylov subspace size: thick-restart
#'   cycles for (block, generalized, and shift-invert callback) Lanczos,
#'   Krylov-Schur restarts for nonsymmetric Arnoldi (including shift-invert
#'   Arnoldi), restart cycles for the reference Arnoldi, LOBPCG iterations,
#'   and Lanczos steps for the unrestarted reference and native-kernel
#'   shift-invert Lanczos routes. Dense direct routes ignore it. The resolved
#'   limit is recorded in `plan$controls$iteration_limit` and
#'   `plan$controls$iteration_limit_kind`. To set the subspace size (the
#'   ARPACK `ncv`), use the method descriptor's `max_subspace`
#'   ([lanczos()], [auto()], [shift_invert()]). A `lanczos(max_restarts =)` or
#'   `lobpcg(maxit =)` that disagrees with `maxit` is an error.
#' @param vectors Whether to compute vectors.
#' @param seed Optional random seed for stochastic solver components. The
#'   global random number stream is restored on exit.
#' @param certify Whether to compute certification diagnostics.
#' @param allow_dense_fallback Dense fallback policy.
#' @param left_vectors Left-eigenvector policy for nonsymmetric problems.
#'   `"auto"` (default) computes and certifies left eigenvectors (and the
#'   biorthogonality) on routes that support them, such as Krylov-Schur
#'   Arnoldi, and records the reason when they are unavailable; `"none"` skips
#'   the left solve and its certificate entirely (roughly halving the cost of a
#'   nonsymmetric Arnoldi solve); `"compute"` is like `"auto"` but an error
#'   when the route cannot return left eigenvectors. Hermitian problems are
#'   unaffected.
#' @param initial_subspace Optional numeric matrix of starting directions
#'   (a warm start). Supported on standard real Hermitian Lanczos paths: the
#'   native paths for explicit dense double or `dgCMatrix` operators, the native
#'   matrix-free callback path selected by `lanczos(block > 1)`, and the scalar
#'   matrix-free reference path selected by `lanczos(block = 1)`;
#'   supplying it on any other planned path (generalized, shift-invert, dense
#'   fallback) is an error. Pass `method = lanczos()` to guarantee a Lanczos
#'   route: with the default `method = auto()`, sparse or `nearest()` problems
#'   may be planned as shift-invert, which does not consume a start and will
#'   reject the argument. The subspace is only a starting hint: projected
#'   quantities, residuals, orthogonality, convergence, and the certificate are
#'   recomputed for the current operator on every solve. The columns are
#'   orthonormalized at the solver boundary and fitted to the method's start
#'   block — when the accepted rank exceeds the block width the block is a
#'   seeded random rotation of the full accepted basis, so every supplied
#'   direction contributes. Because a residual certificate proves eigenpair
#'   accuracy but not target identity, a fully supplied subspace that is already
#'   invariant at `tol` is discarded in favor of a cold start; provenance
#'   records that guard decision. Diagnostics distinguish operator block calls,
#'   operator columns, and certification columns. `NULL` (the default)
#'   preserves the cold random start exactly.
#' @return An `eigencore_eigen_result` containing computed values, optional
#'   vectors, certificate diagnostics, method/plan metadata, and convergence
#'   diagnostics.
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
eig_partial <- function(A, k = NULL, target = largest(), B = NULL, method = auto(),
                        tol = 1e-8, maxit = NULL, vectors = TRUE, seed = NULL,
                        certify = TRUE,
                        allow_dense_fallback = c("auto", "never", "always"),
                        initial_subspace = NULL,
                        left_vectors = c("auto", "none", "compute")) {
  allow_dense_fallback <- match.arg(allow_dense_fallback)
  left_vectors <- match.arg(left_vectors)
  if (missing(target) && inherits(method, "eigencore_method") &&
      identical(method$kind, "shift_invert") && !is.null(method$sigma)) {
    # shift_invert(sigma) computes the eigenvalues nearest sigma.
    target <- nearest(method$sigma)
  }
  if (!is.null(seed)) {
    seed_state <- saved_random_seed()
    on.exit(restore_random_seed(seed_state), add = TRUE)
    set.seed(seed)
  }
  P <- eigen_problem(A, metric = B, target = target,
                     transform = if (is_transform_method(method)) method else NULL)
  solve(P, k = k, method = method, tol = tol, maxit = maxit, vectors = vectors,
        certify = certify, allow_dense_fallback = allow_dense_fallback,
        initial_subspace = initial_subspace, left_vectors = left_vectors)
}

#' Compute a partial singular-value decomposition.
#'
#' @param A Matrix or eigencore operator.
#' @param rank Number of singular values to compute.
#' @param target Eigencore singular-value target descriptor.
#' @param method Solver method descriptor.
#' @param tol Convergence and certification tolerance.
#' @param vectors Which singular-vector sides to compute.
#' @param seed Optional random seed for stochastic solver components. The
#'   global random number stream is restored on exit.
#' @param certify Whether to compute certification diagnostics.
#' @param allow_dense_fallback Dense fallback policy.
#' @return An `eigencore_svd_result` containing singular values, optional left
#'   and right singular vectors, certificate diagnostics, method/plan metadata,
#'   and convergence diagnostics.
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
  dims <- if (inherits(A, "eigencore_operator")) A$dim else dim(A)
  if (length(dims) == 2L) {
    rank <- validate_solution_count(rank, min(dims), "rank")
  }
  fast_started <- proc.time()[["elapsed"]]
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
    work_values <- unclass(fast$work)
    work_values$total_seconds <- proc.time()[["elapsed"]] - fast_started
    fast$work <- new_typed_work_record(work_values)
    return(fast)
  }
  if (!is.null(seed)) {
    seed_state <- saved_random_seed()
    on.exit(restore_random_seed(seed_state), add = TRUE)
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
#'
#' @param a Eigencore eigen problem object.
#' @param b Unused second argument reserved by the base [solve()] generic.
#' @param k Number of eigenpairs to compute.
#' @param method Solver method descriptor.
#' @param tol Convergence and certification tolerance.
#' @param maxit Optional iteration (restart) limit; never the subspace size.
#'   See [eig_partial()].
#' @param vectors Whether to compute vectors.
#' @param certify Whether to compute certification diagnostics.
#' @param allow_dense_fallback Dense fallback policy.
#' @param left_vectors Left-eigenvector policy (`"auto"`, `"none"`,
#'   `"compute"`); see [eig_partial()].
#' @param initial_subspace Optional numeric matrix of starting directions
#'   (a warm start). Supported on standard real Hermitian Lanczos paths —
#'   native dense double / `dgCMatrix`, native matrix-free callbacks for
#'   `lanczos(block > 1)`, and the scalar matrix-free reference path for
#'   `lanczos(block = 1)`; supplying it on any other planned path is an error.
#'   The subspace is only
#'   a starting hint, never a source of reused convergence: every solve
#'   recomputes projected quantities, residuals, orthogonality, convergence,
#'   and a fresh current-operator certificate. An already-invariant supplied
#'   subspace is discarded to a cold start because residual certification alone
#'   cannot establish that it contains the requested extremal eigenpairs.
#'   `NULL` (the default) preserves the cold random start exactly.
#' @param ... Reserved for future solver options.
#' @return An `eigencore_eigen_result`.
#' @export
solve.eigencore_eigen_problem <- function(a, b, k = NULL, method = auto(), tol = 1e-8,
                                          maxit = NULL, vectors = TRUE,
                                          certify = TRUE,
                                          allow_dense_fallback = c("auto", "never", "always"),
                                          initial_subspace = NULL,
                                          left_vectors = c("auto", "none", "compute"),
                                          ...) {
  allow_dense_fallback <- match.arg(allow_dense_fallback)
  left_vectors <- match.arg(left_vectors)
  plan <- plan_solver(
    a,
    k = k,
    method = method,
    tol = tol,
    maxit = maxit,
    vectors = vectors,
    certify = certify,
    allow_dense_fallback = allow_dense_fallback,
    initial_subspace = initial_subspace,
    left_vectors = left_vectors
  )
  solve(plan)
}

#' Solve a planned SVD problem.
#'
#' S3 method that runs the planned solver for an SVD problem built by
#' [svd_problem()]. Most users call [svd_partial()], which constructs the
#' problem and dispatches here; call `solve()` directly when you want to build
#' a problem once and reuse or inspect it. Returns a certified partial
#' singular-value decomposition.
#'
#' @param a Eigencore SVD problem object.
#' @param b Unused second argument reserved by the base [solve()] generic.
#' @param rank Number of singular values to compute.
#' @param method Solver method descriptor.
#' @param tol Convergence and certification tolerance.
#' @param vectors Which singular-vector sides to compute.
#' @param certify Whether to compute certification diagnostics.
#' @param allow_dense_fallback Dense fallback policy.
#' @param ... Reserved for future solver options.
#' @return An `eigencore_svd_result`.
#' @export
solve.eigencore_svd_problem <- function(a, b, rank, method = auto(), tol = 1e-8,
                                        vectors = c("both", "left", "right", "none"),
                                        certify = TRUE,
                                        allow_dense_fallback = c("auto", "never", "always"), ...) {
  vectors <- match.arg(vectors)
  allow_dense_fallback <- match.arg(allow_dense_fallback)
  plan <- plan_solver(
    a,
    rank = rank,
    method = method,
    tol = tol,
    vectors = vectors,
    certify = certify,
    allow_dense_fallback = allow_dense_fallback
  )
  solve(plan)
}

#' Execute a frozen eigencore solver plan.
#'
#' `solve(plan)` validates and runs the problem, method, controls, policy, and
#' execution arguments captured by [plan_solver()]. It never invokes the
#' planner unless `replan = TRUE`. Execution arguments cannot be overridden
#' through `...`.
#'
#' @param a An executable `eigencore_plan`.
#' @param b Unused second argument reserved by the base [solve()] generic.
#' @param restart_state Optional [restart_state()] built from a certified result.
#'   Version 1 accepts the public basis only on admitted standard real Hermitian
#'   Lanczos routes and validates it before applying the destination operator.
#' @param reuse Reuse policy for an optional restart state. `"basis_only"`
#'   always fits the public basis to the receiving method, `"same_operator"`
#'   requires a matching operator identity and compatible retained method state,
#'   and `"auto"` uses eligible same-operator method state or downgrades to the
#'   public basis. No policy reuses convergence, locked vectors, operator
#'   actions, residuals, or a certificate.
#' @param retain_state State-retention request for the returned result. `"basis"`
#'   retains only the certified public basis; `"same_operator"` also requests
#'   an eligible versioned method payload. Unsupported producer routes still
#'   return a basis-only state.
#' @param replan Whether to create and execute a new plan from the embedded
#'   problem under current policy. The original plan is not modified.
#' @param ... Must be empty; execution controls are frozen in the plan.
#' @return An `eigencore_eigen_result` or `eigencore_svd_result` containing the
#'   exact plan that governed execution.
#' @export
solve.eigencore_plan <- function(
    a, b, restart_state = NULL,
    reuse = c("auto", "basis_only", "same_operator"),
    retain_state = c("none", "basis", "same_operator"),
    replan = FALSE, ...) {
  total_started <- proc.time()[["elapsed"]]
  setup_started <- total_started
  reuse <- match.arg(reuse)
  retain_state <- match.arg(retain_state)
  dots <- list(...)
  if (length(dots)) {
    field <- names(dots)[[1L]] %||% "..."
    if (!nzchar(field)) {
      field <- "..."
    }
    plan_error(
      "execution_override", field, "frozen plan value", dots[[1L]],
      paste0("Execution argument '", field,
             "' is frozen in the plan. Create a new plan or use replan = TRUE.")
    )
  }
  if (!is.logical(replan) || length(replan) != 1L || is.na(replan)) {
    plan_error("corrupt_plan", "replan", "TRUE or FALSE", replan)
  }
  validate_eigencore_plan(a)
  if (isTRUE(replan)) {
    a <- replan_eigencore_plan(a)
    validate_eigencore_plan(a)
  }
  restart_preparation <- prepare_restart_state_for_plan(
    restart_state, a, reuse = reuse
  )
  setup_seconds <- proc.time()[["elapsed"]] - setup_started
  context <- new_work_context(a)
  execution_started <- proc.time()[["elapsed"]]
  result <- with_work_context(context, with_execution_policy(a$planner_policy, {
    if (identical(a$problem_type, "eigen")) {
      execute_eigen_plan(a, restart_preparation = restart_preparation)
    } else {
      execute_svd_plan(a)
    }
  }))
  # An interval() target has no fixed k: its count comes from the inertia
  # factorisation and the completeness check verifies the returned set
  # against it, so `requested` (an upper bound) must not withhold it.
  if (!identical(a$problem$target$kind, "interval")) {
    result <- withhold_short_certificate(result, a$requested)
  }
  result <- require_verified_completeness(result)
  finished <- proc.time()[["elapsed"]]
  result$work <- finalize_work_record(
    result,
    context,
    setup_seconds = setup_seconds,
    execution_seconds = finished - execution_started,
    total_seconds = finished - total_started
  )
  result$state_transition <- finalize_restart_transition(
    result, restart_preparation
  )
  result["restart_state"] <- list(
    if (identical(retain_state, "none")) {
      NULL
    } else {
      construct_restart_state(result, retention = retain_state)
    }
  )
  result$memory <- result_memory_record(result)
  result
}

# A certificate proves the returned pairs; it cannot pass for a request it
# did not fill. Krylov routes return fewer than k pairs when the Krylov space
# is exhausted (e.g. shift-invert on a spectrum with repeated eigenvalues:
# 1, 7, 9 for k = 6 of {1,1,7,7,7,9,9,9}); found by the oracle sweep.
#' @keywords internal
# Completeness states that prove (exact, inertia_verified) or give strong
# evidence (probed, repaired) that the returned set is the requested one.
#' @keywords internal
verified_completeness_states <- function() {
  c("exact", "inertia_verified", "probed", "repaired")
}

# `certificate$passed` means "the requested pairs, each accurate": the residual
# certificate passed AND the returned set was verified to be the requested
# one. A result whose residuals pass but whose set was not verified keeps
# residual_passed = TRUE and reports passed = FALSE with a note naming the
# completeness status. Set options(eigencore.require_completeness = FALSE) to
# restore the residual-only meaning of `passed`.
#' @keywords internal
require_verified_completeness <- function(result) {
  cert <- result$certificate
  if (is.null(cert) || !isTRUE(getOption("eigencore.require_completeness", TRUE))) {
    return(result)
  }
  cert$residual_passed <- cert$residual_passed %||% isTRUE(cert$passed)
  status <- cert$target_completeness %||% "not_checked"
  cert$target_completeness <- status
  verified <- status %in% verified_completeness_states()
  cert$target_passed <- cert$target_passed %||%
    (if (verified) TRUE else if (status %in% c("failed", "inertia_failed")) FALSE else NA)
  if (isTRUE(cert$passed) && !verified) {
    cert$passed <- FALSE
    note <- paste0(
      "residuals certified, but the returned set was not verified to be the ",
      "requested one (target_completeness = \"", status, "\"); passed = FALSE"
    )
    cert$notes <- unique(c(cert$notes, note))
  }
  result$certificate <- cert
  result
}

withhold_short_certificate <- function(result, k) {
  cert <- result$certificate
  returned <- length(result$values)
  if (is.null(k) || is.null(cert) || !isTRUE(cert$passed) || returned >= k) {
    return(result)
  }
  note <- sprintf(
    "only %d of the %d requested pairs were returned; certificate withheld",
    returned, as.integer(k)
  )
  cert$passed <- FALSE
  cert$notes <- unique(c(cert$notes, note))
  if (!is.null(cert$target_passed) && isTRUE(cert$target_passed)) {
    cert$target_passed <- FALSE
  }
  result$certificate <- cert
  result$warnings <- unique(c(result$warnings, note))
  result
}

#' @keywords internal
replan_eigencore_plan <- function(plan) {
  execution <- plan$execution
  args <- list(
    problem = plan$problem,
    method = plan$method_descriptor,
    tol = execution$tol,
    vectors = execution$vectors,
    certify = execution$certify,
    allow_dense_fallback = execution$allow_dense_fallback
  )
  if (identical(plan$problem_type, "eigen")) {
    args$k <- plan$requested
    args$maxit <- execution$maxit
    args$initial_subspace <- execution$initial_subspace
    args$left_vectors <- execution$left_vectors %||% "auto"
  } else {
    args$rank <- plan$requested
  }
  do.call(plan_solver, args)
}

#' @keywords internal
execute_eigen_plan <- function(plan, restart_preparation = NULL) {
  # Target completeness (C50): residual certificates prove each returned pair
  # but not that the returned set is the requested one. Eligible Krylov results
  # are probed on the deflated complement after dispatch; the probe needs the
  # Ritz vectors, so they are kept internally and dropped afterwards when the
  # caller asked for values only.
  a <- plan$problem
  if (is_interval_target(a$target)) {
    return(solve_interval_eigen(plan))
  }
  k <- plan$requested
  vectors_requested <- isTRUE(plan$execution$vectors)
  mode <- target_completeness_mode(plan$method_descriptor)
  need_vectors <- !identical(mode, "none") &&
    isTRUE(plan$execution$certify) &&
    identical(target_completeness_route_class(plan, a), "krylov") &&
    (target_completeness_eligible(a, k) ||
       (mode %in% c("auto", "inertia") && inertia_completeness_eligible(a, k)))
  need_vectors <- need_vectors ||
    hermitian_completeness_wants_vectors(plan, a, k, mode)
  started <- proc.time()[["elapsed"]]
  result <- execute_eigen_plan_dispatch(
    plan, restart_preparation = restart_preparation,
    vectors = vectors_requested || need_vectors
  )
  solve_seconds <- proc.time()[["elapsed"]] - started
  seed <- result$transform$inertia_seed %||% NULL
  if (!is.null(result$transform)) {
    result$transform["inertia_seed"] <- NULL
  }
  result <- apply_target_completeness(result, plan, a, k, mode, vectors_requested,
                                      solve_seconds = solve_seconds, seed = seed)
  if (!is.null(result$transform)) {
    result$transform["completeness_solve"] <- NULL
  }
  result
}

#' @keywords internal
execute_eigen_plan_dispatch <- function(plan, restart_preparation = NULL,
                                        vectors = plan$execution$vectors) {
  a <- plan$problem
  execution <- plan$execution
  k <- plan$requested
  method <- plan$method_descriptor
  tol <- execution$tol
  maxit <- execution$maxit
  certify <- execution$certify
  allow_dense_fallback <- execution$allow_dense_fallback
  prepared_restart <- restart_preparation$prepared_start %||% NULL
  initial_subspace <- if (is.null(prepared_restart)) {
    execution$initial_subspace
  } else {
    NULL
  }
  if (is_transform_method(a$transform) &&
      identical(a$transform$kind, "shift_invert")) {
    a <- shift_invert_prepare_tridiagonal(a)
  }
  validate_complex_eigen_plan(a, plan)
  if (!is.null(initial_subspace)) {
    validate_initial_subspace_plan_support(a, plan)
  }
  if (is_transform_method(a$transform)) {
    if (identical(a$transform$kind, "shift_invert") &&
        plan$method %in% shift_invert_arnoldi_labels()) {
      return(solve_shift_invert_general(
        a, k = k, method = a$transform, tol = tol, vectors = vectors,
        certify = certify, plan = plan,
        left_vectors = execution$left_vectors %||% "auto"
      ))
    }
    if (identical(a$transform$kind, "shift_invert")) {
      return(solve_shift_invert_hermitian(
        a, k = k, method = a$transform, tol = tol,
        maxit = maxit, vectors = vectors, certify = certify, plan = plan
      ))
    }
    stop(
      "transform method '", a$transform$kind,
      "' is registered in transform_method_kinds() but has no executable plan dispatch.",
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
    return(solve_eigen_lanczos(
      a, k, method, tol, maxit, vectors, certify, plan,
      initial_subspace = initial_subspace,
      prepared_restart = prepared_restart
    ))
  }
  if (plan_dispatches_sparse_general_pencil_arnoldi(plan)) {
    return(solve_eigen_sparse_general_pencil_arnoldi(
      a, k, method, tol, maxit, vectors, certify, plan
    ))
  }
  if (plan_rejects_sparse_general_pencil(plan)) {
    stop(sparse_general_pencil_unsupported_message(a), call. = FALSE)
  }
  if (plan_dispatches_arnoldi(plan)) {
    return(solve_eigen_arnoldi(a, k, method, tol, maxit, vectors, certify, plan))
  }
  if (plan$method %in% c("native dense Hermitian LAPACK fallback",
                         native_dense_complex_hermitian_label())) {
    return(solve_eigen_native_dense_hermitian(
      a, k, tol, vectors, certify, allow_dense_fallback, plan
    ))
  }
  if (identical(plan$method, native_dense_complex_general_label())) {
    return(solve_eigen_native_dense_general(
      a, k, tol, vectors, certify, allow_dense_fallback, plan
    ))
  }
  solve_eigen_dense_oracle(
    a, k, tol, vectors, certify, allow_dense_fallback, plan
  )
}

#' @keywords internal
execute_svd_plan <- function(plan) {
  a <- plan$problem
  execution <- plan$execution
  rank <- plan$requested
  method <- plan$method_descriptor
  tol <- execution$tol
  vectors <- execution$vectors
  certify <- execution$certify
  allow_dense_fallback <- execution$allow_dense_fallback
  validate_complex_svd_plan(a, plan)
  if (identical(plan$method, "reference randomized SVD prototype") ||
      identical(plan$method, native_dense_randomized_svd_label()) ||
      identical(plan$method, native_csc_randomized_svd_label())) {
    return(solve_svd_randomized(a, rank, method, tol, vectors, certify, plan))
  }
  if (identical(plan$method, "native certified Gram SVD special case")) {
    return(solve_svd_gram(a, rank, tol, vectors, certify, plan))
  }
  if (identical(plan$method, native_implicit_gram_svd_label())) {
    return(solve_svd_implicit_gram(a, rank, tol, vectors, certify, plan))
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
  # Cholesky reduction + zheev yields B-orthonormal vectors even for repeated
  # eigenvalues; zggev vectors are only B-normalized column by column.
  .Call("eigencore_dense_complex_generalized_hpd_eigen",
        eig_full_as_complex_matrix(as.matrix(A)),
        eig_full_as_complex_matrix(as.matrix(B)),
        PACKAGE = "eigencore")
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
    interval = order(Re(x), decreasing = FALSE),
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
plan_dispatches_sparse_general_pencil_arnoldi <- function(plan) {
  identical(plan$method, sparse_general_pencil_diagonal_arnoldi_label())
}

#' @keywords internal
plan_rejects_sparse_general_pencil <- function(plan) {
  identical(plan$method, sparse_general_pencil_unsupported_label())
}

#' @keywords internal
plan_dispatches_lanczos <- function(plan) {
  plan$method %in% c(
    "native scalar thick-restart Hermitian Lanczos",
    "native block Hermitian Lanczos thick-restart candidate",
    "native block Hermitian Lanczos (thick restart, locking)",
    native_matrix_free_block_lanczos_label(),
    native_generalized_lanczos_label(),
    generalized_lanczos_label(),
    native_both_ends_lanczos_label(),
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
    "native block Hermitian Lanczos (thick restart, locking)",
    native_matrix_free_block_lanczos_label()
  )
}

#' @keywords internal
plan_dispatches_reference_hermitian_lanczos <- function(plan) {
  plan$method %in% c(
    "reference Hermitian Lanczos (target unsupported by native path)",
    "reference Hermitian Lanczos (prototype/oracle fallback)"
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
plan_dispatch_available <- function(plan) {
  if (identical(plan$problem_type, "eigen")) {
    problem <- plan$problem
    if (is_interval_target(problem$target)) {
      return(plan$method %in% interval_route_labels())
    }
    if (is_transform_method(problem$transform)) {
      return(identical(problem$transform$kind, "shift_invert"))
    }
    if (plan_dispatches_lobpcg(plan) ||
        plan_dispatches_structured_grid_laplacian_2d(plan) ||
        plan_dispatches_native_tridiagonal_eigen(plan) ||
        plan_dispatches_lanczos(plan) ||
        plan_dispatches_sparse_general_pencil_arnoldi(plan) ||
        plan_rejects_sparse_general_pencil(plan) ||
        plan_dispatches_arnoldi(plan)) {
      return(TRUE)
    }
    if (plan$method %in% c(
      "native dense Hermitian LAPACK fallback",
      native_dense_complex_hermitian_label(),
      native_dense_complex_general_label(),
      "native dense generalized SPD LAPACK fallback"
    )) {
      return(TRUE)
    }
    return(grepl(
      "^(dense LAPACK|dense generalized SPD LAPACK)",
      plan$method
    ))
  }
  if (!identical(plan$problem_type, "svd")) {
    return(FALSE)
  }
  if (plan$method %in% c(
    "reference randomized SVD prototype",
    native_dense_randomized_svd_label(),
    native_csc_randomized_svd_label(),
    "native certified Gram SVD special case",
    native_implicit_gram_svd_label(),
    native_retained_golub_kahan_diagnostic_label(),
    "native dense LAPACK SVD fallback",
    native_dense_complex_svd_label()
  )) {
    return(TRUE)
  }
  plan_dispatches_golub_kahan(plan) ||
    grepl("^dense LAPACK SVD", plan$method)
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
  plan$planned_method <- plan$method
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
  plan$serialization <- new_plan_serialization(
    plan$operator_identity,
    plan$planned_method,
    method_descriptor = plan$method_descriptor,
    controls = plan$controls,
    execution = plan$execution,
    planner_policy = plan$planner_policy
  )
  plan$memory <- new_memory_record(list(
    problem = plan$problem,
    controls = plan$controls,
    policy = plan$planner_policy,
    metadata = list(
      schema_version = plan$schema_version,
      method_descriptor = plan$method_descriptor,
      execution = plan$execution,
      serialization = plan$serialization
    )
  ))
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
  if (!(native_kernel_kind(problem$A) %in% c("csc", "dense"))) {
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
native_dense_symmetric_eigen <- function(A, vectors = TRUE) {
  # dsyevr (MRRR, the driver base eigen() uses), ascending values; dsyev only
  # as an internal-failure fallback (`driver` says which ran). With
  # vectors = FALSE only eigenvalues are computed and `vectors` is NULL.
  .Call("eigencore_dense_symmetric_eigen", as.matrix(A), isTRUE(vectors),
        PACKAGE = "eigencore")
}

#' @keywords internal
native_dense_symmetric_eigen_dsyev <- function(A) {
  .Call("eigencore_dense_symmetric_eigen_dsyev", as.matrix(A), PACKAGE = "eigencore")
}

#' @keywords internal
native_dense_symmetric_eigen_dsyevd <- function(A) {
  .Call("eigencore_dense_symmetric_eigen_dsyevd", as.matrix(A), PACKAGE = "eigencore")
}

#' @keywords internal
native_dense_symmetric_eigen_selected <- function(A, k, target, vectors = TRUE) {
  .Call(
    "eigencore_dense_symmetric_eigen_selected",
    as.matrix(A),
    as.integer(k),
    as.integer(lanczos_target_kind(target)),
    isTRUE(vectors),
    PACKAGE = "eigencore"
  )
}

#' @keywords internal
native_dense_complex_hermitian_eigen <- function(A, vectors = TRUE) {
  # zheev returns ascending real values and unitary vectors, so repeated
  # eigenvalues keep an orthonormal eigenbasis (zgeev does not). R's LAPACK
  # header does not declare zheevr, so the MRRR driver is real-only.
  .Call("eigencore_dense_complex_hermitian_eigen",
        eig_full_as_complex_matrix(as.matrix(A)),
        isTRUE(vectors),
        PACKAGE = "eigencore")
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
  # Thin SVD by dgesdd (divide and conquer, as base svd()); dgesvd only when
  # dgesdd fails to converge. `driver` records which LAPACK routine ran.
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
    native_kernel_kind(problem$A) %in% c("csc", "centered_scaled_csc")
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
  if (!(native_kernel_kind(problem$A) %in% c("csc", "dense"))) {
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

#' @keywords internal
should_use_native_implicit_gram_svd <- function(problem, method, rank = NULL) {
  if (!inherits(method, "eigencore_method") || !identical(method$kind, "auto")) {
    return(FALSE)
  }
  if (!implicit_gram_svd_target_supported(problem$target)) {
    return(FALSE)
  }
  storage <- problem$A$metadata$storage %||% NULL
  source <- source_or_null(problem$A)
  is_csc <- identical(storage, "dgCMatrix")
  is_dense <- is.matrix(source) && is.double(source) && !is.complex(source)
  if (!is_csc && !is_dense) {
    return(FALSE)
  }
  dims <- as.integer(problem$A$dim)
  reduced <- min(dims)
  rank <- as.integer(rank %||% problem$rank %||% 1L)
  if (reduced < getOption("eigencore.implicit_gram_svd_min_dimension", 64L)) {
    # Tiny problems: the explicit Gram or dense LAPACK path is already cheap.
    return(FALSE)
  }
  # Beyond half the reduced dimension the Lanczos subspace approaches the
  # full space; the explicit paths are the honest choice there.
  rank >= 1L && rank <= 0.5 * reduced
}

#' @keywords internal
saved_random_seed <- function() {
  get0(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
}

#' @keywords internal
restore_random_seed <- function(seed_state) {
  if (is.null(seed_state)) {
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  } else {
    assign(".Random.seed", seed_state, envir = .GlobalEnv)
  }
  invisible(NULL)
}
