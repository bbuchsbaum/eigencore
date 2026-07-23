# Warm-start (initial_subspace) plumbing for standard Hermitian Lanczos.
#
# Scope (planning/prd.json `initial_subspace_contract`): expose the Hermitian
# Lanczos start-block capability behind a public `initial_subspace` argument,
# on native dense/CSC paths and on the matrix-free reference Hermitian Lanczos
# path (the spectral-continuation surface for operator-only `A - rho * B`
# sequences). The subspace is only a starting hint; every solve recomputes
# projected quantities, residuals, orthogonality, convergence, and a fresh
# current-operator certificate. Reusable restart objects, recurrence reuse,
# and generalized/transformed promotion are V3.

#' Cold-start provenance record.
#'
#' Returned when no `initial_subspace` is supplied so downstream diagnostics
#' carry a uniform provenance schema for cold and warm solves alike.
#' @return A named list with `start_source = "cold"` and zeroed
#'   supplied/accepted/rejected/augmented/rank counts.
#' @keywords internal
warm_start_cold_provenance <- function() {
  list(
    start_source = "cold",
    supplied = 0L,
    accepted = 0L,
    rejected = 0L,
    augmented = 0L,
    rank = 0L,
    compressed = FALSE,
    invariant_guard_used = FALSE,
    invariant_relative_residual = NA_real_,
    guard_operator_block_calls = 0L,
    guard_operator_columns = 0L
  )
}

#' Whether the resolved plan consumes a user-supplied starting subspace.
#'
#' The boundary is standard (non-generalized, non-transformed) real Hermitian
#' Lanczos: the native dense double / dgCMatrix CSC paths plus the matrix-free
#' reference Hermitian Lanczos path. Every other dispatch is out of scope and
#' must reject `initial_subspace` rather than ignore it, densify, or borrow a
#' production label.
#' @return `TRUE` if the plan consumes a user-supplied start block, else `FALSE`.
#' @keywords internal
warm_start_plan_consumes_start <- function(problem, plan) {
  is.null(problem$metric) &&
    !is_transform_method(problem$transform) &&
    (plan_dispatches_native_lanczos(plan) ||
       plan_dispatches_reference_hermitian_lanczos(plan))
}

#' Guard: error when `initial_subspace` reaches an unsupported plan.
#' @return Invisibly `TRUE` when the plan is supported; otherwise throws an
#'   error.
#' @keywords internal
validate_initial_subspace_plan_support <- function(problem, plan) {
  if (warm_start_plan_consumes_start(problem, plan)) {
    return(invisible(TRUE))
  }
  stop(
    "initial_subspace was supplied but the resolved plan '", plan$method,
    "' does not consume a starting subspace. The warm-start seam is limited ",
    "to standard real Hermitian Lanczos (native dense double / dgCMatrix ",
    "paths and the matrix-free reference path); pass method = lanczos() on ",
    "such an operator, or omit initial_subspace.",
    call. = FALSE
  )
}

#' Validate, orthonormalize, and fit a supplied subspace to a start block.
#'
#' Produces an `n x width` start block for the Lanczos paths together with
#' provenance counts. Orthonormalization happens here at the solver boundary
#' for rank detection and honest reporting; the native kernel
#' re-orthonormalizes and rank-guards the block again before iterating.
#'
#' When the accepted numerical rank exceeds `width` (e.g. a k-column
#' continuation subspace handed to a scalar method), the start block is a
#' seeded random rotation of the full accepted basis rather than a truncation:
#' every supplied direction contributes generic weight to every start column,
#' so a width-1 start still overlaps all k target directions instead of only
#' the first. Provenance records this as `compressed = TRUE`.
#'
#' Seed policy: augmented directions are drawn from the active RNG stream
#' (controlled by the `seed` argument of [eig_partial()] / [solve()]) and
#' orthonormalized against the accepted directions, matching how the cold
#' random start is generated.
#'
#' @param initial_subspace User-supplied numeric matrix (or vector) of start
#'   directions.
#' @param n Operator domain dimension.
#' @param width Number of start-block columns the chosen method consumes
#'   (the Lanczos block size).
#' @param tol Numerical-rank tolerance for boundary orthonormalization.
#' @return A list with `start` (the `n x width` block) and provenance fields
#'   `start_source`, `supplied`, `accepted`, `rejected`, `augmented`, `rank`,
#'   and `compressed`. The internal `accepted_basis` field retains the full
#'   accepted basis for the invariant-subspace safety guard.
#' @keywords internal
prepare_initial_subspace <- function(initial_subspace, n, width,
                                     tol = sqrt(.Machine$double.eps)) {
  width <- as.integer(width)
  if (length(width) != 1L || is.na(width) || width < 1L) {
    stop("Internal error: warm-start block width must be a positive integer.",
         call. = FALSE)
  }
  if (width > n) {
    stop("Warm-start block width cannot exceed the operator dimension.",
         call. = FALSE)
  }

  x <- initial_subspace
  if (is.null(dim(x))) {
    x <- matrix(x, ncol = 1L)
  }
  x <- as.matrix(x)
  if (!is.numeric(x)) {
    stop("initial_subspace must be numeric (a matrix or vector).", call. = FALSE)
  }
  storage.mode(x) <- "double"
  if (nrow(x) != n) {
    stop("initial_subspace has ", nrow(x), " rows but the operator domain has ",
         "dimension ", n, ".", call. = FALSE)
  }
  if (ncol(x) < 1L) {
    stop("initial_subspace must have at least one column.", call. = FALSE)
  }
  if (!all(is.finite(x))) {
    stop("initial_subspace must contain only finite numeric entries.",
         call. = FALSE)
  }

  supplied <- ncol(x)

  # A subspace is invariant to per-column scaling. Normalize robustly before
  # rank detection so tiny but valid directions are not mistaken for zeros and
  # huge directions cannot overflow the norm calculation.
  normalized <- matrix(0, nrow = n, ncol = supplied)
  for (j in seq_len(supplied)) {
    scale <- max(abs(x[, j]))
    if (is.finite(scale) && scale > 0) {
      scaled <- x[, j] / scale
      scaled_norm <- sqrt(sum(scaled^2))
      if (is.finite(scaled_norm) && scaled_norm > 0) {
        normalized[, j] <- scaled / scaled_norm
      }
    }
  }

  # Boundary orthonormalization; MGS2 now sees unit-scale columns, so `tol`
  # expresses a relative numerical-rank threshold.
  ortho <- mgs2(normalized, tol = tol)
  accepted_basis <- ortho$Q
  rank <- ncol(accepted_basis)

  accepted <- rank
  compressed <- rank > width
  start <- if (compressed) {
    # Rotate, don't truncate: mix all accepted directions into the narrower
    # start block via a seeded Gaussian rotation so no supplied direction is
    # wholly discarded. Falls back to truncation on the measure-zero event
    # that the rotated block loses rank.
    rotation <- matrix(stats::rnorm(rank * width), nrow = rank, ncol = width)
    mixed <- mgs2(accepted_basis %*% rotation, tol = tol)$Q
    if (ncol(mixed) == width) {
      mixed
    } else {
      accepted_basis[, seq_len(width), drop = FALSE]
    }
  } else if (rank > 0L) {
    accepted_basis
  } else {
    matrix(numeric(0), nrow = n, ncol = 0L)
  }

  used <- min(rank, width)
  augmented <- width - used
  if (augmented > 0L) {
    # Draw until the fitted block is complete. MGS2 both projects against the
    # accepted block and orthonormalizes the new directions among themselves.
    while (ncol(start) < width) {
      needed <- width - ncol(start)
      raw <- matrix(stats::rnorm(n * needed), nrow = n, ncol = needed)
      filler <- if (ncol(start)) {
        mgs2(raw, against = start, tol = tol)$Q
      } else {
        mgs2(raw, tol = tol)$Q
      }
      if (!ncol(filler)) {
        next
      }
      take <- min(needed, ncol(filler))
      start <- cbind(start, filler[, seq_len(take), drop = FALSE])
    }
  }

  rejected <- supplied - rank
  start_source <- if (rank > 0L) {
    "user_supplied"
  } else {
    # Supplied but numerically unusable (e.g. all-zero / rank-collapsed): the
    # block is fully random, so this is an honest degenerate warm start.
    "user_supplied_degenerate"
  }

  list(
    start = start,
    start_source = start_source,
    supplied = supplied,
    accepted = accepted,
    rejected = rejected,
    augmented = augmented,
    rank = rank,
    compressed = compressed,
    invariant_guard_used = FALSE,
    invariant_relative_residual = NA_real_,
    guard_operator_block_calls = 0L,
    guard_operator_columns = 0L,
    accepted_basis = accepted_basis
  )
}

#' Detect a supplied invariant subspace that cannot establish target identity.
#'
#' Residual certification establishes that returned pairs are eigenpairs; it
#' does not prove that an exact invariant start contains the requested extremal
#' pairs. If a fully supplied basis is already invariant at the requested
#' tolerance, discard it and use the solver's cold start rather than risk
#' certifying a non-target invariant block.
#'
#' @param op Hermitian operator.
#' @param basis Orthonormal accepted user basis.
#' @param tol Solver tolerance.
#' @return Guard decision, relative invariance residual, and exact operator work.
#' @keywords internal
warm_start_invariant_guard <- function(op, basis, tol) {
  basis <- as.matrix(basis)
  if (!ncol(basis)) {
    return(list(
      discard = FALSE,
      relative_residual = NA_real_,
      operator_block_calls = 0L,
      operator_columns = 0L
    ))
  }
  applied <- apply_operator(op, basis)
  residual <- applied - basis %*% crossprod(basis, applied)
  denominator <- max(sqrt(sum(applied^2)), sqrt(sum(basis^2)),
                     .Machine$double.eps)
  relative <- sqrt(sum(residual^2)) / denominator
  threshold <- max(as.numeric(tol), 100 * .Machine$double.eps)
  list(
    discard = is.finite(relative) && relative <= threshold,
    relative_residual = relative,
    operator_block_calls = 1L,
    operator_columns = ncol(basis)
  )
}
