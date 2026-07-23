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

#' Relative size of the random escape blended into a fully user-supplied start.
#'
#' A start block composed entirely of user directions could, by accident, be an
#' exact non-target invariant subspace of the operator. Krylov iteration cannot
#' leave an invariant subspace, so at a minimal `max_subspace` it would converge
#' to and certify the (wrong) non-target eigenpairs it was handed. Blending a
#' small random component (this fraction of each unit column) guarantees the
#' start has generic overlap with every eigendirection while keeping the
#' supplied directions dominant, exactly as a cold random start would.
#'
#' The blend is sized so a trapped non-target subspace cannot satisfy the
#' certificate: its perturbed backward error is on the order of the blend, which
#' must exceed the solver tolerance `tol`. `sqrt(tol)` achieves that for any
#' `tol < 1`, floored at `1e-4` (which clears the empirical trap threshold with
#' margin) so extremely tight tolerances do not shrink the escape below it.
#' @param tol Solver convergence/certification tolerance.
#' @return A positive scalar blend fraction.
#' @keywords internal
warm_start_escape_blend <- function(tol = 1e-8) {
  max(sqrt(tol), 1e-4)
}

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
    escape_blended = FALSE
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
#' @param solver_tol Solver convergence/certification tolerance; sizes the
#'   random escape blended into a fully user-supplied start.
#' @param tol Numerical-rank tolerance for boundary orthonormalization.
#' @return A list with `start` (the `n x width` block) and provenance fields
#'   `start_source`, `supplied`, `accepted`, `rejected`, `augmented`, `rank`,
#'   `compressed`, `escape_blended`.
#' @keywords internal
prepare_initial_subspace <- function(initial_subspace, n, width,
                                     solver_tol = 1e-8,
                                     tol = sqrt(.Machine$double.eps)) {
  width <- as.integer(width)
  if (length(width) != 1L || is.na(width) || width < 1L) {
    stop("Internal error: warm-start block width must be a positive integer.",
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

  # Boundary orthonormalization; MGS2 drops linearly dependent / negligible
  # columns, so ncol(Q) is the accepted numerical rank.
  ortho <- mgs2(x, tol = tol)
  accepted_basis <- ortho$Q
  rank <- ncol(accepted_basis)

  accepted <- min(rank, width)
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
  } else if (accepted > 0L) {
    accepted_basis[, seq_len(accepted), drop = FALSE]
  } else {
    matrix(numeric(0), nrow = n, ncol = 0L)
  }

  augmented <- width - accepted
  if (augmented > 0L) {
    raw <- matrix(stats::rnorm(n * augmented), nrow = n, ncol = augmented)
    filler <- if (accepted > 0L) {
      reorthogonalize_against(raw, start, passes = 2L)
    } else {
      raw
    }
    start <- cbind(start, filler)
  }

  # When the start block is entirely user-supplied (no random augmentation),
  # blend a small deterministic random escape so it can never be an exact
  # non-target invariant subspace that traps Krylov iteration into certifying
  # non-target eigenpairs. The augmented case already carries random content.
  escape_blended <- FALSE
  if (augmented == 0L && accepted >= 1L) {
    escape <- matrix(stats::rnorm(n * width), nrow = n, ncol = width)
    cn <- sqrt(colSums(escape^2))
    cn[cn == 0] <- 1
    escape <- sweep(escape, 2L, cn, `/`)
    start <- start + warm_start_escape_blend(solver_tol) * escape
    escape_blended <- TRUE
  }

  rejected <- supplied - accepted
  start_source <- if (accepted > 0L) {
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
    escape_blended = escape_blended
  )
}
