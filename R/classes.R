#' @noRd
new_target <- function(kind, value = NULL) {
  structure(list(kind = kind, value = value), class = "eigencore_target")
}

#' Target the largest algebraic values.
#'
#' @return An `eigencore_target` descriptor selecting the largest algebraic
#'   eigenvalues or singular values.
#' @export
largest <- function() {
  new_target("largest")
}

#' Target the smallest algebraic values.
#'
#' @return An `eigencore_target` descriptor selecting the smallest algebraic
#'   eigenvalues or singular values.
#' @export
smallest <- function() {
  new_target("smallest")
}

#' Target the largest values by magnitude.
#'
#' @return An `eigencore_target` descriptor selecting values with the largest
#'   modulus (ARPACK `LM`).
#' @export
largest_magnitude <- function() {
  new_target("largest_magnitude")
}

#' Target the smallest values by magnitude.
#'
#' @return An `eigencore_target` descriptor selecting values with the smallest
#'   modulus (ARPACK `SM`).
#' @export
smallest_magnitude <- function() {
  new_target("smallest_magnitude")
}

#' Target values nearest a shift.
#'
#' @param sigma Numeric shift used to rank values by distance `|x - sigma|`.
#' @return An `eigencore_target` descriptor for nearest-to-`sigma` selection.
#' @export
nearest <- function(sigma) {
  new_target("nearest", sigma)
}

#' Target the largest real part.
#'
#' @return An `eigencore_target` descriptor selecting values by largest
#'   `Re(x)` (ARPACK `LR`).
#' @export
largest_real <- function() {
  new_target("largest_real")
}

#' Target the smallest real part.
#'
#' @return An `eigencore_target` descriptor selecting values by smallest
#'   `Re(x)` (ARPACK `SR`).
#' @export
smallest_real <- function() {
  new_target("smallest_real")
}

#' Target the largest imaginary part.
#'
#' @return An `eigencore_target` descriptor selecting values by largest
#'   `Im(x)` (ARPACK `LI`).
#' @export
largest_imaginary <- function() {
  new_target("largest_imaginary")
}

#' Target the smallest imaginary part.
#'
#' @return An `eigencore_target` descriptor selecting values by smallest
#'   `Im(x)` (ARPACK `SI`).
#' @export
smallest_imaginary <- function() {
  new_target("smallest_imaginary")
}

#' Target both algebraic ends.
#'
#' @param k_low Number of values to select from the smallest algebraic end.
#' @param k_high Number of values to select from the largest algebraic end.
#' @return An `eigencore_target` descriptor selecting both algebraic ends
#'   (ARPACK `BE`).
#' @export
both_ends <- function(k_low, k_high) {
  k_low <- as.integer(k_low)
  k_high <- as.integer(k_high)
  if (length(k_low) != 1L || is.na(k_low) || k_low < 0L) {
    stop("k_low must be a single non-negative integer.", call. = FALSE)
  }
  if (length(k_high) != 1L || is.na(k_high) || k_high < 0L) {
    stop("k_high must be a single non-negative integer.", call. = FALSE)
  }
  if ((k_low + k_high) < 1L) {
    stop("k_low + k_high must be at least 1.", call. = FALSE)
  }
  new_target("both_ends", list(k_low = k_low, k_high = k_high))
}

#' Target every eigenvalue in an interval.
#'
#' `interval(a, b)` selects all eigenvalues of a Hermitian problem (or of a
#' symmetric-definite pencil `A x = lambda B x`) that lie in the closed
#' interval `[a, b]`. The number of such eigenvalues is not supplied: the
#' solver counts them first with Sylvester's law of inertia (see
#' [eigen_count()]), so `eig_partial(A, target = interval(a, b))` needs no
#' `k`. A supplied `k` is an upper bound and must be at least the count.
#'
#' **End points.** The interval is closed. Counting factors `A - t B` at
#' `t = a` and `t = b`; when a factorisation there is not reliable (an
#' eigenvalue numerically at the end point) the end point is moved outward
#' by a small recorded perturbation, so an eigenvalue inside that numerical
#' zero band is counted as inside. Eigenvalues whose residual bound
#' straddles an end point are reported in
#' `certificate(fit)$completeness$boundary_ambiguous` rather than dropped.
#'
#' **Routes.** Dense matrices use LAPACK `dsyevr` with `RANGE = "V"`
#' (exact: the tridiagonal Sturm bisection selects the eigenvalues);
#' generalized dense pencils reduce with the Cholesky factor of `B` first.
#' Sparse matrices are counted and solved with \eqn{LDL^T} shift-invert
#' Lanczos at the interval centre; wide intervals are split into slices by
#' inertia counts (spectrum slicing), each slice is solved at its own centre
#' (reusing the symbolic factorisation) and the slices are merged. Every
#' result is certified by residuals and by the inertia count:
#' `certificate(fit)$target_completeness` is `"inertia_verified"` when the
#' returned values provably are all eigenvalues in the interval (`"exact"`
#' for the dense LAPACK route).
#'
#' @param a,b Interval end points, `a < b`. One of them may be infinite
#'   (`interval(-Inf, b)` selects every eigenvalue up to `b`); both infinite
#'   is the full spectrum, use [eig_full()].
#' @return An `eigencore_target` descriptor of kind `"interval"`.
#' @examples
#' A <- diag(c(1, 2, 2, 5, 7, 9))
#' fit <- eig_partial(A, target = interval(1.5, 6))
#' values(fit)
#' certificate(fit)$target_completeness
#' @export
interval <- function(a, b) {
  check <- function(x, name) {
    x <- suppressWarnings(as.numeric(x))
    if (length(x) != 1L || is.na(x)) {
      stop(name, " must be a single number.", call. = FALSE)
    }
    x
  }
  a <- check(a, "a")
  b <- check(b, "b")
  if (!(a < b)) {
    stop("interval(a, b) requires a < b, got a = ", format(a), ", b = ",
         format(b), ".", call. = FALSE)
  }
  if (is.infinite(a) && is.infinite(b)) {
    stop("interval(-Inf, Inf) is the full spectrum; use eig_full().",
         call. = FALSE)
  }
  new_target("interval", list(lower = a, upper = b))
}

#' @noRd
is_interval_target <- function(target) {
  inherits(target, "eigencore_target") && identical(target$kind, "interval")
}

#' @noRd
new_method <- function(kind, ...) {
  structure(c(list(kind = kind), list(...)), class = "eigencore_method")
}

#' @noRd
validate_max_subspace <- function(max_subspace) {
  if (is.null(max_subspace)) {
    return(NULL)
  }
  if (!is.numeric(max_subspace) || length(max_subspace) != 1L ||
      is.na(max_subspace) || max_subspace != round(max_subspace) ||
      max_subspace < 2) {
    stop("max_subspace must be NULL or a single whole number >= 2.",
         call. = FALSE)
  }
  as.integer(max_subspace)
}

#' Automatic solver choice.
#'
#' @param max_subspace Optional maximum Krylov subspace size (the ARPACK
#'   `ncv`). When the planner selects a restarted Krylov route (thick-restart
#'   or block Lanczos, Krylov-Schur Arnoldi, shift-invert Lanczos or Arnoldi,
#'   Golub-Kahan, or implicit-Gram Lanczos for SVD) this caps the active
#'   basis; it is capped at the problem dimension and must leave room for the
#'   wanted pairs (at least `k + 1`, `k + block` for block Lanczos). Routes
#'   without a Krylov basis (dense LAPACK, explicit Gram SVD, LOBPCG) ignore
#'   it. `NULL` (default) uses the route's own default. The solve's `maxit`
#'   argument is the iteration (restart) limit and never changes the
#'   subspace size.
#' @return An `eigencore_method` descriptor that lets the planner choose a
#'   solver based on problem structure.
#' @export
auto <- function(max_subspace = NULL) {
  new_method("auto", max_subspace = validate_max_subspace(max_subspace))
}

#' Hermitian Lanczos method descriptor.
#'
#' @param max_subspace Optional maximum active Krylov subspace size `m` (the
#'   ARPACK `ncv`). Must be at least `k + 1` (`k + block` for block Lanczos).
#'   The native thick-restart path keeps the active basis bounded by this
#'   value across restart cycles; unrestarted reference paths build at most
#'   this many Lanczos vectors. This is the only subspace-size control: the
#'   solve's `maxit` argument is an iteration (restart) limit.
#' @param max_restarts Optional non-negative integer giving the maximum
#'   number of thick-restart cycles allowed before stopping with whatever
#'   has converged. Default `100L`. Equivalent to the solve's `maxit`
#'   argument on thick-restart routes; supplying both with different values
#'   is an error.
#' @param block Native block size. `1L` selects the scalar path; for a
#'   matrix-free operator this remains the reference Hermitian Lanczos
#'   boundary. Values greater than one select the native block Krylov path
#'   where supported, including real Hermitian matrix-free callbacks.
#' @param check_stride Native block thick-restart mid-sweep convergence stride.
#'   `0L` (default) evaluates convergence once per full sweep (legacy). A
#'   positive `N` evaluates convergence every `N` block iterations within a
#'   sweep, letting a warm start that converges after a few blocks stop early
#'   instead of paying a full cold-sized sweep. Mid-sweep checks never consume
#'   extra operator applications and never change results at
#'   `check_stride = 0L`. This control applies only to native block paths.
#' @param reorthogonalize Whether to apply full reorthogonalization. The
#'   native path always reorthogonalizes (DGKS x2) and ignores this flag;
#'   it is preserved for the R reference solver's public API.
#' @param completeness Target-completeness check run after a certified
#'   Hermitian solve with a `largest()`, `smallest()`,
#'   `largest_magnitude()` (or, for the inertia check, `nearest()`) target.
#'   `"inertia"` proves completeness deterministically by counting
#'   eigenvalues with an \eqn{LDL^T} factorisation of `A - t B` (see
#'   [eigen_count()]); it needs an explicit dense or sparse matrix and
#'   reports `"inertia_verified"`, `"inertia_failed"` (repaired when
#'   possible, otherwise `passed = FALSE`) or `"inertia_inconclusive"` (an
#'   eigenvalue cluster straddles the target edge within the residual
#'   bound). `"probe"` runs a short block Lanczos process on the operator
#'   deflated against the returned eigenvectors (from a fixed-seed start; the
#'   global random stream is not touched) and, if it finds a more-preferred
#'   eigenvalue outside the returned set (for example a missed copy of a
#'   repeated eigenvalue), repairs the result with a deflated complement
#'   solve; it is probabilistic (it can prove a set incomplete but not
#'   complete) and is the check for matrix-free operators. `"auto"` uses the
#'   inertia check for edge targets (not `nearest()`) when the matrix is
#'   explicit and its predicted
#'   factorisation time is at most
#'   `max(getOption("eigencore.completeness_inertia_seconds", 0.5),
#'   getOption("eigencore.completeness_inertia_ratio", 1) * solve time)`,
#'   and the probe otherwise. `"none"` skips the check. `NULL` (default) uses
#'   `getOption("eigencore.target_completeness", "auto")`. The outcome is
#'   recorded in `certificate(fit)$target_completeness`; see the
#'   "Certificates" vignette.
#' @return An `eigencore_method` descriptor selecting Lanczos iteration.
#' @export
lanczos <- function(max_subspace = NULL, max_restarts = NULL, block = 1L,
                    check_stride = 0L, reorthogonalize = TRUE,
                    completeness = NULL) {
  completeness <- validate_completeness_mode(completeness)
  block <- as.integer(block)
  if (length(block) != 1L || is.na(block) || block < 1L) {
    stop("block must be a single positive integer.", call. = FALSE)
  }
  check_stride <- as.integer(check_stride)
  if (length(check_stride) != 1L || is.na(check_stride) || check_stride < 0L) {
    stop("check_stride must be a single non-negative integer.", call. = FALSE)
  }
  if (!is.null(max_restarts)) {
    max_restarts <- as.integer(max_restarts)
    if (length(max_restarts) != 1L || is.na(max_restarts) || max_restarts < 0L) {
      stop("max_restarts must be a single non-negative integer.", call. = FALSE)
    }
  }
  method <- new_method(
    "lanczos",
    max_subspace = validate_max_subspace(max_subspace),
    max_restarts = max_restarts,
    block = block,
    check_stride = check_stride,
    reorthogonalize = reorthogonalize
  )
  if (!is.null(completeness)) {
    method$completeness <- completeness
  }
  method
}

#' Golub-Kahan bidiagonalization method descriptor.
#'
#' @param max_subspace Optional maximum Krylov subspace size (the ARPACK
#'   `ncv`); fixes the subspace instead of the default adaptive growth.
#' @param reorthogonalize Whether to apply full two-sided
#'   reorthogonalization. `FALSE` selects the native one-sided small-side
#'   policy where supported, with final acceptance still controlled by the
#'   exact two-sided certificate.
#' @return An `eigencore_method` descriptor selecting Golub-Kahan
#'   bidiagonalization.
#' @export
golub_kahan <- function(max_subspace = NULL, reorthogonalize = TRUE) {
  new_method(
    "golub_kahan",
    max_subspace = validate_max_subspace(max_subspace),
    reorthogonalize = reorthogonalize
  )
}

#' Randomized SVD method descriptor.
#'
#' @param oversample Number of extra samples beyond the requested rank.
#' @param n_iter Number of subspace-iteration refinement passes.
#' @param block Optional block size.
#' @param normalizer Basis normalizer to use (`"qr"`, `"lu"`, or `"none"`).
#' @param refine Whether to refine with a certified Lanczos pass.
#' @return An `eigencore_method` descriptor selecting randomized SVD.
#' @export
randomized <- function(oversample = 10, n_iter = 2, block = NULL,
                       normalizer = c("qr", "lu", "none"), refine = TRUE) {
  normalizer <- match.arg(normalizer)
  new_method(
    "randomized",
    oversample = oversample,
    n_iter = n_iter,
    block = block,
    normalizer = normalizer,
    refine = refine
  )
}

#' LOBPCG method descriptor.
#'
#' @param maxit Maximum LOBPCG iterations. `NULL` (default) uses the solve's
#'   `maxit` argument when given, else the `eigencore.lobpcg_maxit` option
#'   (200). Supplying both this and a different solve-level `maxit` is an
#'   error.
#' @param preconditioner Optional function taking a residual block and
#'   returning a preconditioned block with the same dimensions.
#' @param constraints Optional matrix whose columns span a subspace to deflate.
#'   Iterates are kept orthogonal to this subspace in the Euclidean or
#'   generalized `B` inner product. Native constrained LOBPCG is not promoted
#'   yet; constrained problems use the labelled reference path.
#' @return An `eigencore_method` descriptor selecting LOBPCG. Built-in
#'   standard Hermitian dense/CSC operators may use a native prototype;
#'   unsupported cases route to the reference prototype.
#' @importFrom utils modifyList tail
#' @export
lobpcg <- function(maxit = NULL, preconditioner = NULL, constraints = NULL) {
  if (!is.null(maxit)) {
    maxit <- as.integer(maxit)
    if (length(maxit) != 1L || is.na(maxit) || maxit < 1L) {
      stop("maxit must be NULL or a single positive integer.", call. = FALSE)
    }
  }
  if (!is.null(preconditioner) && !is.function(preconditioner)) {
    stop("preconditioner must be NULL or a function.", call. = FALSE)
  }
  if (!is.null(constraints)) {
    constraints <- as.matrix(constraints)
    storage.mode(constraints) <- "double"
    if (nrow(constraints) < 1L || ncol(constraints) < 1L ||
        any(!is.finite(constraints))) {
      stop("constraints must be a finite numeric matrix with at least one column.", call. = FALSE)
    }
  }
  new_method(
    "lobpcg",
    maxit = maxit,
    preconditioner = preconditioner,
    constraints = constraints
  )
}

#' Shift-invert method descriptor.
#'
#' @param sigma Shift value `sigma`.
#' @param solve Optional user-supplied solve operator for `(A - sigma B)`.
#' @param factorization Optional precomputed factorization handle.
#' @param max_subspace Optional maximum Krylov subspace size for the Lanczos
#'   (Hermitian) or Krylov-Schur Arnoldi (nonsymmetric) iteration on the
#'   inverted operator. `NULL` uses the route default.
#' @return An `eigencore_method` descriptor selecting shift-invert.
#' @export
shift_invert <- function(sigma, solve = NULL, factorization = NULL,
                         max_subspace = NULL) {
  new_method("shift_invert", sigma = sigma, solve = solve,
             factorization = factorization,
             max_subspace = validate_max_subspace(max_subspace))
}

#' General operator structure descriptor.
#'
#' @return An `eigencore_structure` descriptor for general operators.
#' @export
general <- function() {
  structure(list(kind = "general"), class = "eigencore_structure")
}

#' Hermitian/symmetric operator structure descriptor.
#'
#' @return An `eigencore_structure` descriptor marking an operator as
#'   Hermitian / symmetric.
#' @export
hermitian <- function() {
  structure(list(kind = "hermitian"), class = "eigencore_structure")
}

#' Euclidean vector space descriptor.
#'
#' @param dim Dimension of the space.
#' @param dtype Scalar type (currently only `"double"`).
#' @return An `eigencore_space` descriptor for the Euclidean space
#'   `R^dim` or `C^dim`.
#' @export
euclidean <- function(dim, dtype = "double") {
  structure(list(dim = dim, dtype = dtype, metric = NULL), class = "eigencore_space")
}

#' @noRd
target_label <- function(target) {
  if (!inherits(target, "eigencore_target")) {
    return(as.character(target))
  }
  if (identical(target$kind, "nearest")) {
    return(paste0("nearest(", target$value, ")"))
  }
  if (identical(target$kind, "interval")) {
    return(paste0("interval(", format(target$value$lower, digits = 10), ", ",
                  format(target$value$upper, digits = 10), ")"))
  }
  if (identical(target$kind, "both_ends")) {
    return(paste0("both_ends(", target$value$k_low, ", ", target$value$k_high, ")"))
  }
  target$kind
}

#' @noRd
method_label <- function(method) {
  if (is.null(method)) {
    return("auto")
  }
  if (!inherits(method, "eigencore_method")) {
    return(as.character(method))
  }
  method$kind
}

#' Whitelist of method kinds that act as a problem-level spectral transform.
#' New transforms (e.g. polynomial, chebyshev) must be added here so that
#' eig_partial and solve.eigencore_eigen_problem route them through the
#' transform plumbing instead of dropping them silently.
#' @noRd
transform_method_kinds <- function() {
  c("shift_invert")
}

#' @noRd
is_transform_method <- function(method) {
  inherits(method, "eigencore_method") &&
    !is.null(method$kind) &&
    method$kind %in% transform_method_kinds()
}
