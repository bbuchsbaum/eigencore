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

#' @noRd
new_method <- function(kind, ...) {
  structure(c(list(kind = kind), list(...)), class = "eigencore_method")
}

#' Automatic solver choice.
#'
#' @return An `eigencore_method` descriptor that lets the planner choose a
#'   solver based on problem structure.
#' @export
auto <- function() {
  new_method("auto")
}

#' Hermitian Lanczos method descriptor.
#'
#' @param max_subspace Optional maximum active Krylov subspace size `m`. Must
#'   be at least `k + 1`. The native thick-restart path keeps the active
#'   basis bounded by this value across restart cycles.
#' @param max_restarts Optional non-negative integer giving the maximum
#'   number of thick-restart cycles allowed before stopping with whatever
#'   has converged. Default `100L`.
#' @param block Native block size. `1L` selects the scalar thick-restart path;
#'   values greater than one select the native block Krylov prototype where
#'   supported.
#' @param reorthogonalize Whether to apply full reorthogonalization. The
#'   native path always reorthogonalizes (DGKS x2) and ignores this flag;
#'   it is preserved for the R reference solver's public API.
#' @return An `eigencore_method` descriptor selecting Lanczos iteration.
#' @export
lanczos <- function(max_subspace = NULL, max_restarts = NULL, block = 1L,
                    reorthogonalize = TRUE) {
  block <- as.integer(block)
  if (length(block) != 1L || is.na(block) || block < 1L) {
    stop("block must be a single positive integer.", call. = FALSE)
  }
  if (!is.null(max_restarts)) {
    max_restarts <- as.integer(max_restarts)
    if (length(max_restarts) != 1L || is.na(max_restarts) || max_restarts < 0L) {
      stop("max_restarts must be a single non-negative integer.", call. = FALSE)
    }
  }
  new_method(
    "lanczos",
    max_subspace = max_subspace,
    max_restarts = max_restarts,
    block = block,
    reorthogonalize = reorthogonalize
  )
}

#' Golub-Kahan bidiagonalization method descriptor.
#'
#' @param max_subspace Optional maximum Krylov subspace size.
#' @param reorthogonalize Whether to apply full reorthogonalization.
#' @return An `eigencore_method` descriptor selecting Golub-Kahan
#'   bidiagonalization.
#' @export
golub_kahan <- function(max_subspace = NULL, reorthogonalize = TRUE) {
  new_method(
    "golub_kahan",
    max_subspace = max_subspace,
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
#' @param maxit Maximum LOBPCG iterations.
#' @param preconditioner Optional function taking a residual block and
#'   returning a preconditioned block with the same dimensions.
#' @return An `eigencore_method` descriptor selecting LOBPCG. Built-in
#'   standard Hermitian dense/CSC operators may use a native prototype;
#'   unsupported cases route to the reference prototype.
#' @export
lobpcg <- function(maxit = 200L, preconditioner = NULL) {
  maxit <- as.integer(maxit)
  if (length(maxit) != 1L || is.na(maxit) || maxit < 1L) {
    stop("maxit must be a single positive integer.", call. = FALSE)
  }
  if (!is.null(preconditioner) && !is.function(preconditioner)) {
    stop("preconditioner must be NULL or a function.", call. = FALSE)
  }
  new_method("lobpcg", maxit = maxit, preconditioner = preconditioner)
}

#' Shift-invert method descriptor.
#'
#' @param sigma Shift value `sigma`.
#' @param solve Optional user-supplied solve operator for `(A - sigma B)`.
#' @param factorization Optional precomputed factorization handle.
#' @return An `eigencore_method` descriptor selecting shift-invert.
#' @export
shift_invert <- function(sigma, solve = NULL, factorization = NULL) {
  new_method("shift_invert", sigma = sigma, solve = solve, factorization = factorization)
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
