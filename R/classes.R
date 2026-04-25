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
#' @param max_subspace Optional maximum Krylov subspace size.
#' @param reorthogonalize Whether to apply full reorthogonalization.
#' @return An `eigencore_method` descriptor selecting Lanczos iteration.
#' @export
lanczos <- function(max_subspace = NULL, reorthogonalize = TRUE) {
  new_method(
    "lanczos",
    max_subspace = max_subspace,
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
