#' Shifted Cholesky preconditioner.
#'
#' @param A Symmetric positive semidefinite or positive definite matrix.
#' @param shift Non-negative diagonal shift added before factorization.
#' @return A typed preconditioner function mapping residual blocks to
#'   preconditioned blocks.
#' @export
shifted_cholesky_preconditioner <- function(A, shift = 0) {
  if (!inherits(A, "Matrix")) {
    A <- Matrix::Matrix(A, sparse = FALSE)
  }
  shift <- as.numeric(shift)
  if (length(shift) != 1L || is.na(shift) || shift < 0) {
    stop("shift must be a single non-negative number.", call. = FALSE)
  }
  n <- nrow(A)
  if (ncol(A) != n) {
    stop("A must be square.", call. = FALSE)
  }
  factor <- Matrix::Cholesky(A + Matrix::Diagonal(n) * shift, LDL = FALSE)
  apply <- function(R) {
    as.matrix(Matrix::solve(factor, R))
  }
  new_eigencore_preconditioner(
    apply,
    kind = "shifted_cholesky",
    native = FALSE,
    n = n,
    shift = shift,
    factorization = "Matrix::Cholesky"
  )
}

#' Shifted diagonal preconditioner.
#'
#' @param A Diagonal matrix or numeric vector containing the diagonal entries.
#' @param shift Non-negative diagonal shift added before inversion.
#' @return A typed native-backed preconditioner function mapping residual blocks
#'   to preconditioned blocks.
#' @export
shifted_diagonal_preconditioner <- function(A, shift = 0) {
  diag <- if (is.numeric(A) && is.null(dim(A))) {
    as.numeric(A)
  } else {
    if (!inherits(A, "Matrix")) {
      A <- Matrix::Matrix(A, sparse = TRUE)
    }
    n <- nrow(A)
    if (ncol(A) != n) {
      stop("A must be square.", call. = FALSE)
    }
    if (!inherits(A, "diagonalMatrix")) {
      trip <- as.data.frame(Matrix::summary(A))
      if (nrow(trip) && any(trip$i != trip$j)) {
        stop("A must be diagonal.", call. = FALSE)
      }
    }
    as.numeric(Matrix::diag(A))
  }
  if (!length(diag) || any(!is.finite(diag))) {
    stop("A must provide at least one finite diagonal entry.", call. = FALSE)
  }
  shift <- as.numeric(shift)
  if (length(shift) != 1L || is.na(shift) || shift < 0) {
    stop("shift must be a single non-negative number.", call. = FALSE)
  }
  diag <- diag + shift
  if (any(diag <= 0)) {
    stop("shifted diagonal entries must be positive.", call. = FALSE)
  }
  n <- length(diag)
  lower <- numeric(max(n - 1L, 0L))
  upper <- numeric(max(n - 1L, 0L))
  force(lower)
  force(diag)
  force(upper)
  apply <- function(R) {
    R <- as.matrix(R)
    if (nrow(R) != n) {
      stop("residual block has incompatible row count.", call. = FALSE)
    }
    R / diag
  }
  new_eigencore_preconditioner(
    apply,
    kind = "shifted_diagonal",
    native = TRUE,
    n = n,
    shift = shift,
    factorization = "diagonal_scale",
    lower = lower,
    diag = diag,
    upper = upper
  )
}

#' Shifted tridiagonal preconditioner.
#'
#' @param A Real symmetric tridiagonal matrix.
#' @param shift Non-negative diagonal shift added to the tridiagonal system.
#'   `A + shift * I` must be nonsingular: it is factored once at
#'   construction, and a singular system is an error naming this
#'   preconditioner.
#' @return A typed preconditioner function mapping residual blocks to
#'   preconditioned blocks.
#' @export
shifted_tridiagonal_preconditioner <- function(A, shift = 0) {
  if (!inherits(A, "Matrix")) {
    A <- Matrix::Matrix(A, sparse = TRUE)
  }
  shift <- as.numeric(shift)
  if (length(shift) != 1L || is.na(shift) || shift < 0) {
    stop("shift must be a single non-negative number.", call. = FALSE)
  }
  n <- nrow(A)
  if (ncol(A) != n) {
    stop("A must be square.", call. = FALSE)
  }

  diag <- numeric(n)
  lower <- numeric(max(n - 1L, 0L))
  upper <- numeric(max(n - 1L, 0L))
  if (inherits(A, "diagonalMatrix")) {
    diag <- as.numeric(Matrix::diag(A))
  } else {
    # Symmetric/triangular storage keeps one triangle; expand to general CSC
    # (which also sums any duplicate triplets) before reading the bands.
    A <- methods::as(methods::as(A, "generalMatrix"), "CsparseMatrix")
    i_slot <- methods::slot(A, "i") + 1L
    j_slot <- rep.int(seq_len(n), diff(methods::slot(A, "p")))
    x_slot <- as.numeric(methods::slot(A, "x"))
    if (any(abs(i_slot - j_slot) > 1L)) {
      stop("A must be tridiagonal.", call. = FALSE)
    }
    on_diag <- i_slot == j_slot
    diag[j_slot[on_diag]] <- x_slot[on_diag]
    below <- i_slot == j_slot + 1L
    lower[j_slot[below]] <- x_slot[below]
    above <- i_slot == j_slot - 1L
    upper[i_slot[above]] <- x_slot[above]
  }
  if (n > 1L && !isTRUE(all.equal(lower, upper, tolerance = 1e-12))) {
    stop("A must be symmetric tridiagonal.", call. = FALSE)
  }
  diag <- diag + shift
  force(lower)
  force(diag)
  force(upper)
  # Factor once up front (O(n)): a singular A + shift * I used to surface
  # only inside LOBPCG as a bare "status=-5" (F6).
  if (n > 0L) {
    tryCatch(
      .Call("eigencore_tridiagonal_solve", as.numeric(lower), as.numeric(diag),
            as.numeric(upper), matrix(1, n, 1L), PACKAGE = "eigencore"),
      error = function(e) {
        stop(
          "shifted_tridiagonal_preconditioner(): A + shift * I is singular ",
          "to working precision (shift = ", format(shift), "; ",
          conditionMessage(e), "). Use a positive shift that moves the ",
          "system away from the spectrum of A.",
          call. = FALSE
        )
      }
    )
  }
  apply <- function(R) {
    .Call(
      "eigencore_tridiagonal_solve",
      as.numeric(lower),
      as.numeric(diag),
      as.numeric(upper),
      as.matrix(R),
      PACKAGE = "eigencore"
    )
  }
  new_eigencore_preconditioner(
    apply,
    kind = "shifted_tridiagonal",
    native = TRUE,
    n = n,
    shift = shift,
    factorization = "tridiagonal_lu",
    lower = lower,
    diag = diag,
    upper = upper
  )
}

#' @keywords internal
new_eigencore_preconditioner <- function(apply, kind, native, ...) {
  if (!is.function(apply)) {
    stop("preconditioner apply object must be a function.", call. = FALSE)
  }
  metadata <- c(
    list(
      kind = kind,
      native = isTRUE(native),
      typed = TRUE
    ),
    list(...)
  )
  attr(apply, "eigencore_preconditioner") <- metadata
  class(apply) <- unique(c("eigencore_preconditioner", class(apply)))
  apply
}

#' @keywords internal
eigencore_preconditioner_info <- function(preconditioner, include_arrays = TRUE) {
  if (is.null(preconditioner)) {
    return(list(
      supplied = FALSE,
      typed = FALSE,
      kind = "none",
      native = FALSE
    ))
  }
  metadata <- attr(preconditioner, "eigencore_preconditioner", exact = TRUE)
  if (is.null(metadata)) {
    return(list(
      supplied = TRUE,
      typed = FALSE,
      kind = "user_function",
      native = FALSE
    ))
  }
  if (!isTRUE(include_arrays)) {
    metadata[c("lower", "diag", "upper")] <- NULL
  }
  c(list(supplied = TRUE), metadata)
}

#' @keywords internal
preconditioner_plan_reason <- function(preconditioner) {
  info <- eigencore_preconditioner_info(preconditioner, include_arrays = FALSE)
  if (!isTRUE(info$supplied)) {
    return(NULL)
  }
  paste0(
    "preconditioner: ",
    info$kind,
    if (isTRUE(info$native)) " (native-backed)" else " (R-level)"
  )
}
