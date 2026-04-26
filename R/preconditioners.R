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

#' Shifted tridiagonal preconditioner.
#'
#' @param A Real symmetric tridiagonal matrix.
#' @param shift Non-negative diagonal shift added to the tridiagonal system.
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
  } else if (inherits(A, "CsparseMatrix")) {
    i_slot <- methods::slot(A, "i") + 1L
    p_slot <- methods::slot(A, "p")
    x_slot <- methods::slot(A, "x")
    for (j in seq_len(n)) {
      start <- p_slot[[j]] + 1L
      end <- p_slot[[j + 1L]]
      if (start > end) {
        next
      }
      for (pos in start:end) {
        row <- i_slot[[pos]]
        value <- x_slot[[pos]]
        if (abs(row - j) > 1L) {
          stop("A must be tridiagonal.", call. = FALSE)
        }
        if (row == j) {
          diag[j] <- diag[j] + value
        } else if (row == j + 1L) {
          lower[j] <- lower[j] + value
        } else if (row == j - 1L) {
          upper[j - 1L] <- upper[j - 1L] + value
        }
      }
    }
  } else {
    trip <- as.data.frame(Matrix::summary(A))
    if (!all(c("i", "j", "x") %in% names(trip))) {
      stop("A must be tridiagonal.", call. = FALSE)
    }
    if (any(abs(trip$i - trip$j) > 1L)) {
      stop("A must be tridiagonal.", call. = FALSE)
    }
    for (idx in seq_len(nrow(trip))) {
      i <- trip$i[idx]
      j <- trip$j[idx]
      x <- trip$x[idx]
      if (i == j) {
        diag[i] <- x
      } else if (i == j + 1L) {
        lower[j] <- x
      } else if (j == i + 1L) {
        upper[i] <- x
      }
    }
  }
  if (n > 1L && !isTRUE(all.equal(lower, upper, tolerance = 1e-12))) {
    stop("A must be symmetric tridiagonal.", call. = FALSE)
  }
  diag <- diag + shift
  force(lower)
  force(diag)
  force(upper)
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
    factorization = "tridiagonal_thomas",
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
