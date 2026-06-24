#' Create a block-native linear operator.
#'
#' @param dim Integer vector of length two giving row and column dimensions.
#' @param apply Function implementing block multiplication by the operator.
#' @param apply_adjoint Optional function implementing block multiplication by
#'   the adjoint operator.
#' @param dtype Scalar character type label, currently `"double"` or
#'   `"complex"`.
#' @param structure Eigencore structure descriptor such as [general()] or
#'   [hermitian()].
#' @param name Optional operator label used in plans and diagnostics.
#' @param metadata Optional list of implementation metadata.
#' @examples
#' A <- diag(c(3, 2, 1))
#' op <- linear_operator(
#'   dim = dim(A),
#'   apply = function(X, alpha = 1, beta = 0, Y = NULL) {
#'     Z <- alpha * (A %*% X)
#'     if (is.null(Y) || beta == 0) Z else Z + beta * Y
#'   },
#'   structure = hermitian(),
#'   metadata = list(frobenius_norm = sqrt(sum(A^2)))
#' )
#' fit <- eig_partial(op, k = 1, target = largest())
#' values(fit)
linear_operator <- function(dim, apply, apply_adjoint = NULL, dtype = "double",
                            structure = general(), name = NULL,
                            metadata = list()) {
  stopifnot(is.numeric(dim), length(dim) == 2L)
  stopifnot(is.function(apply))
  if (!is.null(apply_adjoint)) {
    stopifnot(is.function(apply_adjoint))
  }

  op <- list(
    dim = as.integer(dim),
    apply = apply,
    apply_adjoint = apply_adjoint,
    dtype = dtype,
    structure = structure,
    name = name %||% "linear_operator",
    metadata = metadata
  )
  class(op) <- "eigencore_operator"
  op
}

#' Convert an object to an eigencore operator.
#'
#' @param x Object to convert.
#' @param ... Additional arguments passed to methods.
#' @examples
#' op <- as_operator(diag(c(3, 2, 1)))
#' op$dim
#' op$structure$kind
as_operator <- function(x, ...) {
  UseMethod("as_operator")
}

#' @export
as_operator.eigencore_operator <- function(x, ...) {
  x
}

#' @export
as_operator.matrix <- function(x, ...) {
  if (is.complex(x)) {
    return(complex_dense_matrix_as_operator(x))
  }
  storage.mode(x) <- "double"
  dim_x <- dim(x)
  linear_operator(
    dim = dim_x,
    apply = function(X, alpha = 1, beta = 0, Y = NULL) {
      dense_block_apply(x, X, alpha = alpha, beta = beta, Y = Y, transpose = FALSE)
    },
    apply_adjoint = function(X, alpha = 1, beta = 0, Y = NULL) {
      dense_block_apply(x, X, alpha = alpha, beta = beta, Y = Y, transpose = TRUE)
    },
    dtype = "double",
    structure = if (is_square_symmetric(x)) hermitian() else general(),
    name = "dense_matrix",
    metadata = list(source = x, native = TRUE)
  )
}

#' @export
as_operator.default <- function(x, ...) {
  stop_if_complex_matrix_input(x)
  if (inherits(x, "ddiMatrix")) {
    return(diagonal_matrix_as_operator(x))
  }
  if (inherits(x, "dsCMatrix")) {
    return(csc_matrix_as_operator(x))
  }
  if (inherits(x, "dgCMatrix")) {
    return(csc_matrix_as_operator(x))
  }
  if (inherits(x, "Matrix")) {
    return(matrix_as_operator(x))
  }
  stop("Cannot convert object of class ", paste(class(x), collapse = "/"), " to an eigencore operator.", call. = FALSE)
}

#' @keywords internal
stop_if_complex_matrix_input <- function(x) {
  complex_input <- if (is.matrix(x)) {
    FALSE
  } else if (inherits(x, "Matrix")) {
    inherits(x, "zMatrix") ||
      inherits(x, "complexMatrix") ||
      ("x" %in% methods::slotNames(x) && is.complex(methods::slot(x, "x")))
  } else {
    FALSE
  }
  if (isTRUE(complex_input)) {
    stop(
      "Complex-valued Matrix inputs are future scope in eigencore's native sparse/operator API. ",
      "Base complex dense matrices use native dense complex LAPACK kernels; ",
      "real-valued matrices may still return complex eigenpairs through eigs().",
      call. = FALSE
    )
  }
  invisible(x)
}

#' Return the adjoint operator.
#'
#' @param x Operator-like object.
#' @param ... Additional arguments passed to methods.
adjoint <- function(x, ...) {
  UseMethod("adjoint")
}

#' @export
adjoint.eigencore_operator <- function(x, ...) {
  if (is.null(x$apply_adjoint)) {
    stop("Operator does not define apply_adjoint().", call. = FALSE)
  }
  source <- x$metadata$source
  if (!is.null(source)) {
    source <- if (identical(x$dtype, "complex")) Conj(t(source)) else t(source)
  }
  matrix <- x$metadata$matrix
  if (!is.null(matrix)) {
    matrix <- if (identical(x$dtype, "complex")) Conj(Matrix::t(matrix)) else Matrix::t(matrix)
  }
  storage <- x$metadata$storage
  if (!is.null(storage)) {
    storage <- paste0("adjoint:", storage)
  }
  linear_operator(
    dim = rev(x$dim),
    apply = x$apply_adjoint,
    apply_adjoint = x$apply,
    dtype = x$dtype,
    structure = x$structure,
    name = paste0("adjoint(", x$name, ")"),
    metadata = list(
      parent = x,
      fused = "adjoint",
      native = isTRUE(x$metadata$native),
      storage = storage,
      source = source,
      matrix = matrix
    )
  )
}

#' @keywords internal
complex_dense_matrix_as_operator <- function(x) {
  dim_x <- dim(x)
  linear_operator(
    dim = dim_x,
    apply = function(X, alpha = 1, beta = 0, Y = NULL) {
      complex_dense_block_apply(x, X, alpha = alpha, beta = beta, Y = Y, adjoint = FALSE)
    },
    apply_adjoint = function(X, alpha = 1, beta = 0, Y = NULL) {
      complex_dense_block_apply(x, X, alpha = alpha, beta = beta, Y = Y, adjoint = TRUE)
    },
    dtype = "complex",
    structure = if (is_square_symmetric(x)) hermitian() else general(),
    name = "complex_dense_matrix",
    metadata = list(
      source = x,
      native = TRUE,
      storage = "complex_dense_matrix",
      native_operator_kernel = "dense_complex_zgemm",
      native_scalar_type = "complex128"
    )
  )
}

#' @export
print.eigencore_operator <- function(x, ...) {
  cat("<eigencore operator>\n")
  cat("  name:", x$name, "\n")
  cat("  dim:", paste(x$dim, collapse = " x "), "\n")
  cat("  dtype:", x$dtype, "\n")
  cat("  structure:", x$structure$kind, "\n")
  invisible(x)
}

#' @keywords internal
apply_operator <- function(op, X, alpha = 1, beta = 0, Y = NULL) {
  op$apply(X, alpha = alpha, beta = beta, Y = Y)
}

#' @keywords internal
apply_adjoint_operator <- function(op, X, alpha = 1, beta = 0, Y = NULL) {
  if (is.null(op$apply_adjoint)) {
    stop("Operator does not define apply_adjoint().", call. = FALSE)
  }
  op$apply_adjoint(X, alpha = alpha, beta = beta, Y = Y)
}

#' @keywords internal
dense_block_apply <- function(A, X, alpha = 1, beta = 0, Y = NULL, transpose = FALSE) {
  X <- as.matrix(X)
  if (is.null(Y)) {
    out_nrow <- if (transpose) ncol(A) else nrow(A)
    Y <- matrix(0, out_nrow, ncol(X))
  } else {
    Y <- as.matrix(Y)
  }
  .Call(
    "eigencore_dense_block_apply",
    A,
    X,
    as.numeric(alpha),
    as.numeric(beta),
    Y,
    isTRUE(transpose),
    PACKAGE = "eigencore"
  )
}

#' @keywords internal
complex_dense_block_apply <- function(A, X, alpha = 1, beta = 0, Y = NULL, adjoint = FALSE) {
  X <- as.matrix(X)
  if (!is.complex(X)) {
    X <- X + 0i
  }
  if (is.null(Y)) {
    out_nrow <- if (adjoint) ncol(A) else nrow(A)
    Y <- matrix(0 + 0i, out_nrow, ncol(X))
  } else {
    Y <- as.matrix(Y)
    if (!is.complex(Y)) {
      Y <- Y + 0i
    }
  }
  .Call(
    "eigencore_dense_complex_block_apply",
    A,
    X,
    as.complex(alpha),
    as.complex(beta),
    Y,
    isTRUE(adjoint),
    PACKAGE = "eigencore"
  )
}

#' @keywords internal
csc_block_apply <- function(A, X, alpha = 1, beta = 0, Y = NULL, transpose = FALSE) {
  X <- as.matrix(X)
  if (is.null(Y)) {
    out_nrow <- if (transpose) ncol(A) else nrow(A)
    Y <- matrix(0, out_nrow, ncol(X))
  } else {
    Y <- as.matrix(Y)
  }
  .Call(
    "eigencore_csc_block_apply",
    methods::slot(A, "i"),
    methods::slot(A, "p"),
    methods::slot(A, "x"),
    methods::slot(A, "Dim"),
    X,
    as.numeric(alpha),
    as.numeric(beta),
    Y,
    isTRUE(transpose),
    PACKAGE = "eigencore"
  )
}

#' @keywords internal
csc_native_matrix <- function(x) {
  if (inherits(x, "dsCMatrix") || inherits(x, "symmetricMatrix")) {
    return(methods::as(methods::as(x, "generalMatrix"), "CsparseMatrix"))
  }
  x
}

#' @keywords internal
csc_matrix_as_operator <- function(x) {
  dim_x <- dim(x)
  input_storage <- class(x)[[1L]]
  symmetric_storage <- inherits(x, "symmetricMatrix")
  frobenius_norm <- if (inherits(x, "dgCMatrix") && !inherits(x, "symmetricMatrix")) {
    sqrt(sum(methods::slot(x, "x")^2))
  } else {
    as.numeric(Matrix::norm(x, type = "F"))
  }
  x <- csc_native_matrix(x)
  linear_operator(
    dim = dim_x,
    apply = function(X, alpha = 1, beta = 0, Y = NULL) {
      csc_block_apply(x, X, alpha = alpha, beta = beta, Y = Y, transpose = FALSE)
    },
    apply_adjoint = function(X, alpha = 1, beta = 0, Y = NULL) {
      csc_block_apply(x, X, alpha = alpha, beta = beta, Y = Y, transpose = TRUE)
    },
    dtype = "double",
    structure = if (is_square_symmetric(x)) hermitian() else general(),
    name = "sparse_csc_matrix",
    metadata = list(
      matrix = x,
      native = TRUE,
      storage = "dgCMatrix",
      input_storage = input_storage,
      symmetric_storage = symmetric_storage,
      frobenius_norm = frobenius_norm
    )
  )
}

#' @keywords internal
diagonal_block_apply <- function(A, X, alpha = 1, beta = 0, Y = NULL) {
  X <- as.matrix(X)
  if (is.null(Y)) {
    Y <- matrix(0, nrow(A), ncol(X))
  } else {
    Y <- as.matrix(Y)
  }
  .Call(
    "eigencore_diagonal_block_apply",
    methods::slot(A, "x"),
    methods::slot(A, "Dim"),
    identical(methods::slot(A, "diag"), "U"),
    X,
    as.numeric(alpha),
    as.numeric(beta),
    Y,
    PACKAGE = "eigencore"
  )
}

#' @keywords internal
diagonal_matrix_as_operator <- function(x) {
  dim_x <- dim(x)
  unit <- identical(methods::slot(x, "diag"), "U")
  values <- if (unit) rep(1, dim_x[1L]) else methods::slot(x, "x")
  linear_operator(
    dim = dim_x,
    apply = function(X, alpha = 1, beta = 0, Y = NULL) {
      diagonal_block_apply(x, X, alpha = alpha, beta = beta, Y = Y)
    },
    apply_adjoint = function(X, alpha = 1, beta = 0, Y = NULL) {
      diagonal_block_apply(x, X, alpha = alpha, beta = beta, Y = Y)
    },
    dtype = "double",
    structure = hermitian(),
    name = "diagonal_matrix",
    metadata = list(
      matrix = x,
      native = TRUE,
      storage = "ddiMatrix",
      frobenius_norm = sqrt(sum(values^2)),
      two_norm_upper = max(abs(values))
    )
  )
}

#' @keywords internal
native_apply_noalloc_check <- function(kind, A, X) {
  Y <- matrix(0, nrow(A), ncol(X))
  .Call(
    "eigencore_native_apply_noalloc_check",
    kind,
    A,
    as.matrix(X),
    Y,
    PACKAGE = "eigencore"
  )
}

#' @keywords internal
matrix_as_operator <- function(x) {
  dim_x <- dim(x)
  linear_operator(
    dim = dim_x,
    apply = function(X, alpha = 1, beta = 0, Y = NULL) {
      out <- alpha * as.matrix(x %*% X)
      if (!is.null(Y) && beta != 0) {
        out <- out + beta * Y
      }
      out
    },
    apply_adjoint = function(X, alpha = 1, beta = 0, Y = NULL) {
      out <- alpha * as.matrix(Matrix::t(x) %*% X)
      if (!is.null(Y) && beta != 0) {
        out <- out + beta * Y
      }
      out
    },
    dtype = "double",
    structure = if (is_square_symmetric(x)) hermitian() else general(),
    name = "matrix_sparse",
    metadata = list(matrix = x, native = FALSE)
  )
}

#' @keywords internal
operator_source_matrix <- function(A) {
  if (inherits(A, "eigencore_operator")) {
    src <- A$metadata$source
    if (is.null(src) && !is.null(A$metadata$matrix)) {
      src <- A$metadata$matrix
    }
    if (is.null(src)) {
      stop("This prototype solver needs an explicit matrix-backed operator.", call. = FALSE)
    }
    return(as.matrix(src))
  }
  as.matrix(A)
}

#' @keywords internal
is_square_symmetric <- function(x, tol = sqrt(.Machine$double.eps)) {
  d <- dim(x)
  if (length(d) != 2L || d[1L] != d[2L]) {
    return(FALSE)
  }
  if (inherits(x, "Matrix")) {
    return(isTRUE(Matrix::isSymmetric(x, tol = tol)))
  }
  if (is.matrix(x) && is.double(x)) {
    return(isTRUE(.Call("eigencore_dense_is_symmetric", x, as.numeric(tol), PACKAGE = "eigencore")))
  }
  isTRUE(isSymmetric.matrix(x, tol = tol))
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
