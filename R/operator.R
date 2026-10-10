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
#' @param metadata Optional list of implementation metadata. Certificates
#'   read two optional entries, trusted as given: `two_norm`, the exact
#'   spectral norm `||A||_2` (makes the certificate scale exact), and
#'   `frobenius_norm`, whose `||A||_F / sqrt(min(dim))` gives a structural
#'   lower bound on `||A||_2`. Without them certificates use lower bounds from
#'   the computed vectors.
#' @param operator_id Optional non-empty character identity for a logical
#'   callback-operator lineage. Supply together with `revision`.
#' @param revision Optional non-empty character revision for the values and
#'   behavior implemented by the callback. Supply together with `operator_id`.
#' @param portable Whether an explicitly identified callback operator may be
#'   restored and validated in another R session. Built-in matrix-backed
#'   operators derive portable deterministic identity automatically.
#' @return An `eigencore_operator` list containing dimensions, apply callbacks,
#'   scalar type, structure metadata, a display name, and implementation
#'   metadata.
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
                            metadata = list(), operator_id = NULL,
                            revision = NULL, portable = FALSE) {
  stopifnot(is.numeric(dim), length(dim) == 2L)
  stopifnot(is.function(apply))
  if (!is.null(apply_adjoint)) {
    stopifnot(is.function(apply_adjoint))
  }
  if (!is.list(metadata)) {
    stop("metadata must be a list.", call. = FALSE)
  }

  identity <- make_operator_identity(
    dim = dim,
    dtype = dtype,
    structure = structure,
    metadata = metadata,
    operator_id = operator_id,
    revision = revision,
    portable = portable,
    defer_builtin = TRUE
  )
  # This evaluation frame, already captured by the apply closures, doubles as
  # the per-operator cache (lazy built-in identity, norm estimates); see
  # operator_identity_cache_resolve() and operator_cache_matches().
  operator_cache <- environment()
  memo_token <- next_operator_memo_token()
  # The work accounting below needs the identity only to recognise metric
  # applies (when the plan has a B), so a deferred built-in identity is
  # resolved lazily there.
  work_identity <- function() {
    if (isTRUE(attr(identity, "deferred", exact = TRUE))) {
      operator_identity_cache_resolve(operator_cache)
    } else {
      identity
    }
  }

  adjoint_view <- identical(metadata$fused %||% NULL, "adjoint")
  raw_apply <- apply
  raw_apply_adjoint <- apply_adjoint
  wrapped_apply <- function(X, alpha = 1, beta = 0, Y = NULL) {
    token <- work_operator_enter(
      work_identity,
      kind = if (adjoint_view) "adjoint" else "operator",
      X = X
    )
    on.exit(work_operator_exit(token), add = TRUE)
    raw_apply(X, alpha = alpha, beta = beta, Y = Y)
  }
  wrapped_apply_adjoint <- if (is.null(raw_apply_adjoint)) {
    NULL
  } else {
    function(X, alpha = 1, beta = 0, Y = NULL) {
      token <- work_operator_enter(
        work_identity,
        kind = if (adjoint_view) "operator" else "adjoint",
        X = X
      )
      on.exit(work_operator_exit(token), add = TRUE)
      raw_apply_adjoint(X, alpha = alpha, beta = beta, Y = Y)
    }
  }

  op <- list(
    dim = as.integer(dim),
    apply = wrapped_apply,
    apply_adjoint = wrapped_apply_adjoint,
    dtype = dtype,
    structure = structure,
    name = name %||% "linear_operator",
    metadata = metadata,
    identity = identity,
    cache = operator_cache
  )
  class(op) <- "eigencore_operator"
  op
}

#' @keywords internal
next_operator_memo_token <- local({
  counter <- 0
  function() {
    counter <<- counter + 1
    paste0(eigencore_session_id(), ":", format(counter, scientific = FALSE))
  }
})

#' Convert an object to an eigencore operator.
#'
#' @param x Object to convert.
#' @param ... Additional arguments passed to methods.
#' @return An `eigencore_operator` representation of `x`.
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
    stop_if_nonfinite_input(x)
    return(complex_dense_matrix_as_operator(x))
  }
  # Only numeric input was ever screened for non-finite entries (logical /
  # character matrices are coerced below); keep that contract.
  screen_finite <- is.numeric(x)
  # Assigning storage.mode unconditionally forces a full copy of an
  # already-double matrix (copy-on-modify), which dominates operator
  # construction for large dense inputs.
  if (storage.mode(x) != "double") {
    storage.mode(x) <- "double"
  }
  # One native pass decides finiteness and (relative-tolerance) symmetry;
  # an R-level all(is.finite(x)) plus a separate symmetry scan used to
  # dominate plan-time cost for large dense input.
  flags <- dense_finite_symmetric(x)
  if (screen_finite && !flags[[1L]]) {
    stop("Matrix input contains NA, NaN, or Inf entries.", call. = FALSE)
  }
  symmetric <- flags[[2L]]
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
    structure = if (symmetric) hermitian() else general(),
    name = "dense_matrix",
    metadata = list(source = x, native = TRUE)
  )
}

#' @export
as_operator.default <- function(x, ...) {
  stop_if_complex_matrix_input(x)
  if (inherits(x, "denseMatrix")) {
    # Dense Matrix classes become a base double matrix whose finiteness and
    # symmetry as_operator.matrix() screens in one native pass; screening
    # the x slot here as well would add another full pass.
    converted <- native_matrix_storage(x)
    if (is.matrix(converted)) {
      return(as_operator.matrix(converted))
    }
  }
  stop_if_nonfinite_input(x)
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
    # Convert other Matrix classes once to a storage with a native kernel
    # instead of dispatching R-level %*% / Matrix::t() on every apply.
    converted <- native_matrix_storage(x)
    if (!is.null(converted)) {
      if (is.matrix(converted)) {
        return(as_operator.matrix(converted))
      }
      return(csc_matrix_as_operator(converted, input_storage = class(x)[[1L]]))
    }
    return(matrix_as_operator(x))
  }
  stop("Cannot convert object of class ", paste(class(x), collapse = "/"), " to an eigencore operator.", call. = FALSE)
}

#' @keywords internal
#' Convert a Matrix-package object once to a storage that has a native block
#' apply: a base double matrix for dense classes (dge/dsy/dtr/dpo/...), or a
#' dgCMatrix (dsCMatrix when stored symmetric) for sparse classes
#' (dgR/dgT/dtC/dsT/...). Returns NULL when no conversion applies.
native_matrix_storage <- function(x) {
  tryCatch({
    if (inherits(x, "denseMatrix")) {
      out <- base::as.matrix(x)
      if (!is.matrix(out) || is.complex(out)) {
        return(NULL)
      }
      if (storage.mode(out) != "double") {
        storage.mode(out) <- "double"
      }
      return(out)
    }
    if (inherits(x, "sparseMatrix")) {
      out <- methods::as(x, "CsparseMatrix")
      if (!methods::is(out, "dMatrix")) {
        out <- methods::as(out, "dMatrix")
      }
      if (!inherits(out, "dsCMatrix")) {
        out <- methods::as(out, "generalMatrix")
      }
      if (inherits(out, "dgCMatrix") || inherits(out, "dsCMatrix")) {
        return(out)
      }
    }
    NULL
  }, error = function(e) NULL)
}

#' @keywords internal
stop_if_nonfinite_input <- function(x) {
  values <- if (inherits(x, "Matrix")) {
    if (methods::.hasSlot(x, "x")) methods::slot(x, "x") else NULL
  } else if (is.numeric(x) || is.complex(x)) {
    x
  } else {
    NULL
  }
  if (!is.null(values) && !all(is.finite(values))) {
    stop("Matrix input contains NA, NaN, or Inf entries.", call. = FALSE)
  }
  invisible(TRUE)
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
      "pass a complex sparse matrix as complex_operator(re, im) (real and imaginary parts); ",
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
#' @return An `eigencore_operator` representing the adjoint map.
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
  if (is.environment(x$metadata$native_composite)) {
    # C52: the adjoint of a native composite is itself a native composite.
    return(native_algebra_operator(
      dim = rev(x$dim),
      apply = x$apply_adjoint,
      apply_adjoint = x$apply,
      dtype = x$dtype,
      structure = x$structure,
      name = paste0("adjoint(", x$name, ")"),
      metadata = list(
        parent = x,
        fused = "adjoint",
        algebra = "adjoint",
        native = isTRUE(x$metadata$native)
      )
    ))
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
#' Native composed-operator kernel (C52): an environment holding the R spec
#' (see native_composite_spec() in operator_algebra.R) and the external
#' pointer built from it. The spec keeps the leaf storage alive and lets a
#' pointer cleared by serialisation be rebuilt on the next apply. Returns NULL
#' when the native build rejects the spec.
new_native_composite_kernel <- function(spec) {
  ptr <- tryCatch(
    .Call("eigencore_composite_operator_build", spec, PACKAGE = "eigencore"),
    error = function(e) NULL
  )
  if (is.null(ptr)) {
    return(NULL)
  }
  kernel <- new.env(parent = emptyenv())
  kernel$spec <- spec
  kernel$ptr <- ptr
  kernel
}

#' @keywords internal
native_composite_block_apply <- function(kernel, X, alpha = 1, beta = 0,
                                         Y = NULL, adjoint = FALSE) {
  X <- as.matrix(X)
  if (!is.double(X)) {
    storage.mode(X) <- "double"
  }
  Y <- block_apply_y(Y, beta)
  if (!is.null(Y) && !is.double(Y)) {
    storage.mode(Y) <- "double"
  }
  alpha <- as.numeric(alpha)
  beta <- as.numeric(beta)
  adjoint <- isTRUE(adjoint)
  out <- .Call("eigencore_composite_block_apply", kernel$ptr, X, alpha, beta,
               Y, adjoint, PACKAGE = "eigencore")
  if (is.null(out)) {
    # Pointer cleared by serialisation: rebuild it from the spec.
    kernel$ptr <- .Call("eigencore_composite_operator_build", kernel$spec,
                        PACKAGE = "eigencore")
    out <- .Call("eigencore_composite_block_apply", kernel$ptr, X, alpha,
                 beta, Y, adjoint, PACKAGE = "eigencore")
  }
  out
}

#' @keywords internal
#' Attach a native composite kernel to an operator's apply closures as
#' attr(, "eigencore_native_kernel") = list(kernel, adjoint, record,
#' namespace). eigencore_r_operator_apply() (the R-callback operator every
#' matrix-free native solver uses) recognises it and applies the composite
#' natively instead of evaluating the closure; record(cols) books that apply
#' in the active work context exactly as linear_operator()'s wrapper would.
attach_native_composite_kernel <- function(op, kernel) {
  frame <- environment(op$apply)
  identity_fn <- frame$work_identity
  adjoint_view <- isTRUE(frame$adjoint_view)
  recorder <- function(kind) {
    force(kind)
    function(cols) {
      token <- work_operator_enter(
        identity_fn, kind = kind, X = matrix(0, 0L, cols)
      )
      work_operator_exit(token)
      invisible(NULL)
    }
  }
  ns <- environment(attach_native_composite_kernel)
  attr(op$apply, "eigencore_native_kernel") <- list(
    kernel, FALSE, recorder(if (adjoint_view) "adjoint" else "operator"), ns
  )
  if (!is.null(op$apply_adjoint)) {
    attr(op$apply_adjoint, "eigencore_native_kernel") <- list(
      kernel, TRUE, recorder(if (adjoint_view) "operator" else "adjoint"), ns
    )
  }
  op
}

#' @keywords internal
#' TRUE when `op` applies through a native composed-operator kernel (C52).
has_native_composite_kernel <- function(op) {
  is.environment(op$metadata$native_composite %||% NULL)
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
#' Output argument for the native block applies (P11). The C entry points
#' allocate the result themselves when Y is NULL (beta is then treated as 0),
#' so no zero matrix is built here; with beta == 0 a supplied Y is never read
#' and is dropped unless it carries dimnames, which the C side copies onto the
#' fresh output.
block_apply_y <- function(Y, beta) {
  if (is.null(Y)) {
    return(NULL)
  }
  if (length(beta) == 1L && isTRUE(beta == 0) && is.null(dimnames(Y))) {
    return(NULL)
  }
  as.matrix(Y)
}

#' @keywords internal
dense_block_apply <- function(A, X, alpha = 1, beta = 0, Y = NULL, transpose = FALSE) {
  X <- as.matrix(X)
  Y <- block_apply_y(Y, beta)
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
  Y <- block_apply_y(Y, beta)
  if (!is.null(Y) && !is.complex(Y)) {
    Y <- Y + 0i
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
  Y <- block_apply_y(Y, beta)
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
csc_column_moments <- function(A) {
  .Call(
    "eigencore_csc_column_moments",
    methods::slot(A, "p"),
    methods::slot(A, "x"),
    methods::slot(A, "Dim"),
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
csc_matrix_as_operator <- function(x, input_storage = class(x)[[1L]]) {
  force(input_storage)
  dim_x <- dim(x)
  symmetric_storage <- inherits(x, "symmetricMatrix")
  moments <- if (inherits(x, "dgCMatrix") &&
                 !inherits(x, "symmetricMatrix")) {
    csc_column_moments(x)
  } else {
    NULL
  }
  frobenius_norm <- if (is.null(moments)) {
    as.numeric(Matrix::norm(x, type = "F"))
  } else {
    sqrt(sum(moments$sum_squares))
  }
  x <- csc_native_matrix(x)
  if (is.null(moments)) {
    moments <- csc_column_moments(x)
  }
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
      frobenius_norm = frobenius_norm,
      column_sums = moments$sum,
      column_sum_squares = moments$sum_squares,
      column_means = moments$mean,
      column_centered_sum_squares = moments$centered_sum_squares
    )
  )
}

#' @keywords internal
diagonal_block_apply <- function(A, X, alpha = 1, beta = 0, Y = NULL) {
  X <- as.matrix(X)
  Y <- block_apply_y(Y, beta)
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
#' Finiteness and symmetry of a double matrix in one native pass: returns
#' c(all_finite, is_symmetric) with the semantics of is_square_symmetric().
dense_finite_symmetric <- function(x, tol = sqrt(.Machine$double.eps)) {
  .Call("eigencore_dense_finite_symmetric", x, as.numeric(tol), PACKAGE = "eigencore")
}

#' @keywords internal
is_square_symmetric <- function(x, tol = sqrt(.Machine$double.eps)) {
  d <- dim(x)
  if (length(d) != 2L || d[1L] != d[2L]) {
    return(FALSE)
  }
  if (inherits(x, "Matrix")) {
    if (inherits(x, "symmetricMatrix")) {
      return(TRUE)
    }
    # Matrix::isSymmetric(tol =) is absolute; scale it by the largest entry
    # so a rescaled nonsymmetric matrix is never classified as symmetric.
    values <- if (methods::.hasSlot(x, "x")) methods::slot(x, "x") else NULL
    scale <- if (length(values)) max(abs(values)) else 1
    if (!is.finite(scale)) {
      return(FALSE)
    }
    return(isTRUE(Matrix::isSymmetric(x, tol = tol * scale)))
  }
  if (is.matrix(x) && is.double(x)) {
    return(isTRUE(.Call("eigencore_dense_is_symmetric", x, as.numeric(tol), PACKAGE = "eigencore")))
  }
  if (is.matrix(x) && is.complex(x)) {
    # isSymmetric.matrix() goes through all.equal(), which switches to an
    # ABSOLUTE comparison when mean(Mod(x)) < tol: a 1e-12-scale general
    # complex matrix was classified Hermitian and solved with zheev (C38 for
    # complex input; found by the oracle sweep). Use a relative test.
    scale <- max(Mod(x))
    if (!is.finite(scale)) {
      return(FALSE)
    }
    if (scale == 0) {
      return(TRUE)
    }
    return(max(Mod(x - Conj(t(x)))) <= tol * scale)
  }
  isTRUE(isSymmetric.matrix(x, tol = tol))
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
