#' Compose two operators.
compose <- function(A, B, name = NULL) {
  A <- as_operator(A)
  B <- as_operator(B)
  if (A$dim[2L] != B$dim[1L]) {
    stop("Cannot compose operators with dimensions ",
         paste(A$dim, collapse = " x "), " and ",
         paste(B$dim, collapse = " x "), ".", call. = FALSE)
  }

  src <- source_or_null(A)
  rhs <- source_or_null(B)
  # NN-3: only fold source when both operands are dense double matrices.
  # A sparse source threaded through metadata$source is later materialized
  # via as.matrix(src) by downstream certification helpers, silently
  # densifying a potentially huge product. A non-dense source must stay
  # NULL here so the operator remains matrix-free.
  source <- if (is_dense_double_matrix(src) && is_dense_double_matrix(rhs)) {
    src %*% rhs
  } else {
    NULL
  }

  linear_operator(
    dim = c(A$dim[1L], B$dim[2L]),
    apply = function(X, alpha = 1, beta = 0, Y = NULL) {
      Z <- apply_operator(B, X)
      apply_operator(A, Z, alpha = alpha, beta = beta, Y = Y)
    },
    apply_adjoint = if (!is.null(A$apply_adjoint) && !is.null(B$apply_adjoint)) {
      function(X, alpha = 1, beta = 0, Y = NULL) {
        Z <- apply_adjoint_operator(A, X)
        apply_adjoint_operator(B, Z, alpha = alpha, beta = beta, Y = Y)
      }
    } else {
      NULL
    },
    dtype = common_dtype(A, B),
    # AB is Hermitian only if A and B are Hermitian AND commute (e.g. when
    # B == A and A is Hermitian). The previous heuristic compared A$name to
    # B$name, but names are user-controlled metadata and must never drive
    # certificate-relevant flags. Default to general() and let the caller
    # explicitly call symmetric() when they can prove the composition is
    # Hermitian (cf. crossprod_operator() which always returns hermitian()).
    structure = general(),
    name = name %||% paste0("compose(", A$name, ",", B$name, ")"),
    metadata = list(left = A, right = B, source = source, native = FALSE)
  )
}

#' Sum compatible operators.
operator_sum <- function(..., name = NULL) {
  ops <- lapply(list(...), as_operator)
  if (!length(ops)) {
    stop("At least one operator is required.", call. = FALSE)
  }
  dims <- vapply(ops, function(op) paste(op$dim, collapse = "x"), character(1))
  if (length(unique(dims)) != 1L) {
    stop("All summed operators must have the same dimensions.", call. = FALSE)
  }

  # NN-3: only fold source when every term carries a dense double-matrix
  # source. A single sparse-source term must collapse the fold to NULL so
  # the summed operator stays matrix-free.
  source <- Reduce(function(a, b) {
    if (is_dense_double_matrix(a) && is_dense_double_matrix(b)) a + b else NULL
  }, lapply(ops, source_or_null))

  linear_operator(
    dim = ops[[1L]]$dim,
    apply = function(X, alpha = 1, beta = 0, Y = NULL) {
      acc <- Reduce(`+`, lapply(ops, function(op) apply_operator(op, X)))
      combine_block(alpha * acc, beta = beta, Y = Y)
    },
    apply_adjoint = if (all(vapply(ops, function(op) !is.null(op$apply_adjoint), logical(1)))) {
      function(X, alpha = 1, beta = 0, Y = NULL) {
        acc <- Reduce(`+`, lapply(ops, function(op) apply_adjoint_operator(op, X)))
        combine_block(alpha * acc, beta = beta, Y = Y)
      }
    } else {
      NULL
    },
    dtype = ops[[1L]]$dtype,
    structure = if (all(vapply(ops, function(op) identical(op$structure$kind, "hermitian"), logical(1)))) hermitian() else general(),
    name = name %||% "operator_sum",
    metadata = list(terms = ops, source = source, native = FALSE)
  )
}

#' Multiply an operator by a scalar.
operator_scale <- function(A, scalar, name = NULL) {
  A <- as_operator(A)
  scalar <- as.numeric(scalar)
  if (length(scalar) != 1L || !is.finite(scalar)) {
    stop("scalar must be one finite numeric value.", call. = FALSE)
  }
  fused <- native_scaled_operator_or_null(A, scalar, axis = "scalar", name = name)
  if (!is.null(fused)) {
    return(fused)
  }
  source <- source_or_null(A)
  # NN-3: only fold a dense-matrix source. A sparse source would later be
  # as.matrix()-ed by certification helpers and silently densify.
  if (is_dense_double_matrix(source)) {
    source <- scalar * source
  } else {
    source <- NULL
  }

  linear_operator(
    dim = A$dim,
    apply = function(X, alpha = 1, beta = 0, Y = NULL) {
      apply_operator(A, X, alpha = alpha * scalar, beta = beta, Y = Y)
    },
    apply_adjoint = if (!is.null(A$apply_adjoint)) {
      function(X, alpha = 1, beta = 0, Y = NULL) {
        apply_adjoint_operator(A, X, alpha = alpha * scalar, beta = beta, Y = Y)
      }
    } else {
      NULL
    },
    dtype = A$dtype,
    structure = A$structure,
    name = name %||% paste0(scalar, "*", A$name),
    metadata = list(parent = A, scalar = scalar, source = source, native = FALSE)
  )
}

#' Scale operator rows.
scale_rows <- function(A, weights, name = NULL) {
  A <- as_operator(A)
  weights <- as.numeric(weights)
  if (length(weights) != A$dim[1L]) {
    stop("row weights length must equal operator row dimension.", call. = FALSE)
  }
  fused <- native_scaled_operator_or_null(A, weights, axis = "rows", name = name)
  if (!is.null(fused)) {
    return(fused)
  }
  source <- source_or_null(A)
  if (is_dense_double_matrix(source)) {
    source <- weights * source
  } else {
    source <- NULL
  }

  linear_operator(
    dim = A$dim,
    apply = function(X, alpha = 1, beta = 0, Y = NULL) {
      out <- weights * apply_operator(A, X)
      combine_block(alpha * out, beta = beta, Y = Y)
    },
    apply_adjoint = if (!is.null(A$apply_adjoint)) {
      function(X, alpha = 1, beta = 0, Y = NULL) {
        apply_adjoint_operator(A, weights * X, alpha = alpha, beta = beta, Y = Y)
      }
    } else {
      NULL
    },
    dtype = A$dtype,
    structure = general(),
    name = name %||% paste0("scale_rows(", A$name, ")"),
    metadata = list(parent = A, weights = weights, axis = "rows", source = source, native = FALSE)
  )
}

#' Scale operator columns.
scale_cols <- function(A, weights, name = NULL) {
  A <- as_operator(A)
  weights <- as.numeric(weights)
  if (length(weights) != A$dim[2L]) {
    stop("column weights length must equal operator column dimension.", call. = FALSE)
  }
  fused <- native_scaled_operator_or_null(A, weights, axis = "cols", name = name)
  if (!is.null(fused)) {
    return(fused)
  }
  source <- source_or_null(A)
  if (is_dense_double_matrix(source)) {
    source <- sweep(source, 2L, weights, `*`)
  } else {
    source <- NULL
  }

  linear_operator(
    dim = A$dim,
    apply = function(X, alpha = 1, beta = 0, Y = NULL) {
      apply_operator(A, weights * X, alpha = alpha, beta = beta, Y = Y)
    },
    apply_adjoint = if (!is.null(A$apply_adjoint)) {
      function(X, alpha = 1, beta = 0, Y = NULL) {
        out <- weights * apply_adjoint_operator(A, X)
        combine_block(alpha * out, beta = beta, Y = Y)
      }
    } else {
      NULL
    },
    dtype = A$dtype,
    structure = general(),
    name = name %||% paste0("scale_cols(", A$name, ")"),
    metadata = list(parent = A, weights = weights, axis = "cols", source = source, native = FALSE)
  )
}

#' Center an operator by rows or columns.
center <- function(A, rows = FALSE, columns = TRUE, row_means = NULL,
                   col_means = NULL, name = NULL) {
  A <- as_operator(A)
  source <- source_or_null(A)
  # NN-3: matrix-free centering must stay matrix-free when the underlying
  # source is sparse. Treat a non-dense source as if no source were
  # available so col_means/row_means must be supplied explicitly rather
  # than computed by colMeans()/rowMeans() on a sparse object that the
  # downstream centered_source path would later as.matrix() into a huge
  # dense fallback.
  if (!is.null(source) && !is_dense_double_matrix(source)) {
    source <- NULL
  }
  if (is.null(col_means) && isTRUE(columns)) {
    if (is.null(source)) {
      stop("col_means must be supplied for matrix-free column centering.", call. = FALSE)
    }
    col_means <- colMeans(source)
  }
  if (is.null(row_means) && isTRUE(rows)) {
    if (is.null(source)) {
      stop("row_means must be supplied for matrix-free row centering.", call. = FALSE)
    }
    row_means <- rowMeans(source)
  }

  if (isTRUE(columns) && length(col_means) != A$dim[2L]) {
    stop("col_means length must equal operator column dimension.", call. = FALSE)
  }
  if (isTRUE(rows) && length(row_means) != A$dim[1L]) {
    stop("row_means length must equal operator row dimension.", call. = FALSE)
  }

  centered_source <- source
  if (!is.null(centered_source) && isTRUE(columns)) {
    centered_source <- sweep(centered_source, 2L, col_means, `-`)
  }
  if (!is.null(centered_source) && isTRUE(rows)) {
    centered_source <- sweep(centered_source, 1L, row_means, `-`)
  }

  linear_operator(
    dim = A$dim,
    apply = function(X, alpha = 1, beta = 0, Y = NULL) {
      out <- apply_operator(A, X)
      if (isTRUE(columns)) {
        out <- out - matrix(1, A$dim[1L], 1L) %*% matrix(crossprod(col_means, X), nrow = 1L)
      }
      if (isTRUE(rows)) {
        out <- out - matrix(row_means, ncol = 1L) %*% matrix(colSums(X), nrow = 1L)
      }
      combine_block(alpha * out, beta = beta, Y = Y)
    },
    apply_adjoint = if (!is.null(A$apply_adjoint)) {
      function(X, alpha = 1, beta = 0, Y = NULL) {
        out <- apply_adjoint_operator(A, X)
        if (isTRUE(columns)) {
          out <- out - matrix(col_means, ncol = 1L) %*% matrix(colSums(X), nrow = 1L)
        }
        if (isTRUE(rows)) {
          out <- out - matrix(1, A$dim[2L], 1L) %*% matrix(crossprod(row_means, X), nrow = 1L)
        }
        combine_block(alpha * out, beta = beta, Y = Y)
      }
    } else {
      NULL
    },
    dtype = A$dtype,
    structure = general(),
    name = name %||% paste0("center(", A$name, ")"),
    metadata = list(
      parent = A,
      rows = rows,
      columns = columns,
      row_means = row_means,
      col_means = col_means,
      source = centered_source,
      native = FALSE
    )
  )
}

#' Mark an operator as symmetric/Hermitian.
symmetric_operator <- function(A, validate = TRUE, tol = 1e-10) {
  A <- as_operator(A)
  if (A$dim[1L] != A$dim[2L]) {
    stop("A symmetric operator must be square.", call. = FALSE)
  }
  if (isTRUE(validate)) {
    check_adjoint(A, trials = 5L, tol = tol)
  }
  A$structure <- hermitian()
  A$name <- paste0("symmetric(", A$name, ")")
  A
}

#' Create A^* A as an operator.
crossprod_operator <- function(A, name = NULL) {
  A <- as_operator(A)
  if (is.null(A$apply_adjoint)) {
    stop("Operator does not define apply_adjoint().", call. = FALSE)
  }
  linear_operator(
    dim = c(A$dim[2L], A$dim[2L]),
    apply = function(X, alpha = 1, beta = 0, Y = NULL) {
      Z <- apply_operator(A, X)
      apply_adjoint_operator(A, Z, alpha = alpha, beta = beta, Y = Y)
    },
    apply_adjoint = function(X, alpha = 1, beta = 0, Y = NULL) {
      Z <- apply_operator(A, X)
      apply_adjoint_operator(A, Z, alpha = alpha, beta = beta, Y = Y)
    },
    dtype = A$dtype,
    structure = hermitian(),
    name = name %||% paste0("crossprod(", A$name, ")"),
    metadata = list(parent = A, fused = "crossprod", native = isTRUE(A$metadata$native))
  )
}

#' Check an operator adjoint identity.
check_adjoint <- function(A, trials = 20, tol = 1e-12, seed = NULL) {
  A <- as_operator(A)
  if (is.null(A$apply_adjoint)) {
    stop("Operator does not define apply_adjoint().", call. = FALSE)
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }

  errors <- numeric(trials)
  for (i in seq_len(trials)) {
    block <- sample.int(3L, 1L)
    x <- matrix(stats::rnorm(A$dim[2L] * block), A$dim[2L], block)
    y <- matrix(stats::rnorm(A$dim[1L] * block), A$dim[1L], block)
    lhs <- sum(apply_operator(A, x) * y)
    rhs <- sum(x * apply_adjoint_operator(A, y))
    denom <- max(1, abs(lhs), abs(rhs))
    errors[[i]] <- abs(lhs - rhs) / denom
  }

  result <- list(
    passed = all(errors <= tol),
    tolerance = tol,
    max_error = max(errors),
    errors = errors,
    trials = trials
  )
  class(result) <- "eigencore_adjoint_check"
  if (!result$passed) {
    stop("Adjoint check failed; max relative error = ",
         format(result$max_error), ".", call. = FALSE)
  }
  result
}

#' @keywords internal
source_or_null <- function(A) {
  A$metadata$source %||% NULL
}

#' @keywords internal
#' Returns TRUE only for plain dense double matrices, the only kind of
#' source that operator_algebra is allowed to fold and rethread through
#' metadata$source. Sparse, integer, or non-matrix inputs must collapse
#' the fold to NULL so downstream certification helpers don't silently
#' as.matrix() a huge product (NN-3 in AGENTS.md).
is_dense_double_matrix <- function(x) {
  !is.null(x) && is.matrix(x) && is.double(x) && !inherits(x, "Matrix")
}

#' @keywords internal
#' Classifies an eigencore_operator into the native kernel kind it can be
#' fed to: "csc" (dgCMatrix metadata), "dense" (dense double source), or
#' NA_character_ when no native kernel is available. Centralizes the
#' `identical(storage, "dgCMatrix") || (is.matrix(source) && is.double(source))`
#' pattern that recurs across solve.R / reference_*.R predicates.
native_kernel_kind <- function(op) {
  if (identical(op$metadata$storage %||% NULL, "dgCMatrix")) {
    return("csc")
  }
  src <- source_or_null(op)
  if (is.matrix(src) && is.double(src)) {
    return("dense")
  }
  NA_character_
}

#' @keywords internal
has_native_kernel <- function(op) {
  !is.na(native_kernel_kind(op))
}

#' @keywords internal
common_dtype <- function(A, B) {
  if (identical(A$dtype, B$dtype)) A$dtype else "double"
}

#' @keywords internal
combine_block <- function(out, beta = 0, Y = NULL) {
  if (!is.null(Y) && beta != 0) {
    out <- out + beta * Y
  }
  out
}

#' @keywords internal
native_scaled_operator_or_null <- function(A, weights, axis, name = NULL) {
  if (!isTRUE(A$metadata$native)) {
    return(NULL)
  }
  scaled <- NULL
  storage <- NULL

  source <- A$metadata$source
  if (is_dense_double_matrix(source)) {
    scaled <- switch(
      axis,
      scalar = weights * source,
      rows = weights * source,
      cols = sweep(source, 2L, weights, `*`),
      NULL
    )
  } else {
    # Non-dense source falls through to the metadata$matrix path below;
    # never thread a sparse object back into metadata$source.
    source <- NULL
    matrix <- A$metadata$matrix
    if (is.null(matrix)) {
      return(NULL)
    }
    scaled <- switch(
      axis,
      scalar = weights * matrix,
      rows = Matrix::Diagonal(x = weights) %*% matrix,
      cols = matrix %*% Matrix::Diagonal(x = weights),
      NULL
    )
    if (inherits(matrix, "dgCMatrix") && inherits(scaled, "sparseMatrix")) {
      scaled <- methods::as(scaled, "dgCMatrix")
    }
    if (!(inherits(scaled, "dgCMatrix") || inherits(scaled, "ddiMatrix"))) {
      return(NULL)
    }
    storage <- class(scaled)[[1L]]
  }

  if (is.null(scaled)) {
    return(NULL)
  }
  op <- as_operator(scaled)
  op$name <- name %||% switch(
    axis,
    scalar = paste0(weights, "*", A$name),
    rows = paste0("scale_rows(", A$name, ")"),
    cols = paste0("scale_cols(", A$name, ")")
  )
  op$metadata$parent <- A
  op$metadata$fused <- switch(
    axis,
    scalar = "scalar_scale",
    rows = "scale_rows",
    cols = "scale_cols"
  )
  op$metadata$axis <- axis
  op$metadata$weights <- weights
  if (!is.null(storage)) {
    op$metadata$storage <- storage
  }
  op
}
