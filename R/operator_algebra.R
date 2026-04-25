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
  source <- if (!is.null(src) && !is.null(rhs)) src %*% rhs else NULL

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
    structure = if (A$dim[1L] == B$dim[2L] && identical(A$name, B$name)) hermitian() else general(),
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

  source <- Reduce(function(a, b) {
    if (is.null(a) || is.null(b)) NULL else a + b
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
  source <- source_or_null(A)
  if (!is.null(source)) {
    source <- scalar * source
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
  source <- source_or_null(A)
  if (!is.null(source)) {
    source <- weights * source
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
  source <- source_or_null(A)
  if (!is.null(source)) {
    source <- sweep(source, 2L, weights, `*`)
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
  compose(adjoint(A), A, name = name %||% paste0("crossprod(", A$name, ")"))
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
