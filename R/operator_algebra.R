#' Compose two operators.
#'
#' @param A Left operator-like object.
#' @param B Right operator-like object.
#' @param name Optional label for the composed operator.
#' @return An `eigencore_operator` representing the composition `A %*% B`.
#' @details When both operands are built-in explicit matrices (dense double
#'   or `Matrix` storage), the product is formed once and wrapped as a native
#'   explicit operator (`metadata$fused == "compose"`) only when that is cheap
#'   and memory-safe: always for a small product (at most 65,536 entries);
#'   for a larger dense product only when it holds no more entries than the
#'   two factors together and costs at most 2^30 multiply-adds; for a larger
#'   sparse product only when a structural bound on its nonzeros stays within
#'   four times the factors' nonzeros and below a quarter of its entries, so
#'   a sparse product is never densified. Otherwise (and for callback or
#'   mixed operands) the result is a lazy composition that applies `B` then
#'   `A` and stores no product.
compose <- function(A, B, name = NULL) {
  A <- as_operator(A)
  B <- as_operator(B)
  if (A$dim[2L] != B$dim[1L]) {
    stop("Cannot compose operators with dimensions ",
         paste(A$dim, collapse = " x "), " and ",
         paste(B$dim, collapse = " x "), ".", call. = FALSE)
  }

  # The product, when it is materialised at all, is formed exactly once.
  # Native operands go through native_compose_operator_or_null(); for
  # non-native dense sources the same size policy decides whether the
  # product is folded into metadata$source. NN-3: never fold a non-dense
  # source (downstream certification helpers would as.matrix() it).
  fused <- native_compose_operator_or_null(A, B, name = name)
  if (!is.null(fused)) {
    return(fused)
  }
  src <- source_or_null(A)
  rhs <- source_or_null(B)
  source <- if (!(isTRUE(A$metadata$native) && isTRUE(B$metadata$native)) &&
                is_dense_double_matrix(src) && is_dense_double_matrix(rhs) &&
                algebra_dense_product_ok(nrow(src), ncol(src), ncol(rhs))) {
    src %*% rhs
  } else {
    NULL
  }

  native_algebra_operator(
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
    metadata = list(left = A, right = B, source = source, native = FALSE,
                    algebra = "compose")
  )
}

#' Sum compatible operators.
#'
#' @param ... Operator-like objects with identical dimensions.
#' @param name Optional label for the summed operator.
#' @keywords internal
operator_sum <- function(..., name = NULL) {
  ops <- lapply(list(...), as_operator)
  if (!length(ops)) {
    stop("At least one operator is required.", call. = FALSE)
  }
  dims <- vapply(ops, function(op) paste(op$dim, collapse = "x"), character(1))
  if (length(unique(dims)) != 1L) {
    stop("All summed operators must have the same dimensions.", call. = FALSE)
  }

  # Native explicit terms are summed once into a native explicit operator
  # (a sum never grows memory beyond one term and never densifies a sparse
  # sum). Only when that does not apply is a dense source folded here.
  fused <- native_sum_operator_or_null(ops, name = name)
  if (!is.null(fused)) {
    return(fused)
  }
  # NN-3: only fold source when every term carries a dense double-matrix
  # source. A single sparse-source term must collapse the fold to NULL so
  # the summed operator stays matrix-free.
  sources <- lapply(ops, source_or_null)
  source <- if (all(vapply(sources, is_dense_double_matrix, logical(1)))) {
    Reduce(`+`, sources)
  } else {
    NULL
  }

  native_algebra_operator(
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
    metadata = list(terms = ops, source = source, native = FALSE,
                    algebra = "sum")
  )
}

#' Multiply an operator by a scalar.
#'
#' @param A Operator-like object.
#' @param scalar Finite numeric scalar multiplier.
#' @param name Optional label for the scaled operator.
#' @keywords internal
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

  native_algebra_operator(
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
    metadata = list(parent = A, scalar = scalar, source = source, native = FALSE,
                    algebra = "scale")
  )
}

#' Scale operator rows.
#'
#' @param A Operator-like object.
#' @param weights Numeric vector of row weights.
#' @param name Optional label for the scaled operator.
#' @return An `eigencore_operator` representing row-wise left scaling of `A`.
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

  native_algebra_operator(
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
    metadata = list(parent = A, weights = weights, axis = "rows", source = source, native = FALSE,
                    algebra = "scale_rows")
  )
}

#' Scale operator columns.
#'
#' @param A Operator-like object.
#' @param weights Numeric vector of column weights.
#' @param name Optional label for the scaled operator.
#' @return An `eigencore_operator` representing column-wise right scaling of
#'   `A`.
scale_cols <- function(A, weights, name = NULL) {
  A <- as_operator(A)
  weights <- as.numeric(weights)
  if (length(weights) != A$dim[2L]) {
    stop("column weights length must equal operator column dimension.", call. = FALSE)
  }
  if (any(!is.finite(weights))) {
    stop("column weights must be finite.", call. = FALSE)
  }
  fused <- native_centered_scaled_csc_operator_or_null(
    A, weights, name = name
  )
  if (!is.null(fused)) {
    return(fused)
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

  native_algebra_operator(
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
    metadata = list(parent = A, weights = weights, axis = "cols", source = source, native = FALSE,
                    algebra = "scale_cols")
  )
}

#' Center an operator by rows or columns.
#'
#' @param A Operator-like object.
#' @param rows Whether to subtract row means.
#' @param columns Whether to subtract column means.
#' @param row_means Optional row means. Required for matrix-free row centering
#'   when they cannot be derived without densifying.
#' @param col_means Optional column means. Required for matrix-free column
#'   centering when they cannot be derived without densifying.
#' @param name Optional label for the centered operator.
#' @details `row_means` and `col_means` are means of the uncentered `A`. When
#'   both `rows` and `columns` are `TRUE` the result is double centered: the
#'   grand mean is added back so every row and column of the result has mean
#'   zero.
#' @return An `eigencore_operator` representing the centered linear map.
center <- function(A, rows = FALSE, columns = TRUE, row_means = NULL,
                   col_means = NULL, name = NULL) {
  A <- as_operator(A)
  source <- source_or_null(A)
  matrix_source <- A$metadata$matrix %||% NULL
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
    if (is_dense_double_matrix(source)) {
      col_means <- colMeans(source)
    } else if (inherits(matrix_source, "Matrix")) {
      col_means <- Matrix::colMeans(matrix_source)
    } else {
      stop("col_means must be supplied for matrix-free column centering.", call. = FALSE)
    }
  }
  if (is.null(row_means) && isTRUE(rows)) {
    if (is_dense_double_matrix(source)) {
      row_means <- rowMeans(source)
    } else if (inherits(matrix_source, "Matrix")) {
      row_means <- Matrix::rowMeans(matrix_source)
    } else {
      stop("row_means must be supplied for matrix-free row centering.", call. = FALSE)
    }
  }

  if (isTRUE(columns) && length(col_means) != A$dim[2L]) {
    stop("col_means length must equal operator column dimension.", call. = FALSE)
  }
  if (isTRUE(rows) && length(row_means) != A$dim[1L]) {
    stop("row_means length must equal operator row dimension.", call. = FALSE)
  }

  if (isTRUE(rows) && isTRUE(columns)) {
    # Double centering is A - 1 c' - r 1' + g 1 1' with g the grand mean.
    # Folding g into the row means keeps every apply path rank-2.
    row_means <- row_means - mean(col_means)
  }

  centered_source <- source
  if (!is.null(centered_source) && isTRUE(columns)) {
    centered_source <- sweep(centered_source, 2L, col_means, `-`)
  }
  if (!is.null(centered_source) && isTRUE(rows)) {
    centered_source <- sweep(centered_source, 1L, row_means, `-`)
  }
  if (is_dense_double_matrix(centered_source)) {
    fused <- as_operator(centered_source)
    fused$name <- name %||% paste0("center(", A$name, ")")
    fused$metadata <- modifyList(
      fused$metadata,
      list(
        parent = A,
        fused = "center",
        rows = rows,
        columns = columns,
        row_means = row_means,
        col_means = col_means
      )
    )
    return(fused)
  }
  fused <- native_centered_sparse_operator_or_null(
    A,
    rows = rows,
    columns = columns,
    row_means = row_means,
    col_means = col_means,
    name = name
  )
  if (!is.null(fused)) {
    return(fused)
  }

  native_algebra_operator(
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
      native = FALSE,
      algebra = "center"
    )
  )
}

#' @keywords internal
csc_centered_block_apply <- function(A, X, alpha = 1, beta = 0, Y = NULL,
                                     transpose = FALSE, rows = FALSE,
                                     columns = TRUE, row_means = NULL,
                                     col_means = NULL) {
  X <- as.matrix(X)
  Y <- block_apply_y(Y, beta)
  .Call(
    "eigencore_csc_centered_block_apply",
    methods::slot(A, "i"),
    methods::slot(A, "p"),
    methods::slot(A, "x"),
    methods::slot(A, "Dim"),
    as.numeric(row_means %||% numeric()),
    as.numeric(col_means %||% numeric()),
    isTRUE(rows),
    isTRUE(columns),
    X,
    as.numeric(alpha),
    as.numeric(beta),
    Y,
    isTRUE(transpose),
    PACKAGE = "eigencore"
  )
}

#' @keywords internal
csc_centered_scaled_block_apply <- function(
    A, col_means, weights, X, alpha = 1, beta = 0, Y = NULL,
    transpose = FALSE) {
  X <- as.matrix(X)
  Y <- block_apply_y(Y, beta)
  .Call(
    "eigencore_csc_centered_scaled_block_apply",
    methods::slot(A, "i"),
    methods::slot(A, "p"),
    methods::slot(A, "x"),
    methods::slot(A, "Dim"),
    as.numeric(col_means),
    as.numeric(weights),
    X,
    as.numeric(alpha),
    as.numeric(beta),
    Y,
    isTRUE(transpose),
    PACKAGE = "eigencore"
  )
}

#' Mark an operator as symmetric/Hermitian.
#'
#' @param A Operator-like object.
#' @param validate Whether to check the adjoint identity before marking the
#'   operator symmetric.
#' @param tol Relative tolerance for the adjoint check.
#' @return An `eigencore_operator` with Hermitian structure metadata.
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
#'
#' @param A Operator-like object with an adjoint implementation.
#' @param name Optional label for the cross-product operator.
#' @return A Hermitian `eigencore_operator` representing `A^* A`.
#' @details For a built-in explicit `A`, `A^* A` is formed once and wrapped as
#'   a native explicit operator (`metadata$materialized_crossprod`) only when
#'   that is cheap and memory-safe: always when it has at most 65,536
#'   entries; for a larger dense `A` only when `A` has no more columns than
#'   rows and forming it costs at most 2^30 multiply-adds; for a larger
#'   sparse `A` only when a structural bound on its nonzeros stays within
#'   four times `nnz(A)` and below a quarter of its entries (a sparse
#'   cross-product is never densified). Otherwise the result is lazy and
#'   applies `A` then its adjoint.
crossprod_operator <- function(A, name = NULL) {
  A <- as_operator(A)
  if (is.null(A$apply_adjoint)) {
    stop("Operator does not define apply_adjoint().", call. = FALSE)
  }
  fused <- native_crossprod_operator_or_null(A, name = name)
  if (!is.null(fused)) {
    return(fused)
  }
  native_algebra_operator(
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
    metadata = list(parent = A, fused = "crossprod", native = isTRUE(A$metadata$native),
                    algebra = "crossprod")
  )
}

#' Check an operator adjoint identity.
#'
#' @param A Operator-like object.
#' @param trials Number of random block trials.
#' @param tol Relative tolerance for each adjoint identity check.
#' @param seed Optional random seed for reproducible trials.
#' @return An `eigencore_adjoint_check` list with pass/fail status, tolerance,
#'   maximum relative error, per-trial errors, and trial count.
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
#' as.matrix() a huge product, preserving the no-silent-densification policy.
is_dense_double_matrix <- function(x) {
  !is.null(x) && is.matrix(x) && is.double(x) && !inherits(x, "Matrix")
}

# Materialisation policy for compose() / crossprod_operator() (P10). A
# product of built-in explicit operators is formed (once) only when that is
# cheap and memory-safe; otherwise the algebra stays lazy.
#' @keywords internal
algebra_small_product_entries <- function() 65536

#' @keywords internal
algebra_dense_product_flop_budget <- function() 2^30

#' @keywords internal
#' Whether the dense (m x k) %*% (k x p) product should be materialised:
#' always when it has at most algebra_small_product_entries() entries;
#' otherwise only when it holds no more entries than its operands (so
#' materialising never grows memory) and costs at most
#' algebra_dense_product_flop_budget() multiply-adds.
algebra_dense_product_ok <- function(m, k, p, operand_entries = NULL) {
  m <- as.double(m)
  k <- as.double(k)
  p <- as.double(p)
  entries <- m * p
  if (entries <= algebra_small_product_entries()) {
    return(TRUE)
  }
  operand_entries <- operand_entries %||% (m * k + k * p)
  entries <= operand_entries && m * k * p <= algebra_dense_product_flop_budget()
}

#' @keywords internal
#' Nonzero pattern counts of a Matrix-package object (per column, per row).
algebra_sparse_pattern_counts <- function(x) {
  g <- methods::as(methods::as(x, "CsparseMatrix"), "generalMatrix")
  list(
    nnz = length(g@i),
    col = diff(g@p),
    row = tabulate(g@i + 1L, nbins = nrow(g))
  )
}

#' @keywords internal
#' Whether the Matrix-package product A %*% B (or crossprod(A) when
#' crossprod = TRUE) should be materialised. Uses the structural bound
#' nnz(AB) <= sum_j nnz(A[, j]) * nnz(B[j, ]) (for A^T A, the sum of squared
#' row counts of A): materialise when the product is small, or when the
#' bound stays within four times the operands' nonzeros and below a quarter
#' of the product's entries, so a sparse product is never densified.
algebra_sparse_product_ok <- function(A, B, crossprod = FALSE) {
  m <- if (crossprod) ncol(A) else nrow(A)
  p <- ncol(B)
  entries <- as.double(m) * as.double(p)
  if (entries <= algebra_small_product_entries()) {
    return(TRUE)
  }
  ok <- tryCatch({
    a <- algebra_sparse_pattern_counts(A)
    if (crossprod) {
      bound <- sum(as.double(a$row)^2)
      operand_nnz <- a$nnz
    } else {
      b <- algebra_sparse_pattern_counts(B)
      bound <- sum(as.double(a$col) * as.double(b$row))
      operand_nnz <- a$nnz + b$nnz
    }
    bound <= 4 * operand_nnz && bound <= 0.25 * entries
  }, error = function(e) FALSE)
  isTRUE(ok)
}

#' @keywords internal
#' Classifies an eigencore_operator into the native kernel kind it can be
#' fed to: "csc" (dgCMatrix metadata), "centered_scaled_csc" (the native
#' sparse-PCA fusion), "dense" (dense double source), or NA_character_ when no
#' native kernel is available. Centralizes the
#' `identical(storage, "dgCMatrix") || (is.matrix(source) && is.double(source))`
#' pattern that recurs across solve.R / reference_*.R predicates.
native_kernel_kind <- function(op) {
  if (identical(op$metadata$storage %||% NULL, "dgCMatrix")) {
    return("csc")
  }
  if (identical(
    op$metadata$storage %||% NULL, "centered_scaled_dgCMatrix"
  )) {
    return("centered_scaled_csc")
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

# Native composed operators (C52). Lazy algebra over operators that all have
# native kernels (dense, CSC, centered / centered-scaled CSC, diagonal, and
# earlier composites) is compiled into one native expression kernel
# (src/native_operators.cpp): products, weighted sums, rank-one centering
# terms and adjoints, with reusable workspace for intermediates. The operator
# keeps its R apply closures for R-level callers (they call the kernel with
# one .Call), and the kernel is attached to those closures so every
# matrix-free native solver (Golub-Kahan, Arnoldi, block Lanczos) applies it
# without crossing back into R. Such operators carry
# metadata$storage == "native_composite" and metadata$native_composite (the
# kernel). native_kernel_kind() stays NA for them: they plan onto the
# matrix-free native paths, which is where the kernel removes the callbacks.
# options(eigencore.native_composite = FALSE) disables the kernel.

#' @keywords internal
native_composite_enabled <- function() {
  !identical(getOption("eigencore.native_composite", TRUE), FALSE)
}

#' @keywords internal
native_composite_csc_spec <- function(x) {
  list(
    type = "csc",
    i = methods::slot(x, "i"),
    p = methods::slot(x, "p"),
    x = methods::slot(x, "x"),
    dim = methods::slot(x, "Dim")
  )
}

#' @keywords internal
#' Native composite spec of an operator, or NULL when some part has no
#' native kernel. Leaves: dense double sources, dgCMatrix, centered and
#' centered-scaled dgCMatrix, ddiMatrix; inner nodes: an existing composite's
#' spec, or the lazy algebra recorded in metadata$algebra.
native_composite_spec <- function(op) {
  if (!inherits(op, "eigencore_operator") || !identical(op$dtype, "double")) {
    return(NULL)
  }
  kernel <- op$metadata$native_composite %||% NULL
  if (is.environment(kernel)) {
    return(kernel$spec)
  }
  kind <- native_kernel_kind(op)
  if (identical(kind, "dense")) {
    source <- source_or_null(op)
    if (!all(dim(source) == op$dim)) {
      return(NULL)
    }
    return(list(type = "dense", x = source))
  }
  storage <- op$metadata$storage %||% NULL
  if (identical(kind, "csc")) {
    matrix <- op$metadata$matrix
    if (!inherits(matrix, "dgCMatrix") || !all(dim(matrix) == op$dim)) {
      return(NULL)
    }
    return(native_composite_csc_spec(matrix))
  }
  if (identical(kind, "centered_scaled_csc")) {
    matrix <- op$metadata$base_matrix
    if (!inherits(matrix, "dgCMatrix")) {
      return(NULL)
    }
    spec <- native_composite_csc_spec(matrix)
    spec$type <- "centered_scaled_csc"
    spec$means <- as.numeric(op$metadata$col_means)
    spec$weights <- as.numeric(op$metadata$weights)
    return(spec)
  }
  if (identical(storage, "ddiMatrix")) {
    matrix <- op$metadata$matrix
    if (!inherits(matrix, "ddiMatrix")) {
      return(NULL)
    }
    values <- if (identical(methods::slot(matrix, "diag"), "U")) {
      rep(1, nrow(matrix))
    } else {
      as.numeric(methods::slot(matrix, "x"))
    }
    return(list(type = "diagonal", x = values))
  }
  if (identical(storage, "centered_dgCMatrix")) {
    matrix <- op$metadata$base_matrix
    if (!inherits(matrix, "dgCMatrix")) {
      return(NULL)
    }
    return(native_composite_centered_spec(
      native_composite_csc_spec(matrix), dim(matrix),
      rows = op$metadata$rows, columns = op$metadata$columns,
      row_means = op$metadata$row_means, col_means = op$metadata$col_means
    ))
  }
  parent <- op$metadata$parent %||% NULL
  if (identical(op$metadata$fused %||% NULL, "adjoint") &&
      inherits(parent, "eigencore_operator") &&
      all(rev(parent$dim) == op$dim)) {
    # adjoint() view of a native leaf (e.g. "adjoint:dgCMatrix").
    spec <- native_composite_spec(parent)
    return(if (is.null(spec)) NULL else native_composite_adjoint_spec(spec))
  }
  NULL
}

#' @keywords internal
#' A - 1 c^T - r 1^T (the centering convention of center(), whose row means
#' already absorb the grand mean when both sides are centered).
native_composite_centered_spec <- function(base, dim, rows, columns,
                                           row_means, col_means) {
  children <- list(base)
  if (isTRUE(columns)) {
    children[[length(children) + 1L]] <- list(
      type = "rank1", u = rep(-1, dim[[1L]]), v = as.numeric(col_means)
    )
  }
  if (isTRUE(rows)) {
    children[[length(children) + 1L]] <- list(
      type = "rank1", u = -as.numeric(row_means), v = rep(1, dim[[2L]])
    )
  }
  if (length(children) == 1L) {
    return(base)
  }
  list(type = "sum", children = children, weights = rep(1, length(children)))
}

#' @keywords internal
native_composite_product_spec <- function(factors) {
  children <- list()
  for (factor in factors) {
    if (identical(factor$type, "product")) {
      children <- c(children, factor$children)
    } else {
      children[[length(children) + 1L]] <- factor
    }
  }
  list(type = "product", children = children)
}

#' @keywords internal
native_composite_adjoint_spec <- function(spec) {
  if (identical(spec$type, "adjoint")) {
    return(spec$child)
  }
  list(type = "adjoint", child = spec)
}

#' @keywords internal
#' Spec of a lazy algebra node from its metadata (NULL when any operand has no
#' native kernel).
native_algebra_spec <- function(metadata, dim) {
  algebra <- metadata$algebra %||% ""
  spec_of <- native_composite_spec
  switch(
    algebra,
    compose = {
      left <- spec_of(metadata$left)
      right <- spec_of(metadata$right)
      if (is.null(left) || is.null(right)) NULL else
        native_composite_product_spec(list(left, right))
    },
    sum = {
      terms <- lapply(metadata$terms, spec_of)
      if (!length(terms) || any(vapply(terms, is.null, logical(1)))) NULL else
        list(type = "sum", children = terms, weights = rep(1, length(terms)))
    },
    scale = {
      parent <- spec_of(metadata$parent)
      if (is.null(parent)) NULL else
        list(type = "sum", children = list(parent),
             weights = as.numeric(metadata$scalar))
    },
    scale_rows = {
      parent <- spec_of(metadata$parent)
      if (is.null(parent)) NULL else native_composite_product_spec(list(
        list(type = "diagonal", x = as.numeric(metadata$weights)), parent
      ))
    },
    scale_cols = {
      parent <- spec_of(metadata$parent)
      if (is.null(parent)) NULL else native_composite_product_spec(list(
        parent, list(type = "diagonal", x = as.numeric(metadata$weights))
      ))
    },
    center = {
      parent <- spec_of(metadata$parent)
      if (is.null(parent)) NULL else native_composite_centered_spec(
        parent, dim, rows = metadata$rows, columns = metadata$columns,
        row_means = metadata$row_means, col_means = metadata$col_means
      )
    },
    crossprod = {
      parent <- spec_of(metadata$parent)
      if (is.null(parent)) NULL else native_composite_product_spec(list(
        native_composite_adjoint_spec(parent), parent
      ))
    },
    adjoint = {
      parent <- spec_of(metadata$parent)
      if (is.null(parent)) NULL else native_composite_adjoint_spec(parent)
    },
    NULL
  )
}

#' @keywords internal
#' linear_operator() for a lazy algebra node. When every operand has a native
#' kernel (and the node did not fold a dense source, which native solvers use
#' directly), the node is compiled into a native composite kernel: the apply
#' closures call it, and it is attached to them for the native solvers. The
#' given R closures remain the semantics otherwise.
native_algebra_operator <- function(dim, apply, apply_adjoint = NULL,
                                    dtype = "double", structure = general(),
                                    name = NULL, metadata = list()) {
  spec <- NULL
  if (identical(dtype, "double") && native_composite_enabled() &&
      !is_dense_double_matrix(metadata$source %||% NULL)) {
    spec <- tryCatch(
      native_algebra_spec(metadata, dim),
      error = function(e) NULL
    )
  }
  kernel <- if (is.null(spec)) NULL else new_native_composite_kernel(spec)
  if (is.null(kernel)) {
    return(linear_operator(
      dim = dim, apply = apply, apply_adjoint = apply_adjoint, dtype = dtype,
      structure = structure, name = name, metadata = metadata
    ))
  }
  native_apply <- function(X, alpha = 1, beta = 0, Y = NULL) {
    native_composite_block_apply(kernel, X, alpha = alpha, beta = beta,
                                 Y = Y, adjoint = FALSE)
  }
  native_apply_adjoint <- if (is.null(apply_adjoint)) {
    NULL
  } else {
    function(X, alpha = 1, beta = 0, Y = NULL) {
      native_composite_block_apply(kernel, X, alpha = alpha, beta = beta,
                                   Y = Y, adjoint = TRUE)
    }
  }
  metadata$storage <- "native_composite"
  metadata$native_composite <- kernel
  op <- linear_operator(
    dim = dim, apply = native_apply, apply_adjoint = native_apply_adjoint,
    dtype = dtype, structure = structure, name = name, metadata = metadata
  )
  attach_native_composite_kernel(op, kernel)
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

#' @keywords internal
native_sum_operator_or_null <- function(ops, name = NULL) {
  if (!all(vapply(ops, function(op) isTRUE(op$metadata$native), logical(1)))) {
    return(NULL)
  }
  dense_sources <- lapply(ops, source_or_null)
  if (all(vapply(dense_sources, is_dense_double_matrix, logical(1)))) {
    summed <- Reduce(`+`, dense_sources)
    op <- native_explicit_operator_or_null(
      summed,
      name = name %||% "operator_sum",
      fused = "sum",
      metadata = list(terms = ops)
    )
    if (!is.null(op)) {
      return(op)
    }
  }

  matrices <- lapply(ops, function(op) op$metadata$matrix %||% NULL)
  if (!all(vapply(matrices, inherits, logical(1), what = "Matrix"))) {
    return(NULL)
  }
  summed <- Reduce(`+`, matrices)
  native_explicit_operator_or_null(
    summed,
    name = name %||% "operator_sum",
    fused = "sum",
    metadata = list(terms = ops)
  )
}

#' @keywords internal
native_compose_operator_or_null <- function(A, B, name = NULL) {
  if (!isTRUE(A$metadata$native) || !isTRUE(B$metadata$native)) {
    return(NULL)
  }
  left_source <- source_or_null(A)
  right_source <- source_or_null(B)
  if (is_dense_double_matrix(left_source) && is_dense_double_matrix(right_source)) {
    if (!algebra_dense_product_ok(nrow(left_source), ncol(left_source),
                                  ncol(right_source))) {
      return(NULL)
    }
    composed <- left_source %*% right_source
    return(native_explicit_operator_or_null(
      composed,
      name = name %||% paste0("compose(", A$name, ",", B$name, ")"),
      fused = "compose",
      metadata = list(left = A, right = B)
    ))
  }

  left_matrix <- A$metadata$matrix
  right_matrix <- B$metadata$matrix
  if (!inherits(left_matrix, "Matrix") || !inherits(right_matrix, "Matrix")) {
    return(NULL)
  }
  if (!algebra_sparse_product_ok(left_matrix, right_matrix)) {
    return(NULL)
  }
  composed <- left_matrix %*% right_matrix
  native_explicit_operator_or_null(
    composed,
    name = name %||% paste0("compose(", A$name, ",", B$name, ")"),
    fused = "compose",
    metadata = list(left = A, right = B)
  )
}

#' @keywords internal
native_crossprod_operator_or_null <- function(A, name = NULL) {
  if (!isTRUE(A$metadata$native)) {
    return(NULL)
  }
  source <- source_or_null(A)
  if (is_dense_double_matrix(source)) {
    # A^T A has ncol^2 entries: materialise it only when small, or when it
    # is no larger than A itself (ncol <= nrow) and cheap to form.
    if (!algebra_dense_product_ok(ncol(source), nrow(source), ncol(source),
                                  operand_entries = length(source))) {
      return(NULL)
    }
    cp <- crossprod(source)
    op <- native_explicit_operator_or_null(
      cp,
      name = name %||% paste0("crossprod(", A$name, ")"),
      fused = "crossprod",
      metadata = list(parent = A, materialized_crossprod = TRUE)
    )
    if (!is.null(op)) {
      op$structure <- hermitian()
    }
    return(op)
  }

  matrix <- A$metadata$matrix
  if (!inherits(matrix, "Matrix")) {
    return(NULL)
  }
  if (!algebra_sparse_product_ok(matrix, matrix, crossprod = TRUE)) {
    return(NULL)
  }
  cp <- Matrix::crossprod(matrix)
  op <- native_explicit_operator_or_null(
    cp,
    name = name %||% paste0("crossprod(", A$name, ")"),
    fused = "crossprod",
    metadata = list(parent = A, materialized_crossprod = TRUE)
  )
  if (!is.null(op)) {
    op$structure <- hermitian()
  }
  op
}

#' @keywords internal
native_centered_sparse_operator_or_null <- function(A, rows, columns,
                                                   row_means = NULL,
                                                   col_means = NULL,
                                                   name = NULL) {
  if (!isTRUE(A$metadata$native)) {
    return(NULL)
  }
  matrix <- A$metadata$matrix %||% NULL
  if (!inherits(matrix, "dgCMatrix")) {
    return(NULL)
  }
  row_means <- if (isTRUE(rows)) as.numeric(row_means) else numeric()
  col_means <- if (isTRUE(columns)) as.numeric(col_means) else numeric()
  frobenius_norm <- if (isTRUE(columns) && !isTRUE(rows)) {
    centered_csc_frobenius_norm(A, matrix, col_means)
  } else {
    NULL
  }
  op <- linear_operator(
    dim = A$dim,
    apply = function(X, alpha = 1, beta = 0, Y = NULL) {
      csc_centered_block_apply(
        matrix, X, alpha = alpha, beta = beta, Y = Y,
        transpose = FALSE, rows = rows, columns = columns,
        row_means = row_means, col_means = col_means
      )
    },
    apply_adjoint = if (!is.null(A$apply_adjoint)) {
      function(X, alpha = 1, beta = 0, Y = NULL) {
        csc_centered_block_apply(
          matrix, X, alpha = alpha, beta = beta, Y = Y,
          transpose = TRUE, rows = rows, columns = columns,
          row_means = row_means, col_means = col_means
        )
      }
    } else {
      NULL
    },
    dtype = A$dtype,
    structure = general(),
    name = name %||% paste0("center(", A$name, ")"),
    metadata = list(
      parent = A,
      base_matrix = matrix,
      fused = "center",
      rows = rows,
      columns = columns,
      row_means = row_means,
      col_means = col_means,
      native = TRUE,
      storage = "centered_dgCMatrix",
      low_rank_correction = TRUE,
      frobenius_norm = frobenius_norm,
      column_sums = A$metadata$column_sums,
      column_sum_squares = A$metadata$column_sum_squares,
      column_means = A$metadata$column_means,
      column_centered_sum_squares = A$metadata$column_centered_sum_squares
    )
  )
  # Matrix-free native solvers (no dedicated centered-CSC entry) apply it
  # through the composite kernel instead of the per-apply R callback (C52).
  spec <- if (native_composite_enabled()) native_composite_spec(op) else NULL
  kernel <- if (is.null(spec)) NULL else new_native_composite_kernel(spec)
  if (!is.null(kernel)) {
    op <- attach_native_composite_kernel(op, kernel)
  }
  op
}

#' @keywords internal
centered_csc_frobenius_norm <- function(A, matrix, col_means) {
  base_means <- A$metadata$column_means %||% NULL
  base_centered_sum_squares <-
    A$metadata$column_centered_sum_squares %||% NULL
  if (is.null(base_means) || is.null(base_centered_sum_squares)) {
    moments <- csc_column_moments(matrix)
    base_means <- moments$mean
    base_centered_sum_squares <- moments$centered_sum_squares
  }
  sum_squares <- pmax(
    base_centered_sum_squares + nrow(matrix) * (base_means - col_means)^2,
    0
  )
  out <- sqrt(sum(sum_squares))
  if (is.finite(out)) out else NULL
}

#' @keywords internal
native_centered_scaled_csc_operator_or_null <- function(A, weights,
                                                        name = NULL) {
  if (!isTRUE(A$metadata$native) ||
      !identical(A$metadata$storage %||% NULL, "centered_dgCMatrix") ||
      !isTRUE(A$metadata$columns) || isTRUE(A$metadata$rows)) {
    return(NULL)
  }
  matrix <- A$metadata$base_matrix %||% NULL
  if (!inherits(matrix, "dgCMatrix")) {
    return(NULL)
  }
  col_means <- as.numeric(A$metadata$col_means %||% numeric())
  if (length(col_means) != ncol(matrix)) {
    return(NULL)
  }
  column_sums <- A$metadata$column_sums %||% NULL
  column_sum_squares <- A$metadata$column_sum_squares %||% NULL
  base_means <- A$metadata$column_means %||% NULL
  base_centered_sum_squares <-
    A$metadata$column_centered_sum_squares %||% NULL
  if (is.null(column_sums) || is.null(column_sum_squares) ||
      is.null(base_means) || is.null(base_centered_sum_squares)) {
    moments <- csc_column_moments(matrix)
    column_sums <- moments$sum
    column_sum_squares <- moments$sum_squares
    base_means <- moments$mean
    base_centered_sum_squares <- moments$centered_sum_squares
  }
  centered_sum_squares <- base_centered_sum_squares +
    nrow(matrix) * (base_means - col_means)^2
  roundoff <- 64 * .Machine$double.eps * pmax(
    abs(column_sum_squares),
    abs(2 * col_means * column_sums),
    abs(nrow(matrix) * col_means^2),
    1
  )
  if (any(centered_sum_squares < -roundoff)) {
    stop(
      "Centered CSC column moments are internally inconsistent.",
      call. = FALSE
    )
  }
  centered_sum_squares <- pmax(centered_sum_squares, 0)
  scaled_sum_squares <- weights^2 * centered_sum_squares
  frobenius_norm <- sqrt(sum(scaled_sum_squares))
  if (!is.finite(frobenius_norm)) {
    frobenius_norm <- NULL
  }

  linear_operator(
    dim = A$dim,
    apply = function(X, alpha = 1, beta = 0, Y = NULL) {
      csc_centered_scaled_block_apply(
        matrix, col_means, weights, X,
        alpha = alpha, beta = beta, Y = Y, transpose = FALSE
      )
    },
    apply_adjoint = function(X, alpha = 1, beta = 0, Y = NULL) {
      csc_centered_scaled_block_apply(
        matrix, col_means, weights, X,
        alpha = alpha, beta = beta, Y = Y, transpose = TRUE
      )
    },
    dtype = A$dtype,
    structure = general(),
    name = name %||% paste0("scale_cols(", A$name, ")"),
    metadata = list(
      parent = A,
      base_matrix = matrix,
      fused = "center_scale_cols",
      rows = FALSE,
      columns = TRUE,
      col_means = col_means,
      weights = weights,
      native = TRUE,
      storage = "centered_scaled_dgCMatrix",
      native_entry_point = "eigencore_golub_kahan_centered_scaled_csc",
      callback_boundary = FALSE,
      materialized_centered_matrix = FALSE,
      low_rank_correction = TRUE,
      column_sums = weights * (column_sums - nrow(matrix) * col_means),
      column_sum_squares = scaled_sum_squares,
      frobenius_norm = frobenius_norm
    )
  )
}

#' @keywords internal
native_explicit_operator_or_null <- function(x, name, fused, metadata = list()) {
  if (is_dense_double_matrix(x)) {
    op <- as_operator(x)
  } else if (inherits(x, "ddiMatrix") || inherits(x, "dgCMatrix")) {
    op <- as_operator(x)
  } else if (inherits(x, "sparseMatrix")) {
    x <- methods::as(methods::as(x, "generalMatrix"), "dgCMatrix")
    op <- as_operator(x)
  } else {
    return(NULL)
  }
  op$name <- name
  op$metadata <- modifyList(op$metadata, c(list(fused = fused), metadata))
  op
}
