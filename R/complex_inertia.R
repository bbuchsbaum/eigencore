# Native complex Hermitian inertia counts (dense).
#
# A dense complex Hermitian count factors A - sigma B by the Hermitian
# Bunch-Kaufman LDL^H (LAPACK zhetrf, src/complex_inertia.cpp) and reads the
# counts off the signs of D: 1x1 blocks are real, 2x2 blocks Hermitian with
# their inertia from determinant and trace. This replaces dsytrf on the real
# symmetric 2n x 2n embedding [Re -Im; Im Re] (8 n^3 / 3 real flops and
# 32 n^2 bytes against 4 n^3 / 3 and 16 n^2 for zhetrf). The tally feeds the
# same reliability test as the real dense count (inertia_tally_reliable): the
# backward bound is p * eps * (||M||_1 + growth) with p = 2n, the dimension
# of the real embedding the count used to run on, which also covers the
# slightly larger rounding constants of complex arithmetic.
#
# Two entry points:
#   * eigen_count() / inertia_context() on a dense complex matrix (counts of
#     the complex problem; inertia_context() calls
#     complex_inertia_fill_dense() instead of building the embedding);
#   * contexts over the real 2n embedding of a realified complex problem
#     (R/complex_hermitian.R): the realified dense operator carries its
#     complex source in `op$complex_embedding`, and inertia_context() asks
#     complex_inertia_embedding_context() first. Those contexts count the
#     real 2n problem, whose eigenvalues are the complex ones twice: counts
#     are doubled and ctx$n is 2n, so every caller of the real count
#     (completeness certificate, metric lower bound) runs unchanged.
#
# getOption("eigencore.complex_inertia_native", TRUE) = FALSE restores the
# embedding count.

#' @keywords internal
complex_inertia_native_enabled <- function() {
  !isFALSE(getOption("eigencore.complex_inertia_native", TRUE))
}

#' @keywords internal
complex_inertia_method_label <- function() {
  "dense_bunch_kaufman_hermitian"
}

# Attach the complex source to a realified dense operator (the embedding of
# `source`, a dense complex Hermitian matrix or a real symmetric one for a
# real metric of a complex pencil).
#' @keywords internal
complex_inertia_mark_embedding <- function(op, source) {
  if (!is.matrix(source) || !complex_inertia_native_enabled()) {
    return(op)
  }
  op$complex_embedding <- list(source = source, n = as.integer(nrow(source)))
  op
}

#' @keywords internal
complex_inertia_embedding_of <- function(x) {
  if (!inherits(x, "eigencore_operator")) {
    return(NULL)
  }
  emb <- x$complex_embedding %||% NULL
  if (is.null(emb) || !is.matrix(emb$source) ||
      !identical(as.integer(x$dim), 2L * c(emb$n, emb$n))) {
    return(NULL)
  }
  emb
}

# Hermitian part of a dense complex (or real) matrix, as the factorisation
# reads it (lower triangle; real diagonal).
#' @keywords internal
complex_inertia_norm1 <- function(x) {
  if (!length(x)) {
    return(0)
  }
  if (is.complex(x)) max(colSums(Mod(x)), 0) else norm(x, "O")
}

# Fill a dense inertia context with complex native factorisation fields.
# `multiplier` = 2 for contexts that count the real 2n embedding.
#' @keywords internal
complex_inertia_fill_dense <- function(ctx, A, B = NULL, multiplier = 1L) {
  if (!is.complex(A)) {
    A <- A + 0i
  }
  if (!is.null(B) && !is.complex(B) && !is.double(B)) {
    storage.mode(B) <- "double"
  }
  ctx$kind <- "dense"
  ctx$A <- A
  ctx$B <- B
  ctx$complex <- TRUE
  ctx$complex_native <- TRUE
  ctx$count_multiplier <- as.integer(multiplier)
  ctx$method <- complex_inertia_method_label()
  ctx
}

# Context over the real 2n embedding of a realified complex problem, or NULL
# when A (or B) is not a marked dense embedding.
#' @keywords internal
complex_inertia_embedding_context <- function(A, B = NULL) {
  if (!complex_inertia_native_enabled()) {
    return(NULL)
  }
  ea <- complex_inertia_embedding_of(A)
  if (is.null(ea)) {
    return(NULL)
  }
  eb <- NULL
  if (!is.null(B)) {
    eb <- complex_inertia_embedding_of(B)
    if (is.null(eb) || !identical(eb$n, ea$n)) {
      return(NULL)
    }
  }
  ctx <- new.env(parent = emptyenv())
  ctx$n <- 2L * ea$n
  ctx$generalized <- !is.null(eb)
  ctx$normA <- complex_inertia_norm1(ea$source)
  ctx$normB <- if (is.null(eb)) 1 else complex_inertia_norm1(eb$source)
  ctx$factor <- NULL
  ctx$factorizations <- 0L
  complex_inertia_fill_dense(ctx, ea$source, if (is.null(eb)) NULL else eb$source,
                             multiplier = 2L)
  if (!is.null(eb)) {
    # B must be positive definite (the real path checks the embedding with
    # chol()): the inertia of its LDL^H has no negative or zero pivot.
    t <- complex_inertia_at_dense(ctx, 0, list(sigma = 0), B_only = TRUE)
    if (!isTRUE(t$ok) || t$neg != 0 || t$zero != 0) {
      stop("B must be symmetric positive definite.", call. = FALSE)
    }
    ctx$factorizations <- 0L
  }
  ctx
}

# One native factorisation of A - s B (or of B alone, `B_only`), filling the
# raw tally `out` the way inertia_at() does for the real dense count.
#' @keywords internal
complex_inertia_at_dense <- function(ctx, s, out, B_only = FALSE) {
  A <- if (B_only) ctx$B else ctx$A
  B <- if (B_only) NULL else ctx$B
  t <- tryCatch(
    .Call("eigencore_dense_hermitian_inertia", A, as.double(s), B,
          PACKAGE = "eigencore"),
    error = function(e) e
  )
  if (inherits(t, "error")) {
    out$ok <- FALSE
    out$error <- conditionMessage(t)
    return(out)
  }
  mult <- ctx$count_multiplier %||% 1L
  out$ok <- TRUE
  out$neg <- mult * t[["neg"]]
  out$zero <- mult * t[["zero"]]
  out$pos <- mult * t[["pos"]]
  out$odd_embedding <- FALSE
  out$min_pivot <- t[["min_abs_pivot"]]
  out$max_pivot <- t[["max_abs_pivot"]]
  out$scale <- t[["norm1"]]
  out$growth <- max(t[["growth"]], t[["norm1"]])
  # Same Bunch-Kaufman bound as the real count, with p = 2n (the real
  # embedding's dimension; complex rounding constants are slightly larger).
  out$backward_bound <- 2 * max(nrow(A), 1) * .Machine$double.eps *
    (out$scale + out$growth)
  out$info <- t[["info"]]
  out
}

# Predicted cost of one native factorisation: zhetrf does n^3 / 3 complex
# multiply-adds, about 4 n^3 / 3 real flops.
#' @keywords internal
complex_inertia_factor_cost <- function(ctx) {
  m <- as.numeric(nrow(ctx$A))
  list(flops = 4 * m^3 / 3, factor_nnz = m^2 / 2, source = "exact")
}
