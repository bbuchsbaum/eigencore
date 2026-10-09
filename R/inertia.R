# Sylvester inertia: eigenvalue counting by symmetric indefinite factorisation
# (tranche 5, capability gap 4).
#
# For a real symmetric (or complex Hermitian) A and a symmetric positive
# definite B, A - sigma B = P L D L' P' is congruent to the block-diagonal D,
# so by Sylvester's law of inertia the number of negative / zero / positive
# eigenvalues of D equals the number of eigenvalues of the pencil (A, B) that
# lie below / at / above sigma. The factorisation used depends on the
# storage:
#
#   diagonal      direct comparison of a_i / b_i with sigma (exact)
#   tridiagonal   Sturm recurrence (unpivoted LDL'), backward stable for
#                 counting (LAPACK dstebz semantics)
#   dense         Bunch-Kaufman LDL' (LAPACK dsytrf), complex Hermitian via
#                 the real symmetric embedding [Re -Im; Im Re] (every
#                 eigenvalue doubled)
#   sparse        CHOLMOD simplicial LDL' with a fill-reducing (AMD) ordering
#                 through Matrix::Cholesky(LDL = TRUE, super = FALSE); the
#                 symbolic factor is reused (Matrix::update) across shifts.
#
# CHOLMOD's LDL' does not pivot, so a tiny pivot can make the computed
# inertia that of a matrix far from A - sigma B. Every count therefore carries
# a backward-error bound ||E|| <= c * eps * || |L| |D| |L'| || and a
# reliability verdict; an unreliable shift is perturbed (and the perturbation
# recorded) instead of being trusted.
#
# Internal API (for spectrum slicing, phase 2):
#   ctx <- inertia_context(A, B)          normalise once, cache norms/factors
#   inertia_at(ctx, s)                    one factorisation at s (raw tally)
#   inertia_count_nudged(ctx, s, dir)     count at s, nudged in direction
#                                         dir = +1 / -1 until reliable
#   inertia_interval_count(ctx, a, b)     eigenvalues in (a, b)
#   inertia_factor_cost(ctx)              predicted factor size / flops

.eigencore_inertia_state <- new.env(parent = emptyenv())

#' Count eigenvalues below, at and above a shift (Sylvester inertia)
#'
#' `eigen_count()` counts the eigenvalues of a real symmetric or complex
#' Hermitian matrix `A` (or of the symmetric-definite pencil `(A, B)`) that
#' lie below, at, and above `sigma`, without computing any eigenvalue. It
#' factors `A - sigma * B` as `P L D L' P'` and reads the counts off the
#' signs of the pivots of `D` (Sylvester's law of inertia).
#'
#' Factorisations by storage:
#' * diagonal `A` (and diagonal `B`): exact comparison;
#' * tridiagonal `A` with identity or diagonal `B`: Sturm sequence
#'   (backward stable for counting);
#' * dense real symmetric: LAPACK Bunch-Kaufman `dsytrf`; complex Hermitian:
#'   the same on the real symmetric embedding `[Re(A) -Im(A); Im(A) Re(A)]`;
#' * sparse (`dsCMatrix`, symmetric `dgCMatrix`, other sparse classes):
#'   CHOLMOD simplicial \eqn{LDL^T} with an AMD fill-reducing ordering via
#'   [Matrix::Cholesky()].
#'
#' **Reliability.** The computed inertia is exact for a nearby matrix
#' `A - sigma B + E`. `eigen_count()` bounds `||E||` by
#' `c * eps * || |L| |D| |L'| ||` (`backward_bound`) and declares a count
#' `reliable` only when every pivot is larger than that bound, the pivot
#' growth `|| |L| |D| |L'| || / ||A - sigma B||` is below `1 / sqrt(eps)` and
#' the factorisation completed. CHOLMOD's \eqn{LDL^T} does not pivot, so a
#' shift near an eigenvalue (or an unlucky ordering) can produce a tiny pivot.
#' With `perturb = TRUE` an unreliable shift is retried at `sigma - delta`
#' and `sigma + delta` for increasing `delta` (relative to the matrix scale);
#' the counts then refer to the half-lines below `sigma - delta` and above
#' `sigma + delta`, and `zero` counts the eigenvalues inside the band
#' `[sigma - delta, sigma + delta]` (`zero_band`). `reliable = FALSE` is
#' reported whenever no attempt was trustworthy; such counts must not be used
#' as a certificate. `reliable = TRUE` means no warning sign was found; it is
#' a numerical diagnosis, not an interval-arithmetic proof.
#'
#' @param A A real symmetric or complex Hermitian matrix: base `matrix`, a
#'   `Matrix` object (dense, sparse or diagonal), or an `eigencore_operator`
#'   with an explicit matrix source. The input must be symmetric (checked);
#'   only one triangle enters the factorisation.
#' @param sigma A single finite shift.
#' @param B Optional symmetric positive definite matrix for the pencil
#'   `A x = lambda B x` (checked). Counts are then of the generalized
#'   eigenvalues.
#' @param perturb Logical; retry an unreliable shift at nearby shifts.
#' @param pivot_tol Optional extra relative pivot floor: a count is
#'   unreliable when the smallest pivot is below `pivot_tol * scale`.
#'   Default `getOption("eigencore.inertia_pivot_tol", 0)` (the backward-error
#'   bound alone decides).
#' @return An object of class `eigencore_inertia`, a list with
#'   * `sigma`: the requested shift; `n`: the dimension;
#'   * `below`, `zero`, `above`: eigenvalue counts below `sigma` (below
#'     `sigma - delta` when perturbed), inside the zero band, and above;
#'   * `zero_band`: `c(sigma - delta, sigma + delta)` (`delta = 0` unless
#'     perturbed); `perturbed`, `perturbation` (`delta`);
#'   * `reliable`: whether the counts are trustworthy (see Details);
#'   * `method`: `"diagonal"`, `"tridiagonal_sturm"`,
#'     `"dense_bunch_kaufman"`, `"dense_bunch_kaufman_hermitian_embedding"` or
#'     `"sparse_cholmod_ldl"`; `generalized`: whether `B` was given;
#'   * `diagnostics`: a list with `scale` (an upper bound of
#'     `||A - sigma B||_2`), `min_pivot` and `max_pivot` (absolute),
#'     `min_pivot_relative` (`min_pivot / scale`), `pivot_growth`,
#'     `backward_bound`, `factor_nnz` (sparse), `attempts` (one row per
#'     factorisation: shift, counts, min pivot, reliable, error) and
#'     `seconds`;
#'   * `notes`: character vector of explanations.
#' @seealso [eig_partial()] uses these counts to certify that a computed set
#'   of extreme eigenvalues is complete (`certificate(fit)$target_completeness
#'   == "inertia_verified"`); see the "Certificates" vignette.
#' @examples
#' A <- diag(c(1, 2, 2, 5, 7))
#' eigen_count(A, 3)
#' S <- Matrix::sparseMatrix(i = c(1:6, 1:5), j = c(1:6, 2:6),
#'                           x = c(rep(2, 6), rep(-1, 5)), symmetric = TRUE)
#' eigen_count(S, 1)
#' sum(eigen(as.matrix(S), only.values = TRUE)$values < 1)
#' @export
eigen_count <- function(A, sigma, B = NULL, perturb = TRUE, pivot_tol = NULL) {
  sigma <- suppressWarnings(as.numeric(sigma))
  if (length(sigma) != 1L || !is.finite(sigma)) {
    stop("sigma must be a single finite number.", call. = FALSE)
  }
  if (!is.logical(perturb) || length(perturb) != 1L || is.na(perturb)) {
    stop("perturb must be TRUE or FALSE.", call. = FALSE)
  }
  pivot_tol <- inertia_pivot_tol(pivot_tol)
  started <- proc.time()[["elapsed"]]
  ctx <- inertia_context(A, B)
  out <- inertia_count_at(ctx, sigma, perturb = perturb, pivot_tol = pivot_tol)
  out$diagnostics$seconds <- proc.time()[["elapsed"]] - started
  out
}

#' @keywords internal
inertia_pivot_tol <- function(pivot_tol = NULL) {
  pivot_tol <- pivot_tol %||% getOption("eigencore.inertia_pivot_tol", 0)
  pivot_tol <- suppressWarnings(as.numeric(pivot_tol))
  if (length(pivot_tol) != 1L || !is.finite(pivot_tol) || pivot_tol < 0) {
    stop("pivot_tol must be a single non-negative number.", call. = FALSE)
  }
  pivot_tol
}

#' @export
print.eigencore_inertia <- function(x, ...) {
  cat("<eigencore_inertia>", x$method,
      if (isTRUE(x$generalized)) "(generalized, B SPD)" else "", "\n")
  if (isTRUE(x$perturbed)) {
    cat("  sigma:", format(x$sigma, digits = 10), " zero band: sigma -/+",
        format(x$perturbation, digits = 3), "\n")
  } else {
    cat("  sigma:", format(x$sigma, digits = 10), "\n")
  }
  cat("  below:", x$below, " zero:", x$zero, " above:", x$above,
      " (n =", x$n, ")\n")
  cat("  reliable:", x$reliable,
      " min pivot / scale:", format(x$diagnostics$min_pivot_relative, digits = 3),
      " pivot growth:", format(x$diagnostics$pivot_growth, digits = 3), "\n")
  if (length(x$notes)) {
    cat(paste0("  note: ", x$notes, collapse = "\n"), "\n")
  }
  invisible(x)
}

# ---------------------------------------------------------------------------
# Input normalisation
# ---------------------------------------------------------------------------

#' @keywords internal
inertia_dense_limit <- function() {
  as.integer(getOption("eigencore.inertia_dense_limit", 20000L))
}

# Extract an explicit matrix from a matrix-like input or an eigencore
# operator; NULL when there is none (matrix-free operator).
#' @keywords internal
inertia_matrix_of <- function(x) {
  if (inherits(x, "eigencore_operator")) {
    src <- source_or_null(x)
    if (is.null(src)) {
      src <- x$metadata$matrix %||% NULL
    }
    return(src)
  }
  x
}

#' @keywords internal
inertia_normalize <- function(x, arg = "A") {
  # Operators were validated at construction (finite entries, C8) and carry
  # a Hermitian structure flag only after a symmetry test.
  trusted <- inherits(x, "eigencore_operator") &&
    identical(x$structure$kind, "hermitian")
  x <- inertia_matrix_of(x)
  if (is.null(x)) {
    stop(arg, " has no explicit matrix source; inertia counting needs a ",
         "dense or sparse matrix.", call. = FALSE)
  }
  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }
  d <- dim(x)
  if (length(d) != 2L || d[[1L]] != d[[2L]]) {
    stop(arg, " must be a square matrix.", call. = FALSE)
  }
  n <- as.integer(d[[1L]])
  if (inherits(x, "diagonalMatrix")) {
    v <- as.numeric(Matrix::diag(x))
    if (any(!is.finite(v))) stop(arg, " has non-finite entries.", call. = FALSE)
    return(list(kind = "diagonal", n = n, d = v))
  }
  if (inherits(x, "sparseMatrix")) {
    if (inherits(x, "dsparseMatrix") || inherits(x, "lsparseMatrix") ||
        inherits(x, "nsparseMatrix") || inherits(x, "isparseMatrix")) {
      x <- methods::as(x, "dMatrix")
    } else {
      stop(arg, ": complex sparse matrices are not supported.", call. = FALSE)
    }
    x <- methods::as(x, "CsparseMatrix")
    xs <- methods::slot(x, "x")
    if (!trusted && length(xs) && any(!is.finite(xs))) {
      stop(arg, " has non-finite entries.", call. = FALSE)
    }
    if (!inherits(x, "symmetricMatrix")) {
      if (!trusted &&
          !isTRUE(Matrix::isSymmetric(x, tol = 100 * .Machine$double.eps))) {
        stop(arg, " must be symmetric.", call. = FALSE)
      }
      x <- Matrix::forceSymmetric(x, uplo = "U")
    }
    x <- methods::as(x, "CsparseMatrix")
    band <- inertia_sparse_bandwidth(x)
    if (band == 0L) {
      return(list(kind = "diagonal", n = n, d = as.numeric(Matrix::diag(x))))
    }
    if (band == 1L) {
      tri <- inertia_sparse_tridiagonal(x)
      return(list(kind = "tridiagonal", n = n, d = tri$d, e = tri$e, sparse = x))
    }
    return(list(kind = "sparse", n = n, matrix = x))
  }
  if (inherits(x, "Matrix")) {
    x <- as.matrix(x)
  }
  if (!is.matrix(x)) {
    stop(arg, " must be a matrix, a Matrix object or an eigencore operator ",
         "with a matrix source.", call. = FALSE)
  }
  if (is.complex(x)) {
    if (any(!is.finite(Re(x))) || any(!is.finite(Im(x)))) {
      stop(arg, " has non-finite entries.", call. = FALSE)
    }
    scale <- max(Mod(x), 0)
    if (max(Mod(x - Conj(t(x))), 0) > 100 * .Machine$double.eps * max(scale, 1e-300) * n) {
      stop(arg, " must be Hermitian.", call. = FALSE)
    }
    return(list(kind = "dense", n = n, matrix = x, complex = TRUE))
  }
  if (!is.numeric(x) && !is.logical(x)) {
    stop(arg, " must be numeric.", call. = FALSE)
  }
  if (!is.double(x)) {
    storage.mode(x) <- "double"
  }
  if (!trusted) {
    check <- .Call("eigencore_dense_finite_symmetric", x,
                   100 * .Machine$double.eps, PACKAGE = "eigencore")
    if (!isTRUE(check[[1L]])) {
      stop(arg, " has non-finite entries.", call. = FALSE)
    }
    if (!isTRUE(check[[2L]])) {
      stop(arg, " must be symmetric.", call. = FALSE)
    }
  }
  list(kind = "dense", n = n, matrix = x, complex = FALSE)
}

#' @keywords internal
inertia_sparse_bandwidth <- function(x) {
  i <- methods::slot(x, "i")
  if (!length(i)) {
    return(0L)
  }
  j <- rep.int(seq_len(ncol(x)) - 1L, diff(methods::slot(x, "p")))
  as.integer(max(abs(i - j)))
}

# Diagonal and off-diagonal of a symmetric CsparseMatrix of bandwidth <= 1.
#' @keywords internal
inertia_sparse_tridiagonal <- function(x) {
  n <- nrow(x)
  i <- methods::slot(x, "i") + 1L
  j <- rep.int(seq_len(n), diff(methods::slot(x, "p")))
  v <- methods::slot(x, "x")
  d <- numeric(n)
  e <- numeric(max(n - 1L, 0L))
  on <- i == j
  d[i[on]] <- v[on]
  off <- abs(i - j) == 1L
  e[pmin(i[off], j[off])] <- v[off]
  list(d = d, e = e)
}

#' @keywords internal
inertia_as_dense <- function(part) {
  switch(
    part$kind,
    dense = part$matrix,
    diagonal = diag(part$d, nrow = part$n),
    tridiagonal = as.matrix(Matrix::forceSymmetric(part$sparse)),
    sparse = as.matrix(part$matrix)
  )
}

#' @keywords internal
inertia_as_sparse <- function(part) {
  switch(
    part$kind,
    diagonal = methods::as(Matrix::forceSymmetric(methods::as(
      Matrix::Diagonal(x = part$d), "CsparseMatrix")), "CsparseMatrix"),
    tridiagonal = part$sparse,
    sparse = part$matrix
  )
}

#' @keywords internal
inertia_norm1 <- function(part) {
  switch(
    part$kind,
    diagonal = max(abs(part$d), 0),
    tridiagonal = {
      e <- abs(part$e)
      max(abs(part$d) + c(0, e) + c(e, 0), 0)
    },
    dense = if (isTRUE(part$complex)) {
      max(colSums(Mod(part$matrix)), 0)
    } else {
      max(colSums(abs(part$matrix)), 0)
    },
    sparse = max(Matrix::colSums(abs(part$matrix)), 0)
  )
}

#' @keywords internal
inertia_check_spd <- function(part) {
  ok <- switch(
    part$kind,
    diagonal = all(part$d > 0),
    tridiagonal = ,
    sparse = isTRUE(tryCatch({
      Matrix::Cholesky(inertia_as_sparse(part), LDL = FALSE, super = NA,
                       perm = TRUE)
      TRUE
    }, error = function(e) FALSE, warning = function(w) FALSE)),
    dense = isTRUE(tryCatch({
      chol(if (isTRUE(part$complex)) inertia_hermitian_embedding(part$matrix) else part$matrix)
      TRUE
    }, error = function(e) FALSE))
  )
  if (!isTRUE(ok)) {
    stop("B must be symmetric positive definite.", call. = FALSE)
  }
  invisible(TRUE)
}

#' @keywords internal
inertia_hermitian_embedding <- function(x) {
  re <- Re(x)
  im <- Im(x)
  out <- rbind(cbind(re, -im), cbind(im, re))
  dimnames(out) <- NULL
  out
}

# Normalise A (and B) once for repeated inertia counts.
#
# Returns an environment holding the storage kind, the normalised matrices,
# the 1-norms (upper bounds of the 2-norms) and, for the sparse route, the
# cached CHOLMOD factor whose symbolic analysis later shifts reuse.
#' @keywords internal
inertia_context <- function(A, B = NULL) {
  a <- inertia_normalize(A, "A")
  b <- if (is.null(B)) NULL else inertia_normalize(B, "B")
  if (!is.null(b) && b$n != a$n) {
    stop("B must have the dimensions of A.", call. = FALSE)
  }
  if (!is.null(b)) {
    inertia_check_spd(b)
  }
  kind <- if (is.null(b)) {
    a$kind
  } else if (a$kind %in% c("diagonal", "tridiagonal") && identical(b$kind, "diagonal")) {
    a$kind
  } else if (identical(a$kind, "dense") || identical(b$kind, "dense")) {
    "dense"
  } else {
    "sparse"
  }
  ctx <- new.env(parent = emptyenv())
  ctx$n <- a$n
  ctx$kind <- kind
  ctx$generalized <- !is.null(b)
  ctx$complex <- FALSE
  ctx$normA <- inertia_norm1(a)
  ctx$normB <- if (is.null(b)) 1 else inertia_norm1(b)
  ctx$factor <- NULL
  ctx$factorizations <- 0L
  if (identical(kind, "diagonal")) {
    ctx$d <- a$d
    ctx$b <- if (is.null(b)) NULL else b$d
    ctx$method <- "diagonal"
  } else if (identical(kind, "tridiagonal")) {
    if (identical(a$kind, "diagonal")) {
      ctx$d <- a$d
      ctx$e <- numeric(max(a$n - 1L, 0L))
    } else {
      ctx$d <- a$d
      ctx$e <- a$e
    }
    ctx$b <- if (is.null(b)) NULL else b$d
    ctx$method <- "tridiagonal_sturm"
  } else if (identical(kind, "dense")) {
    if (a$n > inertia_dense_limit() &&
        (!identical(a$kind, "dense") || (!is.null(b) && !identical(b$kind, "dense")))) {
      stop("eigen_count(): refusing to densify a sparse matrix with n = ", a$n,
           " (> getOption(\"eigencore.inertia_dense_limit\")).", call. = FALSE)
    }
    Ad <- inertia_as_dense(a)
    Bd <- if (is.null(b)) NULL else inertia_as_dense(b)
    cplx <- is.complex(Ad) || is.complex(Bd)
    if (cplx) {
      Ad <- inertia_hermitian_embedding(Ad + 0i)
      Bd <- if (is.null(Bd)) NULL else inertia_hermitian_embedding(Bd + 0i)
    }
    if (!is.double(Ad)) storage.mode(Ad) <- "double"
    if (!is.null(Bd) && !is.double(Bd)) storage.mode(Bd) <- "double"
    ctx$A <- Ad
    ctx$B <- Bd
    ctx$complex <- cplx
    ctx$method <- if (cplx) "dense_bunch_kaufman_hermitian_embedding" else "dense_bunch_kaufman"
  } else {
    ctx$A <- inertia_as_sparse(a)
    ctx$B <- if (is.null(b)) NULL else inertia_as_sparse(b)
    ctx$method <- "sparse_cholmod_ldl"
  }
  ctx
}

# ---------------------------------------------------------------------------
# One factorisation
# ---------------------------------------------------------------------------

#' @keywords internal
inertia_tally_template <- function(sigma, n) {
  list(
    sigma = sigma, ok = FALSE, error = NA_character_,
    neg = NA_real_, zero = NA_real_, pos = NA_real_,
    min_pivot = NA_real_, max_pivot = NA_real_, growth = NA_real_,
    scale = NA_real_, backward_bound = NA_real_, factor_nnz = NA_real_
  )
}

# Factor A - s B once and return the raw inertia tally.
#' @keywords internal
inertia_at <- function(ctx, s) {
  n <- ctx$n
  out <- inertia_tally_template(s, n)
  eps <- .Machine$double.eps
  scale <- ctx$normA + abs(s) * ctx$normB
  out$scale <- scale
  ctx$factorizations <- ctx$factorizations + 1L
  if (identical(ctx$kind, "diagonal")) {
    b <- ctx$b %||% rep(1, n)
    m <- ctx$d - s * b
    lam <- ctx$d / b
    out$ok <- TRUE
    out$neg <- sum(lam < s)
    out$zero <- sum(lam == s)
    out$pos <- n - out$neg - out$zero
    out$min_pivot <- if (n) min(abs(m)) else Inf
    out$max_pivot <- if (n) max(abs(m)) else 0
    out$growth <- out$max_pivot
    out$backward_bound <- 2 * eps * scale
    out$exact <- TRUE
    return(out)
  }
  if (identical(ctx$kind, "tridiagonal")) {
    t <- .Call("eigencore_tridiagonal_inertia", as.double(ctx$d), as.double(ctx$e),
               as.double(s), if (is.null(ctx$b)) NULL else as.double(ctx$b),
               PACKAGE = "eigencore")
    out$ok <- TRUE
    out$neg <- t[["neg"]]
    out$zero <- t[["zero"]]
    out$pos <- t[["pos"]]
    out$min_pivot <- t[["min_abs_pivot"]]
    out$max_pivot <- t[["max_abs_pivot"]]
    out$scale <- max(t[["norm1"]], scale)
    out$growth <- out$scale
    # Sturm counts are exact for a componentwise relative perturbation of
    # the entries of O(eps): ||E|| <= ~4 eps ||T - s B||.
    out$backward_bound <- 4 * eps * out$scale
    return(out)
  }
  if (identical(ctx$kind, "dense")) {
    t <- tryCatch(
      .Call("eigencore_dense_symmetric_inertia", ctx$A, as.double(s), ctx$B,
            PACKAGE = "eigencore"),
      error = function(e) e
    )
    if (inherits(t, "error")) {
      out$error <- conditionMessage(t)
      return(out)
    }
    out$ok <- TRUE
    div <- if (isTRUE(ctx$complex)) 2 else 1
    out$neg <- t[["neg"]] / div
    out$zero <- t[["zero"]] / div
    out$pos <- t[["pos"]] / div
    out$odd_embedding <- isTRUE(ctx$complex) &&
      any(c(t[["neg"]], t[["zero"]], t[["pos"]]) %% 2 != 0)
    out$min_pivot <- t[["min_abs_pivot"]]
    out$max_pivot <- t[["max_abs_pivot"]]
    out$scale <- t[["norm1"]]
    out$growth <- max(t[["growth"]], t[["norm1"]])
    # Bunch-Kaufman: ||E|| <= p(n) eps (||M|| + || |L||D||L'| ||) with p(n)
    # a modest polynomial; n is the usual practical constant.
    out$backward_bound <- max(nrow(ctx$A), 1) * eps * (out$scale + out$growth)
    out$info <- t[["info"]]
    return(out)
  }
  # Sparse CHOLMOD simplicial LDL' (no pivoting).
  F <- tryCatch(
    suppressWarnings(inertia_sparse_factor(ctx, s)),
    error = function(e) e
  )
  if (inherits(F, "error")) {
    out$error <- conditionMessage(F)
    out$zero <- NA_real_
    return(out)
  }
  t <- tryCatch(
    .Call("eigencore_simplicial_ldl_diagnostics", methods::slot(F, "p"),
          methods::slot(F, "i"), methods::slot(F, "x"), methods::slot(F, "nz"),
          PACKAGE = "eigencore"),
    error = function(e) e
  )
  if (inherits(t, "error")) {
    out$error <- conditionMessage(t)
    return(out)
  }
  out$ok <- TRUE
  out$neg <- t[["neg"]]
  out$zero <- t[["zero"]]
  out$pos <- t[["pos"]]
  out$min_pivot <- t[["min_abs_pivot"]]
  out$max_pivot <- t[["max_abs_pivot"]]
  out$growth <- max(t[["growth"]], scale)
  nz <- methods::slot(F, "nz")
  out$factor_nnz <- sum(as.numeric(nz))
  # Componentwise backward error of unpivoted LDL': |E| <= gamma_c |L||D||L'|
  # with c the longest elimination inner product (max column count of L).
  out$backward_bound <- (max(nz, 1) + 3) * eps * out$growth
  out$factor <- F
  out
}

# Numeric CHOLMOD LDL' of A - s B; the symbolic analysis of the first
# successful factorisation is reused for later shifts (identity B only:
# A - s B with a general B can change pattern through cancellation).
#' @keywords internal
inertia_sparse_factor <- function(ctx, s) {
  if (is.null(ctx$B)) {
    F <- if (!is.null(ctx$factor)) {
      Matrix::update(ctx$factor, ctx$A, mult = -s)
    } else {
      Matrix::Cholesky(ctx$A, LDL = TRUE, super = FALSE, perm = TRUE, Imult = -s)
    }
    ctx$factor <- F
    return(F)
  }
  M <- Matrix::forceSymmetric(methods::as(ctx$A - s * ctx$B, "CsparseMatrix"), uplo = "U")
  Matrix::Cholesky(M, LDL = TRUE, super = FALSE, perm = TRUE)
}

#' @keywords internal
inertia_tally_reliable <- function(t, pivot_tol = 0) {
  if (!isTRUE(t$ok)) {
    return(FALSE)
  }
  vals <- c(t$neg, t$zero, t$pos, t$min_pivot, t$growth, t$scale)
  if (any(!is.finite(vals))) {
    return(FALSE)
  }
  if (isTRUE(t$exact)) {
    return(TRUE)
  }
  if (t$zero > 0 || isTRUE(t$odd_embedding)) {
    return(FALSE)
  }
  growth_ratio <- if (t$scale > 0) t$growth / t$scale else 1
  floor <- max(t$backward_bound, pivot_tol * t$scale)
  t$min_pivot > floor && growth_ratio <= 1 / sqrt(.Machine$double.eps)
}

#' @keywords internal
inertia_attempt_row <- function(t, reliable) {
  data.frame(
    sigma = t$sigma,
    below = t$neg, zero = t$zero, above = t$pos,
    min_pivot = t$min_pivot,
    pivot_growth = if (is.finite(t$scale) && t$scale > 0) t$growth / t$scale else NA_real_,
    backward_bound = t$backward_bound,
    reliable = reliable,
    error = t$error,
    stringsAsFactors = FALSE
  )
}

# Relative perturbation sizes tried for an unreliable shift.
#' @keywords internal
inertia_perturbation_steps <- function() {
  c(1e-12, 1e-10, 1e-8, 1e-6)
}

#' @keywords internal
inertia_delta <- function(ctx, s, rel) {
  rel * max(ctx$normA + abs(s) * ctx$normB, .Machine$double.xmin)
}

# Count at sigma for eigen_count(): the shift itself when reliable, otherwise
# a symmetric pair sigma -/+ delta whose two counts bracket a zero band.
#' @keywords internal
inertia_count_at <- function(ctx, sigma, perturb = TRUE, pivot_tol = 0) {
  n <- ctx$n
  t0 <- inertia_at(ctx, sigma)
  r0 <- inertia_tally_reliable(t0, pivot_tol)
  attempts <- inertia_attempt_row(t0, r0)
  chosen <- list(below = t0$neg, zero = t0$zero, above = t0$pos,
                 delta = 0, reliable = r0, lo = t0, hi = t0)
  notes <- character()
  if (!r0 && isTRUE(t0$exact)) {
    r0 <- TRUE
  }
  if (!r0 && isTRUE(perturb)) {
    for (rel in inertia_perturbation_steps()) {
      delta <- inertia_delta(ctx, sigma, rel)
      lo <- inertia_at(ctx, sigma - delta)
      rlo <- inertia_tally_reliable(lo, pivot_tol)
      hi <- inertia_at(ctx, sigma + delta)
      rhi <- inertia_tally_reliable(hi, pivot_tol)
      attempts <- rbind(attempts, inertia_attempt_row(lo, rlo),
                        inertia_attempt_row(hi, rhi))
      if (rlo && rhi && lo$neg <= hi$neg) {
        chosen <- list(below = lo$neg, zero = hi$neg + hi$zero - lo$neg,
                       above = hi$pos, delta = delta, reliable = TRUE,
                       lo = lo, hi = hi)
        notes <- c(notes, paste0(
          "sigma was within the numerical zero band of the factorisation; counts ",
          "use sigma -/+ ", format(delta, digits = 3),
          " and 'zero' counts the eigenvalues inside that band"
        ))
        break
      }
    }
  }
  if (!isTRUE(chosen$reliable) && isTRUE(perturb) &&
      identical(ctx$kind, "sparse") && n <= inertia_dense_fallback_limit()) {
    # Unpivoted sparse LDL' broke down (typically zero or tiny diagonal
    # entries, e.g. saddle-point structure); Bunch-Kaufman pivots.
    dense_ctx <- inertia_dense_fallback_context(ctx)
    fallback <- inertia_count_at(dense_ctx, sigma, perturb = perturb,
                                 pivot_tol = pivot_tol)
    fallback$method <- paste0(dense_ctx$method, "_fallback")
    fallback$notes <- c(
      "sparse LDL' (no pivoting) was unreliable at this shift; counted with dense Bunch-Kaufman instead",
      fallback$notes
    )
    fallback$diagnostics$attempts <- rbind(attempts, fallback$diagnostics$attempts)
    return(fallback)
  }
  if (!isTRUE(chosen$reliable)) {
    if (!isTRUE(t0$ok)) {
      notes <- c(notes, paste0("factorisation failed: ", t0$error))
    }
    notes <- c(notes, "no reliable factorisation was found; counts are not certified")
    if (!isTRUE(t0$ok)) {
      chosen$below <- NA_real_
      chosen$zero <- NA_real_
      chosen$above <- NA_real_
    }
  }
  ref <- if (isTRUE(chosen$reliable) && chosen$delta > 0) chosen$hi else t0
  min_piv <- if (chosen$delta > 0) min(chosen$lo$min_pivot, chosen$hi$min_pivot) else t0$min_pivot
  growth <- if (chosen$delta > 0) max(chosen$lo$growth, chosen$hi$growth) else t0$growth
  scale <- ref$scale
  out <- list(
    sigma = sigma,
    n = n,
    below = as.integer(chosen$below),
    zero = as.integer(chosen$zero),
    above = as.integer(chosen$above),
    zero_band = c(sigma - chosen$delta, sigma + chosen$delta),
    perturbed = chosen$delta > 0,
    perturbation = chosen$delta,
    reliable = isTRUE(chosen$reliable),
    method = ctx$method,
    generalized = isTRUE(ctx$generalized),
    diagnostics = list(
      scale = scale,
      min_pivot = min_piv,
      max_pivot = ref$max_pivot,
      min_pivot_relative = if (is.finite(scale) && scale > 0) min_piv / scale else NA_real_,
      pivot_growth = if (is.finite(scale) && scale > 0) growth / scale else NA_real_,
      backward_bound = ref$backward_bound,
      factor_nnz = ref$factor_nnz,
      attempts = attempts,
      seconds = NA_real_
    ),
    notes = notes
  )
  class(out) <- "eigencore_inertia"
  out
}

#' @keywords internal
inertia_dense_fallback_limit <- function() {
  as.integer(getOption("eigencore.inertia_dense_fallback_limit", 4000L))
}

# Dense Bunch-Kaufman context over the same (A, B) as a sparse context.
#' @keywords internal
inertia_dense_fallback_context <- function(ctx) {
  if (!is.null(ctx$dense_fallback)) {
    return(ctx$dense_fallback)
  }
  d <- new.env(parent = emptyenv())
  d$n <- ctx$n
  d$kind <- "dense"
  d$generalized <- ctx$generalized
  d$complex <- FALSE
  d$normA <- ctx$normA
  d$normB <- ctx$normB
  d$factor <- NULL
  d$factorizations <- 0L
  d$A <- as.matrix(ctx$A)
  d$B <- if (is.null(ctx$B)) NULL else as.matrix(ctx$B)
  storage.mode(d$A) <- "double"
  if (!is.null(d$B)) storage.mode(d$B) <- "double"
  d$method <- "dense_bunch_kaufman"
  ctx$dense_fallback <- d
  d
}

# Count of eigenvalues below s (and above s), nudging s in `direction`
# (+1 / -1) by growing relative steps until the factorisation is reliable.
# Returns list(below, zero, above, s, reliable, tally, attempts).
#' @keywords internal
inertia_count_nudged <- function(ctx, s, direction = 1, pivot_tol = 0,
                                 steps = c(0, inertia_perturbation_steps())) {
  attempts <- NULL
  last <- NULL
  for (rel in steps) {
    si <- s + direction * inertia_delta(ctx, s, rel)
    t <- inertia_at(ctx, si)
    r <- inertia_tally_reliable(t, pivot_tol)
    attempts <- rbind(attempts, inertia_attempt_row(t, r))
    last <- t
    if (r) {
      return(list(below = t$neg, zero = t$zero, above = t$pos, s = si,
                  reliable = TRUE, tally = t, attempts = attempts))
    }
  }
  if (identical(ctx$kind, "sparse") && ctx$n <= inertia_dense_fallback_limit()) {
    fallback <- inertia_count_nudged(inertia_dense_fallback_context(ctx), s,
                                     direction = direction, pivot_tol = pivot_tol,
                                     steps = steps)
    fallback$attempts <- rbind(attempts, fallback$attempts)
    fallback$dense_fallback <- TRUE
    return(fallback)
  }
  list(below = last$neg, zero = last$zero, above = last$pos, s = last$sigma,
       reliable = FALSE, tally = last, attempts = attempts)
}

# Eigenvalues in the open interval (a, b), a < b; the end points are nudged
# outwards when unreliable, so the count covers (a', b') with a' <= a, b' >= b.
#' @keywords internal
inertia_interval_count <- function(ctx, a, b, pivot_tol = 0) {
  lo <- inertia_count_nudged(ctx, a, direction = -1, pivot_tol = pivot_tol)
  hi <- inertia_count_nudged(ctx, b, direction = 1, pivot_tol = pivot_tol)
  list(
    count = hi$below - lo$below - lo$zero,
    interval = c(lo$s, hi$s),
    reliable = lo$reliable && hi$reliable,
    lower = lo,
    upper = hi
  )
}

# ---------------------------------------------------------------------------
# Cost model
# ---------------------------------------------------------------------------

#' @keywords internal
cholmod_bridge_available <- function() {
  cached <- .eigencore_inertia_state$bridge
  if (!is.null(cached)) {
    return(cached)
  }
  ok <- isTRUE(tryCatch({
    abi <- .Call("eigencore_cholmod_abi", PACKAGE = "eigencore")
    v <- Matrix::Matrix.Version()
    ss <- as.integer(strsplit(as.character(v$suitesparse), ".", fixed = TRUE)[[1L]])
    identical(as.integer(as.character(v$abi)), abi[[1L]]) &&
      identical(ss[1:3], abi[2:4])
  }, error = function(e) FALSE))
  .eigencore_inertia_state$bridge <- ok
  ok
}

# Predicted cost of one LDL' factorisation: list(flops, factor_nnz, source).
# Sparse uses CHOLMOD's symbolic analysis (exact fill / flop count of the AMD
# ordering) when the Matrix C API matches the compiled ABI; otherwise NA.
#' @keywords internal
inertia_factor_cost <- function(ctx) {
  n <- as.numeric(ctx$n)
  if (!is.null(ctx$cost)) {
    return(ctx$cost)
  }
  cost <- switch(
    ctx$kind,
    diagonal = list(flops = n, factor_nnz = n, source = "exact"),
    tridiagonal = list(flops = 5 * n, factor_nnz = 2 * n, source = "exact"),
    dense = {
      m <- if (isTRUE(ctx$complex)) 2 * n else n
      list(flops = m^3 / 3, factor_nnz = m^2 / 2, source = "exact")
    },
    sparse = {
      a <- if (cholmod_bridge_available()) {
        M <- if (is.null(ctx$B)) ctx$A else
          Matrix::forceSymmetric(methods::as(ctx$A + ctx$B, "CsparseMatrix"), uplo = "U")
        tryCatch(.Call("eigencore_cholmod_analyze", methods::as(M, "CsparseMatrix"),
                       PACKAGE = "eigencore"),
                 error = function(e) NULL)
      }
      if (is.null(a) || !is.finite(a[["flops"]])) {
        list(flops = NA_real_, factor_nnz = NA_real_, source = "unavailable")
      } else {
        list(flops = a[["flops"]], factor_nnz = a[["lnz"]], source = "cholmod_analyze")
      }
    }
  )
  ctx$cost <- cost
  cost
}
