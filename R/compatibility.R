#' RSpectra-compatible eigen shim.
#'
#' Mirrors `RSpectra::eigs()`. `which` uses ARPACK codes (`"LM"`, `"SM"`,
#' `"LR"`, `"SR"`, `"LI"`, `"SI"`). Unlike ARPACK, `"LI"`/`"SI"` rank by the
#' signed imaginary part (see [largest_imaginary()]), not by its magnitude. A
#' non-`NULL` real `sigma` with `which = "LM"` requests the eigenvalues
#' nearest `sigma`, computed by shift-invert Krylov-Schur Arnoldi on a
#' factorised `A - sigma I` (dense QR or sparse LU).
#'
#' Like `RSpectra::eigs()`, only right eigenvectors are computed by default;
#' `left = TRUE` also computes and certifies left eigenvectors and their
#' biorthogonality with the right ones (the [eig_partial()] default).
#'
#' @param A Matrix, eigencore operator, or a function `f(x, args)` returning
#'   `A %*% x` (then `n` is required).
#' @param k Number of eigenpairs to compute.
#' @param which RSpectra-style target selector.
#' @param sigma Optional real shift; eigenvalues nearest `sigma` are returned.
#' @param opts RSpectra options list. `tol`, `ncv` (Krylov subspace size,
#'   passed as `auto(max_subspace = ncv)` so the planner still chooses the
#'   route), `maxitr` (passed as the iteration limit `maxit`) and `retvec` are
#'   honoured; other keys raise a warning.
#' @param ... Additional arguments passed to [eig_partial()].
#' @param n Dimension of the operator when `A` is a function.
#' @param args Extra argument passed to a function `A`.
#' @param left Whether to compute and certify left eigenvectors
#'   (`left_vectors = "auto"` in [eig_partial()]). `FALSE` (default) skips the
#'   adjoint solve and its certificate entirely.
#' @return A list compatible with `RSpectra::eigs()`, including `values`,
#'   `vectors`, convergence counts, operation counts, certificate diagnostics,
#'   and, with `left = TRUE`, `left_vectors`, `right_vectors`,
#'   `left_certificate` and `biorthogonality`.
#' @examples
#' A <- diag(c(5, 4, 3, 2, 1))
#' A[1, 2] <- 0.5
#' res <- eigs(A, k = 2, which = "LM")
#' res$values
eigs <- function(A, k, which = "LM", sigma = NULL, opts = list(), ...,
                 n = NULL, args = NULL, left = FALSE) {
  if (!is.logical(left) || length(left) != 1L || is.na(left)) {
    stop("left must be TRUE or FALSE.", call. = FALSE)
  }
  A <- compat_function_operator(A, n = n, args = args, structure = general())
  controls <- compat_eigen_opts(opts, allow_initvec = FALSE, dots = list(...))
  retvec <- !isFALSE(opts$retvec)
  target <- compat_eigen_target(which, k = k, sigma = sigma, symmetric = FALSE)
  fit <- do.call(eig_partial, c(
    list(A, k = k, target = target,
         left_vectors = if (isTRUE(left)) "auto" else "none"),
    controls,
    list(...)
  ))
  compat_warn_nconv(fit$nconv, k)
  out <- list(
    values = fit$values,
    vectors = if (retvec) fit$vectors,
    nconv = fit$nconv,
    niter = fit$iterations,
    nops = fit$matvecs,
    certificate = fit$certificate,
    diagnostics = diagnostics(fit)
  )
  if (isTRUE(left)) {
    out$left_vectors <- if (retvec) fit$left_vectors
    out$right_vectors <- if (retvec) right_vectors(fit)
    out$left_certificate <- fit$left_certificate
    out$biorthogonality <- fit$biorthogonality
  }
  out
}

#' RSpectra-compatible symmetric eigen shim.
#'
#' Mirrors `RSpectra::eigs_sym()`: only the `lower` (or upper) triangle of a
#' dense or `Matrix` input is read, values are returned in decreasing order,
#' and `sigma` requests the eigenvalues nearest the shift.
#'
#' @param A Matrix, eigencore operator, or a function `f(x, args)` returning
#'   `A %*% x` (then `n` is required).
#' @param k Number of eigenpairs to compute.
#' @param which RSpectra-style target selector (`"LM"`, `"SM"`, `"LA"`,
#'   `"SA"`, `"BE"`).
#' @param sigma Optional shift; eigenvalues nearest `sigma` are returned.
#' @param opts RSpectra options list. `tol`, `ncv` (Krylov subspace size,
#'   passed as `auto(max_subspace = ncv)`, or `lanczos(max_subspace = ncv)`
#'   when `initvec` forces a Lanczos route), `maxitr` (passed as the
#'   iteration limit `maxit`), `retvec` and `initvec` are honoured; other
#'   keys raise a warning.
#' @param lower Whether to read the lower (`TRUE`) or upper (`FALSE`) triangle
#'   of a matrix input.
#' @param ... Additional arguments passed to [eig_partial()].
#' @param n Dimension of the operator when `A` is a function.
#' @param args Extra argument passed to a function `A`.
#' @return A list compatible with `RSpectra::eigs_sym()`, including `values`,
#'   `vectors`, convergence counts, operation counts, certificate diagnostics,
#'   and eigencore diagnostics.
#' @examples
#' A <- diag(c(5, 4, 3, 2, 1))
#' res <- eigs_sym(A, k = 2, which = "LA")
#' res$values
eigs_sym <- function(A, k, which = "LM", sigma = NULL, opts = list(),
                     lower = TRUE, ..., n = NULL, args = NULL) {
  A <- compat_function_operator(A, n = n, args = args, structure = hermitian())
  A <- compat_symmetric_from_triangle(A, lower = lower)
  controls <- compat_eigen_opts(opts, allow_initvec = is.null(sigma),
                                dots = list(...))
  retvec <- !isFALSE(opts$retvec)
  target <- compat_eigen_target(which, k = k, sigma = sigma, symmetric = TRUE)
  P <- eigen_problem(A, structure = hermitian(), target = target)
  fit <- do.call(solve, c(list(P, k = k), controls, list(...)))
  compat_warn_nconv(fit$nconv, k)
  values <- fit$values
  vectors <- fit$vectors
  ord <- order(values, decreasing = TRUE)
  values <- values[ord]
  vectors <- if (retvec && !is.null(vectors)) vectors[, ord, drop = FALSE]
  # Keep the per-pair certificate fields aligned with the re-sorted values
  # (they were left in solver order; found by the oracle sweep).
  cert <- fit$certificate
  if (!is.null(cert)) {
    for (field in c("residuals", "backward_error", "converged", "scale")) {
      x <- cert[[field]]
      if (is.atomic(x) && length(x) == length(ord)) {
        cert[[field]] <- x[ord]
      }
    }
    if (length(cert$failed_indices) && length(cert$converged) == length(ord)) {
      cert$failed_indices <- which(!cert$converged)
    }
  }
  list(
    values = values,
    vectors = vectors,
    nconv = fit$nconv,
    niter = fit$iterations,
    nops = fit$matvecs,
    certificate = cert,
    diagnostics = diagnostics(fit)
  )
}

#' RSpectra-compatible SVD shim.
#'
#' @param A Matrix, eigencore operator, or a function `f(x, args)` returning
#'   `A %*% x` (then `Atrans` and `dim` are required).
#' @param k Number of singular values to compute.
#' @param nu Number of left singular vectors returned.
#' @param nv Number of right singular vectors returned.
#' @param opts RSpectra options list. `tol`, `center` and `scale` are
#'   honoured (`center`/`scale` are applied as operators, without densifying);
#'   `ncv` is passed as `auto(max_subspace = ncv)` (used by the Golub-Kahan
#'   and implicit-Gram Lanczos routes; the explicit Gram route has no Krylov
#'   subspace). `maxitr` is accepted but not used, since [svd_partial()] has
#'   no iteration limit yet; other keys raise a warning.
#' @param ... Additional arguments passed to [svd_partial()].
#' @param Atrans Function `f(x, args)` returning `t(A) %*% x` when `A` is a
#'   function.
#' @param dim Dimensions of `A` when `A` is a function.
#' @param args Extra argument passed to function inputs.
#' @return A list compatible with `RSpectra::svds()`, including `d`, optional
#'   `u` and `v`, convergence counts, operation counts, certificate
#'   diagnostics, and eigencore diagnostics.
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(60), 10, 6)
#' res <- svds(X, k = 2)
#' res$d
svds <- function(A, k, nu = k, nv = k, opts = list(), ..., Atrans = NULL,
                 dim = NULL, args = NULL) {
  nu <- compat_vector_count(nu, k, "nu")
  nv <- compat_vector_count(nv, k, "nv")
  if (is.function(A)) {
    A <- compat_function_svd_operator(A, Atrans = Atrans, dim = dim, args = args)
  }
  opts <- compat_check_opts(
    opts,
    used = c("tol", "ncv", "center", "scale"),
    ignored = "maxitr"
  )
  if (!is.null(opts$ncv) && !is.null(list(...)$method)) {
    stop("opts$ncv cannot be combined with method =; set max_subspace on ",
         "the method descriptor instead.", call. = FALSE)
  }
  A <- compat_center_scale(A, center = opts$center %||% FALSE,
                           scale = opts$scale %||% FALSE)
  # Like eigs()/eigs_sym() with retvec = FALSE, always solve with both sides
  # so the result can be certified (the certificate needs U and V); nu/nv
  # only trim what is returned. Requesting fewer vectors used to return an
  # uncertified result with an "only 0 eigenvalue(s) converged" warning
  # (oracle sweep), unlike RSpectra.
  vector_mode <- "both"
  controls <- if (is.null(opts$tol)) list() else list(tol = opts$tol)
  if (!is.null(opts$ncv)) {
    controls$method <- auto(max_subspace = as.integer(opts$ncv))
  }
  fit <- do.call(svd_partial, c(
    list(A, rank = k, vectors = vector_mode),
    controls,
    list(...)
  ))
  compat_warn_nconv(fit$nconv, k)
  list(
    d = fit$d,
    u = if (nu > 0 && !is.null(fit$u)) fit$u[, seq_len(min(nu, ncol(fit$u))), drop = FALSE],
    v = if (nv > 0 && !is.null(fit$v)) fit$v[, seq_len(min(nv, ncol(fit$v))), drop = FALSE],
    nconv = fit$nconv,
    niter = fit$iterations,
    nops = fit$matvecs,
    certificate = fit$certificate,
    diagnostics = diagnostics(fit)
  )
}

#' @keywords internal
target_from_which <- function(which, k = NULL) {
  both_ends_from_k <- function(k) {
    k <- as.integer(k)
    if (length(k) != 1L || is.na(k) || k < 1L) {
      stop("ARPACK which = 'BE' requires a positive k.", call. = FALSE)
    }
    k_low <- k %/% 2L
    k_high <- k - k_low
    both_ends(k_low, k_high)
  }
  if (!is.character(which) || length(which) != 1L || is.na(which)) {
    stop("which must be a single ARPACK selector string.", call. = FALSE)
  }
  code <- toupper(which)
  switch(
    code,
    LM = largest_magnitude(),
    SM = smallest_magnitude(),
    LA = largest(),
    SA = smallest(),
    LR = largest_real(),
    SR = smallest_real(),
    LI = largest_imaginary(),
    SI = smallest_imaginary(),
    BE = both_ends_from_k(k),
    stop("Unknown ARPACK selector which = '", which, "'.", call. = FALSE)
  )
}

#' @keywords internal
compat_eigen_target <- function(which, k, sigma, symmetric) {
  if (is.null(sigma)) {
    if (!symmetric && toupper(which) %in% c("LA", "SA", "BE")) {
      stop("which = '", which, "' is only valid for eigs_sym().", call. = FALSE)
    }
    if (symmetric && toupper(which) %in% c("LR", "SR", "LI", "SI")) {
      stop("which = '", which, "' is only valid for eigs().", call. = FALSE)
    }
    return(target_from_which(which, k = k))
  }
  if (!is.numeric(sigma) && !is.complex(sigma) || length(sigma) != 1L ||
      is.na(sigma)) {
    stop("sigma must be a single finite number.", call. = FALSE)
  }
  if (!identical(toupper(which), "LM")) {
    stop(
      "With sigma, only which = 'LM' (eigenvalues nearest sigma) is supported.",
      call. = FALSE
    )
  }
  nearest(sigma)
}

#' @keywords internal
compat_check_opts <- function(opts, used, ignored) {
  if (is.null(opts)) {
    return(list())
  }
  if (!is.list(opts)) {
    stop("opts must be a list.", call. = FALSE)
  }
  keys <- names(opts) %||% rep("", length(opts))
  unknown <- setdiff(keys, c(used, ignored))
  if (length(unknown)) {
    warning("Unknown opts entries ignored: ",
            paste(unknown, collapse = ", "), ".", call. = FALSE)
  }
  skipped <- intersect(keys, ignored)
  if (length(skipped)) {
    warning("opts entries ", paste(skipped, collapse = ", "),
            " are accepted for RSpectra compatibility but not used.",
            call. = FALSE)
  }
  opts[intersect(keys, used)]
}

#' @keywords internal
compat_eigen_opts <- function(opts, allow_initvec, dots = list()) {
  used <- c("tol", "ncv", "maxitr", "retvec", if (allow_initvec) "initvec")
  ignored <- c("mode", if (!allow_initvec) "initvec")
  opts <- compat_check_opts(opts, used = used, ignored = ignored)
  out <- list()
  if (!is.null(opts$tol)) out$tol <- opts$tol
  # ARPACK maxitr is the restart (outer iteration) limit: eigencore's maxit.
  if (!is.null(opts$maxitr) && is.null(dots$maxit)) {
    out$maxit <- as.integer(opts$maxitr)
  }
  ncv <- if (is.null(opts$ncv)) NULL else as.integer(opts$ncv)
  if (!is.null(ncv) && !is.null(dots$method)) {
    stop("opts$ncv cannot be combined with method =; set max_subspace on ",
         "the method descriptor instead.", call. = FALSE)
  }
  # retvec = FALSE still solves with vectors so the result stays certified
  # (as ARPACK does internally); the shims drop them afterwards.
  if (!is.null(opts$initvec)) {
    out$initial_subspace <- matrix(as.numeric(opts$initvec), ncol = 1L)
    out$method <- lanczos(max_subspace = ncv)
  } else if (!is.null(ncv)) {
    # The subspace request rides on auto(), so the planner's route choice is
    # the same as without ncv.
    out$method <- auto(max_subspace = ncv)
  }
  out
}

#' @keywords internal
compat_warn_nconv <- function(nconv, k) {
  if (!is.null(nconv) && length(nconv) == 1L && !is.na(nconv) && nconv < k) {
    warning("only ", nconv, " eigenvalue(s) converged, less than k = ", k,
            call. = FALSE)
  }
  invisible(NULL)
}

#' @keywords internal
compat_vector_count <- function(x, k, name) {
  x <- as.integer(x)
  if (length(x) != 1L || is.na(x) || x < 0L || x > k) {
    stop(name, " must be an integer between 0 and k.", call. = FALSE)
  }
  x
}

#' @keywords internal
compat_symmetric_from_triangle <- function(A, lower) {
  if (is.matrix(A)) {
    if (nrow(A) != ncol(A)) {
      stop("A must be square.", call. = FALSE)
    }
    if (is.double(A) &&
        isTRUE(.Call("eigencore_dense_is_symmetric", A, 0, PACKAGE = "eigencore"))) {
      return(A)
    }
    mirror <- Conj(t(A))
    if (isTRUE(lower)) {
      A[upper.tri(A)] <- mirror[upper.tri(A)]
    } else {
      A[lower.tri(A)] <- mirror[lower.tri(A)]
    }
    return(A)
  }
  if (inherits(A, "Matrix") && !inherits(A, "symmetricMatrix") &&
      !inherits(A, "diagonalMatrix")) {
    return(Matrix::forceSymmetric(A, uplo = if (isTRUE(lower)) "L" else "U"))
  }
  A
}

#' @keywords internal
compat_function_operator <- function(A, n, args, structure) {
  if (!is.function(A)) {
    return(A)
  }
  if (is.null(n)) {
    stop("n must be supplied when A is a function.", call. = FALSE)
  }
  n <- as.integer(n)
  f <- A
  apply_fun <- compat_columnwise(f, n, args)
  linear_operator(
    dim = c(n, n),
    apply = apply_fun,
    apply_adjoint = if (identical(structure$kind, "hermitian")) apply_fun else NULL,
    structure = structure,
    name = "RSpectra function operator"
  )
}

#' @keywords internal
compat_function_svd_operator <- function(A, Atrans, dim, args) {
  if (!is.function(Atrans) || is.null(dim) || length(dim) != 2L) {
    stop("Function A requires Atrans and a length-2 dim.", call. = FALSE)
  }
  dim <- as.integer(dim)
  linear_operator(
    dim = dim,
    apply = compat_columnwise(A, dim[1L], args),
    apply_adjoint = compat_columnwise(Atrans, dim[2L], args),
    name = "RSpectra function operator"
  )
}

#' @keywords internal
compat_columnwise <- function(f, out_rows, args) {
  force(f)
  force(out_rows)
  force(args)
  function(X, alpha = 1, beta = 0, Y = NULL) {
    X <- as.matrix(X)
    Z <- matrix(0, out_rows, ncol(X))
    for (j in seq_len(ncol(X))) {
      Z[, j] <- as.numeric(f(X[, j], args))
    }
    Z <- alpha * Z
    if (!is.null(Y) && beta != 0) Z + beta * Y else Z
  }
}

#' @keywords internal
compat_center_scale <- function(A, center, scale) {
  if (isFALSE(center) && isFALSE(scale)) {
    return(A)
  }
  if (!is.matrix(A) && !inherits(A, "Matrix")) {
    stop("opts$center and opts$scale require a matrix input.", call. = FALSE)
  }
  m <- nrow(A)
  col_sums <- if (inherits(A, "Matrix")) Matrix::colSums(A) else colSums(A)
  col_means <- if (isTRUE(center)) {
    col_sums / m
  } else if (is.numeric(center)) {
    if (length(center) != ncol(A)) {
      stop("opts$center must have length ncol(A).", call. = FALSE)
    }
    as.numeric(center)
  } else {
    NULL
  }
  scale_by <- if (isTRUE(scale)) {
    col_sq <- if (inherits(A, "Matrix")) Matrix::colSums(A^2) else colSums(A^2)
    cm <- col_means %||% 0
    sqrt(pmax(col_sq - 2 * cm * col_sums + m * cm^2, 0) / (m - 1))
  } else if (is.numeric(scale)) {
    if (length(scale) != ncol(A)) {
      stop("opts$scale must have length ncol(A).", call. = FALSE)
    }
    as.numeric(scale)
  } else {
    NULL
  }
  op <- if (is.null(col_means)) {
    A
  } else {
    center(A, rows = FALSE, columns = TRUE, col_means = col_means)
  }
  if (!is.null(scale_by)) {
    if (any(!is.finite(scale_by)) || any(scale_by <= 0)) {
      stop("opts$scale produced zero or non-finite column scales.",
           call. = FALSE)
    }
    op <- scale_cols(op, 1 / scale_by)
  }
  op
}
