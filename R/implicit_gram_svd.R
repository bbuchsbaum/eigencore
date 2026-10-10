#' @keywords internal
implicit_gram_svd_target_supported <- function(target) {
  inherits(target, "eigencore_target") &&
    target$kind %in% c("largest", "largest_magnitude")
}

#' @keywords internal
default_implicit_gram_max_subspace <- function(k, block, sparse = FALSE) {
  k <- as.integer(k)
  block <- as.integer(block)
  # Benchmarks across dense and sparse regimes (4000x1000, 2000x2000,
  # 20000x5000, 100000x2000; k in 10..50) favor a subspace well above the
  # lean 2k+1 sizing: restarts on the normal operator are cheap and the
  # extra room roughly halves the total operator applications. Sparse
  # operators benefit from a larger factor still, since their applies are
  # memory-bound and each avoided restart saves a full pass over the
  # nonzeros.
  factor <- if (isTRUE(sparse)) 6L else 4L
  max(factor * k + 2L * block, 40L)
}

# Column-centred (optionally right-scaled) sparse operator C = (A - 1 mu') D
# from center() / scale_cols(center()) (P20): list(matrix, col_means,
# weights), or NULL when `op` is not one. Row or double centring is not
# covered (its correction is not a column rank-one term).
#' @keywords internal
implicit_gram_centered_csc_parts <- function(op) {
  md <- op$metadata %||% list()
  storage <- md$storage %||% NULL
  if (!isTRUE(md$native) ||
      !(identical(storage, "centered_dgCMatrix") ||
        identical(storage, "centered_scaled_dgCMatrix")) ||
      !isTRUE(md$columns) || isTRUE(md$rows)) {
    return(NULL)
  }
  A <- md$base_matrix %||% NULL
  if (!inherits(A, "dgCMatrix")) {
    return(NULL)
  }
  n <- ncol(A)
  means <- as.numeric(md$col_means %||% numeric())
  weights <- if (identical(storage, "centered_scaled_dgCMatrix")) {
    as.numeric(md$weights %||% numeric())
  } else {
    rep(1, n)
  }
  if (length(means) != n || length(weights) != n ||
      any(!is.finite(means)) || any(!is.finite(weights))) {
    return(NULL)
  }
  list(matrix = A, col_means = means, weights = weights)
}

#' Implicit normal-equations (Gram) partial SVD.
#'
#' Runs the production block thick-restart Lanczos on the smaller-side normal
#' operator (\eqn{A^T A} or \eqn{A A^T}) without materializing the Gram
#' matrix, then recovers the opposite singular factor and certifies the
#' triplets with the exact two-sided residual in original coordinates. This
#' removes the explicit-Gram memory/dimension caps: cost per operator
#' application is one forward and one adjoint apply of \eqn{A}.
#'
#' @keywords internal
native_implicit_gram_svd <- function(op, rank, target = largest(), tol = 1e-8,
                                     vectors = c("both", "left", "right", "none"),
                                     block = NULL,
                                     max_subspace = NULL,
                                     max_restarts = 100L) {
  vectors <- match.arg(vectors)
  op <- as_operator(op)
  source <- source_or_null(op)
  storage <- op$metadata$storage %||% NULL
  is_csc <- identical(storage, "dgCMatrix")
  is_dense <- is.matrix(source) && is.double(source) && !is.complex(source)
  centered <- if (!is_csc && !is_dense) implicit_gram_centered_csc_parts(op) else NULL
  if (!is_csc && !is_dense && is.null(centered)) {
    stop("Native implicit Gram SVD requires a dense double matrix, a dgCMatrix ",
         "or a column-centred dgCMatrix operator.", call. = FALSE)
  }
  sparse <- is_csc || !is.null(centered)
  if (!implicit_gram_svd_target_supported(target)) {
    stop("Native implicit Gram SVD supports largest singular-value targets.",
         call. = FALSE)
  }

  A <- if (is_csc) op$metadata$matrix else if (!is.null(centered)) centered$matrix else source
  m <- as.integer(op$dim[1L])
  n <- as.integer(op$dim[2L])
  limit <- min(m, n)
  rank <- min(as.integer(rank), limit)
  if (rank < 1L) {
    stop("rank must be positive.", call. = FALSE)
  }

  # Operate on the smaller side: side 0 builds the A^T A (right) subspace,
  # side 1 the A A^T (left) subspace.
  side <- if (n <= m) 0L else 1L
  outer <- if (side == 0L) n else m
  # Dense applies amortize the matrix read across block columns; the sparse
  # kernel's per-column working set makes single-vector blocks faster there.
  block <- if (is.null(block)) (if (sparse) 1L else 2L) else as.integer(block)
  m_max <- if (is.null(max_subspace)) {
    min(outer, default_implicit_gram_max_subspace(rank, block, sparse = sparse))
  } else {
    min(outer, as.integer(max_subspace))
  }
  if (m_max < rank + block) {
    m_max <- min(outer, rank + block)
  }
  if (m_max < rank + block) {
    stop("implicit Gram SVD requires max_subspace >= rank + block.", call. = FALSE)
  }
  max_restarts <- as.integer(max_restarts)

  # Locking inside the kernel uses tol * (||A^T A|| + theta) * ||v|| with
  # theta = sigma^2. The recovered triplet then has a right residual of about
  # r_normal / sigma, i.e. an SVD backward error of roughly
  # kernel_tol * (sigma_max / sigma + sigma / sigma_max) against ||A||_2
  # (the certificate's normwise 2-norm definition, C12). That factor is ~2 for
  # a flat top spectrum, so the kernel first runs at tol / 2; if the exact
  # certificate still fails, one tighter run (sized from the observed
  # amplification, warm-started from the top Ritz vectors) follows. The exact
  # original-coordinate certificate is the authoritative pass/fail decision.
  run_kernel <- function(kernel_tol, start) {
    iter <- if (!is.null(centered)) {
      .Call(
        "eigencore_normal_thick_restart_lanczos_centered_scaled_csc",
        methods::slot(A, "i"),
        methods::slot(A, "p"),
        methods::slot(A, "x"),
        methods::slot(A, "Dim"),
        centered$col_means,
        centered$weights,
        as.integer(side),
        as.integer(rank),
        as.integer(m_max),
        as.integer(block),
        1L,  # largest eigenvalues of the normal operator
        as.numeric(kernel_tol),
        max_restarts,
        0.0,
        start,
        PACKAGE = "eigencore"
      )
    } else if (is_csc) {
      .Call(
        "eigencore_normal_thick_restart_lanczos_csc",
        methods::slot(A, "i"),
        methods::slot(A, "p"),
        methods::slot(A, "x"),
        methods::slot(A, "Dim"),
        as.integer(side),
        as.integer(rank),
        as.integer(m_max),
        as.integer(block),
        1L,  # largest eigenvalues of the normal operator
        as.numeric(kernel_tol),
        max_restarts,
        0.0,
        start,
        PACKAGE = "eigencore"
      )
    } else {
      .Call(
        "eigencore_normal_thick_restart_lanczos_dense",
        A,
        as.integer(side),
        as.integer(rank),
        as.integer(m_max),
        as.integer(block),
        1L,
        as.numeric(kernel_tol),
        max_restarts,
        0.0,
        start,
        PACKAGE = "eigencore"
      )
    }
    lambda <- iter$values
    W <- iter$vectors
    # The thick-restart kernel returns pairs in lock order; singular values are
    # reported in decreasing order (C44).
    if (length(lambda) > 1L && is.unsorted(rev(lambda))) {
      perm <- order(lambda, decreasing = TRUE)
      lambda <- lambda[perm]
      W <- W[, perm, drop = FALSE]
      if (length(iter$residuals) == length(perm)) {
        iter$residuals <- iter$residuals[perm]
      }
    }
    sigma <- sqrt(pmax(lambda, 0))
    zero_tol <- gram_svd_zero_tolerance(sigma, tol)
    inv_sigma <- ifelse(sigma > zero_tol, 1 / sigma, 0)

    # Recover the opposite factor, then certify with fresh forward AND adjoint
    # applies in original coordinates (certify_svd_operator); the product used
    # to form u (or v) is not reused as a cached side of the certificate (C13).
    # The 1 / sigma scaling is applied to the short Gram-side block before
    # the product (one pass over an outer x rank block instead of the long
    # side).
    if (side == 0L) {
      v <- W
      vs <- v * rep(inv_sigma, each = nrow(v))
      u <- if (!is.null(centered)) {
        as.matrix(apply_operator(op, vs))
      } else {
        as.matrix(A %*% vs)
      }
    } else {
      u <- W
      us <- u * rep(inv_sigma, each = nrow(u))
      v <- if (!is.null(centered)) {
        as.matrix(apply_adjoint_operator(op, us))
      } else {
        as.matrix(Matrix::crossprod(A, us))
      }
    }
    cert <- certify_svd_operator(op, sigma, u, v, tol = tol)
    list(iter = iter, lambda = lambda, W = W, sigma = sigma,
         zero_tol = zero_tol, u = u, v = v, cert = cert,
         kernel_tol = kernel_tol)
  }

  first_tol <- tol / 2
  run <- run_kernel(
    first_tol,
    matrix(stats::rnorm(outer * block), nrow = outer, ncol = block)
  )
  total_iterations <- run$iter$iterations %||% NA_integer_
  total_matvecs <- run$iter$matvecs %||% NA_integer_
  retried <- FALSE
  be <- run$cert$backward_error
  if (!isTRUE(run$cert$passed) && isTRUE(run$cert$orthogonality_passed) &&
      length(be) && all(is.finite(be)) && max(be) > tol) {
    # Observed amplification max(be) / kernel_tol; aim for half the tolerance.
    retry_tol <- max(first_tol * min(0.5, 0.5 * tol / max(be)), tol * 1e-4)
    retry <- run_kernel(retry_tol, run$W[, seq_len(block), drop = FALSE])
    retried <- TRUE
    total_iterations <- total_iterations + (retry$iter$iterations %||% NA_integer_)
    total_matvecs <- total_matvecs + (retry$iter$matvecs %||% NA_integer_)
    if (isTRUE(retry$cert$passed) ||
        max(retry$cert$backward_error) <= max(be)) {
      run <- retry
    }
  }
  iter <- run$iter
  lambda <- run$lambda
  sigma <- run$sigma
  zero_tol <- run$zero_tol
  u <- run$u
  v <- run$v
  cert <- run$cert

  u_out <- u
  v_out <- v
  if (vectors == "left") {
    v_out <- NULL
  } else if (vectors == "right") {
    u_out <- NULL
  } else if (vectors == "none") {
    u_out <- NULL
    v_out <- NULL
  }

  list(
    d = sigma,
    u = u_out,
    v = v_out,
    values = sigma,
    residuals = cert$residuals,
    backward_error = cert$backward_error,
    orthogonality = cert$orthogonality,
    certificate = cert,
    iterations = total_iterations,
    matvecs = total_matvecs,
    # Each normal-operator application is one forward and one adjoint apply
    # of A, so work() can split the base applies (P20).
    adjoint_matvecs = if (is.na(total_matvecs)) NA_integer_ else as.integer(total_matvecs %/% 2L),
    stage_seconds = iter$stage_seconds %||% numeric(),
    restart = list(
      kind = "implicit_gram_thick_restart_lanczos",
      implemented = TRUE,
      native = TRUE,
      gram_side = if (side == 0L) "right" else "left",
      gram_dimension = outer,
      normal_operator_implicit = TRUE,
      centered_csc = !is.null(centered),
      materialized_gram = FALSE,
      block = block,
      max_subspace = m_max,
      restarts = iter$restarts %||% NA_integer_,
      n_locked = iter$n_locked %||% NA_integer_,
      locking_events = iter$locking_events %||% NA_integer_,
      ortho_passes = iter$ortho_passes %||% NA_integer_,
      normal_lambda = lambda,
      normal_residuals = iter$residuals,
      zero_singular_threshold = zero_tol,
      zero_singular_completion = any(sigma <= zero_tol),
      certified_in_original_coordinates = TRUE,
      kernel_tol = run$kernel_tol,
      tightened_kernel_retry = retried
    )
  )
}
