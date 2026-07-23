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
  if (!is_csc && !is_dense) {
    stop("Native implicit Gram SVD requires a dense double matrix or dgCMatrix operator.",
         call. = FALSE)
  }
  if (!implicit_gram_svd_target_supported(target)) {
    stop("Native implicit Gram SVD supports largest singular-value targets.",
         call. = FALSE)
  }

  A <- if (is_csc) op$metadata$matrix else source
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
  block <- if (is.null(block)) (if (is_csc) 1L else 2L) else as.integer(block)
  m_max <- if (is.null(max_subspace)) {
    min(outer, default_implicit_gram_max_subspace(rank, block, sparse = is_csc))
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

  # Locking inside the kernel uses tol * (1 + theta) * ||v||, matching the
  # package's tol * max(|value|, 1) convention: theta = sigma^2, so the
  # implied singular residual bound is ~ tol * sigma * ||A|| for the dominant
  # triplets. The exact original-coordinate certificate below is the
  # authoritative pass/fail decision.
  start <- matrix(stats::rnorm(outer * block), nrow = outer, ncol = block)
  iter <- if (is_csc) {
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
      as.numeric(tol),
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
      as.numeric(tol),
      max_restarts,
      0.0,
      start,
      PACKAGE = "eigencore"
    )
  }

  lambda <- iter$values
  W <- iter$vectors
  sigma <- sqrt(pmax(lambda, 0))
  zero_tol <- gram_svd_zero_tolerance(sigma, tol)
  inv_sigma <- ifelse(sigma > zero_tol, 1 / sigma, 0)

  Av <- NULL
  if (side == 0L) {
    v <- W
    Av <- as.matrix(A %*% v)
    u <- sweep(Av, 2L, inv_sigma, `*`)
  } else {
    u <- W
    Atu <- as.matrix(Matrix::crossprod(A, u))
    v <- sweep(Atu, 2L, inv_sigma, `*`)
    Av <- as.matrix(A %*% v)
  }

  cert <- certify_svd_operator_cached_av(op, sigma, u, v, Av, tol = tol)

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
    iterations = iter$iterations %||% NA_integer_,
    matvecs = iter$matvecs %||% NA_integer_,
    stage_seconds = iter$stage_seconds %||% numeric(),
    restart = list(
      kind = "implicit_gram_thick_restart_lanczos",
      implemented = TRUE,
      native = TRUE,
      gram_side = if (side == 0L) "right" else "left",
      gram_dimension = outer,
      normal_operator_implicit = TRUE,
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
      certified_in_original_coordinates = TRUE
    )
  )
}
