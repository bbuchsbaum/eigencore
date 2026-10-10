# Complex Hermitian eigenproblems on the iterative routes.
#
# Realification. A complex Hermitian H = A + iB (A real symmetric, B real
# skew-symmetric) acts on C^n as the real symmetric
#
#   R(H) = [ A  -B ]
#          [ B   A ]
#
# acts on R^{2n}: z = x + iy  <->  u = (x; y), and R(H) u = (Re Hz; Im Hz).
# R(H) commutes with J = [0 -I; I 0] (multiplication by i), so every
# eigenvalue of H is an eigenvalue of R(H) of twice the multiplicity, with
# the real eigenvectors (x; y) and J (x; y) = (-y; x). For a Hermitian-
# definite pencil (H, M), R(M) is symmetric positive definite and the same
# holds for (R(H), R(M)). Every norm the certificates use carries over:
# ||R(H) u - lambda R(M) u|| = ||Hz - lambda Mz||, ||R(H)||_2 = ||H||_2.
#
# Route. The complex problem is solved as the real problem
# (R(H), R(M)) for 2k eigenpairs through the ordinary real pipeline (native
# block Lanczos, shift-invert, generalized Lanczos/LOBPCG, interval slicing,
# dense LAPACK), so every native real kernel, the completeness probe and the
# LDL' inertia counts run unchanged. The 2k real Ritz vectors are mapped back
# by a complex Rayleigh-Ritz step: Z = top + i bottom spans (to the solve's
# accuracy) the k-dimensional complex invariant subspace; an M-orthonormal
# basis of span(Z) (rank decided by the Gram spectrum) is projected onto H
# and the k target pairs are selected, then certified from scratch on the
# complex operator. The real route's completeness verdict transfers: the
# returned real set is exactly the doubled complex set (the requested 2k
# real eigenvalues are the k requested complex ones, each twice), and the
# eigenvalue counts of R(H) - sigma R(M) are twice those of H - sigma M.
#
# Completeness of complex results from other routes (reference complex
# Lanczos / LOBPCG) is checked the same way: the returned pairs are
# realified (2k real pairs) and the real Hermitian completeness check
# (inertia count or deflated-complement probe, with repair) runs on the
# realified problem; a repaired set is mapped back by the same Rayleigh-Ritz
# step and re-certified.
#
# Why realify rather than a native complex Krylov kernel: see
# docs/review-2026-10.md ("Complex Hermitian iterative routes"). In short,
# the realified block route reuses the certified real kernels and the whole
# completeness machinery unchanged, is 2.6x / 4x / 6.5x faster than the
# dense zheev route at n = 500 / 1000 / 2000 (k = 10, one thread, reference
# BLAS) with the same eigenvalue accuracy, and costs at most about twice the
# operator work a native complex Lanczos kernel would need (every complex
# direction is found twice).

# ---------------------------------------------------------------------------
# Labels and policy
# ---------------------------------------------------------------------------

#' @keywords internal
complex_hermitian_realified_label <- function() {
  "native realified complex Hermitian route (real 2n embedding)"
}

#' @keywords internal
plan_dispatches_complex_hermitian_realified <- function(plan) {
  identical(plan$method %||% "", complex_hermitian_realified_label())
}

# Does the problem carry a complex operator (A or B)?
#' @keywords internal
complex_hermitian_problem <- function(problem) {
  identical(problem$structure$kind %||% NULL, "hermitian") &&
    (identical(problem$A$dtype %||% "double", "complex") ||
       (!is.null(problem$metric) &&
          identical(problem$metric$dtype %||% "double", "complex")))
}

#' @keywords internal
complex_hermitian_realify_enabled <- function() {
  !isFALSE(getOption("eigencore.complex_hermitian_realify", TRUE))
}

# Largest complex dense n realified as an explicit 2n x 2n real matrix
# (32 n^2 bytes); above it the embedding applies the complex source.
#' @keywords internal
complex_hermitian_dense_realify_max_n <- function() {
  n <- suppressWarnings(as.integer(getOption("eigencore.complex_realify_dense_max_n", 6000L)))
  if (length(n) != 1L || is.na(n) || n < 1L) 6000L else n
}

# Is a standard dense complex Hermitian auto() solve small enough for zheev?
# Same thresholds as the real dense partial-Lanczos switch.
#' @keywords internal
complex_hermitian_dense_lapack_preferred <- function(problem, k) {
  src <- source_or_null(problem$A)
  if (!(is.matrix(src) && is.complex(src)) || !is.null(problem$metric)) {
    return(FALSE)
  }
  n <- as.integer(problem$A$dim[[1L]])
  min_n <- as.integer(getOption("eigencore.dense_partial_lanczos_min_n", 128L))
  max_fraction <- as.numeric(getOption("eigencore.dense_partial_lanczos_max_fraction", 0.25))
  if (length(min_n) != 1L || is.na(min_n) || min_n < 1L) min_n <- 128L
  if (length(max_fraction) != 1L || is.na(max_fraction) || max_fraction <= 0 ||
      max_fraction > 1) {
    max_fraction <- 0.25
  }
  n < min_n || (k / n) > max_fraction
}

# Planner hook (plan_solver): returns the frozen plan for a complex
# Hermitian problem, or NULL to continue with the ordinary planner (small
# dense standard auto() -> zheev, lobpcg() -> reference complex LOBPCG,
# realification disabled -> reference complex Lanczos).
#' @keywords internal
complex_hermitian_plan <- function(problem, k, method, method_descriptor,
                                   execution, planner_policy, maxit,
                                   initial_subspace = NULL) {
  if (!complex_hermitian_problem(problem) || !complex_hermitian_realify_enabled()) {
    return(NULL)
  }
  kind <- if (inherits(method, "eigencore_method")) method$kind else "auto"
  if (identical(kind, "lobpcg")) {
    return(NULL)
  }
  if (!kind %in% c("auto", "lanczos", "shift_invert")) {
    return(NULL)
  }
  transform <- problem$transform
  if (identical(kind, "auto") && is.null(transform) && is.null(initial_subspace) &&
      !is_interval_target(problem$target) &&
      complex_hermitian_dense_lapack_preferred(problem, k)) {
    # Small dense standard problem: zheev returns the whole spectrum, so
    # every target (nearest() included) is selected exactly; no implicit
    # shift-invert (which has no complex factorisation).
    label <- native_dense_complex_hermitian_label()
    return(new_plan(
      problem,
      k = k,
      method = label,
      method_descriptor = method_descriptor,
      reasons = c(
        paste0("structure: ", problem$structure$kind, " (complex)"),
        paste0("target: ", target_label(problem$target)),
        "standard eigenproblem",
        "small dense complex Hermitian matrix: full zheev decomposition",
        operator_kernel_reason(problem$A)
      ),
      fallback = "dense oracle prototype",
      controls = resolve_iteration_limit_controls(list(), problem, label, maxit),
      execution = execution,
      planner_policy = planner_policy
    ))
  }
  if (identical(kind, "shift_invert") || is_transform_method(transform)) {
    sigma <- (if (identical(kind, "shift_invert")) method else transform)$sigma
    if (is.complex(sigma) && isTRUE(Im(sigma) != 0)) {
      stop("shift_invert(sigma) on a complex Hermitian problem needs a real ",
           "sigma: the eigenvalues are real.", call. = FALSE)
    }
  }
  inner <- complex_hermitian_inner_method(problem, method, k)
  n <- as.integer(problem$A$dim[[1L]])
  reasons <- c(
    paste0("structure: ", problem$structure$kind, " (complex)"),
    paste0("target: ", target_label(problem$target)),
    if (is.null(problem$metric)) "standard eigenproblem" else
      "complex Hermitian-definite pencil (metric B supplied)",
    paste0("complex Hermitian operator solved through its real symmetric 2n = ",
           2L * n, " embedding [Re -Im; Im Re] (eigenvalues doubled; k -> 2k)"),
    paste0("inner real route method: ", complex_hermitian_method_summary(inner)),
    operator_kernel_reason(problem$A)
  )
  controls <- list(
    embedding = "real_2n",
    inner_method = complex_hermitian_method_summary(inner),
    inner_k = if (is_interval_target(problem$target)) NA_integer_ else 2L * as.integer(k),
    k_cap = if (is_interval_target(problem$target)) as.integer(k %||% NA_integer_) else NA_integer_,
    iteration_limit = if (is.null(maxit)) NA_integer_ else as.integer(maxit),
    iteration_limit_kind = "inner_real_route"
  )
  new_plan(
    problem,
    k = k %||% n,
    method = complex_hermitian_realified_label(),
    method_descriptor = method_descriptor,
    reasons = reasons,
    fallback = "inner real route decides its own fallback",
    controls = controls,
    execution = execution,
    planner_policy = planner_policy
  )
}

# interval() targets: dense standard complex problems keep the exact dense
# route (zheev on the window); pencils, complex sparse (complex_operator())
# and matrix-free operators are counted and sliced on the real embedding.
#' @keywords internal
complex_hermitian_interval_plan <- function(problem, k, method, tol, maxit, vectors,
                                            certify, allow_dense_fallback,
                                            left_vectors) {
  if (!complex_hermitian_problem(problem) || !complex_hermitian_realify_enabled()) {
    return(NULL)
  }
  src <- source_or_null(problem$A)
  if (is.null(problem$metric) && is.matrix(src) && is.complex(src)) {
    return(NULL)
  }
  if (!is_auto_method(method)) {
    stop("interval() targets choose their own route (dense LAPACK or LDL' ",
         "shift-invert spectrum slicing); use method = auto().", call. = FALSE)
  }
  n <- as.integer(problem$A$dim[[1L]])
  if (!is.null(k)) {
    k <- validate_solution_count(k, n, "k")
  }
  execution <- new_plan_execution(
    "eigen", tol = tol, maxit = maxit, vectors = vectors, certify = certify,
    allow_dense_fallback = allow_dense_fallback, initial_subspace = NULL
  )
  execution$left_vectors <- left_vectors
  complex_hermitian_plan(problem, k, method, method, execution,
                         planner_policy_snapshot(), maxit)
}

# Routes that run complex arithmetic directly (complex matrix-free operators
# are accepted there).
#' @keywords internal
complex_hermitian_reference_labels <- function() {
  c(
    "reference LOBPCG prototype",
    "reference generalized SPD LOBPCG prototype",
    "reference Hermitian Lanczos (prototype/oracle fallback)",
    "reference Hermitian Lanczos (target unsupported by native path)"
  )
}

#' @keywords internal
complex_hermitian_method_summary <- function(method) {
  kind <- method$kind %||% "auto"
  extra <- switch(
    kind,
    lanczos = paste0("(block = ", method$block %||% 1L, ")"),
    shift_invert = paste0("(sigma = ", format(Re(method$sigma), digits = 8), ")"),
    ""
  )
  paste0(kind, extra)
}

# Auto block size for the realified block Lanczos route: every eigenvalue is
# (at least) a double one, so block >= 2. Block 4 suits the BLAS-3 dense
# kernel from n = 500 (dense n = 500 / 1000 / 2000, k = 10, one thread:
# 0.27 / 1.43 / 5.44 s at block 2 vs 0.17 / 1.13 / 4.93 s at block 4); sparse
# and matrix-free operators stay at block 2, which needs fewer operator
# columns (random-phase grid Laplacian n = 3600, k = 10: 876 vs 1005).
#' @keywords internal
complex_hermitian_auto_block <- function(problem) {
  src <- source_or_null(problem$A)
  n <- as.integer(problem$A$dim[[1L]])
  if (is.matrix(src) && n >= 500L) 4L else 2L
}

# The method the real 2n problem is solved with.
#' @keywords internal
complex_hermitian_inner_method <- function(problem, method, k) {
  kind <- if (inherits(method, "eigencore_method")) method$kind else "auto"
  completeness <- if (inherits(method, "eigencore_method")) method$completeness else NULL
  double_or_null <- function(x) if (is.null(x)) NULL else 2L * as.integer(x)
  out <- switch(
    kind,
    lanczos = lanczos(
      max_subspace = double_or_null(method$max_subspace),
      max_restarts = method$max_restarts,
      block = max(2L, 2L * as.integer(method$block %||% 1L)),
      check_stride = method$check_stride %||% 0L,
      reorthogonalize = method$reorthogonalize %||% TRUE
    ),
    shift_invert = complex_hermitian_inner_shift_invert(method),
    {
      target_kind <- if (inherits(problem$target, "eigencore_target")) problem$target$kind else ""
      sparse_split <- identical(problem$A$metadata$storage %||% "", "complex_split_sparse")
      block_targets <- c("largest", "largest_magnitude",
                         if (!sparse_split) "smallest")
      if (is.null(problem$transform) && is.null(problem$metric) &&
          target_kind %in% block_targets) {
        # Block Lanczos resolves the doubled eigenvalues directly; scalar
        # routes would find one copy per eigenvalue and lean on the
        # completeness repair. A sparse smallest() target keeps the real
        # auto() planner, which shift-inverts below the spectrum when the
        # LDL' factor is cheap (C60): random-phase grid Laplacian n = 10000,
        # k = 10: 185 operator columns / 1.6 s vs 1034 / 3.7 s.
        lanczos(block = complex_hermitian_auto_block(problem),
                max_subspace = double_or_null(method$max_subspace))
      } else {
        auto(max_subspace = double_or_null(method$max_subspace))
      }
    }
  )
  if (!is.null(completeness)) {
    out$completeness <- completeness
  }
  out
}

#' @keywords internal
complex_hermitian_inner_shift_invert <- function(method) {
  solve_fn <- method$solve
  if (is.function(solve_fn)) {
    user_solve <- solve_fn
    solve_fn <- function(X) {
      X <- as.matrix(X)
      n <- nrow(X) %/% 2L
      Z <- X[seq_len(n), , drop = FALSE] + 1i * X[n + seq_len(n), , drop = FALSE]
      W <- as.matrix(user_solve(Z))
      rbind(Re(W), Im(W))
    }
  }
  shift_invert(Re(method$sigma), solve = solve_fn,
               factorization = method$factorization,
               max_subspace = if (is.null(method$max_subspace)) NULL else
                 2L * as.integer(method$max_subspace))
}

# ---------------------------------------------------------------------------
# Realification
# ---------------------------------------------------------------------------

#' @keywords internal
complex_hermitian_realify_dense <- function(H, symmetrize = FALSE) {
  H <- as.matrix(H)
  re <- Re(H)
  im <- Im(H)
  if (isTRUE(symmetrize)) {
    # Embedding of the Hermitian part (H + H^H) / 2: exactly symmetric.
    re <- (re + t(re)) / 2
    im <- (im - t(im)) / 2
  }
  n <- nrow(re)
  i1 <- seq_len(n)
  i2 <- n + i1
  out <- matrix(0, 2L * n, 2L * n)
  out[i1, i1] <- re
  out[i2, i2] <- re
  out[i2, i1] <- im
  out[i1, i2] <- -im
  out
}

#' @keywords internal
complex_hermitian_realify_sparse <- function(re, im) {
  re <- methods::as(methods::as(methods::as(re, "CsparseMatrix"), "generalMatrix"), "dMatrix")
  im <- methods::as(methods::as(methods::as(im, "CsparseMatrix"), "generalMatrix"), "dMatrix")
  # Embedding of the Hermitian part: exactly symmetric.
  re <- (re + Matrix::t(re)) / 2
  im <- (im - Matrix::t(im)) / 2
  out <- rbind(cbind(re, -im), cbind(im, re))
  out <- methods::as(methods::as(out, "CsparseMatrix"), "generalMatrix")
  Matrix::drop0(out)
}

# u = (x; y) and J u = (-y; x) for every complex column z = x + iy,
# interleaved: columns 2j - 1, 2j belong to z_j.
#' @keywords internal
complex_hermitian_realify_vectors <- function(V) {
  V <- as.matrix(V)
  k <- ncol(V)
  n <- nrow(V)
  x <- Re(V)
  y <- Im(V)
  out <- matrix(0, 2L * n, 2L * k)
  if (k) {
    odd <- 2L * seq_len(k) - 1L
    out[seq_len(n), odd] <- x
    out[n + seq_len(n), odd] <- y
    out[seq_len(n), odd + 1L] <- -y
    out[n + seq_len(n), odd + 1L] <- x
  }
  out
}

#' @keywords internal
complex_hermitian_complexify <- function(U, n) {
  U <- as.matrix(U)
  U[seq_len(n), , drop = FALSE] + 1i * U[n + seq_len(n), , drop = FALSE]
}

# Apply an operator (complex, or real to a complex block) without densifying.
#' @keywords internal
complex_hermitian_apply <- function(op, X) {
  apply_operator_split_complex(op, X)
}

# A complex-typed view of a real operator (for certifying complex vectors
# against a real metric B of a complex pencil).
#' @keywords internal
complex_hermitian_complex_view <- function(op) {
  if (is.null(op) || identical(op$dtype %||% "double", "complex")) {
    return(op)
  }
  linear_operator(
    dim = op$dim,
    apply = function(X, alpha = 1, beta = 0, Y = NULL) {
      out <- alpha * apply_operator_split_complex(op, as.matrix(X) + 0i)
      if (is.null(Y) || beta == 0) out else out + beta * Y
    },
    dtype = "complex",
    structure = op$structure,
    name = paste0("complex view of ", op$name %||% "operator"),
    metadata = list(
      two_norm = op$metadata$two_norm %||% NULL,
      frobenius_norm = op$metadata$frobenius_norm %||% NULL,
      positive_definite = if (isTRUE(op$metadata$positive_definite)) TRUE else NULL
    )
  )
}

#' @keywords internal
complex_hermitian_norm_metadata <- function(op) {
  meta <- list()
  if (is.numeric(op$metadata$two_norm %||% NULL)) {
    meta$two_norm <- op$metadata$two_norm
  }
  if (is.numeric(op$metadata$frobenius_norm %||% NULL)) {
    # ||R(H)||_F = sqrt(2) ||H||_F
    meta$frobenius_norm <- sqrt(2) * op$metadata$frobenius_norm
  }
  for (flag in c("positive_definite", "symmetric_positive_definite", "spd")) {
    if (isTRUE(op$metadata[[flag]])) {
      meta[[flag]] <- TRUE
    }
  }
  meta
}

# Real symmetric operator of the 2n embedding. Explicit sources become an
# explicit dense or dgCMatrix embedding (native kernels, CHOLMOD inertia);
# matrix-free operators get a realified apply.
#' @keywords internal
complex_hermitian_realify_operator <- function(op) {
  op <- as_operator(op)
  n <- as.integer(op$dim[[1L]])
  src <- source_or_null(op)
  mat <- op$metadata$matrix %||% NULL
  real_op <- !identical(op$dtype %||% "double", "complex")
  if (real_op) {
    # A real operator in a complex problem (e.g. a real SPD metric): its
    # embedding is blockdiag(M, M).
    if (is.matrix(src) && is.double(src)) {
      out <- matrix(0, 2L * n, 2L * n)
      out[seq_len(n), seq_len(n)] <- src
      out[n + seq_len(n), n + seq_len(n)] <- src
      return(as_operator(out))
    }
    if (inherits(mat, "diagonalMatrix")) {
      d <- as.numeric(Matrix::diag(mat))
      return(as_operator(Matrix::Diagonal(x = c(d, d))))
    }
    if (inherits(mat, "sparseMatrix")) {
      return(as_operator(Matrix::bdiag(mat, mat)))
    }
    return(linear_operator(
      dim = c(2L * n, 2L * n),
      apply = function(X, alpha = 1, beta = 0, Y = NULL) {
        X <- as.matrix(X)
        p <- ncol(X)
        out <- as.matrix(apply_operator(op, cbind(X[seq_len(n), , drop = FALSE],
                                                  X[n + seq_len(n), , drop = FALSE])))
        out <- alpha * rbind(out[, seq_len(p), drop = FALSE],
                             out[, p + seq_len(p), drop = FALSE])
        if (is.null(Y) || beta == 0) out else out + beta * Y
      },
      structure = op$structure,
      name = paste0("realified ", op$name %||% "operator"),
      metadata = complex_hermitian_norm_metadata(op)
    ))
  }
  split <- op$metadata$complex_split %||% NULL
  if (!is.null(split)) {
    if (inherits(split$re, "sparseMatrix") || inherits(split$im, "sparseMatrix")) {
      return(as_operator(complex_hermitian_realify_sparse(split$re, split$im)))
    }
    H <- split$re + 1i * split$im
    return(as_operator(complex_hermitian_realify_dense(H, symmetrize = TRUE)))
  }
  if (is.matrix(src) && is.complex(src) && n <= complex_hermitian_dense_realify_max_n()) {
    # Symmetrise first: the embedding of an exactly Hermitian matrix is
    # exactly symmetric (certificates are recomputed on the original).
    return(as_operator(complex_hermitian_realify_dense(src, symmetrize = TRUE)))
  }
  linear_operator(
    dim = c(2L * n, 2L * n),
    apply = function(X, alpha = 1, beta = 0, Y = NULL) {
      X <- as.matrix(X)
      Z <- X[seq_len(n), , drop = FALSE] + 1i * X[n + seq_len(n), , drop = FALSE]
      W <- as.matrix(apply_operator(op, Z))
      out <- alpha * rbind(Re(W), Im(W))
      if (is.null(Y) || beta == 0) out else out + beta * Y
    },
    structure = hermitian(),
    name = paste0("realified ", op$name %||% "complex operator"),
    metadata = complex_hermitian_norm_metadata(op)
  )
}

#' @keywords internal
complex_hermitian_realify_target <- function(target) {
  if (inherits(target, "eigencore_target") && identical(target$kind, "both_ends")) {
    return(both_ends(2L * target$value$k_low, 2L * target$value$k_high))
  }
  target
}

#' @keywords internal
complex_hermitian_realified_problem <- function(problem) {
  A2 <- complex_hermitian_realify_operator(problem$A)
  B2 <- if (is.null(problem$metric)) NULL else
    complex_hermitian_realify_operator(problem$metric)
  transform <- problem$transform
  if (is_transform_method(transform) && identical(transform$kind, "shift_invert")) {
    transform <- complex_hermitian_inner_shift_invert(transform)
  }
  P <- list(
    type = "eigen",
    A = A2,
    metric = B2,
    structure = hermitian(),
    target = complex_hermitian_realify_target(problem$target),
    transform = transform
  )
  class(P) <- "eigencore_eigen_problem"
  P
}

# ---------------------------------------------------------------------------
# Back to complex: Rayleigh-Ritz on span(top + i bottom)
# ---------------------------------------------------------------------------

# Select k of the Rayleigh-Ritz values by the complex problem's target.
#' @keywords internal
complex_hermitian_select <- function(values, target, k) {
  if (inherits(target, "eigencore_target") && identical(target$kind, "interval")) {
    if (length(values) <= k) {
      return(order(values))
    }
    inside <- which(values >= target$value$lower & values <= target$value$upper)
    return(inside[order(values[inside])])
  }
  # Target order (both_ends: low ascending, then high descending).
  utils::head(order_indices(values, target), k)
}

# Map real 2n Ritz vectors U to k complex pairs. Returns values, vectors,
# the complex rank r of span(U) and the Gram gap (r-th / (r+1)-th
# singular value ratio; Inf when r = ncol).
#' @keywords internal
complex_hermitian_from_realified <- function(problem, U, k, expected_rank = k,
                                             rank_tol = 1e-6, min_gap = 1e3) {
  n <- as.integer(problem$A$dim[[1L]])
  complex_hermitian_rr(problem, complex_hermitian_complexify(U, n), k,
                       expected_rank = expected_rank, rank_tol = rank_tol,
                       min_gap = min_gap)
}

# Complex Rayleigh-Ritz of (A, B) on span(Z): an M-orthonormal basis of the
# numerical column space of Z (rank from the Gram spectrum: expected_rank
# when the gap shows it, else every direction above rank_tol), projected
# onto A; the k target pairs are selected.
#' @keywords internal
complex_hermitian_rr <- function(problem, Z, k, expected_rank = NA_integer_,
                                 rank_tol = 1e-6, min_gap = 1e3) {
  Aop <- problem$A
  Bop <- problem$metric
  Z <- as.matrix(Z) + 0i
  BZ <- if (is.null(Bop)) Z else complex_hermitian_apply(Bop, Z)
  G <- Conj(t(Z)) %*% BZ
  G <- (G + Conj(t(G))) / 2
  eg <- eigen(G, symmetric = TRUE)
  lam <- pmax(Re(eg$values), 0)
  top <- max(lam, 0)
  sv_gap <- function(r) {
    if (r < length(lam)) sqrt(lam[[r]] / max(lam[[r + 1L]], .Machine$double.xmin)) else Inf
  }
  # A doubled real set spans a complex subspace of half its size; its
  # singular values are sqrt(2) (x2 copies) and the solve's error level.
  # Use that rank when the gap shows it, else every direction above
  # rank_tol (an inconsistent set: the caller then re-checks completeness).
  r0 <- as.integer(expected_rank)
  r <- if (!is.na(r0) && r0 >= 1L && r0 <= length(lam) && lam[[r0]] > 0 && sv_gap(r0) >= min_gap) {
    r0
  } else {
    sum(lam > rank_tol^2 * top)
  }
  if (r < 1L) {
    stop("complex Rayleigh-Ritz: the realified basis is empty.", call. = FALSE)
  }
  gap <- sv_gap(r)
  Q <- Z %*% (eg$vectors[, seq_len(r), drop = FALSE] %*%
                diag(1 / sqrt(lam[seq_len(r)]), r))
  BQ <- if (is.null(Bop)) Q else complex_hermitian_apply(Bop, Q)
  # Second pass: re-orthonormalise in the B inner product.
  G2 <- Conj(t(Q)) %*% BQ
  G2 <- (G2 + Conj(t(G2))) / 2
  e2 <- eigen(G2, symmetric = TRUE)
  W <- e2$vectors %*% diag(1 / sqrt(pmax(Re(e2$values), .Machine$double.eps)), r)
  Q <- Q %*% W
  AQ <- complex_hermitian_apply(Aop, Q)
  T <- Conj(t(Q)) %*% AQ
  T <- (T + Conj(t(T))) / 2
  et <- eigen(T, symmetric = TRUE)
  values <- Re(et$values)
  sel <- complex_hermitian_select(values, problem$target, k)
  list(
    values = values[sel],
    vectors = Q %*% et$vectors[, sel, drop = FALSE],
    rank = r,
    gap = gap,
    operator_columns = r + if (is.null(Bop)) 0L else 0L,
    metric_columns = if (is.null(Bop)) 0L else ncol(Z) + r
  )
}

# ---------------------------------------------------------------------------
# Solve
# ---------------------------------------------------------------------------

#' @keywords internal
complex_hermitian_certify <- function(problem, values, vectors, tol) {
  Bc <- complex_hermitian_complex_view(problem$metric)
  Ac <- complex_hermitian_complex_view(problem$A)
  if (!length(values)) {
    return(interval_empty_certificate(tol, "no eigenvalues returned"))
  }
  certify_eigen_operator(Ac, values, vectors, Bop = Bc, tol = tol)
}

# ---------------------------------------------------------------------------
# Polishing mapped-back pairs that miss the tolerance
# ---------------------------------------------------------------------------
#
# A completeness repair on the embedding (a deflated complement solve, for
# interior targets on the squared-shift operator) can return real pairs
# that are complete but only accurate to ~1e-6, and the complex
# Rayleigh-Ritz map cannot improve vectors it is given. When the mapped set
# misses the tolerance it is refined on the complex operator: each
# unconverged pair gets one Rayleigh-quotient (shifted inverse) iteration
# step, solved with a dense complex LU or a sparse LU of the real embedding
# where that is cheap, else the span is enlarged by a short block Krylov
# sequence started from the residuals; then a complex Rayleigh-Ritz over the
# span of all returned plus refinement vectors selects the target pairs, and
# the set is re-certified from scratch. The caller re-checks completeness.

# Shifted solve y = (A - theta B)^{-1} rhs, or NULL when no cheap solve
# exists (matrix-free operators, large dense sources, sparse pencils).
#' @keywords internal
complex_hermitian_shift_solver <- function(problem) {
  Aop <- problem$A
  Bop <- problem$metric
  n <- as.integer(Aop$dim[[1L]])
  dense_max <- suppressWarnings(as.integer(
    getOption("eigencore.complex_polish_dense_max_n", 4000L)))
  if (length(dense_max) != 1L || is.na(dense_max)) dense_max <- 4000L
  dense_of <- function(op) {
    if (is.null(op)) {
      return(NULL)
    }
    src <- source_or_null(op)
    if (is.matrix(src)) {
      return(src)
    }
    split <- op$metadata$complex_split %||% NULL
    if (!is.null(split) && !inherits(split$re, "sparseMatrix")) {
      return(split$re + 1i * split$im)
    }
    mat <- op$metadata$matrix %||% NULL
    if (!is.null(mat) && n <= dense_max) {
      return(as.matrix(mat))
    }
    # A small matrix-free operator is materialised by n applies, as the
    # completeness count does (same eigencore.completeness_materialize_limit).
    if (is.null(src) && is.null(mat) &&
        n <= hermitian_completeness_controls()$materialize_limit) {
      m <- hermitian_completeness_materialize(op)
      if (!is.null(m)) {
        return(m$matrix)
      }
    }
    NULL
  }
  split <- Aop$metadata$complex_split %||% NULL
  if (is.null(Bop) && !is.null(split) && inherits(split$re, "sparseMatrix")) {
    R2 <- complex_hermitian_realify_sparse(split$re, split$im)
    I2 <- Matrix::Diagonal(2L * n)
    return(function(theta, rhs) {
      X <- as.matrix(Matrix::solve(R2 - theta * I2, rbind(Re(rhs), Im(rhs))))
      X[seq_len(n), , drop = FALSE] + 1i * X[n + seq_len(n), , drop = FALSE]
    })
  }
  if (n > dense_max) {
    return(NULL)
  }
  A <- dense_of(Aop)
  B <- if (is.null(Bop)) NULL else dense_of(Bop)
  if (is.null(A) || (!is.null(Bop) && is.null(B))) {
    return(NULL)
  }
  function(theta, rhs) {
    M <- if (is.null(B)) A - diag(theta, n) else A - theta * B
    solve(M + 0i, rhs + 0i)
  }
}

#' @keywords internal
complex_hermitian_polish <- function(problem, values, vectors, cert, tol,
                                     max_rounds = 3L) {
  k <- length(values)
  if (is.null(cert) || isTRUE(cert$passed) || !k || is.null(vectors)) {
    return(NULL)
  }
  be <- as.numeric(cert$backward_error)
  if (length(be) != k || !all(is.finite(be)) || max(be) > 1e-3) {
    # Far from converged: not a polishing job.
    return(NULL)
  }
  Aop <- problem$A
  Bc <- complex_hermitian_complex_view(problem$metric)
  n <- as.integer(Aop$dim[[1L]])
  solver <- complex_hermitian_shift_solver(problem)
  scale <- max(abs(values), as.numeric(cert$scale %||% 0), 1)
  apply_B <- function(X) if (is.null(Bc)) X else as.matrix(apply_operator(Bc, X))
  unit <- function(Y) {
    nr <- sqrt(colSums(Mod(Y)^2))
    Y[, nr > 0, drop = FALSE] %*% diag(1 / nr[nr > 0], sum(nr > 0))
  }
  V <- as.matrix(vectors) + 0i
  vals <- as.numeric(values)
  current <- cert
  rounds <- 0L
  method <- if (is.null(solver)) "block_krylov_residual_expansion" else "shifted_inverse_iteration"
  # Inverse iteration converges in one or two rounds; the Krylov expansion
  # (no solve available) gains less per round and gets more of them.
  if (is.null(solver)) {
    max_rounds <- max_rounds + 3L
  }
  for (round in seq_len(max_rounds)) {
    rounds <- round
    bad <- which(!as.logical(current$converged %||% rep(FALSE, k)))
    if (!length(bad)) {
      bad <- seq_len(k)
    }
    Vb <- V[, bad, drop = FALSE]
    Y <- if (!is.null(solver)) {
      BVb <- apply_B(Vb)
      cols <- lapply(seq_along(bad), function(j) {
        for (rel in c(0, 1e-10, -1e-8)) {
          y <- tryCatch(solver(vals[[bad[[j]]]] + rel * scale, BVb[, j, drop = FALSE]),
                        error = function(e) NULL)
          if (!is.null(y) && all(is.finite(Re(y))) && all(is.finite(Im(y)))) {
            return(y)
          }
        }
        NULL
      })
      cols <- Filter(Negate(is.null), cols)
      if (length(cols)) do.call(cbind, cols) else NULL
    } else {
      R <- as.matrix(complex_hermitian_apply(Aop, Vb)) -
        sweep(apply_B(Vb), 2L, vals[bad], `*`)
      # Block Krylov sequence from the residuals, fully reorthogonalised
      # against V and the earlier blocks (only the span matters).
      basis <- qr.Q(qr(V))
      orth_block <- function(W) {
        for (pass in 1:2) {
          W <- W - basis %*% (Conj(t(basis)) %*% W)
        }
        nr <- sqrt(colSums(Mod(W)^2))
        W <- W[, nr > 1e-10 * max(nr, .Machine$double.xmin), drop = FALSE]
        if (ncol(W)) qr.Q(qr(W)) else W
      }
      blocks <- list()
      W <- orth_block(R)
      for (step in seq_len(12L)) {
        if (!ncol(W) || ncol(basis) + ncol(W) >= n) break
        blocks[[step]] <- W
        basis <- cbind(basis, W)
        W <- orth_block(as.matrix(complex_hermitian_apply(Aop, W)))
      }
      if (length(blocks)) do.call(cbind, blocks) else NULL
    }
    if (is.null(Y) || !ncol(Y)) {
      break
    }
    rr <- tryCatch(complex_hermitian_rr(problem, cbind(V, unit(Y)), k,
                                        rank_tol = 1e-10),
                   error = function(e) NULL)
    if (is.null(rr) || length(rr$values) != k) {
      break
    }
    new_cert <- complex_hermitian_certify(problem, rr$values, rr$vectors, tol)
    if (!(max(new_cert$backward_error) < max(current$backward_error))) {
      break
    }
    V <- rr$vectors
    vals <- rr$values
    current <- new_cert
    if (isTRUE(current$passed)) {
      break
    }
  }
  if (identical(current, cert)) {
    return(NULL)
  }
  current$notes <- unique(c(cert$notes, current$notes, sprintf(paste0(
    "mapped-back pairs missed the tolerance (max backward error %.3g); ",
    "polished on the complex operator (%s, %d round(s)) and re-certified"),
    max(be), method, rounds)))
  list(values = vals, vectors = V, certificate = current, rounds = rounds,
       method = method)
}

#' @keywords internal
complex_hermitian_add_work <- function(fit) {
  context <- .eigencore_work_context$current
  w <- fit$work
  if (is.null(context) || !inherits(w, "eigencore_work")) {
    return(invisible(NULL))
  }
  for (field in work_counter_fields()) {
    value <- w[[field]]
    if (is.numeric(value) && length(value) == 1L && is.finite(value)) {
      context$counters[[field]] <- context$counters[[field]] + as.integer(value)
    }
  }
  invisible(NULL)
}

#' @keywords internal
solve_complex_hermitian_realified <- function(plan) {
  started <- proc.time()[["elapsed"]]
  problem <- plan$problem
  execution <- plan$execution
  tol <- execution$tol
  interval <- is_interval_target(problem$target)
  k <- as.integer(plan$requested)
  inner_method <- complex_hermitian_inner_method(problem, plan$method_descriptor, k)
  P2 <- complex_hermitian_realified_problem(problem)
  n <- as.integer(problem$A$dim[[1L]])
  start <- execution$initial_subspace %||% NULL
  if (!is.null(start)) {
    start <- as.matrix(start)
    if (nrow(start) != n) {
      stop("initial_subspace must have ", n, " rows.", call. = FALSE)
    }
    start <- complex_hermitian_realify_vectors(start + 0i)
  }
  k_cap <- plan$controls$k_cap %||% NA_integer_
  k2 <- if (interval) {
    if (is.na(k_cap)) NULL else 2L * k_cap
  } else {
    2L * k
  }
  inner_plan <- plan_solver(
    P2, k = k2, method = inner_method, tol = tol, maxit = execution$maxit,
    vectors = TRUE, certify = execution$certify,
    allow_dense_fallback = execution$allow_dense_fallback,
    initial_subspace = start
  )
  fit <- solve(inner_plan)
  complex_hermitian_add_work(fit)
  inner_cert <- fit$certificate
  inner_status <- inner_cert$target_completeness %||% "not_checked"
  m2 <- length(fit$values)
  expected <- if (interval) m2 %/% 2L else k
  warnings <- character()
  if (is.null(fit$vectors) || !m2) {
    values <- numeric()
    vectors <- matrix(0i, n, 0L)
    mapped <- list(rank = 0L, gap = Inf, operator_columns = 0L, metric_columns = 0L)
  } else {
    mapped <- complex_hermitian_from_realified(problem, fit$vectors, expected)
    values <- mapped$values
    vectors <- mapped$vectors
  }
  cert <- if (isTRUE(execution$certify)) {
    complex_hermitian_certify(problem, values, vectors, tol)
  } else {
    NULL
  }
  polished <- if (length(values) == expected) {
    complex_hermitian_polish(problem, values, vectors, cert, tol)
  } else {
    NULL
  }
  if (!is.null(polished)) {
    values <- polished$values
    vectors <- polished$vectors
    cert <- polished$certificate
  }
  # The real verdict transfers when the real set was the doubled complex
  # set: an even count whose complex span has exactly that rank.
  consistent <- (m2 %% 2L == 0L) && identical(as.integer(mapped$rank), as.integer(m2 %/% 2L)) &&
    length(values) == expected
  record <- inner_cert$completeness %||% list()
  record$realified <- TRUE
  record$embedding <- "real_2n"
  record$inner_status <- inner_status
  record$complex_rank <- as.integer(mapped$rank)
  record$realified_pairs <- m2
  if (!is.null(polished)) {
    record$polished <- polished$method
    record$polish_rounds <- polished$rounds
  }
  status <- if (consistent) inner_status else "not_checked"
  if (!consistent && m2) {
    warnings <- c(warnings, sprintf(paste0(
      "realified solve returned %d real pairs spanning a complex subspace of ",
      "rank %d; its completeness verdict (%s) does not transfer"),
      m2, as.integer(mapped$rank), inner_status))
  }
  if (!is.null(cert)) {
    cert$notes <- unique(c(cert$notes, paste0(
      "complex Hermitian problem solved through its real 2n embedding (",
      inner_plan$method, "); pairs mapped back by complex Rayleigh-Ritz and ",
      "certified on the complex operator")))
    cert$target_completeness <- NULL
    cert <- certificate_with_completeness(cert, status, record)
  }
  extra_cols <- as.integer(mapped$operator_columns %||% 0L) + length(values)
  iter <- list(
    iterations = fit$iterations %||% 1L,
    matvecs = fit$matvecs %||% 0L,
    operator_block_calls = as.integer((fit$work$operator_block_calls %||% 0L) + 2L),
    operator_columns = as.integer((fit$work$operator_columns %||% 0L) + extra_cols),
    certification_operator_columns = as.integer(
      (fit$work$certification_operator_columns %||% 0L) + length(values))
  )
  inner_warnings <- fit$warnings %||% character()
  result <- make_eigen_result(
    values = values,
    vectors = vectors,
    certificate = cert %||% empty_certificate(tol, "certification disabled"),
    iter = iter,
    requested = if (interval) length(values) else k,
    method_label = plan$method,
    target_label_value = target_label(problem$target),
    plan = plan,
    warnings = unique(c(warnings, inner_warnings)),
    extras = list(
      operator_block_calls = iter$operator_block_calls,
      operator_columns = iter$operator_columns,
      certification_operator_columns = iter$certification_operator_columns,
      restart = list(
        kind = "realified_complex_hermitian",
        inner_method = inner_plan$method,
        inner_requested = inner_plan$requested,
        inner_returned = m2,
        inner_completeness = inner_status,
        complex_rank = as.integer(mapped$rank),
        rank_gap = mapped$gap,
        seconds = proc.time()[["elapsed"]] - started
      ),
      locked = if (is.null(cert)) integer() else which(cert$converged)
    )
  )
  if (is.null(cert)) {
    result$certificate <- NULL
  }
  if ((!consistent || !is.null(polished)) && isTRUE(cert$passed) &&
      length(values) == expected &&
      !identical(target_completeness_mode(plan$method_descriptor), "none")) {
    # Re-check the mapped complex set directly.
    checked <- complex_hermitian_target_completeness(
      result, plan, problem, length(values),
      target_completeness_mode(plan$method_descriptor), TRUE,
      solve_seconds = proc.time()[["elapsed"]] - started
    )
    if (!is.null(checked)) {
      result <- checked
    }
  }
  if (!isTRUE(execution$vectors)) {
    result["vectors"] <- list(NULL)
  }
  result
}

# ---------------------------------------------------------------------------
# Completeness of complex results (hooked from hermitian_target_completeness)
# ---------------------------------------------------------------------------

#' @keywords internal
complex_hermitian_target_completeness <- function(result, plan, problem, k, mode,
                                                  vectors_requested,
                                                  solve_seconds = NA_real_,
                                                  depth = 0L) {
  cert <- result$certificate
  values <- result$values
  V <- result$vectors
  if (is.null(cert) || is.null(V) || !length(values) || length(values) != k) {
    return(NULL)
  }
  P2 <- complex_hermitian_realified_problem(problem)
  U <- complex_hermitian_realify_vectors(as.matrix(V) + 0i)
  cert2 <- cert
  cert2$target_completeness <- NULL
  cert2$completeness <- NULL
  for (field in c("residuals", "backward_error", "converged", "scale")) {
    x <- cert[[field]]
    if (is.atomic(x) && length(x) == k) {
      cert2[[field]] <- rep(x, each = 2L)
    }
  }
  r2 <- result
  r2$values <- rep(as.numeric(Re(values)), each = 2L)
  r2$vectors <- U
  r2$certificate <- cert2
  plan2 <- plan
  plan2$problem <- P2
  plan2$requested <- 2L * as.integer(k)
  out2 <- hermitian_target_completeness(
    r2, plan2, P2, 2L * as.integer(k), mode, TRUE,
    solve_seconds = solve_seconds, seed = NULL
  )
  if (is.null(out2) || is.null(out2$certificate)) {
    return(NULL)
  }
  status <- out2$certificate$target_completeness %||% "not_checked"
  record <- out2$certificate$completeness %||% list()
  record$realified <- TRUE
  record$embedding <- "real_2n"
  extra_columns <- as.integer(out2$operator_columns %||% 0L) -
    as.integer(r2$operator_columns %||% 0L)
  changed <- !identical(out2$values, r2$values)
  if (changed && !is.null(out2$vectors)) {
    # The real set was repaired: map it back and re-certify on the complex
    # operator, then verify the mapped set once more.
    mapped <- complex_hermitian_from_realified(problem, out2$vectors, k)
    tol <- cert$tolerance %||% plan$execution$tol
    new_cert <- complex_hermitian_certify(problem, mapped$values, mapped$vectors, tol)
    polished <- complex_hermitian_polish(problem, mapped$values, mapped$vectors,
                                         new_cert, tol)
    if (!is.null(polished)) {
      mapped$values <- polished$values
      mapped$vectors <- polished$vectors
      new_cert <- polished$certificate
    }
    new_cert$notes <- unique(c(cert$notes, new_cert$notes,
      "target completeness check repaired the realified set; result mapped back by complex Rayleigh-Ritz"))
    result$values <- mapped$values
    result$vectors <- mapped$vectors
    result$certificate <- new_cert
    result$residuals <- new_cert$residuals
    result$backward_error <- new_cert$backward_error
    result$orthogonality <- new_cert$orthogonality
    result$nconv <- sum(new_cert$converged)
    result$locked <- which(new_cert$converged)
    result$operator_columns <- as.integer((result$operator_columns %||% 0L) +
                                            max(extra_columns, 0L) + mapped$rank + k)
    if (depth < 1L && isTRUE(new_cert$passed) && length(mapped$values) == k) {
      again <- complex_hermitian_target_completeness(
        result, plan, problem, k, mode, vectors_requested,
        solve_seconds = solve_seconds, depth = depth + 1L
      )
      if (!is.null(again)) {
        again$certificate$notes <- unique(c(again$certificate$notes, new_cert$notes))
        again$warnings <- unique(c(again$warnings,
          "target completeness check found a missing eigenvalue; result repaired by a deflated complement solve on the real embedding"))
        return(again)
      }
    }
    status <- if (isTRUE(new_cert$passed)) "not_checked" else status
    record$repaired_unverified <- TRUE
  } else {
    result$operator_columns <- as.integer((result$operator_columns %||% 0L) +
                                            max(extra_columns, 0L))
  }
  result$certificate <- certificate_with_completeness(result$certificate, status, record)
  if (identical(status, "failed") || identical(status, "inertia_failed")) {
    result$warnings <- unique(c(result$warnings, paste0(
      "target completeness check on the real embedding found a more-preferred ",
      "eigenvalue missing from the returned set; certificate withheld")))
  }
  if (!isTRUE(vectors_requested)) {
    result["vectors"] <- list(NULL)
  }
  result
}

# ---------------------------------------------------------------------------
# Public constructor: complex operators from real and imaginary parts
# ---------------------------------------------------------------------------

#' Complex operator from its real and imaginary parts.
#'
#' Builds an eigencore operator for the complex matrix `re + 1i * im` from
#' two real matrices of the same dimensions, each dense or sparse (any
#' `Matrix` class). This is how complex *sparse* matrices enter eigencore:
#' the `Matrix` package has no complex sparse classes (`zgCMatrix` and
#' friends are not implemented), so a complex sparse Hermitian matrix is
#' passed as a symmetric real part and a skew-symmetric imaginary part.
#'
#' A Hermitian operator (`re` symmetric, `im` skew-symmetric) is solved by
#' [eig_partial()] on the real symmetric `2n x 2n` embedding
#' `[re -im; im re]`, kept sparse: native block Lanczos, CHOLMOD LDL'
#' shift-invert and inertia counts all run on the real sparse embedding, and
#' the eigenpairs are mapped back to complex vectors and certified on the
#' complex operator.
#'
#' @param re Real part: a numeric matrix or a real `Matrix` object.
#' @param im Imaginary part with the dimensions of `re`; `NULL` means zero.
#' @param structure Optional structure descriptor. By default the operator is
#'   [hermitian()] when `re` is symmetric and `im` skew-symmetric (to a
#'   relative tolerance of `sqrt(.Machine$double.eps)`), else [general()].
#' @return An `eigencore_operator` of `dtype = "complex"` whose apply takes
#'   real or complex blocks.
#' @examples
#' re <- Matrix::sparseMatrix(i = c(1, 2, 3, 1), j = c(1, 2, 3, 2),
#'                            x = c(2, 3, 4, 1), dims = c(3, 3), symmetric = TRUE)
#' im <- Matrix::sparseMatrix(i = c(1, 2), j = c(2, 1), x = c(0.5, -0.5),
#'                            dims = c(3, 3))
#' op <- complex_operator(re, im)
#' fit <- eig_partial(op, k = 1, target = largest())
#' values(fit)
#' @export
complex_operator <- function(re, im = NULL, structure = NULL) {
  as_real <- function(x, name) {
    if (inherits(x, "Matrix")) {
      if (inherits(x, "sparseMatrix")) {
        x <- methods::as(methods::as(methods::as(x, "CsparseMatrix"), "generalMatrix"),
                         "dMatrix")
        if (!all(is.finite(methods::slot(x, "x")))) {
          stop(name, " contains NA, NaN, or Inf entries.", call. = FALSE)
        }
        return(x)
      }
      x <- base::as.matrix(x)
    }
    if (!is.matrix(x) || !(is.numeric(x) || is.logical(x)) || is.complex(x)) {
      stop(name, " must be a real numeric matrix or a real Matrix object.", call. = FALSE)
    }
    storage.mode(x) <- "double"
    if (!all(is.finite(x))) {
      stop(name, " contains NA, NaN, or Inf entries.", call. = FALSE)
    }
    x
  }
  re <- as_real(re, "re")
  im <- if (is.null(im)) {
    if (inherits(re, "sparseMatrix")) {
      Matrix::sparseMatrix(i = integer(), j = integer(), x = numeric(), dims = dim(re))
    } else {
      matrix(0, nrow(re), ncol(re))
    }
  } else {
    as_real(im, "im")
  }
  if (!identical(as.integer(dim(re)), as.integer(dim(im)))) {
    stop("re and im must have the same dimensions.", call. = FALSE)
  }
  if (inherits(re, "sparseMatrix") != inherits(im, "sparseMatrix")) {
    if (inherits(re, "sparseMatrix")) {
      im <- methods::as(methods::as(Matrix::Matrix(im, sparse = TRUE), "generalMatrix"), "dMatrix")
    } else {
      re <- methods::as(methods::as(Matrix::Matrix(re, sparse = TRUE), "generalMatrix"), "dMatrix")
    }
  }
  sparse <- inherits(re, "sparseMatrix")
  d <- as.integer(dim(re))
  if (is.null(structure)) {
    hermitian_parts <- d[[1L]] == d[[2L]] && {
      scale <- max(abs(if (sparse) methods::slot(re, "x") else re), 0,
                   abs(if (sparse) methods::slot(im, "x") else im))
      tol <- sqrt(.Machine$double.eps) * max(scale, .Machine$double.xmin)
      max(abs(re - Matrix::t(re)), 0) <= tol && max(abs(im + Matrix::t(im)), 0) <= tol
    }
    structure <- if (isTRUE(hermitian_parts)) hermitian() else general()
  }
  apply_parts <- function(X, alpha, beta, Y, adjoint) {
    X <- as.matrix(X)
    sgn <- if (adjoint) -1 else 1
    mult <- function(M, V) {
      if (adjoint) as.matrix(Matrix::crossprod(M, V)) else as.matrix(M %*% V)
    }
    Xr <- Re(X)
    Xi <- Im(X)
    out <- (mult(re, Xr) - sgn * mult(im, Xi)) +
      1i * (mult(re, Xi) + sgn * mult(im, Xr))
    out <- alpha * out
    if (is.null(Y) || beta == 0) out else out + beta * Y
  }
  fro <- sqrt(sum((if (sparse) methods::slot(re, "x") else re)^2) +
                sum((if (sparse) methods::slot(im, "x") else im)^2))
  linear_operator(
    dim = d,
    apply = function(X, alpha = 1, beta = 0, Y = NULL) {
      apply_parts(X, alpha, beta, Y, adjoint = FALSE)
    },
    apply_adjoint = function(X, alpha = 1, beta = 0, Y = NULL) {
      apply_parts(X, alpha, beta, Y, adjoint = TRUE)
    },
    dtype = "complex",
    structure = structure,
    name = if (sparse) "complex_sparse_split" else "complex_dense_split",
    metadata = list(
      complex_split = list(re = re, im = im),
      storage = if (sparse) "complex_split_sparse" else "complex_split_dense",
      frobenius_norm = fro
    )
  )
}
