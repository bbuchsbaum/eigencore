# Target completeness for partial SVD results.
#
# A residual certificate proves each returned triplet (sigma, u, v) is a
# singular triplet to the stated backward error. `certificate$passed` also
# requires the returned SET to be the requested one (R/solve.R,
# require_verified_completeness()). This file establishes that for SVD
# results, route by route:
#
# * Dense LAPACK routes (full or partial SVD of an explicit matrix) compute
#   every singular value and select the requested ones: "exact".
# * Gram routes whose Gram matrix was formed and decomposed by a dense
#   LAPACK eigensolver (all eigenvalues, or a selected index range by
#   bisection) are exact modulo that eigensolver: "exact". A Gram route that
#   ran a Krylov / subspace eigensolver on the Gram matrix is treated like a
#   Krylov route.
# * Krylov routes (Golub-Kahan / IRLBA, retained and block Golub-Kahan,
#   randomized, implicit Gram, matrix-free) and Gram routes that fell back to
#   Golub-Kahan are checked by
#   1. an inertia (eigenvalue-counting) proof on the augmented symmetric
#      matrix K = [0 A; A' 0] when A has an explicit matrix source and the
#      factorisation is affordable (completeness mode "auto"/"inertia");
#      the eigenvalues of K are +/- sigma_i plus |m - n| zeros, so the
#      inertia of K - s I counts the singular values above s, and K is
#      sparse when A is (CHOLMOD LDL', dense Bunch-Kaufman otherwise); and
#   2. otherwise (matrix-free, too expensive, or an inconclusive count) the
#      deflated-complement probe of R/target_completeness.R applied to the
#      Gram operator A'A (or AA', always on the smaller side, so it has
#      exactly min(m, n) eigenvalues sigma_i^2 and no spurious zeros) with
#      the returned right (or left) singular vectors. A probe that finds an
#      intruder repairs the set by the same deflated complement solve and
#      Rayleigh-Ritz step, followed by a Rayleigh-Ritz step for A on the
#      repaired right (left) subspace; the repaired triplets are
#      re-certified from scratch in original coordinates.
#
# The probe is evidence, not proof (see the header of
# R/target_completeness.R). Because it works on the Gram operator its
# resolution in sigma is about margin_gram / (2 sigma_edge): sharp for
# largest targets, weak for tiny smallest singular values. nearest() SVD
# targets have no probe (the Gram compression does not localise interior
# singular values) and are verified only by dense routes or inertia.

#' @keywords internal
svd_completeness_kind <- function(target) {
  kind <- if (inherits(target, "eigencore_target")) target$kind else "largest"
  switch(
    kind,
    largest = ,
    largest_magnitude = "largest",
    smallest = ,
    smallest_magnitude = "smallest",
    nearest = if (is.numeric(target$value) && length(target$value) == 1L &&
                  is.finite(target$value)) "nearest" else NULL,
    NULL
  )
}

# Route class of an SVD result: "exact" (complete by construction), "krylov"
# (inertia / probe), or "unsupported".
#' @keywords internal
svd_completeness_route_class <- function(result, plan) {
  restart <- result$restart %||% list()
  if (svd_plan_dispatches_dense(plan)) {
    return("exact")
  }
  if (identical(restart$kind %||% "", "gram_svd_special_case") &&
      !isTRUE(restart$fallback_used) &&
      grepl("^lapack_", restart$native_gram_eigensolver %||% "")) {
    return("exact")
  }
  if (!identical(plan$problem$A$dtype %||% "double", "double")) {
    return("unsupported")
  }
  "krylov"
}

# Mirrors execute_svd_plan_dispatch(): plans that run solve_svd_dense(), a
# dense LAPACK SVD computing every singular value.
#' @keywords internal
svd_plan_dispatches_dense <- function(plan) {
  method <- plan$method %||% ""
  !(method %in% c(
    "reference randomized SVD prototype",
    native_dense_randomized_svd_label(),
    native_csc_randomized_svd_label(),
    "native certified Gram SVD special case",
    native_implicit_gram_svd_label(),
    native_retained_golub_kahan_diagnostic_label()
  ) || plan_dispatches_golub_kahan(plan))
}

# Lazy Gram operator of `op` on one side: "right" is A'A (n x n), "left" is
# AA' (m x m). Never materialised.
#' @keywords internal
svd_completeness_gram_operator <- function(op, side) {
  n <- if (identical(side, "right")) op$dim[[2L]] else op$dim[[1L]]
  gram_apply <- function(X, alpha = 1, beta = 0, Y = NULL) {
    X <- as.matrix(X)
    # Force the inner apply before entering the outer one (a lazily
    # evaluated argument would run nested and go uncounted).
    Z <- if (identical(side, "right")) {
      AX <- as.matrix(apply_operator(op, X))
      apply_adjoint_operator(op, AX)
    } else {
      AtX <- as.matrix(apply_adjoint_operator(op, X))
      apply_operator(op, AtX)
    }
    Z <- alpha * as.matrix(Z)
    if (is.null(Y) || beta == 0) Z else Z + beta * Y
  }
  G <- linear_operator(
    dim = c(n, n),
    apply = gram_apply,
    structure = hermitian(),
    name = "eigencore SVD completeness Gram operator"
  )
  # Unwrapped apply: the work record then sees the forward and adjoint
  # applies of A themselves (a nested apply is not counted).
  G$apply <- gram_apply
  G
}

#' @keywords internal
svd_completeness_target <- function(kind) {
  if (identical(kind, "largest")) largest() else smallest()
}

# Rayleigh-Ritz for A on span(W) (W orthonormal, on `side`): the SVD of A W
# (right) or A'W (left) gives triplets whose singular values are the best
# approximations from that subspace.
#' @keywords internal
svd_completeness_rayleigh_ritz <- function(op, W, side, target) {
  W <- as.matrix(W)
  AW <- if (identical(side, "right")) {
    as.matrix(apply_operator(op, W))
  } else {
    as.matrix(apply_adjoint_operator(op, W))
  }
  s <- svd(AW)
  idx <- order_indices(s$d, target)
  idx <- idx[seq_len(min(ncol(W), length(idx)))]
  d <- s$d[idx]
  P <- s$u[, idx, drop = FALSE]
  Q <- W %*% s$v[, idx, drop = FALSE]
  if (identical(side, "right")) {
    list(d = d, u = P, v = Q)
  } else {
    list(d = d, u = Q, v = P)
  }
}

#' @keywords internal
svd_completeness_norm <- function(cert, d) {
  scale <- suppressWarnings(max(c(abs(as.numeric(cert$scale)), abs(d), 0),
                                na.rm = TRUE))
  if (!is.finite(scale) || scale <= 0) 1 else scale
}

# Gram probe (and repair) of a certified SVD result. Returns a check object:
# list(status, record, repaired, d, u, v).
#' @keywords internal
svd_completeness_probe_check <- function(op, d, u, v, target, cert, tol,
                                         controls = target_completeness_controls()) {
  kind <- svd_completeness_kind(target)
  # The smaller side's Gram has exactly min(m, n) eigenvalues sigma^2; the
  # larger side adds |m - n| zeros, harmless only for largest targets.
  sw <- svd_completeness_side(op, u, v, kind)
  if (is.null(sw) || ncol(sw$W) >= nrow(sw$W)) {
    return(NULL)
  }
  side <- sw$side
  W <- sw$W
  normA <- svd_completeness_norm(cert, d)
  left <- as.numeric(cert$residuals$left %||% rep(0, length(d)))
  right <- as.numeric(cert$residuals$right %||% rep(0, length(d)))
  # ||A'A v - s^2 v|| <= ||A|| ||A v - s u|| + s ||A'u - s v|| (right side);
  # sides swap for AA'.
  gram_residuals <- if (identical(side, "right")) {
    normA * left + abs(d) * right
  } else {
    normA * right + abs(d) * left
  }
  edge <- if (identical(kind, "largest")) min(d) else max(d)
  # Relative tolerance on the Gram scale ||A||^2 matching tol * ||A|| on
  # the singular values at the edge.
  tol_gram <- max(tol * (2 * edge / normA + tol), 64 * .Machine$double.eps)
  G <- svd_completeness_gram_operator(op, side)
  gram_target <- svd_completeness_target(kind)
  check <- target_completeness_check(
    G, d^2, W, gram_target, tol = tol_gram,
    residuals = gram_residuals, norm_scale = normA^2, controls = controls
  )
  record <- check$record
  record$method <- "gram_deflated_complement_probe"
  record$gram_side <- side
  record$edge_singular_value <- edge
  # Every Gram column is one apply of A and one of A'.
  record$operator_columns <- as.integer(record$operator_columns %||% 0L)
  record$operator_block_calls <- as.integer(record$operator_block_calls %||% 0L)
  record$adjoint_columns <- record$operator_columns
  record$adjoint_block_calls <- record$operator_block_calls
  out <- list(status = check$status, record = record, repaired = FALSE,
              d = d, u = u, v = v)
  if (isTRUE(check$repaired)) {
    W2 <- completeness_orthonormalize(check$vectors)
    if (ncol(W2) < length(d)) {
      out$status <- "failed"
      return(out)
    }
    rr <- svd_completeness_rayleigh_ritz(op, W2, side, target)
    out$d <- rr$d
    out$u <- rr$u
    out$v <- rr$v
    out$repaired <- TRUE
    out$record$operator_columns <- out$record$operator_columns +
      if (identical(side, "right")) ncol(W2) else 0L
    out$record$adjoint_columns <- out$record$adjoint_columns +
      if (identical(side, "left")) ncol(W2) else 0L
  }
  out
}

# Side and basis for a Gram computation on a returned SVD: the smaller side
# when its vectors are present, else the larger side for largest targets.
#' @keywords internal
svd_completeness_side <- function(op, u, v, kind) {
  m <- op$dim[[1L]]
  n <- op$dim[[2L]]
  side <- if (n <= m) "right" else "left"
  if (identical(side, "right") && is.null(v) || identical(side, "left") && is.null(u)) {
    if (!identical(kind, "largest")) {
      return(NULL)
    }
    side <- if (identical(side, "right")) "left" else "right"
  }
  W <- if (identical(side, "right")) v else u
  if (is.null(W)) NULL else list(side = side, W = as.matrix(W))
}

# A short result (fewer triplets than requested, typically a Krylov space
# exhausted by repeated singular values) is completed from the deflated
# complement of the Gram operator: the `deficit` most preferred eigenvectors
# of P G P (P = I - W W'; the returned directions moved to the
# non-preferred side), then a Rayleigh-Ritz step for A on span(W, X_c). Small
# complements are decomposed densely, larger ones by block Lanczos. Returns
# list(d, u, v, columns) or NULL.
#' @keywords internal
svd_completeness_fill <- function(op, d, u, v, target, k_req, tol) {
  kind <- svd_completeness_kind(target)
  if (is.null(kind) || identical(kind, "nearest")) {
    return(NULL)
  }
  sw <- svd_completeness_side(op, u, v, kind)
  if (is.null(sw)) {
    return(NULL)
  }
  W <- completeness_orthonormalize(sw$W)
  if (ncol(W) < length(d)) {
    return(NULL)
  }
  N <- nrow(W)
  deficit <- as.integer(k_req) - ncol(W)
  nc <- N - ncol(W)
  if (deficit < 1L || nc < deficit) {
    return(NULL)
  }
  G <- svd_completeness_gram_operator(op, sw$side)
  gram_target <- svd_completeness_target(kind)
  columns <- 0L
  if (nc <= 64L) {
    Qc <- completeness_orthonormalize(diag(N), W)
    Qc <- Qc[, seq_len(min(nc, ncol(Qc))), drop = FALSE]
    GQ <- as.matrix(apply_operator(G, Qc))
    columns <- ncol(Qc)
    H <- crossprod(Qc, GQ)
    eig <- eigen((H + t(H)) / 2, symmetric = TRUE)
    sel <- utils::head(order_indices(eig$values, gram_target), deficit)
    Xc <- Qc %*% eig$vectors[, sel, drop = FALSE]
  } else {
    scale <- max(abs(d), 1)^2
    shift <- completeness_shift(d^2, kind, scale)
    deflated_apply <- function(X, alpha = 1, beta = 0, Y = NULL) {
      X <- as.matrix(X)
      C <- crossprod(W, X)
      PX <- X - W %*% C
      GPX <- as.matrix(apply_operator(G, PX))
      Z <- GPX - W %*% crossprod(W, GPX) + shift * (W %*% C)
      Z <- alpha * Z
      if (is.null(Y) || beta == 0) Z else Z + beta * Y
    }
    Gop <- linear_operator(dim = c(N, N), apply = deflated_apply,
                           structure = hermitian(),
                           name = "eigencore SVD completeness deflated Gram operator")
    Gop$apply <- deflated_apply
    block <- min(max(2L, deficit + 1L), 8L, nc)
    m_max <- min(N, default_block_lanczos_max_subspace(deficit, block))
    start <- completeness_probe_start(N, block, 1201L)
    comp <- native_block_lanczos_hermitian(
      Gop, k = deficit, target = gram_target, tol = tol, maxit = m_max,
      block = block, max_restarts = 100L, vectors = TRUE,
      full_subspace = FALSE, certificate_fallback = FALSE, start = start
    )
    Xc <- comp$vectors
    columns <- as.integer(comp$operator_columns %||% comp$matvecs %||% 0L)
  }
  Z <- cbind(W, completeness_orthonormalize(Xc, W))
  if (ncol(Z) < k_req) {
    return(NULL)
  }
  rr <- svd_completeness_rayleigh_ritz(op, Z, sw$side, target)
  sel <- seq_len(k_req)
  list(d = rr$d[sel], u = rr$u[, sel, drop = FALSE], v = rr$v[, sel, drop = FALSE],
       columns = columns + ncol(Z), side = sw$side, deficit = deficit)
}

# ---------------------------------------------------------------------------
# Inertia proof on the augmented matrix K = [0 A; A' 0]
# ---------------------------------------------------------------------------

# Explicit (real) matrix behind an SVD operator, or NULL.
#' @keywords internal
svd_inertia_matrix <- function(op) {
  if (!inherits(op, "eigencore_operator") ||
      !identical(op$dtype %||% "double", "double")) {
    return(NULL)
  }
  src <- inertia_matrix_of(op)
  if (is.null(src)) {
    return(NULL)
  }
  if (inherits(src, "sparseMatrix") || inherits(src, "diagonalMatrix")) {
    if (!inherits(src, "dsparseMatrix") && !inherits(src, "diagonalMatrix")) {
      return(NULL)
    }
    return(src)
  }
  if (is.matrix(src) && (is.double(src) || is.integer(src) || is.logical(src))) {
    return(src)
  }
  NULL
}

# Build the augmented symmetric matrix (upper triangle for sparse input).
#' @keywords internal
svd_augmented_matrix <- function(A) {
  m <- nrow(A)
  n <- ncol(A)
  N <- m + n
  if (is.matrix(A)) {
    K <- matrix(0, N, N)
    K[seq_len(m), m + seq_len(n)] <- A
    K[m + seq_len(n), seq_len(m)] <- t(A)
    return(K)
  }
  T <- methods::as(methods::as(A, "CsparseMatrix"), "TsparseMatrix")
  i <- methods::slot(T, "i") + 1L
  j <- methods::slot(T, "j") + 1L
  x <- as.numeric(methods::slot(T, "x"))
  Matrix::sparseMatrix(i = i, j = m + j, x = x, dims = c(N, N),
                       symmetric = TRUE, repr = "C")
}

# Cost gate for the augmented inertia count in auto mode. The Gram probe
# that it replaces costs a few dozen applies of A, so the factorisation must
# stay small next to the solve: it runs when its predicted time is below
# eigencore.svd_completeness_inertia_seconds (default 0.25 s) or below
# eigencore.svd_completeness_inertia_ratio (default 0.1) times the solve.
# The prediction counts two factorisations per threshold (the margin
# re-count is common).
#' @keywords internal
svd_inertia_completeness_controls <- function() {
  controls <- inertia_completeness_controls()
  num <- function(value, default) {
    value <- suppressWarnings(as.numeric(value %||% default))
    if (length(value) != 1L || !is.finite(value) || value < 0) default else value
  }
  controls$seconds <- num(getOption("eigencore.svd_completeness_inertia_seconds"), 0.25)
  controls$ratio <- num(getOption("eigencore.svd_completeness_inertia_ratio"), 0.1)
  controls
}

# Predicted factorisation cost before building a dense context.
#' @keywords internal
svd_inertia_gate <- function(op, kind, mode, solve_seconds,
                             controls = svd_inertia_completeness_controls()) {
  no <- function(reason) list(use = FALSE, reason = reason, ctx = NULL,
                              predicted_seconds = NA_real_)
  if (!mode %in% c("auto", "inertia")) {
    return(no("mode"))
  }
  A <- svd_inertia_matrix(op)
  if (is.null(A)) {
    return(no("no explicit matrix source"))
  }
  N <- as.numeric(sum(dim(A)))
  per_threshold <- if (identical(kind, "nearest")) 4 else 2
  dense <- is.matrix(A)
  if (dense) {
    predicted <- per_threshold * (N^3 / 3) / controls$dense_rate
    if (N > inertia_dense_limit()) {
      return(no("augmented matrix too large to factor densely"))
    }
  }
  budget <- max(controls$seconds,
                controls$ratio * (if (is.finite(solve_seconds)) solve_seconds else 0))
  if (dense && identical(mode, "auto") && !(predicted <= budget)) {
    return(list(use = FALSE, reason = "predicted factorisation cost exceeds the gate",
                ctx = NULL, predicted_seconds = predicted))
  }
  ctx <- tryCatch(inertia_context(svd_augmented_matrix(A)), error = function(e) e)
  if (inherits(ctx, "error")) {
    return(no(paste0("inertia context: ", conditionMessage(ctx))))
  }
  if (!dense) {
    cost <- operator_memoised_value(op, "svd_inertia_factor_cost",
                                    inertia_factor_cost(ctx))
    rate <- if (identical(ctx$kind, "sparse")) controls$sparse_rate else controls$dense_rate
    predicted <- per_threshold * cost$flops / rate
  }
  if (identical(mode, "inertia")) {
    return(list(use = TRUE, reason = "requested", ctx = ctx,
                predicted_seconds = predicted))
  }
  if (!is.finite(predicted)) {
    return(list(use = FALSE, reason = "factorisation cost unknown", ctx = NULL,
                predicted_seconds = predicted))
  }
  if (predicted <= budget) {
    list(use = TRUE, reason = "cost gate passed", ctx = ctx,
         predicted_seconds = predicted)
  } else {
    list(use = FALSE, reason = "predicted factorisation cost exceeds the gate",
         ctx = NULL, predicted_seconds = predicted)
  }
}

# Number of singular values with preference distance f < t, from inertia of
# K - s I (s > 0: the eigenvalues of K above s are exactly the singular
# values above s). `outward = TRUE` may only grow the counted region,
# `outward = FALSE` only shrink it. Returns list(count, t, reliable,
# factorizations, backward_bound).
#' @keywords internal
svd_inertia_region_count <- function(ctx, kind, t, p, center = 0, outward = TRUE) {
  dir <- if (outward) 1 else -1
  parts <- list()
  res <- function(count, reliable, t_eff) {
    list(
      count = if (isTRUE(reliable)) count else NA_real_,
      t = t_eff,
      reliable = isTRUE(reliable),
      factorizations = sum(vapply(parts, function(r) NROW(r$attempts), numeric(1L))),
      backward_bound = max(c(0, vapply(parts, function(r) {
        b <- r$tally$backward_bound %||% NA_real_
        if (is.finite(b)) b else 0
      }, numeric(1L))))
    )
  }
  # Count of singular values strictly above s (s nudged in `direction`).
  above <- function(s, direction) {
    r <- inertia_count_nudged(ctx, s, direction = direction)
    parts[[length(parts) + 1L]] <<- r
    r
  }
  switch(
    kind,
    largest = {
      a <- -t
      if (!(a > 0)) {
        return(res(p, TRUE, t))
      }
      r <- above(a, -dir)
      res(r$above, r$reliable && r$s > 0, -r$s)
    },
    smallest = {
      if (!(t > 0)) {
        return(res(0, TRUE, t))
      }
      r <- above(t, dir)
      res(p - r$above - r$zero, r$reliable && r$s > 0, r$s)
    },
    nearest = {
      if (!(t > 0)) {
        return(res(0, TRUE, t))
      }
      lo <- center - t
      hi <- above(center + t, dir)
      if (lo > 0) {
        lo_r <- above(lo, -dir)
        ok <- hi$reliable && lo_r$reliable && lo_r$s > 0
        t_eff <- if (outward) max(hi$s - center, center - lo_r$s) else
          min(hi$s - center, center - lo_r$s)
        # Singular values in (lo', hi'): above(lo') - above(hi') - zero(hi').
        res(lo_r$above - hi$above - hi$zero, ok, t_eff)
      } else {
        # The region reaches 0: it is {sigma < hi'}.
        t_eff <- if (outward) hi$s - center else min(t, hi$s - center)
        res(p - hi$above - hi$zero, hi$reliable && hi$s > 0, t_eff)
      }
    }
  )
}

# The counting argument of R/completeness_inertia.R on K with
# x_i = [u_i; v_i] / sqrt(2): ||K x_i - sigma_i x_i|| = combined_i / sqrt(2)
# and max |X'X - I| <= (omega_U + omega_V) / 2. Kahan's theorem then matches
# each sigma_i with a distinct eigenvalue of K within rho; when every
# sigma_i - rho > 0 the matched eigenvalues are distinct singular values.
#' @keywords internal
svd_inertia_completeness_check <- function(ctx, d, cert, target, p) {
  started <- proc.time()[["elapsed"]]
  kind <- svd_completeness_kind(target)
  center <- if (identical(kind, "nearest")) as.numeric(target$value) else 0
  d <- as.numeric(d)
  k <- length(d)
  eps <- .Machine$double.eps
  record <- list(
    method = "augmented_inertia_count", kind = kind, k = k,
    edge = NA_real_, rho = NA_real_, margin = NA_real_,
    threshold_upper = NA_real_, count_upper = NA_real_,
    threshold_lower = NA_real_, count_lower = NA_real_,
    factorizations = 0, inertia_method = ctx$method,
    augmented_dimension = ctx$n,
    reason = NA_character_, seconds = NA_real_
  )
  finish <- function(status, reason = NA_character_) {
    record$reason <- reason
    record$seconds <- proc.time()[["elapsed"]] - started
    list(status = status, record = record)
  }
  if (is.null(kind) || !k || any(!is.finite(d))) {
    return(finish("inertia_inconclusive", "target or values not supported"))
  }
  if (k >= p) {
    return(finish("inertia_verified", "all singular values returned"))
  }
  omega <- suppressWarnings(max(abs(as.numeric(cert$orthogonality)), 0, na.rm = TRUE))
  gram_floor <- 1 - k * omega
  if (!is.finite(gram_floor) || gram_floor < 0.5) {
    return(finish("inertia_inconclusive",
                  "returned vectors are too far from orthonormal for the residual bound"))
  }
  combined <- as.numeric(cert$residuals$combined %||% numeric())
  if (length(combined) != k || any(!is.finite(combined))) {
    return(finish("inertia_inconclusive", "certificate residuals unavailable"))
  }
  dmax <- max(abs(d))
  rounding <- sqrt(k) * 32 * eps * (ctx$normA + dmax) * sqrt(1 + omega)
  rho <- (sqrt(sum(combined^2)) / sqrt(2) + rounding) / sqrt(gram_floor)
  scale <- ctx$normA + abs(center) + dmax
  margin <- max(rho, 1e-10 * scale, 64 * eps * scale)
  record$rho <- rho
  if (min(d) - rho - margin <= 0) {
    return(finish("inertia_inconclusive",
                  "a returned singular value is within the residual bound of zero"))
  }
  f <- inertia_preference(d, kind, center)
  E <- max(f)
  record$edge <- E
  upper <- svd_inertia_region_count(ctx, kind, E + rho + margin, p, center, outward = TRUE)
  record$factorizations <- record$factorizations + upper$factorizations
  if (isTRUE(upper$reliable) && upper$backward_bound >= margin) {
    margin <- 2 * upper$backward_bound
    upper <- svd_inertia_region_count(ctx, kind, E + rho + margin, p, center, outward = TRUE)
    record$factorizations <- record$factorizations + upper$factorizations
  }
  record$margin <- margin
  record$threshold_upper <- upper$t
  record$count_upper <- upper$count
  if (!isTRUE(upper$reliable)) {
    return(finish("inertia_inconclusive", "inertia count above the edge was not reliable"))
  }
  if (upper$count == k) {
    return(finish("inertia_verified"))
  }
  if (upper$count < k) {
    return(finish("inertia_inconclusive",
                  "fewer singular values than returned values inside the residual bound"))
  }
  lower <- svd_inertia_region_count(ctx, kind, E - rho - margin, p, center, outward = FALSE)
  record$factorizations <- record$factorizations + lower$factorizations
  record$threshold_lower <- lower$t
  record$count_lower <- lower$count
  if (!isTRUE(lower$reliable)) {
    return(finish("inertia_inconclusive", "inertia count below the edge was not reliable"))
  }
  if (lower$count >= k) {
    return(finish("inertia_failed",
                  "at least k singular values are strictly more preferred than the returned edge"))
  }
  could_match <- sum(f - rho < lower$t)
  if (lower$count > could_match) {
    return(finish("inertia_failed",
                  "a singular value among the k most preferred has no returned value within the residual bound"))
  }
  finish("inertia_inconclusive",
         "the k-th and (k+1)-th singular values are not separated by more than the residual bound")
}

# ---------------------------------------------------------------------------
# Post-dispatch hook
# ---------------------------------------------------------------------------

# Re-certify repaired triplets and copy them into the result.
#' @keywords internal
svd_result_with_repair <- function(result, op, fixed, tol, vectors_mode) {
  cert_old <- result$certificate
  cert <- certify_svd_operator(op, fixed$d, fixed$u, fixed$v, tol = tol)
  cert$notes <- unique(c(cert_old$notes, cert$notes))
  result$d <- fixed$d
  result$values <- fixed$d
  result$u <- if (vectors_mode %in% c("both", "left")) fixed$u else NULL
  result$v <- if (vectors_mode %in% c("both", "right")) fixed$v else NULL
  result$residuals <- cert$residuals
  result$backward_error <- cert$backward_error
  result$orthogonality <- cert$orthogonality
  result$nconv <- sum(cert$converged)
  result$certificate <- cert
  result
}

#' @keywords internal
svd_result_add_work <- function(result, record, recert_columns = 0L) {
  op_cols <- as.integer(record$operator_columns %||% 0L)
  op_blocks <- as.integer(record$operator_block_calls %||% 0L)
  adj_cols <- as.integer(record$adjoint_columns %||% 0L)
  adj_blocks <- as.integer(record$adjoint_block_calls %||% 0L)
  if (inherits(result$work, "eigencore_work")) {
    work <- unclass(result$work)
    add <- function(field, amount) {
      if (!is.na(work[[field]] %||% NA)) {
        work[[field]] <<- as.integer(work[[field]] + amount)
      }
    }
    add("operator_block_calls", op_blocks)
    add("operator_columns", op_cols)
    add("adjoint_block_calls", adj_blocks)
    add("adjoint_columns", adj_cols)
    if (recert_columns > 0L) {
      add("certification_operator_block_calls", 1L)
      add("certification_operator_columns", recert_columns)
      add("certification_adjoint_block_calls", 1L)
      add("certification_adjoint_columns", recert_columns)
    }
    result$work <- new_typed_work_record(work)
  }
  result
}

#' @keywords internal
svd_completeness_warning <- function(status, record) {
  switch(
    status,
    failed = paste0(
      "SVD target completeness probe found a singular value (",
      format(sqrt(max(record$most_preferred_complement %||% NA_real_, 0)), digits = 8),
      ") more preferred than the returned edge outside the returned set; certificate withheld"
    ),
    inertia_failed = paste0(
      "inertia count proves a singular value more preferred than the returned edge is ",
      "missing (", record$reason %||% "count mismatch", "); certificate withheld"
    ),
    repaired = "SVD target completeness check found a missing singular value copy; result repaired by a deflated complement solve",
    NULL
  )
}

# Attach a target-completeness verdict to an SVD result (C50 for SVD).
#' @keywords internal
apply_svd_completeness <- function(result, plan, solve_seconds = NA_real_) {
  cert <- result$certificate
  if (is.null(cert) || !is.null(cert$target_completeness) ||
      !inherits(result, "eigencore_svd_result")) {
    return(result)
  }
  problem <- plan$problem
  op <- problem$A
  route <- svd_completeness_route_class(result, plan)
  if (identical(route, "exact")) {
    result$certificate <- certificate_with_completeness(cert, "exact")
    return(result)
  }
  mode <- target_completeness_mode(plan$method_descriptor)
  kind <- svd_completeness_kind(problem$target)
  d <- result$d
  k <- length(d)
  p <- if (inherits(op, "eigencore_operator")) min(op$dim) else NA_integer_
  eligible <- identical(route, "krylov") && !identical(mode, "none") &&
    isTRUE(plan$execution$certify) && isTRUE(cert$passed) && !is.null(kind) &&
    inherits(op, "eigencore_operator") &&
    is.numeric(d) && !is.complex(d) && k >= 1L && all(is.finite(d))
  if (!eligible) {
    result$certificate <- certificate_with_completeness(cert, "not_checked")
    return(result)
  }
  tol <- cert$tolerance %||% plan$execution$tol
  vectors_mode <- plan$execution$vectors %||% "both"
  started <- proc.time()[["elapsed"]]
  k_req <- as.integer(plan$requested %||% k)
  filled <- 0L
  if (k < k_req && k_req <= p) {
    # Short result: complete it from the deflated Gram complement; the
    # completed set is kept only when it certifies.
    fill <- tryCatch(svd_completeness_fill(op, d, result$u, result$v,
                                           problem$target, k_req, tol),
                     error = function(e) NULL)
    if (!is.null(fill)) {
      candidate <- svd_result_with_repair(result, op, fill, tol, vectors_mode)
      fill_cols <- as.integer(fill$columns)
      candidate <- svd_result_add_work(
        candidate,
        list(operator_columns = fill_cols, adjoint_columns = fill_cols),
        recert_columns = length(fill$d)
      )
      if (isTRUE(candidate$certificate$passed)) {
        result <- candidate
        result$warnings <- c(
          result$warnings,
          sprintf(paste0("the solver returned %d of %d singular triplets; the ",
                         "set was completed by a deflated Gram complement solve"),
                  k, k_req)
        )
        result$certificate$notes <- unique(c(
          result$certificate$notes,
          "short result completed by a deflated Gram complement solve and re-certified"
        ))
        cert <- result$certificate
        filled <- as.integer(fill$deficit)
        d <- result$d
        k <- length(d)
      }
    }
  }
  if (k >= p) {
    # Every singular value was returned (certified, orthonormal triplets).
    result$certificate <- certificate_with_completeness(
      cert, "exact", list(method = "full_set", k = k, filled = filled))
    return(result)
  }
  status <- NULL
  record <- NULL
  gate <- NULL
  inertia_record <- NULL
  if (mode %in% c("auto", "inertia")) {
    gate <- tryCatch(svd_inertia_gate(op, kind, mode, solve_seconds),
                     error = function(e) list(use = FALSE, reason = conditionMessage(e)))
    if (isTRUE(gate$use)) {
      ic <- tryCatch(svd_inertia_completeness_check(gate$ctx, d, cert, problem$target, p),
                     error = function(e) list(status = "inertia_inconclusive",
                                              record = list(reason = conditionMessage(e))))
      inertia_record <- ic$record
      inertia_record$gate <- gate$reason
      inertia_record$predicted_seconds <- gate$predicted_seconds
      if (identical(ic$status, "inertia_verified")) {
        status <- "inertia_verified"
        record <- inertia_record
      } else if (identical(ic$status, "inertia_failed")) {
        status <- "inertia_failed"
        record <- inertia_record
      }
    }
  }
  if (is.null(status) || identical(status, "inertia_failed")) {
    probe <- if (identical(kind, "nearest")) NULL else
      tryCatch(svd_completeness_probe_check(op, d, result$u, result$v,
                                            problem$target, cert, tol),
               error = function(e) NULL)
    if (!is.null(probe)) {
      if (isTRUE(probe$repaired)) {
        result <- svd_result_with_repair(result, op, probe, tol, vectors_mode)
        cert <- result$certificate
        result <- svd_result_add_work(result, list(), recert_columns = length(probe$d))
      }
      result <- svd_result_add_work(result, probe$record)
      if (identical(status, "inertia_failed") && isTRUE(probe$repaired) &&
          isTRUE(gate$use)) {
        again <- tryCatch(
          svd_inertia_completeness_check(gate$ctx, result$d, cert, problem$target, p),
          error = function(e) list(status = "inertia_inconclusive", record = list()))
        status <- if (identical(again$status, "inertia_verified")) "inertia_verified" else probe$status
        record <- again$record
        record$repaired <- TRUE
        record$probe <- probe$record
      } else if (identical(status, "inertia_failed") && !identical(probe$status, "failed") &&
                 !isTRUE(probe$repaired)) {
        # The count proves a value is missing; the probe did not see it.
        record <- inertia_record
        record$probe <- probe$record
      } else {
        status <- probe$status
        record <- probe$record
        if (!is.null(inertia_record)) {
          record$inertia <- inertia_record
        }
      }
    } else if (is.null(status)) {
      if (is.null(inertia_record)) {
        status <- "not_checked"
        record <- list(
          reason = if (identical(kind, "nearest"))
            "nearest SVD targets are verified only by inertia counting" else
              "no singular vectors available for the probe",
          inertia_gate = gate$reason %||% NA_character_)
      } else {
        # The count ran but could not decide, and no probe applies.
        status <- "inertia_inconclusive"
        record <- inertia_record
      }
    }
  }
  record$filled <- filled
  record$seconds <- proc.time()[["elapsed"]] - started
  result$certificate <- certificate_with_completeness(result$certificate, status, record)
  msg <- svd_completeness_warning(status, record)
  if (!is.null(msg)) {
    result$warnings <- c(result$warnings, msg)
  }
  result
}
