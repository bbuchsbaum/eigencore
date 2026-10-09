# interval(a, b) eigenvalue targets and spectrum slicing (tranche 5, phase 2
# of capability gap 4).
#
# A Hermitian problem (or a pencil with SPD B) with an explicit matrix is
# solved for every eigenvalue in the closed interval [a, b]:
#
# 1. Count. m = #{lambda in [a, b]} by Sylvester inertia at both end points
#    (R/inertia.R). An unreliable end point is moved outward, so the counted
#    interval [a', b'] contains [a, b] and eigenvalues inside the numerical
#    zero band of an end point count as inside. A supplied k is a cap.
# 2. Solve.
#    * Dense: LAPACK dsyevr with RANGE = "V" on (a' - eta, b'] (generalized:
#      after the Cholesky reduction). The tridiagonal Sturm bisection selects
#      the eigenvalues, so the set is complete by construction ("exact").
#    * Sparse: spectrum slicing. [a', b'] is split by inertia counts at slice
#      mid points until every slice holds at most `slice_size` eigenvalues
#      (adjacent small slices are merged again). Every slice is solved by
#      shift-invert Lanczos at its centre: the c_j eigenvalues of a slice are
#      exactly the c_j nearest its centre. The centre's LDL' factor comes
#      from the shared inertia context (one symbolic analysis, Matrix::update
#      for every shift) and its inertia is a free count at the centre that is
#      checked against the returned values. A few buffer pairs per slice
#      cover eigenvalues that sit numerically on a slice boundary: pairs near
#      a boundary from both neighbours are merged by a Rayleigh-Ritz step on
#      the union of their vectors, which removes duplicates and fills gaps.
# 3. Certify. Residual certificate on the original problem, then the count
#    certificate: by Kahan's theorem each returned value has its own
#    eigenvalue within rho (R/completeness_inertia.R). If the m returned
#    values lie inside [a', b'] by more than rho (plus the factorisation's
#    backward error) the returned set is exactly the eigenvalues in the
#    interval: "inertia_verified". A value within rho of an end point is
#    boundary-ambiguous; one recount on the interval widened by rho decides:
#    the same count proves the set, a larger count leaves it
#    "inertia_inconclusive" (an eigenvalue sits within the residual bound
#    just outside). A returned count different from m is "inertia_failed"
#    (passed = FALSE).

#' @keywords internal
interval_dense_label <- function() {
  "interval: native dense LAPACK dsyevr (RANGE = V)"
}

#' @keywords internal
interval_dense_generalized_label <- function() {
  "interval: native dense generalized SPD (Cholesky + dsyevr RANGE = V)"
}

#' @keywords internal
interval_slicing_label <- function() {
  "interval: inertia spectrum slicing (sparse LDL' shift-invert Lanczos)"
}

#' @keywords internal
interval_route_labels <- function() {
  c(interval_dense_label(), interval_dense_generalized_label(),
    interval_slicing_label())
}

#' @keywords internal
interval_dense_limit <- function() {
  as.integer(getOption("eigencore.interval_dense_limit", 800L))
}

# ---------------------------------------------------------------------------
# Planning
# ---------------------------------------------------------------------------

#' @keywords internal
interval_operator_matrix <- function(op) {
  if (is.null(op)) {
    return(NULL)
  }
  source_or_null(op) %||% op$metadata$matrix %||% NULL
}

#' @keywords internal
interval_route <- function(problem, allow_dense_fallback = "auto") {
  A <- interval_operator_matrix(problem$A)
  B <- interval_operator_matrix(problem$metric)
  if (is.null(A) || (!is.null(problem$metric) && is.null(B))) {
    stop("interval() targets need an explicit dense or sparse matrix (A and B): ",
         "the eigenvalues are counted by LDL' inertia.", call. = FALSE)
  }
  n <- problem$A$dim[[1L]]
  dense_A <- is.matrix(A)
  dense_B <- is.null(B) || is.matrix(B) || inherits(B, "diagonalMatrix")
  sparse_ok <- (inherits(A, "sparseMatrix") || inherits(A, "diagonalMatrix")) &&
    (is.null(B) || inherits(B, "sparseMatrix") || inherits(B, "diagonalMatrix"))
  if (is.complex(A) || is.complex(B)) {
    if (!dense_A || !dense_B || !is.null(B)) {
      stop("interval() targets support complex Hermitian matrices only for dense ",
           "standard problems.", call. = FALSE)
    }
    return(interval_dense_label())
  }
  small <- n <= interval_dense_limit() && !identical(allow_dense_fallback, "never")
  if ((dense_A || small) && (dense_B || small)) {
    return(if (is.null(B)) interval_dense_label() else interval_dense_generalized_label())
  }
  if (!sparse_ok) {
    stop("interval() targets need dense matrices or sparse A with sparse or ",
         "diagonal B.", call. = FALSE)
  }
  interval_slicing_label()
}

#' @keywords internal
plan_interval_eigen <- function(problem, k, method, tol, maxit, vectors,
                                certify, allow_dense_fallback,
                                initial_subspace, left_vectors) {
  n <- as.integer(problem$A$dim[[1L]])
  if (!identical(problem$structure$kind, "hermitian")) {
    stop("interval() targets require a Hermitian (real symmetric or complex ",
         "Hermitian) problem: the eigenvalues must be real to lie in an interval.",
         call. = FALSE)
  }
  if (!is_auto_method(method)) {
    stop("interval() targets choose their own route (dense LAPACK or LDL' ",
         "shift-invert spectrum slicing); use method = auto().", call. = FALSE)
  }
  if (!is.null(initial_subspace)) {
    stop("initial_subspace is not supported for interval() targets.", call. = FALSE)
  }
  if (!is.null(k)) {
    k <- validate_solution_count(k, n, "k")
  }
  if (!is.null(problem$metric) && !isTRUE(tryCatch(
    generalized_spd_metric_known(problem$metric), error = function(e) FALSE))) {
    stop("interval() targets with B require a symmetric positive definite B.",
         call. = FALSE)
  }
  chosen <- interval_route(problem, allow_dense_fallback)
  planner_policy <- planner_policy_snapshot()
  execution <- new_plan_execution(
    "eigen", tol = tol, maxit = maxit, vectors = vectors, certify = certify,
    allow_dense_fallback = allow_dense_fallback, initial_subspace = NULL
  )
  execution$left_vectors <- left_vectors
  value <- problem$target$value
  controls <- list(
    interval = c(value$lower, value$upper),
    k_cap = if (is.null(k)) NA_integer_ else k,
    dense_limit = interval_dense_limit(),
    slice_size = getOption("eigencore.interval_slice_size", NULL),
    max_subspace = method$max_subspace %||% NULL,
    iteration_limit = maxit %||% 100L,
    iteration_limit_kind = "restarts_per_slice"
  )
  reasons <- c(
    paste0("structure: ", problem$structure$kind),
    paste0("target: ", target_label(problem$target)),
    if (is.null(problem$metric)) "standard eigenproblem" else "metric/operator B supplied",
    if (is.null(k)) "k inferred from the inertia count of the interval" else
      paste0("k = ", k, " is an upper bound on the interval count"),
    if (identical(chosen, interval_slicing_label())) {
      "sparse matrix above eigencore.interval_dense_limit: LDL' inertia counts and shift-invert spectrum slicing"
    } else {
      "dense (or small) matrix: LAPACK dsyevr value range selects the interval exactly"
    }
  )
  new_plan(
    problem,
    k = if (is.null(k)) n else k,
    method = chosen,
    method_descriptor = method,
    reasons = reasons,
    fallback = "none; interval routes count and certify their own result",
    controls = controls,
    execution = execution,
    planner_policy = planner_policy
  )
}

# ---------------------------------------------------------------------------
# Execution
# ---------------------------------------------------------------------------

#' @keywords internal
solve_interval_eigen <- function(plan) {
  problem <- plan$problem
  if (plan$method %in% c(interval_dense_label(), interval_dense_generalized_label())) {
    return(solve_interval_dense(problem, plan))
  }
  solve_interval_sparse(problem, plan)
}

#' @keywords internal
interval_k_cap_check <- function(plan, m) {
  cap <- plan$controls$k_cap %||% NA_integer_
  if (!is.na(cap) && m > cap) {
    stop("interval(", format(plan$controls$interval[[1L]]), ", ",
         format(plan$controls$interval[[2L]]), ") holds ", m,
         " eigenvalues, more than k = ", cap,
         "; omit k (it is inferred from the count) or increase it.",
         call. = FALSE)
  }
  invisible(TRUE)
}

#' @keywords internal
interval_result <- function(plan, values, vectors, cert, record, status,
                            iter, extras = list(), warnings = character()) {
  problem <- plan$problem
  ord <- order(values)
  values <- values[ord]
  if (!is.null(vectors)) {
    vectors <- vectors[, ord, drop = FALSE]
  }
  if (!is.null(cert) && length(cert$residuals) == length(ord)) {
    for (field in c("residuals", "backward_error", "converged")) {
      if (length(cert[[field]]) == length(ord)) {
        cert[[field]] <- cert[[field]][ord]
      }
    }
    cert$failed_indices <- which(!cert$converged)
  }
  if (!is.null(record$boundary_ambiguous) && length(record$boundary_ambiguous)) {
    record$boundary_ambiguous <- sort(match(record$boundary_ambiguous, ord))
  }
  cert <- certificate_with_completeness(cert, status, record)
  if (identical(status, "inertia_failed")) {
    warnings <- c(warnings, paste0(
      "interval count certificate failed (", record$reason %||% "count mismatch",
      "); certificate withheld"))
  }
  if (length(record$boundary_ambiguous %||% integer())) {
    warnings <- c(warnings, paste0(
      length(record$boundary_ambiguous), " returned eigenvalue(s) lie within the ",
      "residual bound of an interval end point (boundary-ambiguous)"))
  }
  make_eigen_result(
    values = values,
    vectors = if (isTRUE(plan$execution$vectors)) vectors else NULL,
    certificate = cert,
    iter = iter,
    requested = length(values),
    method_label = plan$method,
    target_label_value = target_label(problem$target),
    plan = plan,
    warnings = warnings,
    extras = extras
  )
}

#' @keywords internal
interval_empty_certificate <- function(tol, note) {
  new_certificate(
    tol = tol, residuals = numeric(), backward_error = numeric(),
    orthogonality = numeric(), converged = logical(), scale = NA_real_,
    notes = note, certificate_type = "interval_count", norm_bound_type = "none"
  )
}

# Dense route: dsyevr RANGE = "V" (Cholesky reduction for a pencil).
#' @keywords internal
solve_interval_dense <- function(problem, plan) {
  started <- proc.time()[["elapsed"]]
  tol <- plan$execution$tol
  certify <- isTRUE(plan$execution$certify)
  lower <- plan$controls$interval[[1L]]
  upper <- plan$controls$interval[[2L]]
  # Dense sources as they are; sparse sources only below
  # eigencore.interval_dense_limit (the planner's choice).
  densify <- function(op) {
    x <- interval_operator_matrix(op)
    if (is.null(x)) NULL else if (is.matrix(x)) x else as.matrix(x)
  }
  A <- densify(problem$A)
  B <- densify(problem$metric)
  if (!is.complex(A)) storage.mode(A) <- "double"
  n <- nrow(A)
  R <- NULL
  C <- A
  if (!is.null(B)) {
    B <- (B + t(B)) / 2
    R <- tryCatch(chol(B), error = function(e) NULL)
    if (is.null(R)) {
      stop("interval() targets with B require a symmetric positive definite B.",
           call. = FALSE)
    }
    C <- backsolve(R, t(backsolve(R, A, transpose = TRUE)), transpose = TRUE)
    C <- (C + t(C)) / 2
  }
  norm1 <- if (is.complex(C)) max(colSums(Mod(C)), 0) else max(colSums(abs(C)), 0)
  eps <- .Machine$double.eps
  # Numerical zero band of the reduction: eigenvalues within eta of an end
  # point are counted as inside (closed interval).
  eta <- 4 * max(n, 1) * eps * max(norm1, .Machine$double.xmin)
  lo <- if (is.finite(lower)) lower - eta else -(norm1 + 1)
  hi <- if (is.finite(upper)) upper + eta else norm1 + 1
  if (is.complex(C)) {
    eig <- native_dense_complex_hermitian_eigen(C, vectors = TRUE)
    keep <- which(eig$values > lo & eig$values <= hi)
    values <- eig$values[keep]
    Y <- eig$vectors[, keep, drop = FALSE]
    driver <- "zheev_filter"
  } else {
    eig <- .Call("eigencore_dense_symmetric_eigen_value_range", C, as.double(lo),
                 as.double(hi), TRUE, PACKAGE = "eigencore")
    values <- eig$values
    Y <- eig$vectors
    driver <- "dsyevr_range_v"
  }
  m <- length(values)
  interval_k_cap_check(plan, m)
  X <- if (is.null(R)) Y else backsolve(R, Y)
  cert <- if (!m) {
    interval_empty_certificate(tol, "no eigenvalue in the interval")
  } else if (certify) {
    certify_eigen_operator(problem$A, values, X, Bop = problem$metric, tol = tol)
  } else {
    empty_certificate(tol, note = "certification disabled by caller")
  }
  ambiguous <- integer()
  if (m) {
    # Residuals of the (reduced) standard problem C y = lambda y, whose
    # vectors are orthonormal: Kahan's bound needs no lambda_min(B).
    CY <- C %*% Y
    res <- sqrt(colSums(Mod(CY - sweep(Y, 2L, values, `*`))^2))
    rho <- sqrt(sum(res^2)) + eta
    ambiguous <- which((is.finite(lower) & values - rho <= lower) |
                         (is.finite(upper) & values + rho >= upper))
  }
  record <- list(
    method = "lapack_value_range",
    driver = driver,
    lower = lower, upper = upper,
    counted_interval = c(lo, hi),
    zero_band = eta,
    count = m,
    returned = m,
    boundary_ambiguous = ambiguous,
    reason = "LAPACK tridiagonal bisection selects every eigenvalue in the range",
    seconds = proc.time()[["elapsed"]] - started
  )
  status <- if (!certify && m) "not_checked" else "exact"
  interval_result(
    plan, values, X, cert, record, status,
    iter = list(iterations = 1L, matvecs = 0L),
    extras = list(
      interval = record,
      restart = list(kind = "dense_interval_lapack", implemented = TRUE,
                     native = TRUE, eigensolver = driver,
                     materialized_dense_operator = TRUE)
    )
  )
}

# ---------------------------------------------------------------------------
# Sparse: counts, slicing, per-slice shift-invert, merge
# ---------------------------------------------------------------------------

# A "sparse" inertia context over the problem (counts on a tridiagonal or
# diagonal input stay exact through the general context; solves always need
# the sparse factor).
#' @keywords internal
interval_contexts <- function(problem) {
  ctx <- inertia_context(problem$A, problem$metric)
  if (identical(ctx$kind, "sparse")) {
    return(list(count = ctx, solve = ctx))
  }
  sctx <- new.env(parent = emptyenv())
  A <- interval_operator_matrix(problem$A)
  B <- interval_operator_matrix(problem$metric)
  sym <- function(x) {
    x <- methods::as(x, "CsparseMatrix")
    if (!inherits(x, "symmetricMatrix")) {
      x <- Matrix::forceSymmetric(x, uplo = "U")
    }
    methods::as(x, "CsparseMatrix")
  }
  sctx$n <- ctx$n
  sctx$kind <- "sparse"
  sctx$generalized <- ctx$generalized
  sctx$complex <- FALSE
  sctx$normA <- ctx$normA
  sctx$normB <- ctx$normB
  sctx$factor <- NULL
  sctx$factorizations <- 0L
  sctx$A <- sym(A)
  sctx$B <- if (is.null(B)) NULL else sym(B)
  sctx$method <- "sparse_cholmod_ldl"
  list(count = ctx, solve = sctx)
}

# Finite bounds of the spectrum: |lambda| <= ||A||_1 / lambda_min(B).
#' @keywords internal
interval_spectrum_bound <- function(ctx, problem) {
  beta <- 1
  if (!is.null(problem$metric)) {
    beta <- inertia_metric_lower_bound(problem$metric)
    if (!is.finite(beta) || beta <= 0) {
      stop("could not bound lambda_min(B) for an infinite interval end point.",
           call. = FALSE)
    }
  }
  1.01 * ctx$normA / beta + 1e-300
}

# Count below t (nudged in `direction` until reliable).
#' @keywords internal
interval_count_point <- function(ctx, t, direction) {
  r <- inertia_count_nudged(ctx, t, direction = direction)
  list(s = r$s, below = r$below, zero = r$zero, reliable = isTRUE(r$reliable),
       factorizations = NROW(r$attempts),
       backward_bound = r$tally$backward_bound %||% NA_real_)
}

#' @keywords internal
interval_slice_size <- function(ctx, plan, m) {
  forced <- plan$controls$slice_size %||% NULL
  if (!is.null(forced)) {
    forced <- suppressWarnings(as.integer(forced))
    if (length(forced) == 1L && !is.na(forced) && forced >= 1L) {
      return(forced)
    }
  }
  # Cost per slice ~ F (numeric factorisation) + per-slice overhead, plus
  # ~40 n s^2 for the Lanczos basis work of s pairs; minimising
  # (m / s) * (F + overhead) + 40 n m s gives s = sqrt((F + overhead) / 40 n).
  cost <- tryCatch(inertia_factor_cost(ctx), error = function(e) NULL)
  flops <- cost$flops %||% NA_real_
  if (!is.finite(flops)) {
    flops <- 50 * ctx$n
  }
  overhead <- 5e7
  s <- sqrt((flops + overhead) / (40 * ctx$n))
  as.integer(min(200L, max(40L, ceiling(s))))
}

# Split [lo, hi] (counts below_lo, below_hi) by mid-point counts until every
# slice holds at most `size` eigenvalues. Returns a data frame of slices.
#' @keywords internal
interval_partition <- function(ctx, lo, hi, below_lo, below_hi, size, scale,
                               max_depth = 60L, max_stall = 6L) {
  slices <- list()
  factorizations <- 0
  # `stall` counts consecutive halvings that left every eigenvalue on one
  # side: a cluster narrower than the slice. After max_stall of them (the
  # slice shrank 64-fold around the cluster) the cluster becomes one slice;
  # block Lanczos resolves its multiplicity.
  stack <- list(list(lo = lo, hi = hi, blo = below_lo, bhi = below_hi, depth = 0L,
                     stall = 0L))
  while (length(stack)) {
    cur <- stack[[length(stack)]]
    stack[[length(stack)]] <- NULL
    count <- cur$bhi - cur$blo
    width <- cur$hi - cur$lo
    if (count <= size || cur$depth >= max_depth || cur$stall >= max_stall ||
        width <= 1e-8 * scale) {
      slices[[length(slices) + 1L]] <- data.frame(
        lower = cur$lo, upper = cur$hi, below_lower = cur$blo,
        below_upper = cur$bhi, count = count, reliable = TRUE
      )
      next
    }
    mid <- cur$lo + width / 2
    p <- interval_count_point(ctx, mid, direction = 1)
    factorizations <- factorizations + p$factorizations
    if (!isTRUE(p$reliable) || !is.finite(p$below) ||
        p$s <= cur$lo || p$s >= cur$hi || p$below < cur$blo || p$below > cur$bhi) {
      slices[[length(slices) + 1L]] <- data.frame(
        lower = cur$lo, upper = cur$hi, below_lower = cur$blo,
        below_upper = cur$bhi, count = count, reliable = isTRUE(p$reliable)
      )
      next
    }
    stall <- if (p$below == cur$blo || p$below == cur$bhi) cur$stall + 1L else 0L
    stack[[length(stack) + 1L]] <- list(lo = p$s, hi = cur$hi, blo = p$below,
                                        bhi = cur$bhi, depth = cur$depth + 1L,
                                        stall = stall)
    stack[[length(stack) + 1L]] <- list(lo = cur$lo, hi = p$s, blo = cur$blo,
                                        bhi = p$below, depth = cur$depth + 1L,
                                        stall = stall)
  }
  out <- do.call(rbind, slices)
  out <- out[order(out$lower), , drop = FALSE]
  # Merge adjacent slices while the merged count stays within `size` (fewer
  # factorisations), and drop empty slices.
  merged <- list()
  for (i in seq_len(nrow(out))) {
    row <- out[i, ]
    if (length(merged)) {
      last <- merged[[length(merged)]]
      if (last$count + row$count <= size) {
        last$upper <- row$upper
        last$below_upper <- row$below_upper
        last$count <- last$count + row$count
        last$reliable <- last$reliable && row$reliable
        merged[[length(merged)]] <- last
        next
      }
    }
    merged[[length(merged) + 1L]] <- row
  }
  out <- do.call(rbind, merged)
  rownames(out) <- NULL
  attr(out, "factorizations") <- factorizations
  out
}

# Shift-invert Lanczos for the `want` eigenvalues nearest `center`, of
# which `count` lie in the slice [lower, upper]. A single-vector Krylov
# method sees one copy of an exactly repeated eigenvalue per start vector, so
# when fewer than `count` converged values land in the slice the solve is
# repeated with block Lanczos, the block sized by the deficit.
#' @keywords internal
interval_slice_solve <- function(problem, sctx, center, want, count, lower,
                                 upper, tol, plan, metric_factor,
                                 restart_limit) {
  started <- proc.time()[["elapsed"]]
  n <- sctx$n
  # The shift must not sit on an eigenvalue: move it inside the slice until
  # the LDL' factor is reliable (its inertia is then a valid count too).
  half <- (upper - lower) / 2
  tally <- NULL
  for (offset in c(0, 0.25, -0.25, 0.4, -0.4, 0.1, -0.1)) {
    shift <- center + offset * half
    tally <- tryCatch(suppressWarnings(inertia_at(sctx, shift)), error = function(e) NULL)
    if (!is.null(tally) && isTRUE(inertia_tally_reliable(tally))) {
      center <- shift
      break
    }
  }
  factor <- if (!is.null(tally) && isTRUE(tally$ok) &&
                isTRUE(all.equal(tally$sigma, center))) tally$factor else NULL
  prep <- prepare_shift_invert_operator(problem, center, tol = tol,
                                        ldl_factor = factor,
                                        metric_factor = metric_factor,
                                        try_spd = FALSE)
  M <- prep$operator
  center_count <- if (!is.null(tally) && isTRUE(inertia_tally_reliable(tally))) {
    tally$neg
  } else {
    NA_real_
  }
  k <- min(as.integer(want), n - 1L)
  subspace <- as.integer(plan$controls$max_subspace %||%
                           default_shift_invert_max_subspace(n, k))
  subspace <- min(n, max(subspace, k + 2L))
  slack <- 1e-8 * (sctx$normA + abs(center) * sctx$normB)
  run <- function(block, subspace, inner_tol) {
    iter <- if (block <= 1L) {
      shift_invert_transformed_lanczos(M, k = k, tol = inner_tol, maxit = subspace,
                                       max_restarts = restart_limit)
    } else {
      native_block_lanczos_hermitian(
        M, k = k, target = largest_magnitude(), tol = inner_tol,
        maxit = min(n, max(subspace, k + 3L * block)), block = block,
        max_restarts = restart_limit, vectors = TRUE, full_subspace = FALSE,
        certificate_fallback = FALSE
      )
    }
    mu <- as.numeric(iter$values)
    ok <- is.finite(mu) & abs(mu) > .Machine$double.eps
    lambda <- center + 1 / mu[ok]
    X <- prep$recover_vectors(iter$vectors[, ok, drop = FALSE])
    AX <- as.matrix(apply_operator(problem$A, X))
    BX <- if (is.null(problem$metric)) X else as.matrix(apply_operator(problem$metric, X))
    res <- sqrt(colSums((AX - sweep(BX, 2L, lambda, `*`))^2))
    scale <- sctx$normA + abs(lambda) * sctx$normB
    xn <- sqrt(colSums(X * BX))
    converged <- res <= tol * scale * pmax(xn, .Machine$double.xmin)
    inside <- converged & lambda >= lower - slack & lambda <= upper + slack
    # Radius around the shift that the converged values cover: the slice is
    # covered once it reaches both slice ends.
    radius <- if (any(converged)) max(abs(lambda[converged] - center)) else 0
    list(values = lambda, vectors = X, residuals = res, converged = converged,
         found = sum(inside), block = block, k = k,
         covered = radius >= max(center - lower, upper - center) - slack,
         matvecs = as.integer(iter$matvecs %||% iter$operator_columns %||% 0L),
         iterations = as.integer(iter$iterations %||% 0L))
  }
  best <- NULL
  matvecs <- 0L
  iterations <- 0L
  block <- 1L
  inner_tol <- tol
  for (attempt in 1:4) {
    cand <- run(block, subspace, inner_tol)
    matvecs <- matvecs + cand$matvecs
    iterations <- iterations + cand$iterations
    if (is.null(best) || cand$found > best$found ||
        (cand$found == best$found && sum(cand$converged) > sum(best$converged))) {
      best <- cand
    }
    if (best$found >= count) {
      break
    }
    deficit <- count - best$found
    if (all(cand$converged) && !isTRUE(cand$covered) && k < n - 1L) {
      # The shift is off the slice centre (it had to avoid an eigenvalue), so
      # the `want` nearest values do not reach both slice ends: ask for more.
      k <- as.integer(min(n - 1L, k + deficit + max(4L, deficit)))
      subspace <- min(n, max(subspace, default_shift_invert_max_subspace(n, k)))
    } else if (all(cand$converged)) {
      # Everything converged but copies are missing: repeated eigenvalues.
      block <- as.integer(min(count + 1L, max(2L * block, deficit + 1L), n %/% 4L + 1L))
    } else {
      inner_tol <- max(inner_tol * 1e-3, 10 * .Machine$double.eps)
      subspace <- min(n, 2L * subspace)
      if (attempt >= 2L) {
        block <- as.integer(min(count + 1L, max(2L * block, deficit + 1L), n %/% 4L + 1L))
      }
    }
  }
  best$center <- center
  best$center_count <- center_count
  best$matvecs <- matvecs
  best$iterations <- iterations
  best$label_kind <- prep$label_kind
  best$native_solve <- isTRUE(M$metadata$native_solve)
  best$ldl_fallback <- isTRUE(prep$ldl_fallback)
  best$seconds <- proc.time()[["elapsed"]] - started
  best
}

# Rayleigh-Ritz of (A, B) on span(X): B-orthonormal basis with near-dependent
# directions (duplicate copies of one eigenvector) dropped.
#' @keywords internal
interval_rayleigh_ritz <- function(problem, X) {
  X <- as.matrix(X)
  BX <- if (is.null(problem$metric)) X else as.matrix(apply_operator(problem$metric, X))
  G <- crossprod(X, BX)
  G <- (G + t(G)) / 2
  eg <- eigen(G, symmetric = TRUE)
  keep <- eg$values > 1e-8 * max(eg$values[[1L]], .Machine$double.xmin)
  Q <- X %*% sweep(eg$vectors[, keep, drop = FALSE], 2L, sqrt(eg$values[keep]), `/`)
  AQ <- as.matrix(apply_operator(problem$A, Q))
  H <- crossprod(Q, AQ)
  H <- (H + t(H)) / 2
  eh <- eigen(H, symmetric = TRUE)
  list(values = eh$values, vectors = Q %*% eh$vectors)
}

# Merge per-slice solutions: per slice its `count` picks nearest the centre,
# plus its buffer pairs near any slice boundary; groups of pairs from several
# slices near a common boundary are replaced by their Rayleigh-Ritz pairs.
#' @keywords internal
interval_merge_slices <- function(problem, slices, solved, width) {
  vals <- numeric()
  vecs <- list()
  owner <- integer()
  pick <- logical()
  for (j in seq_along(solved)) {
    s <- solved[[j]]
    if (is.null(s) || !length(s$values)) next
    # The slice's eigenvalues are the `count` nearest its mid point (the
    # shift itself may sit off centre).
    mid <- (slices$lower[[j]] + slices$upper[[j]]) / 2
    ord <- order(abs(s$values - mid))
    take <- utils::head(ord, slices$count[[j]])
    rest <- setdiff(ord, take)
    edges <- c(slices$lower[[j]], slices$upper[[j]])
    near <- rest[vapply(rest, function(i) any(abs(s$values[[i]] - edges) <= width),
                        logical(1L))]
    idx <- c(take, near)
    vals <- c(vals, s$values[idx])
    vecs[[length(vecs) + 1L]] <- s$vectors[, idx, drop = FALSE]
    owner <- c(owner, rep(j, length(idx)))
    pick <- c(pick, rep(c(TRUE, FALSE), c(length(take), length(near))))
  }
  X <- if (length(vecs)) do.call(cbind, vecs) else matrix(0, problem$A$dim[[1L]], 0L)
  merged_groups <- 0L
  internal <- if (nrow(slices) > 1L) slices$upper[-nrow(slices)] else numeric()
  drop <- rep(FALSE, length(vals))
  add_vals <- numeric()
  add_vecs <- list()
  for (s in internal) {
    g <- which(abs(vals - s) <= width & !drop)
    if (length(unique(owner[g])) < 2L) next
    rr <- interval_rayleigh_ritz(problem, X[, g, drop = FALSE])
    drop[g] <- TRUE
    add_vals <- c(add_vals, rr$values)
    add_vecs[[length(add_vecs) + 1L]] <- rr$vectors
    merged_groups <- merged_groups + 1L
  }
  keep <- which(!drop & pick)
  out_vals <- c(vals[keep], add_vals)
  out_vecs <- do.call(cbind, c(list(X[, keep, drop = FALSE]), add_vecs))
  # Buffer pairs that were not merged are kept as candidates too: when the
  # count assigns a boundary eigenvalue to the other slice than the one that
  # computed it, the buffer copy is the only one.
  extra <- which(!drop & !pick)
  if (length(extra)) {
    out_vals <- c(out_vals, vals[extra])
    out_vecs <- cbind(out_vecs, X[, extra, drop = FALSE])
  }
  list(values = out_vals, vectors = out_vecs, merged_groups = merged_groups,
       n_buffer = length(extra))
}

# Select the m candidates nearest the interval centre, dropping duplicates
# (near-parallel vectors with equal values) first.
#' @keywords internal
interval_select <- function(problem, values, vectors, m, center, width) {
  ord <- order(abs(values - center))
  chosen <- integer()
  for (i in ord) {
    if (length(chosen) >= m) break
    dup <- FALSE
    close <- chosen[abs(values[chosen] - values[[i]]) <= width]
    if (length(close)) {
      v <- vectors[, i]
      Bv <- if (is.null(problem$metric)) v else as.numeric(apply_operator(problem$metric, matrix(v)))
      nv <- sqrt(sum(v * Bv))
      for (j in close) {
        w <- vectors[, j]
        Bw <- if (is.null(problem$metric)) w else as.numeric(apply_operator(problem$metric, matrix(w)))
        if (abs(sum(v * Bw)) > 0.5 * nv * sqrt(sum(w * Bw))) {
          dup <- TRUE
          break
        }
      }
    }
    if (!dup) chosen <- c(chosen, i)
  }
  chosen
}

# Count certificate for an interval result (see the file header).
#' @keywords internal
interval_completeness_check <- function(ctx, problem, values, residuals,
                                        orthogonality, lower_eff, upper_eff,
                                        m, count_bound) {
  started <- proc.time()[["elapsed"]]
  k <- length(values)
  eps <- .Machine$double.eps
  record <- list(
    method = "inertia_count_interval",
    count = m, returned = k,
    counted_interval = c(lower_eff, upper_eff),
    rho = NA_real_, margin = NA_real_,
    boundary_ambiguous = integer(),
    widened_interval = c(NA_real_, NA_real_), widened_count = NA_real_,
    factorizations = 0, metric_lower_bound = NA_real_,
    reason = NA_character_, seconds = NA_real_
  )
  finish <- function(status, reason = NA_character_) {
    record$reason <- reason
    record$seconds <- proc.time()[["elapsed"]] - started
    list(status = status, record = record)
  }
  if (k != m) {
    return(finish("inertia_failed",
                  paste0(k, " values returned but the interval holds ", m)))
  }
  if (!k) {
    return(finish("inertia_verified", "the interval holds no eigenvalue"))
  }
  omega <- suppressWarnings(max(abs(as.numeric(orthogonality)), 0, na.rm = TRUE))
  gram_floor <- 1 - k * omega
  if (!is.finite(gram_floor) || gram_floor < 0.5) {
    return(finish("inertia_inconclusive",
                  "returned vectors are too far from orthonormal for the residual bound"))
  }
  beta <- 1
  if (!is.null(problem$metric)) {
    beta <- inertia_metric_lower_bound(problem$metric)
    record$metric_lower_bound <- beta
    if (!is.finite(beta) || beta <= 0) {
      return(finish("inertia_inconclusive", "no positive lower bound of lambda_min(B)"))
    }
  }
  residuals <- abs(as.numeric(residuals))
  if (length(residuals) != k || any(!is.finite(residuals))) {
    return(finish("inertia_inconclusive", "certificate residuals unavailable"))
  }
  vmax <- max(abs(values))
  xnorm <- sqrt((1 + omega) / beta)
  rounding <- sqrt(k) * 32 * eps * (ctx$normA + vmax * ctx$normB) * xnorm
  rho <- (sqrt(sum(residuals^2)) + rounding) / sqrt(beta * gram_floor)
  scale <- ctx$normA + vmax * ctx$normB
  margin <- max(if (is.finite(count_bound)) count_bound else 0, 64 * eps * scale)
  record$rho <- rho
  record$margin <- margin
  outside <- values + rho + margin < lower_eff | values - rho - margin > upper_eff
  if (any(outside)) {
    return(finish("inertia_failed",
                  "a returned value lies outside the interval by more than its residual bound"))
  }
  ambiguous <- which(values - rho - margin <= lower_eff | values + rho + margin >= upper_eff)
  record$boundary_ambiguous <- ambiguous
  if (!length(ambiguous)) {
    # m distinct eigenvalues matched strictly inside an interval holding m.
    return(finish("inertia_verified"))
  }
  # Recount on the interval widened by rho + margin: if it still holds m,
  # no eigenvalue lies in the ambiguous bands just outside, so every matched
  # eigenvalue is one of the m inside.
  wl <- if (is.finite(lower_eff)) lower_eff - rho - 2 * margin else -Inf
  wu <- if (is.finite(upper_eff)) upper_eff + rho + 2 * margin else Inf
  lo <- if (is.finite(wl)) interval_count_point(ctx, wl, direction = -1) else
    list(below = 0, zero = 0, reliable = TRUE, factorizations = 0, s = -Inf)
  hi <- if (is.finite(wu)) interval_count_point(ctx, wu, direction = 1) else
    list(below = ctx$n, zero = 0, reliable = TRUE, factorizations = 0, s = Inf)
  record$factorizations <- lo$factorizations + hi$factorizations
  record$widened_interval <- c(lo$s, hi$s)
  if (!isTRUE(lo$reliable) || !isTRUE(hi$reliable)) {
    return(finish("inertia_inconclusive",
                  "the widened recount at the end points was not reliable"))
  }
  wide <- hi$below + hi$zero - lo$below
  record$widened_count <- wide
  if (wide == m) {
    return(finish("inertia_verified",
                  "boundary-ambiguous values resolved by the widened recount"))
  }
  finish("inertia_inconclusive",
         "an eigenvalue lies within the residual bound just outside an end point")
}

#' @keywords internal
solve_interval_sparse <- function(problem, plan) {
  started <- proc.time()[["elapsed"]]
  tol <- plan$execution$tol
  certify <- isTRUE(plan$execution$certify)
  lower <- plan$controls$interval[[1L]]
  upper <- plan$controls$interval[[2L]]
  contexts <- interval_contexts(problem)
  ctx <- contexts$count
  sctx <- contexts$solve
  n <- ctx$n
  bound <- if (is.finite(lower) && is.finite(upper)) NA_real_ else
    interval_spectrum_bound(ctx, problem)
  # 1. Count at the end points (moved outward while unreliable).
  lo <- if (is.finite(lower)) interval_count_point(ctx, lower, direction = -1) else
    list(s = -bound, below = 0, zero = 0, reliable = TRUE, factorizations = 0,
         backward_bound = 0)
  hi <- if (is.finite(upper)) interval_count_point(ctx, upper, direction = 1) else
    list(s = bound, below = n, zero = 0, reliable = TRUE, factorizations = 0,
         backward_bound = 0)
  counts_reliable <- isTRUE(lo$reliable) && isTRUE(hi$reliable)
  below_lo <- lo$below
  below_hi <- hi$below + hi$zero
  if (!is.finite(below_lo) || !is.finite(below_hi)) {
    stop("interval(): the inertia count at an end point failed (no usable ",
         "factorisation of A - t B); the interval cannot be counted.", call. = FALSE)
  }
  m <- as.integer(below_hi - below_lo)
  interval_k_cap_check(plan, m)
  count_bound <- max(c(lo$backward_bound, hi$backward_bound, 0), na.rm = TRUE)
  record_base <- list(
    lower = lower, upper = upper,
    counted_interval = c(lo$s, hi$s),
    count = m,
    count_reliable = counts_reliable,
    count_method = ctx$method
  )
  if (m == 0L) {
    record <- c(record_base, list(
      method = "inertia_count_interval", returned = 0L, slices = 0L,
      boundary_ambiguous = integer(),
      reason = "the interval holds no eigenvalue",
      seconds = proc.time()[["elapsed"]] - started
    ))
    status <- if (counts_reliable) "inertia_verified" else "inertia_inconclusive"
    return(interval_result(
      plan, numeric(), matrix(0, n, 0L),
      interval_empty_certificate(tol, "no eigenvalue in the interval"),
      record, status, iter = list(iterations = 0L, matvecs = 0L),
      extras = list(interval = record,
                    restart = list(kind = "interval_spectrum_slicing", native = TRUE))
    ))
  }
  # 2. Slices.
  scale <- ctx$normA + max(abs(c(lo$s, hi$s))) * ctx$normB
  size <- interval_slice_size(sctx, plan, m)
  slices <- interval_partition(ctx, lo$s, hi$s, below_lo, below_hi, size, scale)
  partition_factorizations <- attr(slices, "factorizations")
  metric_factor <- if (is.null(problem$metric)) NULL else
    shift_invert_metric_factor(problem$metric)
  restart_limit <- as.integer(plan$controls$iteration_limit %||% 100L)
  buffer <- function(c) as.integer(min(20L, max(4L, ceiling(0.1 * c))))
  solved <- vector("list", nrow(slices))
  for (j in seq_len(nrow(slices))) {
    if (slices$count[[j]] < 1L) next
    center <- (slices$lower[[j]] + slices$upper[[j]]) / 2
    solved[[j]] <- interval_slice_solve(
      problem, sctx, center, slices$count[[j]] + buffer(slices$count[[j]]),
      slices$count[[j]], slices$lower[[j]], slices$upper[[j]],
      tol, plan, metric_factor, restart_limit
    )
  }
  max_res <- max(c(0, unlist(lapply(solved, function(s) s$residuals))), na.rm = TRUE)
  width <- max(100 * max_res, 1e-10 * scale)
  merged <- interval_merge_slices(problem, slices, solved, width)
  # 3. Select, certify, count.
  center <- (lo$s + hi$s) / 2
  sel <- interval_select(problem, merged$values, merged$vectors, m, center, width)
  values <- merged$values[sel]
  X <- merged$vectors[, sel, drop = FALSE]
  cert <- if (certify) {
    certify_eigen_operator(problem$A, values, X, Bop = problem$metric, tol = tol)
  } else {
    empty_certificate(tol, note = "certification disabled by caller")
  }
  check <- if (certify) {
    interval_completeness_check(ctx, problem, values, cert$residuals,
                                cert$orthogonality, lo$s, hi$s, m, count_bound)
  } else {
    list(status = "not_checked", record = list(reason = "certify = FALSE"))
  }
  status <- check$status
  if (identical(status, "inertia_verified") && !counts_reliable) {
    status <- "inertia_inconclusive"
    check$record$reason <- "an end-point count came from an unreliable factorisation"
  }
  slice_table <- data.frame(
    lower = slices$lower, upper = slices$upper, count = slices$count,
    center = vapply(solved, function(s) s$center %||% NA_real_, numeric(1L)),
    center_count = vapply(solved, function(s) s$center_count %||% NA_real_, numeric(1L)),
    returned_below_center = vapply(seq_along(solved), function(j) {
      s <- solved[[j]]
      if (is.null(s)) return(NA_real_)
      sum(values < s$center & values >= slices$lower[[j]])
    }, numeric(1L)),
    matvecs = vapply(solved, function(s) as.numeric(s$matvecs %||% 0), numeric(1L)),
    native_solve = vapply(solved, function(s) isTRUE(s$native_solve), logical(1L)),
    ldl_fallback = vapply(solved, function(s) isTRUE(s$ldl_fallback), logical(1L)),
    block = vapply(solved, function(s) as.numeric(s$block %||% NA_real_), numeric(1L)),
    seconds = vapply(solved, function(s) as.numeric(s$seconds %||% 0), numeric(1L))
  )
  # Free check from each slice's own factorisation at its centre: the number
  # of returned values in [lower_j, centre_j) must equal the count there.
  centre_mismatch <- which(is.finite(slice_table$center_count) &
    slice_table$center_count - slices$below_lower != slice_table$returned_below_center)
  record <- c(record_base, check$record[setdiff(names(check$record), names(record_base))])
  record$slices <- nrow(slices)
  record$slice_size <- size
  record$slice_table <- slice_table
  record$centre_count_mismatch <- centre_mismatch
  record$merged_boundary_groups <- merged$merged_groups
  record$factorizations <- partition_factorizations +
    sum(c(lo$factorizations, hi$factorizations), na.rm = TRUE) +
    sum(slices$count > 0) + (check$record$factorizations %||% 0)
  record$seconds <- proc.time()[["elapsed"]] - started
  matvecs <- as.integer(sum(slice_table$matvecs))
  warnings <- character()
  if (any(slice_table$ldl_fallback)) {
    warnings <- c(warnings, "a slice's LDL' factor was unreliable; that slice used a sparse LU solve")
  }
  interval_result(
    plan, values, X, cert, record, status,
    iter = list(iterations = as.integer(sum(vapply(solved, function(s)
      as.numeric(s$iterations %||% 0), numeric(1L)))), matvecs = matvecs),
    extras = list(
      interval = record,
      restart = list(kind = "interval_spectrum_slicing", native = TRUE,
                     native_solve = all(slice_table$native_solve[slice_table$count > 0]),
                     slices = nrow(slices))
    ),
    warnings = warnings
  )
}
