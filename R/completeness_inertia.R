# Deterministic target-completeness certificate by inertia counting
# (tranche 5; the exact answer to C50).
#
# A residual certificate proves each returned pair is an eigenpair; it does
# not prove the returned SET is the requested one. For a Hermitian problem
# with an explicit matrix (and SPD B) the set can be proved complete by
# counting eigenvalues with Sylvester's law of inertia:
#
# 1. Kahan's theorem: for Hermitian C, any full-rank X and diagonal Theta,
#    there are k DISTINCT eigenvalues lambda_{j_i} of C with
#    |theta_i - lambda_{j_i}| <= ||C X - X Theta||_2 / sigma_min(X).
#    With the certificate's residual norms r_i, ||R||_2 <= ||R||_F and
#    sigma_min(X)^2 >= 1 - k * omega (omega = max |X'X - I|), so every
#    returned value has its own eigenvalue within
#        rho = ||r||_2 / sqrt(beta (1 - k omega)) + rounding,
#    where beta = 1 for B = I and beta <= lambda_min(B) otherwise (the pencil
#    is C = B^{-1/2} A B^{-1/2} with Y = B^{1/2} X).
# 2. Write the target as "the k eigenvalues with the smallest preference
#    distance f": f(l) = l (smallest), -l (largest), -|l| (largest
#    magnitude), |l - c| (nearest c). f is 1-Lipschitz, so the matched
#    eigenvalues have f <= E + rho, E = max_i f(theta_i).
# 3. Count N(t) = #{l : f(l) < t} by inertia (one or two LDL'
#    factorisations). With a margin m exceeding the factorisation's
#    backward error:
#      N(E + rho + m) == k  =>  the k matched eigenvalues are exactly the k
#                               most preferred ones: "inertia_verified";
#      N(E - rho - m) >= k  =>  k eigenvalues are strictly more preferred
#                               than the edge's matched eigenvalue: a value
#                               is missing ("inertia_failed");
#      also failed when N(E - rho - m) exceeds the number of returned values
#      that could match below that threshold;
#      otherwise the k-th and (k+1)-th eigenvalues are not separated by more
#      than the residual uncertainty (a cluster straddles the target edge):
#      "inertia_inconclusive" (no claim either way).
# Every count must come from a reliable factorisation; an unreliable count
# makes the verdict inconclusive, never verified.

#' @keywords internal
completeness_modes <- function() {
  c("auto", "inertia", "probe", "none")
}

#' @keywords internal
inertia_completeness_controls <- function() {
  num <- function(value, default) {
    value <- suppressWarnings(as.numeric(value %||% default))
    if (length(value) != 1L || !is.finite(value) || value < 0) default else value
  }
  list(
    # Wall-clock budget (seconds) below which auto mode always runs the
    # inertia certificate ...
    seconds = num(getOption("eigencore.completeness_inertia_seconds"), 0.5),
    # ... or, above it, the allowed ratio of predicted factorisation time to
    # the solve's own time.
    ratio = num(getOption("eigencore.completeness_inertia_ratio"), 1),
    # Assumed factorisation rates (flop/s) for the cost prediction.
    dense_rate = num(getOption("eigencore.completeness_dense_flop_rate"), 1e10),
    sparse_rate = num(getOption("eigencore.completeness_sparse_flop_rate"), 1.5e9)
  )
}

# Preference distance used by the counting argument.
#' @keywords internal
inertia_completeness_kind <- function(target) {
  kind <- if (inherits(target, "eigencore_target")) target$kind else NA_character_
  switch(
    kind,
    largest = ,
    largest_real = "largest",
    smallest = ,
    smallest_real = "smallest",
    largest_magnitude = "largest_magnitude",
    nearest = if (is.numeric(target$value) && length(target$value) == 1L &&
                  is.finite(target$value)) "nearest" else NULL,
    NULL
  )
}

#' @keywords internal
inertia_preference <- function(values, kind, center = 0) {
  switch(
    kind,
    smallest = values,
    largest = -values,
    largest_magnitude = -abs(values),
    nearest = abs(values - center)
  )
}

# Number of eigenvalues with preference distance f < t. `outward = TRUE`
# nudges unreliable boundaries so the counted region grows (t' >= t);
# `outward = FALSE` shrinks it (t' <= t). Returns list(count, t, reliable,
# factorizations, attempts).
#' @keywords internal
inertia_region_count <- function(ctx, kind, t, center = 0, outward = TRUE) {
  out <- function(count, reliable, parts, t_eff) {
    list(
      count = if (isTRUE(reliable)) count else NA_real_,
      t = t_eff,
      reliable = isTRUE(reliable),
      dense_fallback = any(vapply(parts, function(p) isTRUE(p$dense_fallback), logical(1L))),
      factorizations = sum(vapply(parts, function(p) NROW(p$attempts), numeric(1L))),
      backward_bound = max(vapply(parts, function(p) {
        b <- p$tally$backward_bound %||% NA_real_
        if (is.finite(b)) b else 0
      }, numeric(1L)))
    )
  }
  dir <- if (outward) 1 else -1
  switch(
    kind,
    smallest = {
      r <- inertia_count_nudged(ctx, t, direction = dir)
      out(r$below, r$reliable, list(r), r$s)
    },
    largest = {
      r <- inertia_count_nudged(ctx, -t, direction = -dir)
      out(r$above, r$reliable, list(r), -r$s)
    },
    largest_magnitude = {
      a <- -t
      if (!(a > 0)) {
        return(list(count = NA_real_, t = t, reliable = FALSE,
                    dense_fallback = FALSE, factorizations = 0, backward_bound = 0))
      }
      up <- inertia_count_nudged(ctx, a, direction = -dir)
      lo <- inertia_count_nudged(ctx, -a, direction = dir)
      a_eff <- min(up$s, -lo$s)
      out(up$above + lo$below, up$reliable && lo$reliable, list(up, lo),
          if (outward) -a_eff else -max(up$s, -lo$s))
    },
    nearest = {
      if (!(t > 0)) {
        return(list(count = 0, t = t, reliable = TRUE,
                    dense_fallback = FALSE, factorizations = 0, backward_bound = 0))
      }
      hi <- inertia_count_nudged(ctx, center + t, direction = dir)
      lo <- inertia_count_nudged(ctx, center - t, direction = -dir)
      t_eff <- if (outward) max(hi$s - center, center - lo$s) else
        min(hi$s - center, center - lo$s)
      out(hi$below - lo$below - lo$zero, hi$reliable && lo$reliable,
          list(hi, lo), t_eff)
    }
  )
}

# Lower bound of lambda_min(B) for SPD B by inertia: the largest
# t = d_min * 2^-j (d_min = min diag(B) >= lambda_min) with no eigenvalue of
# B below t.
#' @keywords internal
inertia_metric_lower_bound <- function(Bop) {
  bctx <- tryCatch(inertia_context(Bop), error = function(e) NULL)
  if (is.null(bctx)) {
    return(NA_real_)
  }
  if (identical(bctx$kind, "diagonal")) {
    return(min(bctx$d))
  }
  dmin <- switch(
    bctx$kind,
    tridiagonal = min(bctx$d),
    dense = min(diag(bctx$A)),
    sparse = min(Matrix::diag(bctx$A))
  )
  if (!is.finite(dmin) || dmin <= 0) {
    return(NA_real_)
  }
  t <- dmin
  for (j in 0:12) {
    r <- inertia_count_nudged(bctx, t, direction = -1)
    if (isTRUE(r$reliable) && r$below == 0) {
      return(r$s)
    }
    t <- t / 4
  }
  NA_real_
}

# Run the counting argument on a certified result. `residuals` and
# `orthogonality` come from the residual certificate (absolute residual
# norms ||A x - theta B x|| and max |X'BX - I|).
#' @keywords internal
inertia_completeness_check <- function(Aop, values, residuals, orthogonality,
                                       target, Bop = NULL, ctx = NULL) {
  started <- proc.time()[["elapsed"]]
  kind <- inertia_completeness_kind(target)
  center <- if (identical(kind, "nearest")) as.numeric(target$value) else 0
  values <- as.numeric(values)
  k <- length(values)
  eps <- .Machine$double.eps
  record <- list(
    method = "inertia_count",
    kind = kind,
    k = k,
    edge = NA_real_,
    rho = NA_real_,
    margin = NA_real_,
    threshold_upper = NA_real_,
    count_upper = NA_real_,
    threshold_lower = NA_real_,
    count_lower = NA_real_,
    factorizations = 0,
    inertia_method = NA_character_,
    metric_lower_bound = NA_real_,
    reason = NA_character_,
    seconds = NA_real_
  )
  finish <- function(status, reason = NA_character_) {
    record$reason <- reason
    record$seconds <- proc.time()[["elapsed"]] - started
    list(status = status, record = record)
  }
  if (is.null(kind) || !k || any(!is.finite(values))) {
    return(finish("inertia_inconclusive", "target or values not supported"))
  }
  if (is.null(ctx)) {
    ctx <- tryCatch(inertia_context(Aop, Bop), error = function(e) e)
    if (inherits(ctx, "error")) {
      return(finish("inertia_inconclusive", conditionMessage(ctx)))
    }
  }
  record$inertia_method <- ctx$method
  n <- ctx$n
  if (k >= n) {
    return(finish("inertia_verified", "full spectrum returned"))
  }
  omega <- suppressWarnings(max(abs(as.numeric(orthogonality)), 0, na.rm = TRUE))
  gram_floor <- 1 - k * omega
  if (!is.finite(gram_floor) || gram_floor < 0.5) {
    return(finish("inertia_inconclusive",
                  "returned vectors are too far from orthonormal for the residual bound"))
  }
  beta <- 1
  if (!is.null(Bop)) {
    beta <- inertia_metric_lower_bound(Bop)
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
  f <- inertia_preference(values, kind, center)
  E <- max(f)
  record$edge <- E
  record$rho <- rho
  scale <- ctx$normA + (abs(center) + vmax) * ctx$normB
  # A relative floor of 1e-10 covers the typical factorisation backward
  # error, so the widening recount below is rarely needed; gaps smaller than
  # that are inside any practical residual bound anyway.
  margin <- max(rho, 1e-10 * scale, 64 * eps * scale)
  # The count is exact for A - t B + E with ||E|| <= backward_bound; the
  # matched eigenvalues sit at least `margin` from the threshold, so the
  # margin must exceed that bound (widen once if it does not).
  upper <- inertia_region_count(ctx, kind, E + rho + margin, center, outward = TRUE)
  record$factorizations <- record$factorizations + upper$factorizations
  if (isTRUE(upper$reliable) && is.finite(upper$backward_bound) &&
      upper$backward_bound >= margin) {
    margin <- 2 * upper$backward_bound
    upper <- inertia_region_count(ctx, kind, E + rho + margin, center, outward = TRUE)
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
    # Contradicts Kahan's bound: the residual certificate or the count is
    # wrong. Make no claim.
    return(finish("inertia_inconclusive",
                  "fewer eigenvalues than returned values inside the residual bound"))
  }
  lower <- inertia_region_count(ctx, kind, E - rho - margin, center, outward = FALSE)
  record$factorizations <- record$factorizations + lower$factorizations
  record$threshold_lower <- lower$t
  record$count_lower <- lower$count
  if (!isTRUE(lower$reliable)) {
    return(finish("inertia_inconclusive", "inertia count below the edge was not reliable"))
  }
  if (lower$count >= k) {
    return(finish("inertia_failed",
                  "at least k eigenvalues are strictly more preferred than the returned edge"))
  }
  could_match <- sum(f - rho < lower$t)
  if (lower$count > could_match) {
    return(finish("inertia_failed",
                  "an eigenvalue among the k most preferred has no returned value within the residual bound"))
  }
  finish("inertia_inconclusive",
         "the k-th and (k+1)-th eigenvalues are not separated by more than the residual bound")
}

# Should this solve use the inertia certificate? Returns list(use, reason,
# predicted_seconds, ctx).
#' @keywords internal
inertia_completeness_gate <- function(problem, mode, k, solve_seconds = NA_real_,
                                      controls = inertia_completeness_controls()) {
  no <- function(reason) list(use = FALSE, reason = reason, ctx = NULL,
                              predicted_seconds = NA_real_)
  if (!mode %in% c("auto", "inertia")) {
    return(no("mode"))
  }
  op <- problem$A
  if (!inherits(op, "eigencore_operator") ||
      !identical(op$structure$kind, "hermitian")) {
    return(no("not Hermitian"))
  }
  if (is.null(inertia_completeness_kind(problem$target))) {
    return(no("target"))
  }
  if (is.null(inertia_matrix_of(op)) ||
      (!is.null(problem$metric) && is.null(inertia_matrix_of(problem$metric)))) {
    return(no("no explicit matrix source"))
  }
  ctx <- tryCatch(inertia_context(op, problem$metric), error = function(e) e)
  if (inherits(ctx, "error")) {
    return(no(paste0("inertia context: ", conditionMessage(ctx))))
  }
  cost <- if (is.null(problem$metric)) {
    operator_memoised_value(op, "inertia_factor_cost", inertia_factor_cost(ctx))
  } else {
    inertia_factor_cost(ctx)
  }
  rate <- if (identical(ctx$kind, "sparse")) controls$sparse_rate else controls$dense_rate
  # Two factorisations in the common verified case (largest_magnitude and
  # nearest need two per threshold).
  per_threshold <- if (inertia_completeness_kind(problem$target) %in%
                       c("largest_magnitude", "nearest")) 2 else 1
  predicted <- per_threshold * cost$flops / rate
  if (identical(mode, "inertia")) {
    return(list(use = TRUE, reason = "requested", ctx = ctx,
                predicted_seconds = predicted))
  }
  if (!is.finite(predicted)) {
    return(list(use = FALSE, reason = "factorisation cost unknown", ctx = NULL,
                predicted_seconds = predicted))
  }
  budget <- max(controls$seconds,
                controls$ratio * (if (is.finite(solve_seconds)) solve_seconds else 0))
  if (predicted <= budget) {
    list(use = TRUE, reason = "cost gate passed", ctx = ctx,
         predicted_seconds = predicted)
  } else {
    list(use = FALSE, reason = "predicted factorisation cost exceeds the gate",
         ctx = NULL, predicted_seconds = predicted)
  }
}

# Inertia certificate with repair for a standard Hermitian result: on
# inertia_failed, the deflated-complement probe/repair runs (and, when the
# probe sees no intruder, a deflated complement solve from a deterministic
# start), then the count is repeated. Returns a check object for
# result_with_completeness().
#' @keywords internal
inertia_completeness_run <- function(problem, values, vectors, cert, tol, ctx,
                                     gate = NULL,
                                     controls = target_completeness_controls()) {
  started <- proc.time()[["elapsed"]]
  target <- problem$target
  op <- problem$A
  check <- inertia_completeness_check(
    op, values, cert$residuals, cert$orthogonality, target,
    Bop = problem$metric, ctx = ctx
  )
  record <- check$record
  record$repaired <- FALSE
  record$rounds <- 0L
  record$operator_columns <- 0L
  record$operator_block_calls <- 0L
  record$gate <- gate$reason %||% NA_character_
  record$predicted_seconds <- gate$predicted_seconds %||% NA_real_
  status <- check$status
  V <- vectors
  repairable <- identical(status, "inertia_failed") && is.null(problem$metric) &&
    !is.null(V) && !is.null(completeness_target_kind(target))
  if (repairable) {
    kind <- completeness_target_kind(target)
    residuals <- cert$residuals
    round <- 0L
    while (identical(status, "inertia_failed") && round < controls$max_rounds) {
      round <- round + 1L
      margin <- completeness_margin(values, residuals, numeric(), tol, cert$scale)
      probe <- completeness_probe(op, values, V, kind, margin,
                                  block = controls$block, steps = controls$steps,
                                  stream = round)
      record$operator_columns <- record$operator_columns + probe$columns
      record$operator_block_calls <- record$operator_block_calls + probe$block_calls
      if (!length(probe$theta) ||
          !any(completeness_beyond(probe$theta, completeness_edge(values, kind),
                                   margin, kind))) {
        # The count proves a value is missing but the short probe did not
        # see it: solve the deflated complement from a deterministic start.
        probe <- list(theta = numeric(), Y = matrix(0, nrow(V), 0L), exhausted = FALSE)
      }
      fixed <- completeness_repair_round(op, values, V, target, kind, tol, probe, margin)
      record$operator_columns <- record$operator_columns + fixed$columns
      record$operator_block_calls <- record$operator_block_calls + fixed$block_calls
      values <- fixed$values
      V <- fixed$vectors
      AV <- as.matrix(apply_operator(op, V))
      residuals <- sqrt(colSums((AV - sweep(V, 2L, values, `*`))^2))
      record$operator_columns <- record$operator_columns + ncol(V)
      record$operator_block_calls <- record$operator_block_calls + 1L
      orth <- max(abs(crossprod(V) - diag(length(values))))
      again <- inertia_completeness_check(op, values, residuals, orth, target, ctx = ctx)
      record$factorizations <- record$factorizations + again$record$factorizations
      keep <- c("edge", "rho", "margin", "threshold_upper", "count_upper",
                "threshold_lower", "count_lower", "reason")
      record[keep] <- again$record[keep]
      status <- again$status
      record$repaired <- TRUE
    }
    record$rounds <- round
  }
  record$seconds <- proc.time()[["elapsed"]] - started
  list(
    status = status,
    record = record,
    repaired = isTRUE(record$repaired),
    values = values,
    vectors = V
  )
}
