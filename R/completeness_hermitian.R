# Target completeness for every Hermitian eigen route (standard, generalized
# SPD, complex Hermitian, matrix-free, constrained LOBPCG).
#
# `certificate$passed` means "the requested pairs, each accurate": the
# residual certificate proves each pair, and the returned SET must also be
# verified (R/solve.R, require_verified_completeness()). The original checks
# (R/target_completeness.R, R/completeness_inertia.R) cover standard real
# problems with an extremal target (largest, smallest, largest_magnitude).
# This file extends them to every other Hermitian target and route:
#
# Targets. A target is split into "parts", each the k_p most preferred
# eigenvalues for one preference function:
#   largest / smallest / largest_magnitude  one part, as before;
#   smallest_magnitude                       one part, nearest(0);
#   nearest(c)                               one part, |lambda - c|;
#   both_ends(k_low, k_high)                 two parts: the k_low smallest
#                                            returned values are checked as
#                                            "smallest", the rest as "largest".
#
# Proofs (status "inertia_verified" / "inertia_failed"). With an explicit
# matrix (and explicit SPD B) every part is counted by Sylvester inertia
# exactly as in R/completeness_inertia.R; both_ends parts are counted
# separately (each count is a proof for its own end, and the two ends are
# disjoint index sets because k_low + k_high < n). nearest() and
# smallest_magnitude() are now counted under completeness = "auto" whenever
# the cost gate passes. A matrix-free operator with n below
# getOption("eigencore.completeness_materialize_limit", 2000) is materialised
# by n operator applies (M = A I) when the gate allows the cost; the count is
# then exact for the symmetrised M, and the residuals used in the count are
# inflated by 2 delta, delta bounding ||M_sym - A|| (asymmetry of M plus an
# apply-rounding floor), so the verdict transfers to A by Weyl's inequality.
#
# Evidence (status "probed" / "repaired" / "failed"). Otherwise the
# deflated-complement block Lanczos probe of R/target_completeness.R runs in a
# "standard space" in which the problem is a Euclidean Hermitian one:
#   B = NULL: the operator itself;
#   explicit SPD B = R'R: C = R^{-T} A R^{-1} with Y = R X (B-orthonormal X
#     gives orthonormal Y), so Euclidean deflation of C is B-deflation of the
#     pencil;
#   LOBPCG constraints Q: the complement is taken against [Q, Y], i.e. the
#     probe runs on A compressed to the constraint complement, which is the
#     problem LOBPCG solved.
# Extremal parts probe the standard operator directly. Interior parts
# (nearest(c), smallest_magnitude) probe a transformed operator G whose
# extremal end is the target: the shift-invert operator (A - cI)^{-1} (kind
# largest_magnitude) when the shift-invert route handed its factorisation
# over, else the squared shift (A - cI)^2 (kind smallest). Every complement
# Ritz value of G is a Rayleigh quotient of the compression, so a value
# beyond the edge proves a more-preferred eigenvalue is missing; a clean
# probe is evidence, not proof. The squared probe converges more slowly
# (interior gaps are squared and divided by the squared spread), so it gets
# a longer step budget (eigencore.completeness_interior_probe_steps). A
# repair solves the deflated complement of G with native block Lanczos,
# merges with a Rayleigh-Ritz step on the standard operator selected by the
# problem's own target, maps back and re-certifies from scratch.
#
# Complex Hermitian results are checked on the real 2n embedding
# [Re -Im; Im Re] (R/complex_hermitian.R): the returned pairs are realified
# (every eigenvalue doubled) and counted or probed there, and a repaired set
# is mapped back by a complex Rayleigh-Ritz step.

#' @keywords internal
hermitian_completeness_parts <- function(target, values = NULL) {
  kind <- if (inherits(target, "eigencore_target")) target$kind else NA_character_
  all_idx <- seq_along(values)
  one <- function(t, probe, center = 0) {
    list(list(idx = all_idx, target = t, probe = probe, center = center))
  }
  switch(
    kind,
    largest = ,
    largest_real = one(largest(), "largest"),
    smallest = ,
    smallest_real = one(smallest(), "smallest"),
    largest_magnitude = one(largest_magnitude(), "largest_magnitude"),
    smallest_magnitude = one(nearest(0), "nearest", 0),
    nearest = {
      c0 <- target$value
      if (is.numeric(c0) && length(c0) == 1L && is.finite(c0)) {
        one(nearest(as.numeric(c0)), "nearest", as.numeric(c0))
      } else {
        NULL
      }
    },
    both_ends = {
      kl <- as.integer(target$value$k_low %||% NA_integer_)
      kh <- as.integer(target$value$k_high %||% NA_integer_)
      if (is.na(kl) || is.na(kh) || (!is.null(values) && length(values) != kl + kh)) {
        return(NULL)
      }
      low <- if (kl > 0L && length(values)) utils::head(order(Re(values)), kl) else integer()
      high <- setdiff(all_idx, low)
      parts <- list()
      if (kl > 0L) {
        parts[[length(parts) + 1L]] <- list(idx = low, target = smallest(),
                                            probe = "smallest", center = 0)
      }
      if (kh > 0L) {
        parts[[length(parts) + 1L]] <- list(idx = high, target = largest(),
                                            probe = "largest", center = 0)
      }
      parts
    },
    NULL
  )
}

#' @keywords internal
hermitian_completeness_applicable <- function(problem, k) {
  op <- problem$A
  inherits(op, "eigencore_operator") &&
    identical(op$structure$kind, "hermitian") &&
    (op$dtype %||% "double") %in% c("double", "complex") &&
    !is.null(hermitian_completeness_parts(problem$target)) &&
    k <= op$dim[[1L]]
}

#' @keywords internal
hermitian_completeness_constraints <- function(plan) {
  desc <- plan$method_descriptor
  if (inherits(desc, "eigencore_method") && identical(desc$kind, "lobpcg")) {
    desc$constraints
  } else {
    NULL
  }
}

# Does a route need its Ritz vectors kept internally for the check?
#' @keywords internal
hermitian_completeness_wants_vectors <- function(plan, problem, k, mode) {
  !identical(mode, "none") &&
    isTRUE(plan$execution$certify) &&
    hermitian_completeness_applicable(problem, k) &&
    !identical(target_completeness_route_class(plan, problem), "exact")
}

# Does this file handle the result? Every non-exact Hermitian route with a
# supported target and no verdict yet (or a route-internal "not_checked");
# full-spectrum routes keep the "exact" label from apply_target_completeness().
#' @keywords internal
hermitian_completeness_handles <- function(result, plan, problem, k, mode) {
  cert <- result$certificate
  status0 <- cert$target_completeness %||% NULL
  if (is.null(cert) || identical(mode, "none") ||
      !(is.null(status0) || identical(status0, "not_checked")) ||
      !hermitian_completeness_applicable(problem, k) ||
      identical(target_completeness_route_class(plan, problem), "exact")) {
    return(FALSE)
  }
  TRUE
}

#' @keywords internal
hermitian_completeness_controls <- function() {
  base <- target_completeness_controls()
  steps <- suppressWarnings(as.integer(
    getOption("eigencore.completeness_interior_probe_steps", 40L)))
  base$interior_steps <- if (length(steps) == 1L && !is.na(steps) && steps >= 1L) {
    steps
  } else {
    40L
  }
  limit <- suppressWarnings(as.integer(
    getOption("eigencore.completeness_materialize_limit", 2000L)))
  base$materialize_limit <- if (length(limit) == 1L && !is.na(limit)) limit else 2000L
  base
}

# ---------------------------------------------------------------------------
# Inertia (proof)
# ---------------------------------------------------------------------------

#' @keywords internal
hermitian_completeness_gate_target <- function(parts) {
  probes <- vapply(parts, function(p) p$probe, character(1L))
  if (any(probes == "nearest")) {
    return(parts[[which(probes == "nearest")[[1L]]]]$target)
  }
  if (length(parts) > 1L || any(probes == "largest_magnitude")) {
    # Two thresholds per count: the gate's cost for magnitude targets.
    return(largest_magnitude())
  }
  parts[[1L]]$target
}

#' @keywords internal
hermitian_completeness_threshold_count <- function(parts) {
  sum(vapply(parts, function(p) {
    if (p$probe %in% c("nearest", "largest_magnitude")) 2 else 1
  }, numeric(1L)))
}

# Materialise an operator by n applies (M = A I) for an exact count of a
# small matrix-free problem. Returns list(matrix, delta) with delta an upper
# bound of ||M_sym - A||_2 up to the apply's own rounding (the asymmetry of
# M plus a 4 n eps ||M||_1 floor), or NULL.
#' @keywords internal
hermitian_completeness_materialize <- function(op) {
  n <- as.integer(op$dim[[1L]])
  M <- tryCatch(as.matrix(apply_operator(op, diag(n))), error = function(e) NULL)
  if (is.null(M) || !identical(dim(M), c(n, n)) ||
      any(!is.finite(if (is.complex(M)) Mod(M) else M))) {
    return(NULL)
  }
  if (!is.complex(M)) {
    storage.mode(M) <- "double"
  }
  Mh <- if (is.complex(M)) Conj(t(M)) else t(M)
  asym <- 0.5 * sqrt(sum(Mod(M - Mh)^2))
  Msym <- (M + Mh) / 2
  norm1 <- max(colSums(Mod(Msym)), 0)
  list(matrix = Msym, delta = asym + 4 * n * .Machine$double.eps * norm1)
}

# Gate and context for the inertia certificate. Explicit sources use
# inertia_completeness_gate(); a matrix-free A (or B) below the
# materialisation limit is materialised when the predicted cost (n applies,
# at the solve's own per-column rate, plus the dense LDL' counts) fits the
# same budget. Returns list(use, reason, ctx, predicted_seconds, delta_A,
# delta_B, metric, columns).
#' @keywords internal
hermitian_completeness_inertia_gate <- function(problem, mode, k, parts,
                                                solve_seconds = NA_real_,
                                                seed = NULL, result = NULL,
                                                controls = hermitian_completeness_controls()) {
  no <- function(reason, predicted = NA_real_) {
    list(use = FALSE, reason = reason, ctx = NULL, predicted_seconds = predicted)
  }
  if (!mode %in% c("auto", "inertia")) {
    return(no("mode"))
  }
  metric <- problem$metric
  if (!is.null(metric) && !isTRUE(tryCatch(generalized_spd_metric_known(metric),
                                           error = function(e) FALSE))) {
    return(no("B not known to be positive definite"))
  }
  gate_problem <- problem
  gate_problem$target <- hermitian_completeness_gate_target(parts)
  A_explicit <- !is.null(inertia_matrix_of(problem$A))
  B_explicit <- is.null(metric) || !is.null(inertia_matrix_of(metric))
  if (A_explicit && B_explicit) {
    gate <- inertia_completeness_gate(gate_problem, mode, k,
                                      solve_seconds = solve_seconds, seed = seed)
    gate$delta_A <- 0
    gate$delta_B <- 0
    gate$metric <- metric
    gate$columns <- 0L
    return(gate)
  }
  n <- as.integer(problem$A$dim[[1L]])
  if (n > controls$materialize_limit) {
    return(no("no explicit matrix source (above the materialisation limit)"))
  }
  icontrols <- inertia_completeness_controls()
  cols_used <- as.numeric(result$operator_columns %||% result$matvecs %||% NA_real_)
  per_column <- if (is.finite(solve_seconds) && is.finite(cols_used) && cols_used > 0) {
    solve_seconds / cols_used
  } else {
    0
  }
  materialized <- (!A_explicit) + (!B_explicit)
  predicted <- materialized * n * per_column +
    hermitian_completeness_threshold_count(parts) * (n^3 / 3) / icontrols$dense_rate
  budget <- max(icontrols$seconds,
                icontrols$ratio * (if (is.finite(solve_seconds)) solve_seconds else 0))
  if (identical(mode, "auto") && !(predicted <= budget)) {
    return(no("predicted materialisation and factorisation cost exceeds the gate",
              predicted))
  }
  Am <- if (A_explicit) list(matrix = problem$A, delta = 0) else
    hermitian_completeness_materialize(problem$A)
  Bm <- if (is.null(metric)) NULL else if (B_explicit) list(matrix = metric, delta = 0) else
    hermitian_completeness_materialize(metric)
  if (is.null(Am) || (!is.null(metric) && is.null(Bm))) {
    return(no("operator could not be materialised", predicted))
  }
  ctx <- tryCatch(inertia_context(Am$matrix, if (is.null(Bm)) NULL else Bm$matrix),
                  error = function(e) e)
  if (inherits(ctx, "error")) {
    return(no(paste0("inertia context: ", conditionMessage(ctx)), predicted))
  }
  ctx$method <- paste0(ctx$method, " (materialised operator)")
  list(use = TRUE,
       reason = if (identical(mode, "inertia")) "requested (materialised operator)" else
         "cost gate passed (materialised operator)",
       ctx = ctx, predicted_seconds = predicted,
       delta_A = Am$delta, delta_B = if (is.null(Bm)) 0 else Bm$delta,
       metric = if (is.null(Bm)) NULL else Bm$matrix,
       columns = as.integer(materialized * n))
}

# Ties at the target edge. inertia_completeness_check() is inconclusive when
# more than k eigenvalues lie within its threshold window around the edge,
# which is always the case when the k-th and (k+1)-th eigenvalues are equal
# (a repeated eigenvalue straddling the edge), although then every choice of
# copies is a correct answer. Its window is wide (margin >= 1e-10 ||A||, to
# cover any factorisation backward error), so the counts are first repeated
# with the tight margin m = max(rho, 2 * the factorisations' own backward
# bound, 64 eps scale): this alone separates near-ties wider than ~4 rho,
# giving "inertia_verified" (N(E + rho + m) = k) or "inertia_failed"
# (a more-preferred eigenvalue is missing) exactly as in the original check.
# If the tight window still holds more than k eigenvalues, with f the
# preference distance, t_lo the lower threshold (N(t_lo) = l < k
# eigenvalues strictly more preferred than the edge cluster) and t_up the
# upper one (N(t_up) > k): a returned value with f(theta) + rho < t_lo is
# matched (Kahan: distinct eigenvalues within rho; f is 1-Lipschitz) to an
# eigenvalue counted in N(t_lo). If at least l returned values are that far
# inside, all l strictly-preferred eigenvalues are returned, and the other
# k - l returned values, like the true (l+1)-th..k-th preferred eigenvalues,
# are matched to eigenvalues in the window [t_lo, t_up] of width
# 2 (rho + m), the residual accuracy of the returned values. The returned set
# is then a correct top-k set up to the choice of eigenvalues inside that
# window: exact for a repeated eigenvalue, and within the residual accuracy
# otherwise. The verdict is "inertia_verified" with record$tie = TRUE and the
# window recorded.
#' @keywords internal
hermitian_completeness_tie <- function(check, values, part, ctx) {
  rec <- check$record
  if (!identical(check$status, "inertia_inconclusive") ||
      !grepl("not separated", rec$reason %||% "", fixed = TRUE) ||
      is.null(ctx) || !is.finite(rec$rho %||% NA_real_) ||
      !is.finite(rec$edge %||% NA_real_)) {
    return(check)
  }
  kind <- inertia_completeness_kind(part$target)
  center <- part$center
  k <- length(values)
  rho <- rec$rho
  E <- rec$edge
  eps <- .Machine$double.eps
  scale <- ctx$normA + (abs(center) + max(abs(values))) * ctx$normB
  margin <- max(rho, 64 * eps * scale)
  upper <- inertia_region_count(ctx, kind, E + rho + margin, center, outward = TRUE)
  factorizations <- upper$factorizations
  if (isTRUE(upper$reliable) && upper$backward_bound >= margin / 2) {
    margin <- 2 * upper$backward_bound
    upper <- inertia_region_count(ctx, kind, E + rho + margin, center, outward = TRUE)
    factorizations <- factorizations + upper$factorizations
  }
  lower <- if (isTRUE(upper$reliable)) {
    inertia_region_count(ctx, kind, E - rho - margin, center, outward = FALSE)
  } else {
    NULL
  }
  factorizations <- factorizations + (lower$factorizations %||% 0)
  rec$factorizations <- (rec$factorizations %||% 0) + factorizations
  rec$tight_margin <- margin
  if (!isTRUE(upper$reliable) || !isTRUE(lower$reliable) ||
      isTRUE(lower$backward_bound >= margin)) {
    return(list(status = check$status, record = rec))
  }
  f <- inertia_preference(values, kind, center)
  set_record <- function(rec) {
    rec$margin <- margin
    rec$threshold_upper <- upper$t
    rec$count_upper <- upper$count
    rec$threshold_lower <- lower$t
    rec$count_lower <- lower$count
    rec
  }
  if (upper$count == k) {
    rec <- set_record(rec)
    rec$reason <- "separated by the tight-margin recount"
    return(list(status = "inertia_verified", record = rec))
  }
  if (upper$count < k) {
    return(list(status = check$status, record = rec))
  }
  if (lower$count >= k || lower$count > sum(f - rho < lower$t)) {
    rec <- set_record(rec)
    rec$reason <- "the tight-margin recount shows a more-preferred eigenvalue is missing"
    return(list(status = "inertia_failed", record = rec))
  }
  inside <- sum(f + rho < lower$t)
  if (inside < lower$count) {
    return(list(status = check$status, record = rec))
  }
  rec <- set_record(rec)
  rec$tie <- TRUE
  rec$tie_width <- upper$t - lower$t
  rec$tie_eigenvalues <- upper$count - lower$count
  rec$reason <- sprintf(paste0(
    "complete up to a tie at the target edge: %d eigenvalue(s) lie in a ",
    "window of width %.3g (the residual accuracy) and %d of them were returned"),
    as.integer(rec$tie_eigenvalues), rec$tie_width,
    as.integer(k - lower$count))
  list(status = "inertia_verified", record = rec)
}

# Count every part. Returns list(status, record, failed_part).
#' @keywords internal
hermitian_completeness_count <- function(problem, gate, values, residuals,
                                         orthogonality, parts) {
  inflate <- 2 * (gate$delta_A %||% 0) + 2 * (gate$delta_B %||% 0) * abs(values)
  residuals <- abs(as.numeric(residuals)) + inflate
  checks <- lapply(parts, function(p) {
    check <- inertia_completeness_check(
      problem$A, values[p$idx], residuals[p$idx], orthogonality, p$target,
      Bop = gate$metric, ctx = gate$ctx
    )
    hermitian_completeness_tie(check, values[p$idx], p, gate$ctx)
  })
  statuses <- vapply(checks, function(ch) ch$status, character(1L))
  status <- if (all(statuses == "inertia_verified")) {
    "inertia_verified"
  } else if (any(statuses == "inertia_failed")) {
    "inertia_failed"
  } else {
    "inertia_inconclusive"
  }
  worst <- which(statuses == "inertia_failed")
  if (!length(worst)) worst <- which(statuses == "inertia_inconclusive")
  worst <- if (length(worst)) worst[[1L]] else 1L
  record <- checks[[worst]]$record
  record$factorizations <- sum(vapply(checks, function(ch) {
    as.numeric(ch$record$factorizations %||% 0)
  }, numeric(1L)))
  if (length(parts) > 1L) {
    record$parts <- lapply(seq_along(parts), function(i) {
      c(list(part = parts[[i]]$probe, status = statuses[[i]]), checks[[i]]$record)
    })
  }
  if (any(inflate > 0)) {
    record$materialized <- TRUE
    record$materialize_delta <- max(gate$delta_A %||% 0, gate$delta_B %||% 0)
  }
  list(status = status, record = record,
       failed_part = if (identical(statuses[[worst]], "inertia_failed")) worst else NA_integer_)
}

# ---------------------------------------------------------------------------
# Probe (evidence) in a Euclidean standard space
# ---------------------------------------------------------------------------

# B = R'R for an explicit (or materialisable) SPD metric: list(apply_R,
# solve_R, solve_Rt) or NULL.
#' @keywords internal
hermitian_completeness_metric_factor <- function(Bop, controls) {
  src <- inertia_matrix_of(Bop)
  if (is.null(src)) {
    if (as.integer(Bop$dim[[1L]]) > controls$materialize_limit) {
      return(NULL)
    }
    m <- hermitian_completeness_materialize(Bop)
    if (is.null(m) || is.complex(m$matrix)) {
      return(NULL)
    }
    src <- m$matrix
  }
  if (inherits(src, "diagonalMatrix")) {
    d <- as.numeric(Matrix::diag(src))
    if (any(!is.finite(d)) || any(d <= 0)) {
      return(NULL)
    }
    s <- sqrt(d)
    return(list(kind = "diagonal",
                apply_R = function(X) s * X,
                solve_R = function(Y) Y / s,
                solve_Rt = function(Y) Y / s))
  }
  if (inherits(src, "sparseMatrix")) {
    Bs <- tryCatch(methods::as(Matrix::forceSymmetric(
      methods::as(methods::as(src, "dMatrix"), "CsparseMatrix"), uplo = "U"),
      "CsparseMatrix"), error = function(e) NULL)
    FB <- if (is.null(Bs)) NULL else tryCatch(
      Matrix::Cholesky(Bs, LDL = FALSE, super = FALSE, perm = TRUE),
      error = function(e) NULL, warning = function(w) NULL)
    if (is.null(FB)) {
      return(NULL)
    }
    # B = P1' L L' P1, R = L' P1.
    L <- methods::as(Matrix::expand1(FB, "L"), "CsparseMatrix")
    P1 <- Matrix::expand1(FB, "P1")
    Lt <- Matrix::t(L)
    return(list(kind = "sparse",
                apply_R = function(X) as.matrix(Lt %*% (P1 %*% X)),
                solve_R = function(Y) as.matrix(Matrix::crossprod(P1, Matrix::solve(Lt, Y))),
                solve_Rt = function(Y) as.matrix(Matrix::solve(L, P1 %*% Y))))
  }
  B <- as.matrix(src)
  if (is.complex(B) || !is.numeric(B)) {
    return(NULL)
  }
  storage.mode(B) <- "double"
  R <- tryCatch(chol((B + t(B)) / 2), error = function(e) NULL)
  if (is.null(R)) {
    return(NULL)
  }
  list(kind = "dense",
       apply_R = function(X) R %*% X,
       solve_R = function(Y) backsolve(R, Y),
       solve_Rt = function(Y) backsolve(R, Y, transpose = TRUE))
}

# The Euclidean standard problem the probe runs on. Returns list(op, Y,
# fixed, to_original, space) or NULL when no standard space is available
# (complex operators, matrix-free B above the materialisation limit).
#' @keywords internal
hermitian_completeness_standard_space <- function(problem, vectors,
                                                  constraints = NULL,
                                                  controls = hermitian_completeness_controls()) {
  op <- problem$A
  if (!identical(op$dtype %||% "double", "double") || is.null(vectors) ||
      is.complex(vectors)) {
    return(NULL)
  }
  X <- as.matrix(vectors)
  storage.mode(X) <- "double"
  n <- as.integer(op$dim[[1L]])
  orth <- function(Q) {
    if (is.null(Q) || !ncol(Q)) {
      return(NULL)
    }
    Q <- qr.Q(qr(Q))
    Q
  }
  if (is.null(problem$metric)) {
    return(list(
      op = op, Y = X,
      fixed = orth(if (is.null(constraints)) NULL else as.matrix(constraints)),
      to_original = function(Y) Y,
      space = "standard"
    ))
  }
  fac <- hermitian_completeness_metric_factor(problem$metric, controls)
  if (is.null(fac)) {
    return(NULL)
  }
  C <- linear_operator(
    dim = c(n, n),
    apply = function(X, alpha = 1, beta = 0, Y = NULL) {
      Z <- as.matrix(X)
      out <- fac$solve_Rt(as.matrix(apply_operator(op, fac$solve_R(Z))))
      out <- alpha * as.matrix(out)
      if (is.null(Y) || beta == 0) out else out + beta * Y
    },
    structure = hermitian(),
    name = "eigencore completeness transformed pencil R^-T A R^-1"
  )
  list(
    op = C,
    Y = as.matrix(fac$apply_R(X)),
    fixed = orth(if (is.null(constraints)) NULL else
      as.matrix(fac$apply_R(as.matrix(constraints)))),
    to_original = function(Y) as.matrix(fac$solve_R(Y)),
    space = paste0("transformed_standard_problem (", fac$kind, " B = R'R)")
  )
}

# A shift-invert solve for repairing a nearest/smallest_magnitude part after
# an inertia count proved a value missing (standard problems with an explicit
# matrix, whose factorisation the count already showed affordable). The
# squared-shift operator converges too slowly for a reliable repair. A
# singular shift is perturbed; the recount decides the verdict either way.
#' @keywords internal
hermitian_completeness_repair_solve <- function(problem, parts) {
  if (!is.null(problem$metric) || is.null(inertia_matrix_of(problem$A))) {
    return(NULL)
  }
  near <- Filter(function(p) identical(p$probe, "nearest"), parts)
  if (!length(near)) {
    return(NULL)
  }
  center <- near[[1L]]$center
  scale <- tryCatch(operator_norm_for_certificate_info(problem$A)$value,
                    error = function(e) NA_real_)
  if (!is.finite(scale) || scale <= 0) {
    scale <- 1
  }
  sub <- list(A = problem$A, metric = NULL, structure = problem$A$structure)
  for (rel in c(0, 1e-9, -1e-7)) {
    shift <- center + rel * scale
    prep <- tryCatch(suppressWarnings(prepare_shift_invert_operator(sub, shift)),
                     error = function(e) NULL)
    if (!is.null(prep)) {
      return(list(sigma = center, shift = shift, operator = prep$operator))
    }
  }
  NULL
}

#' @keywords internal
hermitian_completeness_kind_target <- function(kind) {
  switch(kind,
         largest = largest(),
         smallest = smallest(),
         largest_magnitude = largest_magnitude())
}

# Operator whose extremal end is the part's target: list(G, map, kind, space).
#' @keywords internal
hermitian_completeness_part_operator <- function(std, part, values, solve_T = NULL) {
  if (!identical(part$probe, "nearest")) {
    return(list(G = std$op, map = function(v) v, kind = part$probe,
                space = "direct"))
  }
  center <- part$center
  # solve_T$shift is the shift actually factored; it differs from the
  # target's center only for a repair helper built after a count failure
  # (sigma perturbed off an exact eigenvalue), never for a probe verdict.
  shift <- solve_T$shift %||% solve_T$sigma
  if (!is.null(solve_T) && identical(std$space, "standard") &&
      isTRUE(abs(solve_T$sigma - center) <= 0) &&
      all(abs(values - shift) > 0)) {
    return(list(G = solve_T$operator, map = function(v) 1 / (v - shift),
                kind = "largest_magnitude", space = "shift_invert"))
  }
  base <- std$op
  n <- as.integer(base$dim[[1L]])
  G <- linear_operator(
    dim = c(n, n),
    apply = function(X, alpha = 1, beta = 0, Y = NULL) {
      Z <- as.matrix(X)
      S1 <- as.matrix(apply_operator(base, Z)) - center * Z
      S2 <- as.matrix(apply_operator(base, S1)) - center * S1
      S2 <- alpha * S2
      if (is.null(Y) || beta == 0) S2 else S2 + beta * Y
    },
    structure = hermitian(),
    name = "eigencore completeness squared shift (A - cI)^2"
  )
  list(G = G, map = function(v) (v - center)^2, kind = "smallest",
       space = "squared_shift")
}

# Probe margin at the residual accuracy of the returned pairs. Every
# complement Ritz value is a Rayleigh quotient of the compression of A to the
# complement of V, whose spectrum is the rest of A's spectrum up to the
# residual norm, so a value beyond the edge by more than twice the residual
# mass (plus rounding) is a missing eigenvalue. Unlike completeness_margin()
# there is no tolerance-times-norm cushion: a near-tie wider than the
# residual accuracy is repaired (the Rayleigh-Ritz step then selects the more
# preferred eigenvalue), consistent with the tight inertia count.
#' @keywords internal
hermitian_completeness_margin <- function(values, residuals, theta, tol,
                                          scale = NULL) {
  finite_abs <- function(x) {
    x <- abs(as.numeric(x))
    x[is.finite(x)]
  }
  scale <- max(c(finite_abs(values), finite_abs(theta), finite_abs(scale), 0))
  residual_mass <- sqrt(sum(finite_abs(residuals)^2))
  max(2 * residual_mass, 64 * .Machine$double.eps * scale)
}

# Probe margin in a transformed space G: the eigenvalue-space margin m
# (hermitian_completeness_margin() of the eigenvalue residuals)
# mapped through G's monotone transform of |lambda - c| at the edge distance
# d, so an interior probe is no looser than an extremal one, plus G's own
# residual mass.
#' @keywords internal
hermitian_completeness_mapped_margin <- function(pop, values, vG, rG, m) {
  center_dist <- if (identical(pop$space, "shift_invert")) {
    1 / abs(vG)
  } else {
    sqrt(abs(vG))
  }
  d <- max(center_dist)
  mapped <- if (identical(pop$space, "shift_invert")) {
    if (m < d) 1 / (d - m) - 1 / d else Inf
  } else {
    if (m < d) d^2 - (d - m)^2 else d^2
  }
  max(mapped, 2 * sqrt(sum(rG^2)), 64 * .Machine$double.eps * max(abs(vG), 0))
}

# One repair round for a part with an intruder (or a count failure): solve
# the deflated complement of G, then Rayleigh-Ritz on the standard operator
# over span(Y, X_c), selected by the problem's own target.
#' @keywords internal
hermitian_completeness_repair_round <- function(std, Y, pop, vG, probe, margin,
                                                tol, target, k) {
  n <- nrow(Y)
  basis <- cbind(std$fixed, Y)
  nc <- n - ncol(basis)
  if (nc < 1L) {
    return(NULL)
  }
  kind <- pop$kind
  ktarget <- hermitian_completeness_kind_target(kind)
  columns <- 0L
  block_calls <- 0L
  if (isTRUE(probe$exhausted) && length(probe$theta)) {
    Xc <- probe$Y
  } else {
    edge <- completeness_edge(vG, kind)
    intruding <- if (length(probe$theta)) {
      which(completeness_beyond(probe$theta, edge, margin, kind))
    } else {
      integer()
    }
    kc <- min(length(vG), nc)
    block <- min(max(2L, length(intruding) + 1L), 8L, nc)
    m_max <- min(n, default_block_lanczos_max_subspace(kc, block))
    if (m_max < kc + block) {
      kc <- max(1L, m_max - block)
    }
    if (m_max < kc + block) {
      return(NULL)
    }
    start_cols <- if (length(probe$theta)) {
      probe$Y[, utils::head(order_indices(probe$theta, ktarget), block), drop = FALSE]
    } else {
      matrix(0, n, 0L)
    }
    if (ncol(start_cols) < block) {
      start_cols <- cbind(start_cols,
                          completeness_probe_start(n, block - ncol(start_cols), 977L))
    }
    scale <- max(abs(c(vG, probe$theta)), 0)
    Gdef <- completeness_deflated_operator(pop$G, basis,
                                           completeness_shift(vG, kind, scale))
    comp <- tryCatch(native_block_lanczos_hermitian(
      Gdef, k = kc, target = ktarget, tol = tol, maxit = m_max, block = block,
      max_restarts = 100L, vectors = TRUE, full_subspace = FALSE,
      certificate_fallback = FALSE, start = start_cols
    ), error = function(e) NULL)
    if (is.null(comp) || is.null(comp$vectors)) {
      return(NULL)
    }
    Xc <- comp$vectors
    columns <- as.integer(comp$operator_columns %||% comp$matvecs %||% 0L)
    block_calls <- as.integer(comp$operator_block_calls %||% comp$matvecs %||% 0L)
  }
  Qc <- completeness_orthonormalize(Xc, basis)
  if (!ncol(Qc)) {
    return(NULL)
  }
  Z <- cbind(Y, Qc)
  AZ <- as.matrix(apply_operator(std$op, Z))
  H <- crossprod(Z, AZ)
  H <- (H + t(H)) / 2
  eig <- eigen(H, symmetric = TRUE)
  sel <- utils::head(order_indices(eig$values, target), k)
  list(
    values = eig$values[sel],
    Y = Z %*% eig$vectors[, sel, drop = FALSE],
    columns = columns + ncol(Z),
    block_calls = block_calls + 1L
  )
}

# Probe every part (repairing up to controls$max_rounds). `force_part`
# forces one repair round of that part even when its probe sees no intruder
# (an inertia count already proved a value missing).
#' @keywords internal
hermitian_completeness_probe_check <- function(std, values, target, tol,
                                               solve_T = NULL, norm_scale = NULL,
                                               force_part = NA_integer_,
                                               controls = hermitian_completeness_controls()) {
  Y <- std$Y
  k <- length(values)
  record <- list(
    method = "deflated_complement_probe",
    space = std$space,
    block = controls$block,
    max_steps = controls$steps,
    steps = 0L,
    operator_columns = 0L,
    operator_block_calls = 0L,
    rounds = 0L,
    edge = NA_real_,
    most_preferred_complement = NA_real_,
    margin = NA_real_,
    intruder_found = FALSE,
    exhausted = FALSE,
    probe_operators = character()
  )
  repaired <- FALSE
  round <- 0L
  status <- "probed"
  repeat {
    parts <- hermitian_completeness_parts(target, values)
    basis <- cbind(std$fixed, Y)
    basis <- qr.Q(qr(basis))
    hit <- NULL
    for (p in seq_along(parts)) {
      part <- parts[[p]]
      pop <- hermitian_completeness_part_operator(std, part, values[part$idx], solve_T)
      record$probe_operators <- unique(c(record$probe_operators, pop$space))
      vG <- pop$map(values[part$idx])
      Yp <- Y[, part$idx, drop = FALSE]
      GY <- as.matrix(apply_operator(pop$G, Yp))
      rG <- sqrt(colSums((GY - sweep(Yp, 2L, vG, `*`))^2))
      record$operator_columns <- record$operator_columns + ncol(Yp)
      record$operator_block_calls <- record$operator_block_calls + 1L
      scale <- if (identical(pop$space, "direct")) norm_scale else NULL
      g_margin <- if (identical(pop$space, "direct")) {
        NULL
      } else {
        # The eigenvalue-space margin of the direct probe, mapped into G.
        AY <- as.matrix(apply_operator(std$op, Yp))
        r_lambda <- sqrt(colSums((AY - sweep(Yp, 2L, values[part$idx], `*`))^2))
        record$operator_columns <- record$operator_columns + ncol(Yp)
        record$operator_block_calls <- record$operator_block_calls + 1L
        hermitian_completeness_mapped_margin(
          pop, values[part$idx], vG, rG,
          hermitian_completeness_margin(values[part$idx], r_lambda, numeric(), tol, norm_scale)
        )
      }
      margin <- g_margin %||% hermitian_completeness_margin(vG, rG, numeric(), tol, scale)
      steps <- if (identical(pop$space, "squared_shift")) {
        max(controls$steps, controls$interior_steps)
      } else {
        controls$steps
      }
      probe <- completeness_probe(pop$G, vG, basis, pop$kind, margin,
                                  block = controls$block, steps = steps,
                                  stream = round + 1000L * (p - 1L))
      margin <- g_margin %||% hermitian_completeness_margin(vG, rG, probe$theta, tol, scale)
      record$steps <- record$steps + probe$steps
      record$operator_columns <- record$operator_columns + probe$columns
      record$operator_block_calls <- record$operator_block_calls + probe$block_calls
      edge <- completeness_edge(vG, pop$kind)
      record$edge <- edge
      record$margin <- margin
      record$most_preferred_complement <- completeness_most_preferred(probe$theta, pop$kind)
      record$exhausted <- isTRUE(probe$exhausted)
      intruder <- length(probe$theta) &&
        any(completeness_beyond(probe$theta, edge, margin, pop$kind))
      forced <- !intruder && round == 0L && !is.na(force_part) && p == force_part
      if ((intruder || forced) && is.null(hit)) {
        if (intruder) {
          record$intruder_found <- TRUE
        }
        hit <- list(pop = pop, vG = vG, margin = margin,
                    probe = if (intruder) probe else
                      list(theta = numeric(), Y = matrix(0, nrow(Y), 0L),
                           exhausted = FALSE))
      }
    }
    if (is.null(hit)) {
      status <- if (repaired) "repaired" else "probed"
      break
    }
    if (round >= controls$max_rounds) {
      status <- "failed"
      break
    }
    round <- round + 1L
    fixed <- hermitian_completeness_repair_round(std, Y, hit$pop, hit$vG, hit$probe,
                                                 hit$margin, tol, target, k)
    if (is.null(fixed)) {
      status <- "failed"
      break
    }
    record$operator_columns <- record$operator_columns + fixed$columns
    record$operator_block_calls <- record$operator_block_calls + fixed$block_calls
    values <- fixed$values
    Y <- fixed$Y
    repaired <- TRUE
  }
  record$rounds <- round
  list(status = status, record = record, repaired = repaired,
       values = values, Y = Y)
}

# ---------------------------------------------------------------------------
# Driver (hooked from apply_target_completeness())
# ---------------------------------------------------------------------------

# Install a (possibly repaired) set and its completeness verdict.
#' @keywords internal
hermitian_completeness_finish <- function(result, problem, plan, status, record,
                                          values = NULL, vectors = NULL,
                                          repaired = FALSE) {
  record$operator_columns <- as.integer(record$operator_columns %||% 0L)
  record$operator_block_calls <- as.integer(record$operator_block_calls %||% 0L)
  if (isTRUE(repaired)) {
    cert <- result$certificate
    tol <- cert$tolerance %||% plan$execution$tol
    new_cert <- certify_eigen_operator(problem$A, values, vectors,
                                       Bop = problem$metric, tol = tol)
    new_cert$notes <- unique(c(cert$notes, new_cert$notes))
    result$certificate <- new_cert
    result$values <- values
    result$vectors <- vectors
    result$residuals <- new_cert$residuals
    result$backward_error <- new_cert$backward_error
    result$orthogonality <- new_cert$orthogonality
    result$nconv <- sum(new_cert$converged)
    result$locked <- which(new_cert$converged)
    record$operator_columns <- record$operator_columns + length(values)
    record$operator_block_calls <- record$operator_block_calls + 1L
    result$certification_operator_columns <- as.integer(
      (result$certification_operator_columns %||% 0L) + length(values)
    )
  }
  result_with_completeness(
    result, list(status = status, record = record, repaired = FALSE), problem, plan
  )
}

#' @keywords internal
hermitian_target_completeness <- function(result, plan, problem, k, mode,
                                          vectors_requested,
                                          solve_seconds = NA_real_, seed = NULL) {
  if (!hermitian_completeness_handles(result, plan, problem, k, mode)) {
    return(NULL)
  }
  cert <- result$certificate
  values <- result$values
  if (!isTRUE(plan$execution$certify) || !isTRUE(cert$passed) ||
      !is.numeric(values) || is.complex(values) || length(values) != k ||
      any(!is.finite(values))) {
    return(NULL)
  }
  if (k >= problem$A$dim[[1L]]) {
    # n certified pairs with (near-)orthonormal vectors match n distinct
    # eigenvalues (Kahan): the returned set is the whole spectrum.
    result$certificate <- certificate_with_completeness(
      cert, "exact", list(method = "full_spectrum",
                          reason = "all n eigenpairs returned and certified")
    )
    if (!isTRUE(vectors_requested)) {
      result["vectors"] <- list(NULL)
    }
    return(result)
  }
  if (is.null(result$vectors)) {
    return(NULL)
  }
  if (complex_hermitian_problem(problem)) {
    # Complex results are checked on the real 2n embedding
    # (R/complex_hermitian.R).
    return(complex_hermitian_target_completeness(
      result, plan, problem, k, mode, vectors_requested,
      solve_seconds = solve_seconds
    ))
  }
  drop_vectors <- function(res) {
    if (!isTRUE(vectors_requested)) {
      res["vectors"] <- list(NULL)
    }
    res
  }
  controls <- hermitian_completeness_controls()
  target <- problem$target
  parts <- hermitian_completeness_parts(target, values)
  constraints <- hermitian_completeness_constraints(plan)
  tol <- cert$tolerance %||% plan$execution$tol
  solve_T <- result$transform$completeness_solve %||% NULL
  started <- proc.time()[["elapsed"]]

  gate <- NULL
  if (is.null(constraints)) {
    gate <- hermitian_completeness_inertia_gate(
      problem, mode, k, parts, solve_seconds = solve_seconds, seed = seed,
      result = result, controls = controls
    )
  }
  if (isTRUE(gate$use)) {
    count <- hermitian_completeness_count(problem, gate, values, cert$residuals,
                                          cert$orthogonality, parts)
    record <- count$record
    record$gate <- gate$reason
    record$predicted_seconds <- gate$predicted_seconds
    record$reused_symbolic <- isTRUE(gate$ctx$seeded_symbolic)
    record$repaired <- FALSE
    record$operator_columns <- gate$columns %||% 0L
    record$operator_block_calls <- as.integer((gate$columns %||% 0L) > 0L)
    status <- count$status
    out_values <- values
    out_vectors <- result$vectors
    repaired <- FALSE
    if (identical(status, "inertia_failed")) {
      std <- hermitian_completeness_standard_space(problem, result$vectors,
                                                   controls = controls)
      if (!is.null(std)) {
        repair_T <- solve_T %||% hermitian_completeness_repair_solve(problem, parts)
        pc <- hermitian_completeness_probe_check(
          std, values, target, tol, solve_T = repair_T, norm_scale = cert$scale,
          force_part = count$failed_part, controls = controls
        )
        record$probe <- pc$record
        record$operator_columns <- record$operator_columns + pc$record$operator_columns
        record$operator_block_calls <- record$operator_block_calls +
          pc$record$operator_block_calls
        if (isTRUE(pc$repaired)) {
          out_values <- pc$values
          out_vectors <- std$to_original(pc$Y)
          new_cert <- certify_eigen_operator(problem$A, out_values, out_vectors,
                                             Bop = problem$metric, tol = tol)
          again <- hermitian_completeness_count(
            problem, gate, out_values, new_cert$residuals, new_cert$orthogonality,
            hermitian_completeness_parts(target, out_values)
          )
          keep <- c("edge", "rho", "margin", "threshold_upper", "count_upper",
                    "threshold_lower", "count_lower", "reason", "parts",
                    "tie", "tie_width", "tie_eigenvalues")
          for (field in keep) {
            record[field] <- list(again$record[[field]])
          }
          record <- record[!vapply(record, is.null, logical(1L))]
          record$factorizations <- (record$factorizations %||% 0) +
            (again$record$factorizations %||% 0)
          status <- again$status
          record$repaired <- TRUE
          repaired <- TRUE
        }
      }
    }
    record$seconds <- proc.time()[["elapsed"]] - started
    result <- hermitian_completeness_finish(result, problem, plan, status, record,
                                            out_values, out_vectors, repaired)
    return(drop_vectors(result))
  }

  std <- if (mode %in% c("auto", "inertia", "probe")) {
    hermitian_completeness_standard_space(problem, result$vectors, constraints,
                                          controls = controls)
  } else {
    NULL
  }
  if (is.null(std)) {
    record <- list(reason = "no inertia count and no probe space for this route")
    if (!is.null(gate)) {
      record$inertia_gate <- gate$reason
      record$inertia_predicted_seconds <- gate$predicted_seconds
    }
    result$certificate <- certificate_with_completeness(cert, "not_checked", record)
    return(drop_vectors(result))
  }
  pc <- hermitian_completeness_probe_check(std, values, target, tol,
                                           solve_T = solve_T,
                                           norm_scale = cert$scale,
                                           controls = controls)
  record <- pc$record
  record$seconds <- proc.time()[["elapsed"]] - started
  record$constrained <- !is.null(constraints)
  if (!is.null(gate)) {
    record$inertia_gate <- gate$reason
    record$inertia_predicted_seconds <- gate$predicted_seconds
  }
  result <- hermitian_completeness_finish(
    result, problem, plan, pc$status, record,
    values = pc$values,
    vectors = if (isTRUE(pc$repaired)) std$to_original(pc$Y) else NULL,
    repaired = isTRUE(pc$repaired)
  )
  drop_vectors(result)
}
